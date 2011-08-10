/*
 * $Id: template.c,v 1.4 2001/07/23 15:28:29 lindahl Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Gyas ROwers Mature At Cryogenic Speed
 */

/* This line is only for CVS version info */
static char *SRCID_template_c = "$Id: template.c,v 1.4 2001/07/23 15:28:29 lindahl Exp $";

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "smalloc.h"
#include "macros.h"
#include "math.h"
#include "statutil.h"
#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "copyrite.h"
#include "statutil.h"
#include "tpxio.h"
#include "index.h"
#include "futil.h"
#include "pbc.h"
#include "rmpbc.h"
   double distPBC(double,double,double);
/* runs voronoi program outputs vertices */
   void readVoronoi(char*);
/* keeps track of number of vertices*/
   void typeIncrement(char);
/* truncates a line by shifting vertex onto the boundary of the surface*/
   void fitVertex(int,int,int**,double**,double**);
/* similar to fitVertex*/
   void genVertex(int,int,int**,double**,double**);

   void freeMD(int** marray, int rows);
   void freeMD_d(double** marray, int rows);
   void sfreeMD(int** marray, int rows,int cols);

   int findVertex(double,double, double**);

/*delete element from array*/
   int sliceVerArray(int**,int,int,int);

   int pointInPoly(double,double,double,double);


   int** tri_shapes;
   int* tri_shape_sizes;
   double** vertices_good;
   int vert_good_memsize = 0;
   int** triangles;
   
   int **edges=NULL;
   double **vertices=NULL;
   double **lines=NULL;
   double **pts=NULL;

   int halfGroupSize = 0;   
   int nTriangles = -1;
   int badVertList[200];
   int lastBadVer=-1;
  int lastEdge = -1;
   int lastVer = -1;
   int lastVerAll = -1;
   int lastLine = -1;
   int lastPt = -1;
   int lastArea = -1;
   int triShapeCount = 0;
   
   int k = 0;
   
   int j = 0;

   double min_x = 0;
   double min_y = 0;
   double max_x = 7;
   double max_y = 7;
   double* areas;
   double cutoff = 0.001;
   int* verAssoc;

char* tr_loc = NULL;
char* voro_loc = NULL;

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "This program does Voronoi analysis on a lipid bilayer.  ",
    "It can calculate area per lipid and distances from the  ",
    "surface of the bilayer. The index file should contain  ",
    "a group containing one atom from each lipid residue -  ",
    "a phosphate generally works well.  "
  };
  
  static int n=1;

  double* z_pts;
  /* Extra arguments - but note how you always get the begin/end
   * options when running the program, without mentioning them here!
   */
  
    real g_w = 1;
    real g_h = 1;   /* parameters for grid sizing */

  t_pargs pa[] = {
    { "-n1", FALSE, etINT, {&n},
      "Plot data for atom number n (starting on 1)"
    },
    { "-xgrid", FALSE, etREAL, {&g_w},
      "Grid cell width (nm)"
    },
    { "-ygrid", FALSE, etREAL, {&g_h},
      "Grid cell height (nm)"
    }
  };
  #define NPA asize(pa) 
  t_topology* top=NULL;
  char       title[STRLEN];
  t_trxframe fr;
  rvec       *x=NULL;
  matrix     box;
  t_pbc pbc;
  int        status;
  int        flags = TRX_READ_X;
  /* index stuff */
  int     ngrps;     /* the number of index groups */
  atom_id **index,max;   /* the index for the atom numbers */
  int     *isize;    /* the size of each group */
  char    **grpname; /* the name of each group */
  rvec    *com;
  int natoms;
  real t;

/* look for voronoi and triangle programs */
if ((tr_loc=getenv("TRIANG"))==NULL)
    tr_loc="~/triangle";
if (!fexist(tr_loc))
    gmx_fatal(FARGS,"Triangle not found in %s (use setenv TRIANG)",tr_loc);
if ((voro_loc=getenv("VORON"))==NULL)
    voro_loc="~/voronoi";
if (!fexist(voro_loc))
    gmx_fatal(FARGS,"Voronoi (Fortune) not found in %s (use setenv VORON)",voro_loc);
   verAssoc = malloc(sizeof(int)*500);
/*   edges = malloc(500*sizeof(int));
   vertices = malloc(500*sizeof(double));
   lines = malloc(500*sizeof(double));
   pts = malloc(500*sizeof(double));
*/

  t_filenm fnm[] = {
    { efTPX,  "-s",  NULL, ffREAD },   /* this is for the topology */
    { efTRX, "-f", NULL, ffREAD },      /* and this for the trajectory */
    { efNDX, "-n", NULL, ffOPTRD }
  };
  
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);

  /* This is the routine responsible for adding default options,
   * calling the X/motif interface, etc. */
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_BE_NICE,
		    NFILE,fnm,NPA,pa,asize(desc),desc,0,NULL);

  /* read topology */
  /*top=read_top(ftp2fn(efTPS,NFILE,fnm));
*/
  /* read index files */
  /* may need to change # groups */
  /* We don't need any topology information to write the coordinates,
   * but to show how it works we start by writing the name and
   * charge of the selected atom. It returns a boolean telling us
   * whether the topology was found and could be read
   */
 /*
  natoms=read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
 t_topology *top=NULL;
 */ 
 /* read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,top,&xtop,NULL,box,TRUE);
  sfree(xtop);
 */
  top=read_top(ftp2fn(efTPX,NFILE,fnm));
  /* 2 groups is ok here... want protein as #2, and lipid/whatever as #1 */
  ngrps = 2;
  snew(com,ngrps);
  snew(grpname,ngrps);
  snew(index,ngrps);
  snew(isize,ngrps);
  /* read in index file! */
  get_index(&top->atoms,ftp2fn(efNDX,NFILE,fnm),ngrps,isize,index,grpname);
  n=n-1; 
/* Our enumeration started on 1, but C starts from 0 */
  /* check that this atom exists */
  if(n<0 || n>(top->atoms.nr)) 
  {
    printf("Error: Atom number %d is out of range.\n",n);
    exit(1);
  }
  
/*  printf("Atom name: %s\n",*(top.atoms.atomname[n]));
  printf("Atom charge: %f\n",top.atoms.atom[n].q);
 0 and 1 group names  */
  /* third group is SOLVENT */
  FILE* fareas;
  fareas = fopen("areas.out","w");
  /* The first time we read data is a little special */
  /* read_first_frame(&status,ftp2fn(efTRX,NFILE,fnm),&fr,flags); */
  natoms = read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
 int f = 0;
 real mbox = 0;
 /* half of the group is on either side, with 10% "overdraft" */
 halfGroupSize = (int)floor((isize[0]/2)*1.10);
 double a_z = 0;
 for (f=0;(f<isize[0]);f++)  {
   /* calculate center of mass in z */
   a_z = x[index[0][f]][ZZ];
   if (a_z<0)
	a_z+=box[ZZ][ZZ];
   mbox+=a_z; /* x[index[0][f]][ZZ]; */
} 
 mbox = mbox / isize[0]; 
  printf("Center of mass of group 0 is %1.3f, half Z is %1.3f\n",mbox,0.5*box[ZZ][ZZ]);
  /* This is the main loop over frames */
  rvec dx;
  real dist2;
   FILE* mout1=NULL;
   FILE* mout2=NULL;
   FILE* ngrpcoords=NULL;
   ngrpcoords = fopen("grp_coords.out","w");
   char filename1[50]; /* = "coords_top.txt"; */
   char filename2[50]; /* = "coords_top.txt"; */
    int sel_atoms[isize[0]];
    int sel_atoms2[isize[0]];
    int n_selatoms = 0;
    int n_selatoms2 = 0;
    int i2 = 0;
    double* z_dists=NULL;
/*     z_dists = malloc(sizeof(double)*isize[1]);
 */
   snew(z_dists,isize[1]); 
   for (i2=0;i2<isize[1];i2++)
	{
		z_dists[i2]=0;
	}
   int fShape = 0;
   int fSh = 0;
   double d_xz = 0;
   double nBox = 0;
   double mDist = 0;
   do {
    /* coordinates are available in the vector fr.x
     * you can find this and all other structures in
     * the types directory under the gromacs include dir.
     * Note how flags determines wheter to read x/v/f!
     */
    /* let group 0 be the "slab" - heavy atoms only*/
    /* initialisation for correct distance calculations */
    set_pbc(&pbc,box);
    /* make molecules whole again */
    rm_pbc(&top->idef,natoms,box,x,x); 
/* recalculate COM every 100 steps for posterity */
if ((((int)floor(t))%100)==0)
{
 for (f=0;(f<isize[0]);f++)  {
   /* calculate center of mass in z */
   a_z = x[index[0][f]][ZZ];
   if (a_z<0)
	a_z+=box[ZZ][ZZ];
   mbox+=a_z; /* x[index[0][f]][ZZ]; */
} 
 mbox = mbox / isize[0]; 
  printf("Center of mass of group 0 is %1.3f, half Z is %1.3f\n",mbox,0.5*box[ZZ][ZZ]);
}
 /*  sprintf(filename1,"coords_%4.0f_top.txt",fr.time);
   printf("%s\n",*filename1);
 */ 
   sprintf(filename1, "%.0f_top.txt",t);
   sprintf(filename2, "%.0f_bottom.txt",t);
  mout1 = fopen(filename1,"w");
  mout2 = fopen(filename2,"w");
/* for now, simple - can improve */
/*    mbox = fr.box[ZZ][ZZ] * 0.5;
 */
	n_selatoms = 0;	
	n_selatoms2 = 0;
    real f_nx;
    real f_ny;
    int isTop = -1;
    /* group index 2, 3rd one - solvent */
    real tmp_dx = 0;
    real tmp_z = 0;
    int gridX = 0;  /* parameters for grid cells for atom */
    int gridY = 0;
    
    double** minX;
    int** minX_id;
    double** maxX;
    int** maxX_id;
 
    int n_cells_w = floor(box[XX][XX]/g_w)+1;
    int n_cells_h = floor(box[YY][YY]/g_h)+1;
    minX = malloc(n_cells_w*sizeof(double));
    minX_id = malloc(n_cells_w*sizeof(int));
    maxX = malloc(n_cells_w*sizeof(double));
    maxX_id = malloc(n_cells_w*sizeof(int));
    /* allocate the grid */ 
    for (j=0;j<n_cells_w;j++)
    {
	    minX[j] = malloc(n_cells_h*sizeof(double));
	    minX_id[j] = malloc(n_cells_h*sizeof(int));
	    maxX[j] = malloc(n_cells_h*sizeof(double));
	    maxX_id[j] = malloc(n_cells_h*sizeof(int));
	for (k=0;k<n_cells_h;k++)
	{
		maxX[j][k] = 0;
		minX_id[j][k] = 0;
		minX[j][k] = box[ZZ][ZZ];
		maxX_id[j][k] = 0;
	}
    }
    for (i2=0;(i2<isize[0]);i2++)  {
    tmp_z = x[index[0][i2]][ZZ];
 /*   if (tmp_z<0)
	tmp_z += box[ZZ][ZZ];
    if (tmp_z>box[ZZ][ZZ])
	tmp_z = tmp_z - box[ZZ][ZZ]; 
*/ 
     /* test if tmp_z is mbox */ 
    /* which grid cell is this part of */
    double tmpX = x[index[0][i2]][XX];
    double tmpY = x[index[0][i2]][YY];
    if (tmpX<0)	
	tmpX+=box[XX][XX];
    if (tmpY<0)	
	tmpY+=box[YY][YY];
    if (tmpX>box[XX][XX])
	tmpX=tmpX-box[XX][XX];
    if (tmpY>box[YY][YY])
	tmpY=tmpY-box[YY][YY];
    gridX = floor(tmpX/g_w);
    gridY = floor(tmpY/g_h);
	if (tmp_z<minX[gridX][gridY])
	{
		minX_id[gridX][gridY] = index[0][i2];
		minX[gridX][gridY] = tmp_z;
	}
	if (tmp_z>maxX[gridX][gridY])
	{
		maxX_id[gridX][gridY] = index[0][i2];
		maxX[gridX][gridY] = tmp_z;
	}
}
    for (j=0;j<n_cells_w;j++)
	{
		for (k=0;k<n_cells_h;k++)
		{
			if (minX_id[j][k]!=0)
			{
			sel_atoms[n_selatoms] = minX_id[j][k];
  		 	 f_nx = x[sel_atoms[n_selatoms]][XX];
    if (f_nx<0)
	f_nx=f_nx+box[XX][XX];
    if (f_nx>box[XX][XX])
	f_nx=f_nx-box[XX][XX];
    f_ny = x[sel_atoms[n_selatoms]][YY];
    if (f_ny<0)
        f_ny=f_ny+box[YY][YY];
    if (f_ny>box[YY][YY])
        f_ny=f_ny-box[YY][YY];
    fprintf(mout1,"%2.3f %2.3f\n",f_nx,f_ny);
			n_selatoms++;			
			}
			if (maxX_id[j][k]!=0)
			{
			sel_atoms2[n_selatoms2] = maxX_id[j][k];
    f_nx = x[sel_atoms2[n_selatoms2]][XX];
    if (f_nx<0)
	f_nx=f_nx+box[XX][XX];
    if (f_nx>box[XX][XX])
	f_nx=f_nx-box[XX][XX];
    f_ny = x[sel_atoms2[n_selatoms2]][YY];
    if (f_ny<0)
        f_ny=f_ny+box[YY][YY];
    if (f_ny>box[YY][YY])
        f_ny=f_ny-box[YY][YY];
    fprintf(mout2,"%8.3f %8.3f\n",f_nx, f_ny);
			n_selatoms2++;	
		}		
	} 
}
freeMD_d(minX,n_cells_w);
freeMD(minX_id,n_cells_w);
freeMD_d(maxX,n_cells_w);
freeMD(maxX_id,n_cells_w);
fclose(mout1);
fclose(mout2);
/* printf("\n"); */
   /* tesselation goes here */
   /* produce output file with 2 columns of XY coords */
   /* separate top and bottom  - two arrays */
   /* run voronoi */
/* printf("%d %d\n",n_selatoms,n_selatoms2);  */
   
   char picNewName[100];
   /*sprintf(commandString,"/Users/ananikolic/Desktop/fortune/parse_matrix_3.pl mout1.tmp %2.2f %2.2f",box[XX][XX],box[YY][YY]);
   system(commandString);
   */
     max_x = box[XX][XX];
     max_y = box[YY][YY];
   /*  free(vertices);
     free(lines); */
	 /*areas = malloc(sizeof(double)*(n_selatoms2+n_selatoms));*/
	 snew(areas,n_selatoms2+n_selatoms);
	 lastArea = -1;
     readVoronoi(filename1);
/* point locate the 128 residues to a "shape" id */
z_pts = malloc(sizeof(double)*triShapeCount);
for (i2=0;i2<n_selatoms;i2++)
{
/*	printf("%f %f\n",x[sel_atoms[i2]][XX],x[sel_atoms[i2]][YY]); */
      fShape = pointInPoly(x[sel_atoms[i2]][XX],x[sel_atoms[i2]][YY],box[XX][XX],box[YY][YY]);
	z_pts[fShape] = x[sel_atoms[i2]][ZZ];
}
/* analyse stuff */
 for (i2=0;(i2<isize[1]);i2++)  {
isTop = -1;
tmp_z = x[index[1][i2]][ZZ];
if (tmp_z<0)
	tmp_z += box[ZZ][ZZ];
if (tmp_z>box[ZZ][ZZ])
	tmp_z = tmp_z - box[ZZ][ZZ];
/* tmp_dx, tmp_z  */
  /*  d_xz = x[index[1][i2]][ZZ] - mbox;
    nBox = fabs(d_xz);
    if (((nBox<(0.5*box[ZZ][ZZ]))&&(x[index[1][i2]][ZZ]<mbox))||((nBox>(0.5*box[ZZ][ZZ]))&&(x[index[1][i2]][ZZ]>mbox)))
    {*/
	fSh = pointInPoly(x[index[1][i2]][XX],x[index[1][i2]][YY],box[XX][XX],box[YY][YY]);
tmp_dx = fabs(mbox - tmp_z); /* test with COM of layer */
	if ((tmp_dx<(0.5*box[ZZ][ZZ]))&&(tmp_z < mbox))
	{
		isTop = 1;
	}
	else if (tmp_z>mbox)
	{
	 if ((tmp_dx>(0.5*box[ZZ][ZZ])))
		isTop = 1;
	}
	/* difference in distance is ... */
if (isTop == 1)
{
/*	fprintf(stderr,"top %f %f %f\n",tmp_z,mbox,0.5*box[ZZ][ZZ]); */
	mDist = distPBC(z_pts[fSh],x[index[1][i2]][ZZ],box[ZZ][ZZ]);
        if (x[index[1][i2]][ZZ]>z_pts[fSh])
	{
		mDist = mDist * -1;
	}
}
/*         if (nBox>(0.5*box[ZZ][ZZ]))
                mDist = box[ZZ][ZZ]-mDist;
*/
	z_dists[i2] = mDist;
  /*  } */ 
}
/* free(tri_shapes);
        free(tri_shape_sizes); */
	    freeMD(tri_shapes,triShapeCount); 
		free(tri_shape_sizes); 
  /* PROBLEM: CRASH IN SOME SYSTEMS HERE!! */
		freeMD_d(vertices_good,vert_good_memsize);  
        freeMD(triangles,nTriangles);   
/*free(triangles); */
	/* read all in index 2 */
	triShapeCount = 0;
	free(z_pts);

	/* bottom edge */
	 readVoronoi(filename2);
	z_pts = malloc(sizeof(double)*triShapeCount);

for (i2=0;i2<n_selatoms2;i2++)
{
/*      printf("%f %f\n",x[sel_atoms[i2]][XX],x[sel_atoms[i2]][YY]); */
      fShape = pointInPoly(x[sel_atoms2[i2]][XX],x[sel_atoms2[i2]][YY],box[XX][XX],box[YY][YY]);
        z_pts[fShape] = x[sel_atoms2[i2]][ZZ];
}
/* analyse stuff */

 for (i2=0;(i2<isize[1]);i2++)  {
  /*  if (x[index[1][i2]][ZZ]>mbox)
    { */
/* d_xz = x[index[1][i2]][ZZ] - mbox;
    nBox = fabs(d_xz);
*/
   /* if (((nBox<(0.5*box[ZZ][ZZ]))&&(x[index[1][i2]][ZZ]<mbox))||((nBox>(0.5*box[ZZ][ZZ]))&&(x[index[1][i2]][ZZ]>mbox))) */
/*    if (((nBox<(0.5*box[ZZ][ZZ]))&&(x[index[1][i2]][ZZ]>mbox))||((nBox>(0.5*box[ZZ][ZZ]))&&(x[index[1][i2]][ZZ]<mbox)))
    {   
*/
   isTop = -1;
   tmp_z = x[index[1][i2]][ZZ];
    if (tmp_z<0)
        tmp_z += box[ZZ][ZZ];
    if (tmp_z>box[ZZ][ZZ])
        tmp_z = tmp_z - box[ZZ][ZZ];
    tmp_dx = fabs(tmp_z-mbox);
    if (tmp_z<mbox)
    {
        if (tmp_dx<(0.5*box[ZZ][ZZ]))
        {
                isTop = 1;
        }      
    }
    else
    {
        if (tmp_dx>(0.5*box[ZZ][ZZ]))
        {
                isTop = 1;
        }
    }

   fSh = pointInPoly(x[index[1][i2]][XX],x[index[1][i2]][YY],box[XX][XX],box[YY][YY]);
        /* difference in distance is ... */ /*
        mDist = x[index[1][i2]][ZZ] - z_pts[fSh]; */
  /* TODO: make this cleaner */
if (isTop==-1)
{
   /*fprintf(stderr,"bottom is %f %f %f\n", tmp_z,mbox,box[ZZ][ZZ]); */
   mDist = distPBC(z_pts[fSh],x[index[1][i2]][ZZ],box[ZZ][ZZ]);
        if (x[index[1][i2]][ZZ]<z_pts[fSh])
        {
                mDist = mDist * -1;
        }

  /*      if (nBox>(0.5*box[ZZ][ZZ]))
		mDist = box[ZZ][ZZ]-mDist;
 */
	z_dists[i2] = mDist;	
}
}
 free(tri_shapes);
        free(tri_shape_sizes); 
   /*     freeMD(tri_shapes,triShapeCount);
		free(tri_shape_sizes); */
        freeMD_d(vertices_good,vert_good_memsize);  
        freeMD(triangles,nTriangles); 
/*        free(vertices_good); */
/* free(triangles); */
free(z_pts);
/*
	free(tri_shapes);
	free(tri_shape_sizes);
	free(vertices_good);
	free(triangles);
	free(edges);
	free(pts);
	free(z_pts); */
	for (k=0;k<=lastArea;k++)
	{
	fprintf(fareas,"%f\n",areas[k]);
    }
	sfree(areas);
	fprintf(ngrpcoords,"%f ",t);
	for (k=0;k<isize[1];k++)
	{
		/* output z_dist values */
		fprintf(ngrpcoords,"%f ",z_dists[k]);
	}
	fprintf(ngrpcoords,"\n");
/*    sprintf(picNewName,"mv test3.png frm%4.0f_b.png",time);
	system(picNewName);
  */ 
 /*  system("/Users/ananikolic/Desktop/fortune/parse_matrix_3.pl mout1.tmp");
  */
/*   sprintf(commandString,"/Users/ananikolic/Desktop/fortune/parse_matrix_3.pl mout1.tmp %2.2f %2.2f",box[XX][XX],box[YY][YY]);
   system(commandString);
   system("cat areas.txt >> all_areas.txt");
   system(picNewName);
 */
     /* now somehow match up the chunks of diagram with corresponding points */
   /* match up based on s a b lines -- input XY coords  */ 
   /* then work with group 1 Z coords, match to segment, extract and compare with Z */
 /*   printf("Coordinates at t=%8.3f : %8.5f %8.5f %8.5f\n",time,x[n][XX],fr.x[n][YY],[n][ZZ]);
  */

 /* clean up and remove files */
  } while(read_next_x(status,&t,natoms,x,box));
  close_trj(status);
	/* clean up vertex and stuff */
	freeMD_d(vertices,500);
	freeMD_d(lines,500);
	freeMD(edges,500);
	freeMD_d(pts,500);
  system("rm *top.txt");
  system("rm *bottom.txt");
  sfree(z_dists);
  fclose(fareas);
  fclose(ngrpcoords);
  sfree(grpname);
  sfree(index);
  sfree(isize);
  free(verAssoc);
  thanx(stderr);
  
  return 0;
}
double r = 0;
double nBoxes = 0;
double distPBC(double d1, double d2, double box_d)
{
        r = fabs(d1-d2);
        nBoxes = floor(r/box_d + 0.5) * box_d;
        return r - nBoxes;
}

	int tpoints[3];
	double ab_v[2]; /* vector ab */
	double ab_n[2]; /*clockwise perp of vector ab */
	double p_b[2]; /* vector of point of interest through point b */
	int pt1 = 0;
	int pt2 = 0;
	int pt3 = 0;
        int triID = -1;
	double dotP[3];
int pointInPoly(double ptX_k, double ptY_k,double maxX, double maxY)
{
	/* adjust XY for PBC */
	double ptX = ptX_k;
	if (ptX<0)
	{
		ptX = ptX + maxX;
	}
	if (ptX>maxX)
	{
		ptX = ptX - maxX;
	}
	double ptY = ptY_k;
	if (ptY<0)
	{
		ptY = ptY + maxY;
	}
	else if (ptY > maxY)
	{
		ptY = ptY - maxY;
	}
	/* find which shape the point is in */
	for (k=0;k<triShapeCount;k++)
	{
		for (j=0;j<=tri_shape_sizes[k];j++)
		{
			pt1 = triangles[tri_shapes[k][j]][0];
			pt2 = triangles[tri_shapes[k][j]][1];
			pt3 = triangles[tri_shapes[k][j]][2];
			/* take the points and test intersections */
			/* first */
			ab_n[1] = (vertices_good[pt2][0] - vertices_good[pt1][0])*-1;
			ab_n[0] = vertices_good[pt2][1] - vertices_good[pt1][1];
			p_b[0] = ptX - vertices_good[pt1][0];
			p_b[1] = ptY - vertices_good[pt1][1];
			dotP[0] = p_b[0] * ab_n[0] + p_b[1] * ab_n[1];

			/* repeat for second point */
			ab_n[1] = (vertices_good[pt3][0] - vertices_good[pt2][0])*-1;
			ab_n[0] = vertices_good[pt3][1] - vertices_good[pt2][1];
			p_b[0] = ptX - vertices_good[pt2][0];
			p_b[1] = ptY - vertices_good[pt2][1];
			dotP[1] = p_b[0] * ab_n[0] + p_b[1] * ab_n[1];

			ab_n[1] = (vertices_good[pt1][0] - vertices_good[pt3][0])*-1;
			ab_n[0] = vertices_good[pt1][1] - vertices_good[pt3][1];
			p_b[0] = ptX - vertices_good[pt3][0];
			p_b[1] = ptY - vertices_good[pt3][1];
			dotP[2] = p_b[0] * ab_n[0] + p_b[1] * ab_n[1];

			/* test signedness */
			if (((dotP[0]>0)&&(dotP[1]>0))&&(dotP[2]>0))
			{
				triID = k;
				break;
			}
			else if (((dotP[0]<0)&&(dotP[1]<0))&&(dotP[2]<0))
			{
				triID = k;
				break;
			}
		}
	}
	return triID;
}
int tmpV = 0;
char commandString[200];
char* fname = "mout1.tmp";
FILE* fvor;
int curEdge=0;
FILE* ftri;
FILE* ftri2;

int i = 0;
void readVoronoi(char* filename_x)
{
   sprintf(commandString,"cat %s | %s > mout1.tmp",filename_x,voro_loc); /* filename_x); */
   system(commandString);
	lastBadVer = -1;
	lastEdge = -1;
    lastVer = -1;
    lastVerAll = -1;
    lastLine = -1;
    lastPt = -1;
   
  fname = "mout1.tmp";
  curEdge=0;
  if (edges==NULL)
  {
   edges = malloc(500*sizeof(int));
   vertices = malloc(500*sizeof(double));
   lines = malloc(500*sizeof(double));
   pts = malloc(500*sizeof(double));
   int i =0;
   if (edges!=NULL)
   {
       for (i=0;i<500;i++)
       {
            edges[i] = malloc(2*sizeof(int));
            vertices[i] = malloc(2*sizeof(double));
            lines[i] = malloc(3*sizeof(double));
            pts[i] = malloc(2*sizeof(double));
       }
   }
}

   char lastType = 0;
   fvor = fopen(fname,"r");
   int curField = 0;
   while (!feof(fvor))
   {
       char line[30];
       double mtmp = -55;
       fscanf(fvor,"%s",line);
       char tmp = 0;
       sscanf(line,"%lf%c",&mtmp,&tmp); 
	if ((mtmp==-55)&&(tmp==0))
       {
	    lastType = line[0];
            curField = 0;
            typeIncrement(lastType);
       }
       else
       {
             /* assign value */
             switch (lastType)
             {
                 case 'e':
                 {
			/* need to replace vertex index if necessary */
			tmpV = (int)mtmp;
			sscanf(line,"%d",&tmpV); /* try this */
		/*	printf("%s %d %f\n",line,tmpV,mtmp);
		*/
			if (curField==0)
			{
				curEdge = tmpV;
			}
			else
			{
				if (tmpV!=-1)
				{
					tmpV = verAssoc[tmpV];
				}
				edges[curEdge][curField-1] = tmpV;
			}
                 }
		 break;
		 case 'v':
		 {
			vertices[lastVer][curField] = mtmp;
		 }
	         break;
                 case 's':
		 {
			pts[lastPt][curField] = mtmp;
		 }
		 break;
		 case 'l':
		 {
			lines[lastLine][curField] = mtmp;
		 }
	         break;
                 default:
                 {
			/* do nothing */
                 }
             }
	     curField++;
	     if ((curField==2)&&(lastType=='v')) 
	     {
                 int vIndex = -1;
		 /* vertex control */
		 for (i=0;i<lastVer;i++)
		 {
			double diffA = fabs(vertices[i][0] - vertices[lastVer][0]);
			double diffB = fabs(vertices[i][1] - vertices[lastVer][1]);
/*			if ((vertices[i][0]==vertices[lastVer][0])&&(vertices[i][1]==vertices[lastVer][1]))
*/
			if ((diffA<cutoff)&&(diffB<cutoff))
			{
				/* found vertex */
				vIndex = i;
				break;
			}
		 }
		if (vIndex!=-1)
		{
/*			printf("%d %d %d\n",vIndex,lastVerAll,lastVer);
*/	
			verAssoc[lastVerAll] = vIndex;
			lastVer--;
		}
		else
		{
/*			printf("%d %d %d\n",vIndex,lastVerAll,lastVer);
*/
			verAssoc[lastVerAll] = lastVer;
		}
             }
       }
      if (tmp!=0)
      {
          lastType = tmp;
	  curField = 0;
	  curEdge = 0;
          typeIncrement(lastType);
      }
   }
   fclose(fvor);
/* edge truncation */


/* malloc edgesNew with the number of rows based on lastEdge+1 */
/*int** edges_new;
edges_new = malloc((lastEdge+1)*sizeof(int));
*/
for (j=0;j<=lastEdge;j++)
{
/*edges_new[j] = malloc(2 * sizeof(int));
*/
if (edges[j][0]!=-1)
{
	fitVertex(j,0,edges,vertices,lines);
}
else
{
	genVertex(j,0,edges,vertices,lines);
}
if (edges[j][1]!=-1)
{
	fitVertex(j,1,edges,vertices,lines);
}
else
{
	genVertex(j,1,edges,vertices,lines);
}
}

   /* vertices are all "good" now, except for bad ones */

/* append edge vertices */
int foundindex = -1;
foundindex=-1;
foundindex = findVertex(min_x,min_y,vertices);
if (foundindex==-1)
{
     lastVer++;
     vertices[lastVer][0] = min_x;
     vertices[lastVer][1] = min_y;
}
foundindex=-1;
foundindex = findVertex(max_x,min_y,vertices);
if (foundindex==-1)
{
     lastVer++;
     vertices[lastVer][0] = max_x;
     vertices[lastVer][1] = min_y;
}
foundindex=-1;
foundindex = findVertex(min_x,max_y,vertices);
if (foundindex==-1)
{
     lastVer++;
     vertices[lastVer][0] = min_x;
     vertices[lastVer][1] = max_y;
}
foundindex=-1;
foundindex = findVertex(max_x,max_y,vertices);
if (foundindex==-1)
{
     lastVer++;
     vertices[lastVer][0] = max_x;
     vertices[lastVer][1] = max_y;
}
/* edges are all good */

/* free some space */
/* free(lines); */

/* time for output to triangle */

ftri = fopen("tmp_pts.poly1","w");
fprintf(ftri,"%d 2 0 0\n",lastVer-lastBadVer);
j = 0;
int isbad = -1;
int cur_vref = 0;
/* 10 % overdraft for this */
vertices_good = malloc(((int)floor(1.1*(lastVer-lastBadVer)))*sizeof(double));
vert_good_memsize = lastVer - lastBadVer + 10;
/*[lastVer+10];
*/
int* vert_n_index;
vert_n_index = malloc((lastVer+1)*sizeof(int));

for (k=0;k<=lastVer;k++)
{
      isbad = -1;
      for (j=0;j<=lastBadVer;j++)
      {
	if (k==badVertList[j])
	{
		isbad = 1;
		break;
	}
      }
     if (isbad==-1)
     {     
      /* make sure it's not a "bad" vertex */
      fprintf(ftri,"%d %f %f\n",cur_vref+1, vertices[k][0], vertices[k][1]);
      /* adjust indices accordingly */
      vertices_good[cur_vref] = malloc(2*sizeof(double));
      vertices_good[cur_vref][0] = vertices[k][0];
      vertices_good[cur_vref][1] = vertices[k][1];
      vert_n_index[k] = cur_vref; 
/*      printf("test %d %d %f,%f\n", k,cur_vref,vertices_good[cur_vref][0], vertices_good[cur_vref][1]);  */
      cur_vref++;
     }
}
/*printf("stuff done %d",cur_vref); */
int edge_count_r = lastEdge+1;
ftri2 = fopen("tmp.q","w");
for (k=0;k<=lastEdge;k++)
{
    if (edges[k][0]!=edges[k][1])
    {
         fprintf(ftri2,"%d %d %d\n",k+1,vert_n_index[edges[k][0]]+1,vert_n_index[edges[k][1]]+1);
    } 
    else
    {
	edge_count_r--;
    }
}
fprintf(ftri2,"0 0\n");
fclose(ftri2);
fprintf(ftri,"%d 2\n",edge_count_r);
fclose(ftri);
system("cat tmp_pts.poly1 tmp.q > tmp_pts.poly");
char tmp_instr[200];
sprintf(tmp_instr,"%s -c tmp_pts.poly > tmp.q",tr_loc);
/*system("/Users/ananikolic/Desktop/fortune/triangle -c tmp_pts.poly > tmp.q"); */
system(tmp_instr);
system("sed '$d' tmp_pts.1.node > tmp_pts.1.node.nohead");


/* now load up vertex file and additional vertices, if any */
ftri = fopen("tmp_pts.1.node.nohead","r");
int newVert = -1;
while (!feof(ftri))
{
     if (newVert==-1)
     {
	int tmp1 = 0;
         fscanf(ftri, "%d",&newVert);
         fscanf(ftri, "%d",&tmp1);
         fscanf(ftri, "%d",&tmp1);
         fscanf(ftri, "%d",&tmp1);
   
 	if (cur_vref==newVert)
	{
		break;
	}        
     }
     else
     {
	int verInd = 0;
	double ver1 = 0;
	double ver2 = 0;
	int m_fourth = 0;
          fscanf(ftri,"%d",&verInd);
          fscanf(ftri,"%lf",&ver1);
          fscanf(ftri,"%lf",&ver2);
          fscanf(ftri,"%d",&m_fourth);
	  if (verInd>cur_vref)
	  {
		vertices_good[cur_vref] = malloc(2*sizeof(double));
           	vertices_good[cur_vref][0] = ver1;     
           	vertices_good[cur_vref][1] = ver2;     
		cur_vref++;
	  }
        /*printf("%d | %f %f\n",verInd,ver1,ver2); */
          if (verInd==newVert)
	  {
		break;
	  }
     }
}
vert_good_memsize = cur_vref;
fclose(ftri);
system("sed '$d' tmp_pts.1.ele > tmp_pts.1.ele.nohead");
/* then the triangle file with triangles */
ftri = fopen("tmp_pts.1.ele.nohead","r");
nTriangles = -1;
int rtmp = 0;
int trInd = 0;
int** edge_tr;
edge_tr = malloc(cur_vref*sizeof(int));
for (k=0;k<cur_vref;k++)
{
	edge_tr[k] = malloc(cur_vref*sizeof(int)); 
	for (j=0;j<cur_vref;j++)
	{
		edge_tr[k][j] = 0;
	}
}
int* triMatched;
while (!feof(ftri))
{
	/* read in the header */
	if (nTriangles==-1)
	{
		fscanf(ftri,"%d",&nTriangles);
		/* allocate memory */
		triangles = malloc(nTriangles*sizeof(int));
		triMatched = malloc(nTriangles*sizeof(int));
		fscanf(ftri,"%d",&rtmp);
		fscanf(ftri,"%d",&rtmp);
	}
	else
	{
		/* read a line, 4 items per line */
		 trInd = 0;
		fscanf(ftri,"%d",&trInd);
		triangles[trInd-1] = malloc(3*sizeof(int));
		triMatched[trInd-1] = 0;
		fscanf(ftri,"%d",&triangles[trInd-1][0]);
		fscanf(ftri,"%d",&triangles[trInd-1][1]);
		fscanf(ftri,"%d",&triangles[trInd-1][2]);
		triangles[trInd-1][0]--;
		triangles[trInd-1][1]--;
		triangles[trInd-1][2]--;
		edge_tr[triangles[trInd-1][0]][triangles[trInd-1][1]]++;
		edge_tr[triangles[trInd-1][1]][triangles[trInd-1][2]]++;
		edge_tr[triangles[trInd-1][2]][triangles[trInd-1][0]]++;
		edge_tr[triangles[trInd-1][1]][triangles[trInd-1][0]]++;
		edge_tr[triangles[trInd-1][2]][triangles[trInd-1][1]]++;
		edge_tr[triangles[trInd-1][0]][triangles[trInd-1][2]]++;
/*		printf("triangle %d is (%d,%d,%d)\n",trInd,triangles[trInd-1][0],triangles[trInd-1][1],triangles[trInd-1][2]); */
		if (trInd==nTriangles)
		{
			break;
		}
	}
}
fclose(ftri);

   int** non_nat; /* non-native contacts -- i.e. triangle junctions not corresponding to shapes in the voronoi diagram */
non_nat = malloc((cur_vref*cur_vref)*sizeof(int));

/*snew(non_nat,cur_vref*cur_vref); */
int lastNonNat = -1;
/* vert_n_index[k] = cur_vref;
*/
int l = 0;
int match = -1;
int diff_a = 0;
int diff_a2 = 0;
int diff_b = 0;
int diff_b2 = 0;
for (j=0;j<cur_vref-1;j++)
{
    for (k=j+1;k<cur_vref;k++)
    {
         if (edge_tr[k][j]>=2)
	{
   	 match = -1;
	for (l=0;l<=lastEdge;l++)
	{
		diff_a = vert_n_index[edges[l][0]] - k;	
		diff_a2 = vert_n_index[edges[l][1]] - k;	
		diff_b = vert_n_index[edges[l][1]] - j;	
		diff_b2 = vert_n_index[edges[l][0]] - j;
	/*for debug 	printf("%d %d %d %d\n",k,j,vert_n_index[edges[l][0]],vert_n_index[edges[l][1]]); */
		if (((diff_a==0)&&(diff_b==0))||((diff_a2==0)&&(diff_b2==0)))
		{
			match = 1;
			break;
		}
	}	
 	/*   printf("(%d,%d) match=%d",k,j,match); */
  	  if (match==-1)
    	{
		lastNonNat++;
		non_nat[lastNonNat] = malloc(2*sizeof(int));
		non_nat[lastNonNat][0] = j;
		non_nat[lastNonNat][1] = k;
		/* for debug only printf("%d %d is joined\n",j,k); */
    	}
       }
    }
}
/* done with edge_tr, free that mess */
freeMD(edge_tr,cur_vref);

/* generate a triangle connectivity map */
int** tri_conn;
tri_conn = malloc(500*sizeof(int)); /* make this allocation cleaner */
int lastTriConn = -1;
int match_tr[2];
int matched = 0;
int a = 0;
for (k=0;k<=lastNonNat;k++)
{
	int num_edge = 0;
	for (j=0;j<nTriangles;j++)
	{
		matched = 0;
		for (a=0;a<3;a++)
		{
			if ((triangles[j][a]==non_nat[k][0])||(triangles[j][a]==non_nat[k][1]))
			{
				matched++;
				if (matched==2)
					break;
			}
		}
		if (matched==2)
		{
			match_tr[num_edge] = j;
			num_edge++;
			if (num_edge==2)
				break;
		}
	}
	lastTriConn++;	
	tri_conn[lastTriConn] = malloc(2*sizeof(int));
	tri_conn[lastTriConn][0] = match_tr[0];
	tri_conn[lastTriConn][1] = match_tr[1];
}
/* printf("%d as count\n", lastTriConn); */
freeMD(non_nat,lastNonNat+1);
/* allocate new vertices */
/* leave some space for extras ? */
int init_tr = 0;
int init_tr2 = 0;
int* all_tr;
all_tr = malloc(nTriangles*sizeof(int));
int nShapeSize;
triShapeCount=0;
/* need more intelligent method to allocate these */
tri_shapes = malloc(halfGroupSize*sizeof(int));
int oldLastTriConn = lastTriConn;
tri_shape_sizes = malloc(halfGroupSize*sizeof(int));
while (lastTriConn >= 0)
{
     init_tr = tri_conn[0][0];
     init_tr2 = tri_conn[0][1];
     all_tr[0] = init_tr;   
     all_tr[1] = init_tr2;
     triMatched[init_tr]++;
     triMatched[init_tr2]++;
     nShapeSize = 2;
     int match = 0;
     int nmatches = -1;
     while (nmatches!=0)
     {
         nmatches=0;
         int b = -1;
	/* DONKEY - may need to change */
         while (b<lastTriConn)
         {
		b++;
		match = 0;
		int matched = -1;
		for (k=0;k<nShapeSize;k++)
		{
			if (tri_conn[b][0]==all_tr[k])
			{
				matched=1;
				match++;
			}
			else if (tri_conn[b][1]==all_tr[k])
			{
				matched=0;
				match++;
			}
		}
		if (match==1)
		{
			all_tr[nShapeSize] = tri_conn[b][matched];
			triMatched[tri_conn[b][matched]]++;
			nShapeSize++;
			nmatches++;
		}
		if (match==2)
		{
			nmatches++;
			int tmplast = sliceVerArray(tri_conn, b, lastTriConn,2);
			/*printf("old is %d new is %d\n",lastTriConn,tmplast); */
			b--;
			lastTriConn = tmplast;
		}
         }
     }
	/* continue with stuff */
	tri_shape_sizes[triShapeCount] = nShapeSize-1;
	int n =0;
	tri_shapes[triShapeCount] = malloc(nShapeSize*sizeof(int));
	for (n=0;n<nShapeSize;n++)
	{
		tri_shapes[triShapeCount][n] = all_tr[n];
		/*fprintf(stderr,"adding triangle %d\n",all_tr[n]); */
	}
/*	fprintf(stderr,"\n"); */
	triShapeCount++;
}
/* done with the shapes ! */
/* free(vertices); */
freeMD(tri_conn,oldLastTriConn+1);

free(all_tr);
/* account for singlets */
for (k=0;k<nTriangles;k++)
{
     if (triMatched[k]==0)
     {
		tri_shapes[triShapeCount] = malloc(1*sizeof(int));
		tri_shapes[triShapeCount][0] = k;
		tri_shape_sizes[triShapeCount] = 0;
		triShapeCount++;
     }
}
free(triMatched);
/* now map is complete */

/* may put output for picture */
/* FILE *fpic = fopen("forpic.txt","w");   */
/* output files with p for points */ 
/*for (k=0;k<=lastPt;k++) 
{
	fprintf(fpic,"p%f %f\n",pts[k][0],pts[k][1]);
} */
/* output files with v for vertices */ 
/*for (k=0;k<cur_vref;k++)
{
	fprintf(fpic,"v%f %f\n", vertices_good[k][0],vertices_good[k][1]);
} */
/* output files with e for edges */ 
/*for (k=0;k<=lastEdge;k++)
{
	fprintf(fpic,"e%d %d\n",vert_n_index[edges[k][0]],vert_n_index[edges[k][1]]);
} */
/* output files with t for triangles -- include color */
/*for (k=0;k<triShapeCount;k++)
{
	for (j=0;j<=tri_shape_sizes[k];j++)
	{
		fprintf(fpic,"t%d %d %d %d\n",triangles[tri_shapes[k][j]][0],triangles[tri_shapes[k][j]][1],triangles[tri_shapes[k][j]][2],k);
	}
}
fclose(fpic);
char mCmd1[100];
sprintf(mCmd1,"/Users/ananikolic/Desktop/fortune/draw_only.pl forpic.txt %f %f",max_x,max_y);
system(mCmd1);
sprintf(mCmd1, "mv test3.png %s.png",filename_x);
system(mCmd1);
*/



/* area calculation */
free(vert_n_index);

/* do areas */
int v1 = 0;
int v2 = 0;
int v3 = 0;
double ab_x = 0;
double ab_y = 0;
double ac_x = 0;
double ac_y = 0;
double det = 0;
double tmp_area = 0;
for (k=0;k<triShapeCount;k++)
{
	tmp_area = 0;
     for (j=0;j<=tri_shape_sizes[k];j++)
     {
			v1 = triangles[tri_shapes[k][j]][0];
			v2 = triangles[tri_shapes[k][j]][1];
			v3 = triangles[tri_shapes[k][j]][2];
			ab_x = vertices_good[v2][0] - vertices_good[v1][0];
			ab_y = vertices_good[v2][1] - vertices_good[v1][1];
			ac_x = vertices_good[v3][0] - vertices_good[v1][0];
		    ac_y = vertices_good[v3][1] - vertices_good[v1][1];									
			det = ab_x*ac_y - ac_x*ab_y;
			tmp_area += 0.5*fabs(det);
	 }
	 lastArea++;
	 areas[lastArea] = tmp_area;
}
}
int sliceVerArray(int** srcArray, int remIndex, int rows, int cols)
{
	int i1 = 0;
	int i2 = 0;
        int offset = 0;
	for (i1=0;i1<=rows;i1++)
	{
		if (i1!=remIndex)
		{
		for (i2=0;i2<cols;i2++)
		{
			/* copy columns */
			srcArray[i1-offset][i2] = srcArray[i1][i2];
		}
		}
		else
		{
			offset++;
		}
	}
	return rows-1;
}


void genVertex(int edge_id, int edge_pt,int** edges, double** vertices, double** lines)
{
	double mx = min_x;
	double my = min_y;
	if (edge_pt==1)
	{
		/* second point*/
		if (lines[edge_id][1] != 0)
		{
			mx = max_x;
			my = max_y;
		}
	}
	else if (edge_pt==0)
	{
		if (lines[edge_id][1] == 0)
		{
			mx = max_x;
			my = max_y;
		}
	}
	double tmpX = 0;
	double tmpY = 0;
	if (lines[edge_id][1]!=0)
	{
		int outofbounds = -1;
		tmpY = -1*mx*(lines[edge_id][0]/lines[edge_id][1]) + lines[edge_id][2]/lines[edge_id][1];
		if (tmpY < min_y)
		{
			outofbounds=1;
			tmpY = min_y;
			tmpX = (lines[edge_id][2] - lines[edge_id][1]*min_y)/lines[edge_id][0];
			mx = tmpX;
		}
		else if (tmpY > max_y)
		{
			outofbounds=1;
			tmpX = (lines[edge_id][2] - lines[edge_id][1]*max_y)/lines[edge_id][0];
			tmpY = max_y;
			mx = tmpX;
		}
		if (outofbounds==1)
		{
			if (mx<min_x)
			{
				double oldY = tmpY;
				tmpY = (-1*min_x*lines[edge_id][0]+lines[edge_id][2])/lines[edge_id][1];
				mx = min_x;
			}
			else if (mx>max_x)
			{
				double oldY = tmpY;
				tmpY = (-1*max_x*lines[edge_id][0]+lines[edge_id][2])/lines[edge_id][1];
				mx = max_x;
			}
		}
		/* seek the vertex mx,tmp_y in the vertex */
               int foundIndex = -1;
		foundIndex = findVertex(mx,tmpY,vertices);
	/*	printf("foundindex is %d\n",foundIndex); */
                if (foundIndex==-1)
                {
                        /* append vertex and set thing */
                        lastVer++;
                        vertices[lastVer][0] = mx;
                        vertices[lastVer][1] = tmpY;
                        edges[edge_id][edge_pt] = lastVer;
                }
                else
                {
                        edges[edge_id][edge_pt] = foundIndex;
                }
	}
	else
	{
		/* vertical line */
		tmpY = lines[edge_id][2]/lines[edge_id][0];
		/* seek stuff */
		int foundIndex = -1;
		foundIndex = findVertex(tmpY,my,vertices);
		/*printf("foundindex is %d\n",foundIndex);*/
		if (foundIndex==-1)
		{
			/* append vertex and set thing */
			lastVer++;
			vertices[lastVer][0] = tmpY;
			vertices[lastVer][1] = my;
			edges[edge_id][edge_pt] = lastVer;
		}
		else
		{
			edges[edge_id][edge_pt] = foundIndex;
		}
		/* fill point tmp_y,my if necessary */
	}
}

void fitVertex(int edge_id,int edge_pt,int** edges, double** vertices, double** lines)
{
    double newX1 = -1;
    double newY1 = -1;
    int altered = -1;
    double tmpX = vertices[edges[edge_id][edge_pt]][0];
    double tmpY = vertices[edges[edge_id][edge_pt]][1];
    int refStep = 0;
    double ax = 0;
    double ay = 0;

    REFINE:
    if (refStep > 10)
        goto DONE_REFINE;
    if (((tmpX <= min_x)||(tmpX >= max_x))&&(refStep < 3))
    {
        altered = 1;
        ax = min_x;
        if (ax < tmpX)
             ax = max_x;
        double new_y =  -1*ax*lines[edge_id][0]/lines[edge_id][1] + lines[edge_id][2]/lines[edge_id][1];  
	tmpY = new_y;
        tmpX = ax;
        newX1 = ax;
        newY1 = new_y;
    }
    if (((refStep==0)&&((tmpY <= min_y)||(tmpY >= max_y)))||((tmpY < min_y)||(tmpY > max_y)))
    {
	altered = 1;
        ay = min_y;
        if (ay < tmpY)
	    ay = max_y;
       double new_x = -1*ay*lines[edge_id][1]/lines[edge_id][0] + lines[edge_id][2]/lines[edge_id][0];
       newY1 = ay;
       newX1 = new_x;
       tmpX = newX1;
       tmpY = newY1;
       refStep++;
       goto REFINE;
    }

    DONE_REFINE:
    if (altered==1)
    {
	/* bad vert list */
        int k = 0;
	int found = -1;
	for (k=0;k<=lastBadVer;k++)
	{
		if (badVertList[k]==edges[edge_id][edge_pt])
		{
			found = k;
		}
	}
 	if (found==-1)	
	{
		lastBadVer++;
		badVertList[lastBadVer] = edges[edge_id][edge_pt];
	}
	/*printf("for edge [%d,%d] %f %f modified\n",edge_id, edge_pt, newX1,newY1); */
       /* seek the vertex mx,tmp_y in the vertex */
       int foundIndex = -1;
	/*foundIndex = findVertex(newX1,newY1); */
	foundIndex = findVertex(newX1,newY1,vertices);
	/*printf("foundindex is %d\n",foundIndex);*/
        if (foundIndex==-1)
        {
              /* append vertex and set thing */
              lastVer++;
              vertices[lastVer][0] = newX1;
              vertices[lastVer][1] = newY1;
              edges[edge_id][edge_pt] = lastVer;
        }
        else
        {
              edges[edge_id][edge_pt] = foundIndex;
        }
    }
}
int t = 0;
int findVertex(double v_x, double v_y,double** vertices)
{
	float diffA = 5;
	float diffB = 5;
	int isfound = -1;
	float cutoff_a = 0.0001;
	/* find a vertex */
	for (t=0;t<=lastVer;t++)
	{
		diffA = fabs(vertices[t][0]-v_x);
		diffB = fabs(vertices[t][1]-v_y);
		if ((diffA<cutoff_a)&&(diffB<cutoff_a))
		{
			isfound = t;
			break;
		}
	}
	return isfound;	
}

void typeIncrement(char mType)
{
	   switch (mType)
            {
                case 'e':
                {
                    lastEdge++;
                }
		break;
                case 'v':
                 {
                        lastVer++;
                        lastVerAll++;
                 }
	         break;
                 case 's':
                 {
                        lastPt++;
                 }
		 break;
                 case 'l':
                 {
                        lastLine++;
                 }
		break;
                 default:
                 {
                 }
            }

}

void freeMD(int** marray, int rows)
{
    for (j=0;j<rows;j++)
	{
		if (marray[j]!=NULL)
			free(marray[j]);
	}
	free(marray);
}
void freeMD_d(double** marray, int rows)
{
    for (j=0;j<rows;j++)
	{
		if (marray[j]!=NULL)
			free(marray[j]);
	}
    free(marray);
}
