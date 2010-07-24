/*
 * $Id: gmx_mdmat.c,v 1.7 2008/03/07 11:21:17 hess Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

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
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <string.h>

#include "macros.h"
#include "vec.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "filenm.h"
#include "statutil.h"
#include "copyrite.h"
#include "futil.h"
#include "gmx_fatal.h"
#include "smalloc.h"
#include "matio.h"
#include "xvgr.h"
#include "index.h"
#include "tpxio.h"
#include "rmpbc.h"
#include "pbc.h"

#define FARAWAY 10000

int *res_ndx(t_atoms *atoms)
{
  int *rndx;
  int i,r0;
  
  if (atoms->nr <= 0)
    return NULL;
  snew(rndx,atoms->nr);
  r0=atoms->atom[0].resnr;
  for(i=0; (i<atoms->nr); i++)
    rndx[i]=atoms->atom[i].resnr-r0;
  
  return rndx;
}

int *res_natm(t_atoms *atoms)
{
  int *natm;
  int i,j,r0;
  
  if (atoms->nr <= 0)
    return NULL;
  snew(natm,atoms->nres);
  r0=atoms->atom[0].resnr;
  j=0;
  for(i=0; (i<atoms->nres); i++) {
    while ((atoms->atom[j].resnr)-r0 == i) {
      natm[i]++;
      j++;
    }
  }
  
  return natm;
}

static void calc_mat(int nres, int natoms, int rndx[],
		     rvec x[], atom_id *index,
		     real trunc, real **mdmat, int **mdmat_contact, 
		     int **nmat,int ePBC,matrix box)
{
  int i,j,resi,resj;
  real trunc2,r,r2;
  t_pbc pbc;
  rvec ddx;

  set_pbc(&pbc,ePBC,box);
  trunc2=sqr(trunc);
  for(resi=0; (resi<nres); resi++) {
    for(resj=0; (resj<nres); resj++) {
      mdmat[resi][resj]=FARAWAY;
      mdmat_contact[resi][resj]=0;
    }
  }

  for(i=0; (i<natoms); i++) {
    resi=rndx[i];
    for(j=i+1; (j<natoms); j++) {
      resj=rndx[j];
      pbc_dx(&pbc,x[index[i]],x[index[j]],ddx);
      r2 = norm2(ddx);
      if ( r2 < trunc2 ) {
		nmat[resi][j]++;
		nmat[resj][i]++;

 		//GRACE: calculate number of contact based matrix 
		if(resi!=resj){
			mdmat_contact[resi][resj]++;
		}
      }
      mdmat[resi][resj]=min(r2,mdmat[resi][resj]);
    }
  } 
  
  for(resi=0; (resi<nres); resi++) {
    mdmat[resi][resi]=0;
    for(resj=resi+1; (resj<nres); resj++) {
      r=sqrt(mdmat[resi][resj]);
      mdmat[resi][resj]=r;
      mdmat[resj][resi]=r;
    }
  }


}

static void tot_nmat(int nres, int natoms, int nframes, int **nmat, 
		     int *tot_n, real *mean_n)
{
  int i,j;
  
  for (i=0; (i<nres); i++) {
    for (j=0; (j<natoms); j++)
      if (nmat[i][j] != 0) {
	tot_n[i]++;
	mean_n[i]+=nmat[i][j];
      }
    mean_n[i]/=nframes;
  }
}

void write_txt_real_mat(FILE* fp, real** mat, int nres){
	int i,j;
	for (i=0; (i<nres); i++) {
      	for (j=0; (j<nres); j++) {
         		fprintf(fp, "%.4f ", mat[i][j]);
      	}
      	fprintf(fp, "\n");
 	}
 	fprintf(fp, "\n");
}

void write_txt_int_mat(FILE* fp, int** mat, int nres){
   int i,j;
   for (i=0; (i<nres); i++) {
        for (j=0; (j<nres); j++) {
                fprintf(fp, "%d ", mat[i][j]);
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
}

int calc_frac_native(int** mdmat_tpr, int** mdmat_contact, int NUM_NATIVE, int nres) {
	/* this function just does a matrix compare and totals the number of entries in common*/
	int i,j;
	int curr_frame_num_native=0;

    /*compare only upper right side of the matrix (matrix is symmetric)*/
	for(i=0;i<nres; i++){
		for(j=i+1;j<nres; j++){
			if((j-i)>1) {
				if(mdmat_tpr[i][j] > 0 && (mdmat_tpr[i][j] == mdmat_contact[i][j])){
					curr_frame_num_native++;		
				}			
			}
		}
	} 

	return curr_frame_num_native;
}

int sum_mat(int** mat, int nres){
	int total=0, i,j;
	for(i=0; i<nres; i++){
		for(j=i+1;j<nres;j++){
 			if((j-i)>1){
				total+=mat[i][j];
			}
		}
	}
	return total;
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_mdmat makes distance matrices consisting of the smallest distance",
    "between residue pairs. With -frames these distance matrices can be",
    "stored as a function",
    "of time, to be able to see differences in tertiary structure as a",
    "funcion of time. If you choose your options unwise, this may generate",
    "a large output file. Default only an averaged matrix over the whole",
    "trajectory is output.",
    "Also a count of the number of different atomic contacts between",
    "residues over the whole trajectory can be made.",
    "The output can be processed with xpm2ps to make a PostScript (tm) plot."
  };

  static real truncate=1.5;
  static bool bAtom=FALSE;
  static int  nlevels=40;
  static char* distmat_text_file="dist_matrix.txt";
  static char* contactmat_text_file="contact_matrix.txt";
  static char* native_text_file="fraction_native.txt";
  t_pargs pa[] = { 
    { "-t",   FALSE, etREAL, {&truncate},
      "trunc distance" },
    { "-nlevels",   FALSE, etINT,  {&nlevels},
      "Discretize distance in # levels" },
    { "-txt-dist", FALSE, etSTR, {&distmat_text_file}, 
      "Write the average distance matrix to a text file"},
	{ "-txt-contact", FALSE, etSTR, {&contactmat_text_file},
	  "Write the residue-residue contact matrix to a text file"},
	{ "-txt-native", FALSE, etSTR, {&native_text_file},
	  "Outputs the fraction of native contacts; makes sense for Ca-Ca"}
  };

  t_filenm   fnm[] = {
    { efTRX, "-f",  NULL, ffREAD },
    { efTPS, NULL,  NULL, ffREAD },
    { efNDX, NULL,  NULL, ffOPTRD },
    { efXPM, "-mean", "dm", ffWRITE },
    { efXPM, "-frames", "dmf", ffOPTWR },
    { efXVG, "-no", "num",ffOPTWR },
  };
#define NFILE asize(fnm)

  FILE       *out=NULL,*fp;
 
  t_topology top;
  int        ePBC;
  t_atoms    useatoms;
  int        isize;
  atom_id    *index;
  char       *grpname;
  int        *rndx,*natm,prevres,newres;
  int        i,j,status,nres,natoms,nframes,it,trxnat;
  int        nr0;
  bool       bCalcN,bFrames;
  real       t,ratio;
  char       title[256],label[234];
  t_rgb      rlo,rhi;
  rvec       *x;
  rvec       *xtpr;  /* GRACE */
  real       **mdmat,*resnr,**totmdmat;
  int        **mdmat_contact;
  int        **totmdmat_contact;
  int        **mdmat_tpr;
  int        **nmat,**totnmat;
  real       *mean_n;
  int        *tot_n;
  matrix     box;
 
  int total_native_contacts=0;
 
  CopyRight(stderr,argv[0]);

  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_BE_NICE,NFILE,fnm,
		    asize(pa),pa,asize(desc),desc,0,NULL);
  
  fprintf(stderr,"Will truncate at %f nm\n",truncate);
  bCalcN = opt2bSet("-no",NFILE,fnm);
  bFrames= opt2bSet("-frames",NFILE,fnm);
  if ( bCalcN ) 
    fprintf(stderr,"Will calculate number of different contacts\n");
  

  /*GRACE: read in the structure in the tpr file and save it in a separate array*/ 
  read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&ePBC,&xtpr,NULL,box,FALSE);
  
  fprintf(stderr,"Select group for analysis\n");
  get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&isize,&index,&grpname);
 
  /*  fprintf(stderr, "grpname=%s", grpname); */

  natoms=isize;
  snew(useatoms.atom,natoms);
  snew(useatoms.atomname,natoms);
    
  useatoms.nres = 0;
  snew(useatoms.resname,natoms);
  
  prevres = top.atoms.atom[index[0]].resnr;
  newres  = 0;
  for(i=0;(i<isize);i++) {
    int ii = index[i];
    useatoms.atomname[i]=top.atoms.atomname[ii];
    if (top.atoms.atom[ii].resnr != prevres) {
      prevres = top.atoms.atom[ii].resnr;
      newres++;
      useatoms.resname[i] = top.atoms.resname[prevres];
      if (debug) {
	fprintf(debug,"New residue: atom %5s %5s %6d, index entry %5d, newres %5d\n",
		*(top.atoms.resname[top.atoms.atom[ii].resnr]),
		*(top.atoms.atomname[ii]),
		ii,i,newres);
      }
    }
    useatoms.atom[i].resnr = newres;
  }
  useatoms.nres = newres+1;
  useatoms.nr = isize;
    
  rndx=res_ndx(&(useatoms));
  natm=res_natm(&(useatoms));
  nres=useatoms.nres;
  fprintf(stderr,"There are %d residues with %d atoms\n",nres,natoms);
    
  snew(resnr,nres);
  snew(mdmat,nres);
  snew(mdmat_contact,nres);
  snew(mdmat_tpr,nres);
  snew(nmat,nres);
  snew(totnmat,nres);
  snew(mean_n,nres);
  snew(tot_n,nres);
  for(i=0; (i<nres); i++) {
    snew(mdmat[i],nres);
    snew(mdmat_contact[i],nres);
    snew(mdmat_tpr[i],nres);
    snew(nmat[i],natoms);
    snew(totnmat[i],natoms);
    resnr[i]=i+1;
  }

  snew(totmdmat,nres);
  snew(totmdmat_contact,nres);
  for(i=0; (i<nres); i++) {
    snew(totmdmat[i],nres);
    snew(totmdmat_contact[i],nres);
  }
  
  trxnat=read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
  
  nframes=0;
  rlo.r=1.0, rlo.g=1.0, rlo.b=1.0;
  rhi.r=0.0, rhi.g=0.0, rhi.b=0.0;
  if (bFrames)
    out=opt2FILE("-frames",NFILE,fnm,"w");


  /* calculate the contact mat for the tpr file -- assuming that this is the native protein*/
  calc_mat(nres,natoms,rndx,xtpr,index,truncate,mdmat,mdmat_tpr,nmat,ePBC,box);
  total_native_contacts=sum_mat(mdmat_tpr,nres);

  FILE* fp_native = fopen(native_text_file, "w");
  fprintf(fp_native, "# native_frame, total_contacts, total_native_contacts, Q\n");
  do {
    rm_pbc(&top.idef,ePBC,trxnat,box,x,x);
    nframes++;

    calc_mat(nres,natoms,rndx,x,index,truncate,mdmat,mdmat_contact,nmat,ePBC,box);

	int total_contact = sum_mat(mdmat_contact,nres);
    int qnum = calc_frac_native(mdmat_tpr, mdmat_contact, total_native_contacts, nres);

	fprintf(fp_native, "%d %d %d %f\n", qnum, total_contact, total_native_contacts, (double)qnum/total_native_contacts);

    for (i=0; (i<nres); i++)
      for (j=0; (j<natoms); j++)
		if (nmat[i][j]) 
		  totnmat[i][j]++;

    for (i=0; (i<nres); i++){
      for (j=0; (j<nres); j++){
		totmdmat[i][j] += mdmat[i][j];
		totmdmat_contact[i][j] += mdmat_contact[i][j];
      }
    }

    if (bFrames) {
      sprintf(label,"t=%.0f ps",t);
      write_xpm(out,0,label,"Distance (nm)","Residue Index","Residue Index",
		nres,nres,resnr,resnr,mdmat,0,truncate,rlo,rhi,&nlevels);
    }
  } while (read_next_x(status,&t,trxnat,x,box));

  fclose(fp_native);
  fprintf(stderr,"\n");
  close_trj(status);

  if (bFrames)
    fclose(out);
  
  fprintf(stderr,"Processed %d frames\n",nframes);
    
  for (i=0; (i<nres); i++)
    for (j=0; (j<nres); j++)
      totmdmat[i][j] /= nframes;

  write_xpm(opt2FILE("-mean",NFILE,fnm,"w"),0,"Mean smallest distance",
	    "Distance (nm)","Residue Index","Residue Index",
	    nres,nres,resnr,resnr,totmdmat,0,truncate,rlo,rhi,&nlevels);
 
   /*GRACE: writes out the contact map as a matrix in text format*/
   FILE* fp_test = fopen("native_contact_map.txt", "w");
   if(fp_test){
		write_txt_int_mat(fp_test, mdmat_tpr, nres);
   }else{
		printf("GRACE: file open failed");
   }

   FILE* fp_dist_map_text = fopen(distmat_text_file, "w");
   if(fp_dist_map_text){
		write_txt_real_mat(fp_dist_map_text, totmdmat, nres);
   }else{
		printf("GRACE: -txt-dist file could not be opened\n");
   }

	FILE* fp_contact_map_text = fopen(contactmat_text_file, "w");
	if(fp_contact_map_text) {
		write_txt_int_mat(fp_contact_map_text, totmdmat_contact, nres);
	}else{
		printf("GRACE: -txt-contact file could not be opened\n");
	}	


  if ( bCalcN ) {
    tot_nmat(nres,natoms,nframes,totnmat,tot_n,mean_n);
    fp=xvgropen(ftp2fn(efXVG,NFILE,fnm),
		"Increase in number of contacts","Residue","Ratio");
    fprintf(fp,"@ legend on\n");
    fprintf(fp,"@ legend box on\n");
    fprintf(fp,"@ legend loctype view\n");
    fprintf(fp,"@ legend 0.75, 0.8\n");
    fprintf(fp,"@ legend string 0 \"Total/mean\"\n");
    fprintf(fp,"@ legend string 1 \"Total\"\n");
    fprintf(fp,"@ legend string 2 \"Mean\"\n");
    fprintf(fp,"@ legend string 3 \"# atoms\"\n");
    fprintf(fp,"@ legend string 4 \"Mean/# atoms\"\n");
    fprintf(fp,"#%3s %8s  %3s  %8s  %3s  %8s\n",
	    "res","ratio","tot","mean","natm","mean/atm");
    for (i=0; (i<nres); i++) {
      if (mean_n[i]==0)
	ratio=1;
      else
	ratio=tot_n[i]/mean_n[i];
      fprintf(fp,"%3d  %8.3f  %3d  %8.3f  %3d  %8.3f\n",
	      i+1,ratio,tot_n[i],mean_n[i],natm[i],mean_n[i]/natm[i]);
    }
    fclose(fp);
  }
    
  thanx(stderr);
    
  return 0;
}
