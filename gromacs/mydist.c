#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <typedefs.h>

#include "smalloc.h"
#include "macros.h"
#include "math.h"
#include "xvgr.h"
#include "copyrite.h"
#include "statutil.h"
//#include "string2.h"
#include "vec.h"
#include "index.h"
#include "fatal.h"
#include "futil.h"
#include "gstat.h"
#include "pbc.h"

int main(int argc,char *argv[])
{
  static char *desc[] = {
	"g_mydist is g_dist with periodic boundary conditions disabled completely.",
	"g_mydist does not output the dx,dy,dz components of the coordinate displacement vector.",
	"see g_dist -h for other usage"
  };
  
  t_topology *top=NULL;
  real t,cut2,dist2;
  rvec *x=NULL,*v=NULL,dx;
  matrix box;
  int status;
  int natoms;

  int g,d,i,j,res,teller=0;
  atom_id aid;

  int     ngrps;     /* the number of index groups */
  atom_id **index,max;   /* the index for the atom numbers */
  int     *isize;    /* the size of each group */
  char    **grpname; /* the name of each group */
  rvec    *com;
  real    *mass;
  FILE    *fp=NULL;
  bool    bCutoff;
  t_pbc   pbc;

  char    *leg[4] = { "|d|","d\\sx\\N","d\\sy\\N","d\\sz\\N" };

  static real cut=0;

  static t_pargs pa[] = {
    { "-dist",      FALSE, etREAL, {&cut},
      "Print all atoms in group 2 closer than dist to the center of mass of group 1" },
  };
#define NPA asize(pa)

  t_filenm fnm[] = {
    { efTRX, "-f", NULL, ffREAD },
    { efTPX, NULL, NULL, ffREAD },
    { efNDX, NULL, NULL, ffOPTRD },
    { efXVG, NULL, "dist", ffOPTWR },
  };
#define NFILE asize(fnm)


  CopyRight(stderr,argv[0]);

  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_BE_NICE,
		    NFILE,fnm,NPA,pa,asize(desc),desc,0,NULL);
  
  bCutoff=opt2parg_bSet("-dist",NPA,pa);
  cut2=cut*cut;
  
  top=read_top(ftp2fn(efTPX,NFILE,fnm));
  
  /* read index files */
  ngrps = 2;
  snew(com,ngrps);
  snew(grpname,ngrps);
  snew(index,ngrps);
  snew(isize,ngrps);
  get_index(&top->atoms,ftp2fn(efNDX,NFILE,fnm),ngrps,isize,index,grpname);
  
  /* calculate mass */
  max=0;
  snew(mass,ngrps);
  for(g=0;(g<ngrps);g++) {
    mass[g]=0;
    for(i=0;(i<isize[g]);i++) {
      if (index[g][i]>max)
	max=index[g][i];
      if (index[g][i] >= top->atoms.nr)
	gmx_fatal(FARGS,"Atom number %d, item %d of group %d, is larger than number of atoms in the topolgy (%d)\n",index[g][i]+1,i+1,g+1,top->atoms.nr+1);
      mass[g]+=top->atoms.atom[index[g][i]].m;
    }
  }

  natoms=read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);

  if (max>=natoms)
    gmx_fatal(FARGS,"Atom number %d in an index group is larger than number of atoms in the trajectory (%d)\n",(int)max+1,natoms);

  if (!bCutoff) {
    /* open output file */
    fp = xvgropen(ftp2fn(efXVG,NFILE,fnm),
		  "Distance","Time (ps)","Distance (nm)");
    xvgr_legend(fp,4,leg);
  } else
    ngrps=1;
  
  do {
    /* initialisation for correct distance calculations */

   /* GL : set box diagonal to zero to avoid pbc distance calculations 
 	* note that set_pbc do not use pbc to calc. distances if one of 
 	* the diagonal entries for box is 0	
 	*/

 /* GL: comment out this line to get PBC back */ 

    box[XX][XX]=box[YY][YY]=box[ZZ][ZZ]=0;
/* end comment */

    set_pbc(&pbc,box);

    /* make molecules whole again */
	/*GL: rm_pbc transforms the molecules' coordinates to make them whole.  
 	* Not sure if rm_pbc needs to be called here..but I think it was put here 
 	* in case a user used an input trajectory with non-whole molecule(s)
 	*/

    rm_pbc(&top->idef,natoms,box,x,x);

    /* calculate center of masses */
    for(g=0;(g<ngrps);g++) {
      for(d=0;(d<DIM);d++) {
	com[g][d]=0;
	for(i=0;(i<isize[g]);i++) {
	  com[g][d] += x[index[g][i]][d] * top->atoms.atom[index[g][i]].m;
	}
	com[g][d] /= mass[g];
      }
    }
    
    if (!bCutoff) {
      /* write to output */
      fprintf(fp,"%12.7f ",t);
      for(g=0;(g<ngrps/2);g++) {
	pbc_dx(&pbc,com[2*g],com[2*g+1],dx);
	/*fprintf(fp,"%12.7f %12.7f %12.7f %12.7f",
		norm(dx),dx[XX],dx[YY],dx[ZZ]);*/
	/*GL: prints only the distance*/
	fprintf(fp,"%12.7f",norm(dx));
      }
      fprintf(fp,"\n");
    } else {
      for(i=0;(i<isize[1]);i++) { 
	j=index[1][i];
	pbc_dx(&pbc,x[j],com[0],dx);
	dist2 = norm2(dx);
	if (dist2<cut2) {
	  res=top->atoms.atom[j].resnr;
	  fprintf(stdout,"\rt: %g  %d %s %d %s  %g (nm)\n",
		  t,res+1,*top->atoms.resname[res],
		  j+1,*top->atoms.atomname[j],sqrt(dist2));     
	} 
      }
    }
    
    teller++;
  } while (read_next_x(status,&t,natoms,x,box));

  if (!bCutoff)
    fclose(fp);

  close_trj(status);
  
  thanx(stderr);
  return 0;
}





  
