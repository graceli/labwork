#define CPLUSPLUS

using namespace std;

//c++ libraries
#include <map>
#include <set>
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>

extern "C" {
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <ctype.h>

#include "string2.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "macros.h"
#include "vec.h"
#include "pbc.h"
#include "rmpbc.h"
#include "statutil.h"
#include "xvgr.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "tpxio.h"
#include "index.h"
#include "smalloc.h"
#include "calcgrid.h"
#include "nrnb.h"
#include "physics.h"
#include "coulomb.h"
#include "pme.h"
#include "gstat.h"
#include "matio.h"
#include "mtop_util.h"
}

void output_cluster_info(gmx_mtop_t *mtop, const string &output, int* clust_index, int nindex, real t) {
	// Transform clust_index into a data structure that can be outputted to identify the inositol molecules in a particular
	// cluster

	t_topology top = gmx_mtop_t_to_t_topology(mtop);

	typedef map<int, set<int> >::iterator MapIter;
	typedef set<int>::iterator SetIter;

	MapIter map_iter;
	SetIter set_iter;

	// This vector maps cluster id to a set containing the residue ids of inositol
	map<int, set<int> > cluster_info;
	for(int i=0; (i < nindex); i++) {
		int inos_residue_id = top.atoms.atom[i].resnr;
		int cluster_id = clust_index[i];
		cluster_info[cluster_id].insert(inos_residue_id);
	}

	// Find all the clusters with only one molecule
	set<int> unclustered_residue_ids;
	for(map_iter = cluster_info.begin(); (map_iter != cluster_info.end()); ++map_iter) {
		set<int> cluster = map_iter->second;
		if(cluster.size() == 1) {
			// Since the set has only one element, the first element pointed to by the iterator
			// should be that element.
			unclustered_residue_ids.insert(*(cluster.begin()));
		}
	}

	ofstream f_cluster_info(output.c_str());
	// Output the cluster info data
	for(map_iter = cluster_info.begin(); map_iter != cluster_info.end(); ++map_iter) {
		set<int> a_cluster = map_iter->second;
		if(a_cluster.size() > 1) {
			f_cluster_info << t << " yes";
			for(set_iter = a_cluster.begin(); set_iter != a_cluster.end(); ++set_iter) {
				f_cluster_info << " " << *set_iter;
			}
			f_cluster_info << endl;
		}
	}

	// Output the molecules that are not in a cluster on a single row
	f_cluster_info << t << " no";
	for(set_iter = unclustered_residue_ids.begin(); set_iter != unclustered_residue_ids.end(); ++set_iter) {
		f_cluster_info << " " << *set_iter;
	}
	f_cluster_info << endl;
}

static void clust_size(char *ndx, char *trx, char *xpm,
		char *xpmw, char *ncl, char *acl,
		char *mcl, char *histo, char *tempf,
		char *mcn, bool bMol, bool bPBC, char *tpr,
		real cut, int nskip, int nlevels,
		t_rgb rmid, t_rgb rhi, int ndf)
{
	FILE    *fp,*gp,*hp,*tp;
	atom_id *index=NULL;
	int     nindex,natoms,status;
	rvec    *x=NULL,*v=NULL,dx;
	t_pbc   pbc;
	char    *gname;
	char    timebuf[32];
	bool    bSame,bTPRwarn=TRUE;
	/* Topology stuff */
	t_trxframe  fr;
	t_tpxheader tpxh;
	gmx_mtop_t *mtop=NULL;
	int     ePBC=-1;
	t_block *mols=NULL;
	t_atom  *atom;
	int     version,generation,sss,ii,jj,nsame;
	real    ttt,lll,temp,tfac;
	/* Cluster size distribution (matrix) */
	real    **cs_dist=NULL;
	real    tf,dx2,cut2,*t_x=NULL,*t_y,cmid,cmax,cav,ekin;
	int     i,j,k,ai,aj,ak,ci,cj,nframe,nclust,n_x,n_y,max_size=0;
	int     *clust_index,*clust_size,max_clust_size,max_clust_ind,nav,nhisto;
	t_rgb   rlo = { 1.0, 1.0, 1.0 };

	clear_trxframe(&fr,TRUE);
	sprintf(timebuf,"Time (%s)",time_unit());
	tf     = time_factor();
	fp     = xvgropen(ncl,"Number of clusters",timebuf,"N");
	gp     = xvgropen(acl,"Average cluster size",timebuf,"#molecules");
	hp     = xvgropen(mcl,"Max cluster size",timebuf,"#molecules");
	tp     = xvgropen(tempf,"Temperature of largest cluster",timebuf,"T (K)");

	if (!read_first_frame(&status,trx,&fr,TRX_NEED_X | TRX_READ_V))
		gmx_file(trx);

	natoms = fr.natoms;
	x      = fr.x;

	if (tpr) {
		snew(mtop,1);
		read_tpxheader(tpr, &tpxh, TRUE, &version, &generation);
		if (tpxh.natoms != natoms)
			gmx_fatal(FARGS,"tpr (%d atoms) and xtc (%d atoms) do not match!",
					tpxh.natoms,natoms);
		ePBC = read_tpx(tpr,&sss,&ttt,&lll,NULL,NULL,&natoms,NULL,NULL,NULL,mtop);
	}

	if (ndf <= -1)
		tfac = 1;
	else
		tfac = ndf/(3.0*natoms);

	if (bMol) {
		if (ndx)
			printf("Using molecules rather than atoms. Not reading index file %s\n",
					ndx);
		mols = &(mtop->mols);

		/* Make dummy index */
		nindex = mols->nr;
		snew(index,nindex);
		for(i=0; (i<nindex); i++)
			index[i] = i;
		gname = strdup("mols");
	}
	else
		rd_index(ndx,1,&nindex,&index,&gname);

	snew(clust_index, nindex);
	snew(clust_size, nindex);
	cut2   = cut*cut;
	nframe = 0;
	n_x    = 0;
	snew(t_y, nindex);

	for(i=0; (i < nindex); i++)
		t_y[i] = i+1;

	max_clust_size = 1;
	max_clust_ind  = -1;

	do {
		if ((nskip == 0) || ((nskip > 0) && ((nframe % nskip) == 0))) {
			if (bPBC)
				set_pbc(&pbc,ePBC,fr.box);

			max_clust_size = 1;
			max_clust_ind  = -1;

			/* Put all atoms/molecules in their own cluster, with size 1 */
			for(i=0; (i < nindex); i++) {
				/* Cluster index is indexed with atom index number */
				clust_index[i] = i;
				/* Cluster size is indexed with cluster number */
				clust_size[i]  = 1;
			}

			/* Loop over atoms */
			for(i=0; (i < nindex); i++) {
				ai = index[i];
				ci = clust_index[i];

				/* Loop over atoms (only half a matrix) */
				for(j=i+1; (j < nindex); j++) {
					cj = clust_index[j];

					/* If they are not in the same cluster already */
					if (ci != cj) {
						aj = index[j];

						/* Compute distance */
						if (bMol) {
							bSame = FALSE;
							for(ii=mols->index[ai]; !bSame && (ii<mols->index[ai+1]); ii++) {
								for(jj=mols->index[aj]; !bSame && (jj<mols->index[aj+1]); jj++) {
									if (bPBC)
										pbc_dx(&pbc,x[ii],x[jj],dx);
									else
										rvec_sub(x[ii],x[jj],dx);
									dx2   = iprod(dx,dx);
									bSame = (dx2 < cut2);
								}
							}
						}
						else {
							if (bPBC)
								pbc_dx(&pbc,x[ai],x[aj],dx);
							else
								rvec_sub(x[ai],x[aj],dx);
							dx2 = iprod(dx,dx);
							bSame = (dx2 < cut2);
						}
						/* If distance less than cut-off */
						if (bSame) {
							/* Merge clusters: check for all atoms whether they are in
							 * cluster cj and if so, put them in ci
							 */
							for(k=0; (k<nindex); k++) {
								if ((clust_index[k] == cj)) {
									if (clust_size[cj] <= 0)
										gmx_fatal(FARGS,"negative cluster size %d for element %d",
												clust_size[cj],cj);
									clust_size[cj]--;
									clust_index[k] = ci;
									clust_size[ci]++;
								}
							}
						}
					}
				}
			}
			n_x++;
			srenew(t_x,n_x);
			t_x[n_x-1] = fr.time*tf;
			srenew(cs_dist,n_x);
			snew(cs_dist[n_x-1],nindex);
			nclust = 0;
			cav    = 0;
			nav    = 0;
			for(i=0; (i<nindex); i++) {
				ci = clust_size[i];
				if (ci > max_clust_size) {
					max_clust_size = ci;
					max_clust_ind  = i;
				}
				if (ci > 0) {
					nclust++;
					cs_dist[n_x-1][ci-1] += 1.0;
					max_size = max(max_size,ci);
					if (ci > 1) {
						cav += ci;
						nav++;
					}
				}
			}
			fprintf(fp,"%14.6e  %10d\n",fr.time,nclust);
			if (nav > 0)
				fprintf(gp,"%14.6e  %10.3f\n",fr.time,cav/nav);
			fprintf(hp, "%14.6e  %10d\n",fr.time,max_clust_size);

			output_cluster_info(mtop, "cluster_test.dat", clust_index, nindex, fr.time);
		}
		/* Analyse velocities, if present */
		if (fr.bV) {
			if (!tpr) {
				if (bTPRwarn) {
					printf("You need a tpr file to analyse temperatures\n");
					bTPRwarn = FALSE;
				}
			}
			else {
				v = fr.v;
				/* Loop over clusters and for each cluster compute 1/2 m v^2 */
				if (max_clust_ind >= 0) {
					ekin = 0;
					for(i=0; (i<nindex); i++)
						if (clust_index[i] == max_clust_ind) {
							ai    = index[i];
							gmx_mtop_atomnr_to_atom(mtop,ai,&atom);
							ekin += 0.5*atom->m*iprod(v[ai],v[ai]);
						}
					temp = (ekin*2.0)/(3.0*tfac*max_clust_size*BOLTZ);
					fprintf(tp,"%10.3f  %10.3f\n",fr.time,temp);
				}
			}
		}
		nframe++;
	} while (read_next_frame(status,&fr));
	close_trx(status);
	fclose(fp);
	fclose(gp);
	fclose(hp);
	fclose(tp);
	if (max_clust_ind >= 0) {
		fp = fopen(mcn,"w");
		fprintf(fp,"[ max_clust ]\n");
		for(i=0; (i<nindex); i++)
			if (clust_index[i] == max_clust_ind) {
				if (bMol) {
					for(j=mols->index[i]; (j<mols->index[i+1]); j++)
						fprintf(fp,"%d\n",j+1);
				}
				else {
					fprintf(fp,"%d\n",index[i]+1);
				}
			}
		fclose(fp);
	}

	/* Look for the smallest entry that is not zero
	 * This will make that zero is white, and not zero is coloured.
	 */
	cmid = 100.0;
	cmax = 0.0;
	for(i=0; (i<n_x); i++)
		for(j=0; (j<max_size); j++) {
			if ((cs_dist[i][j] > 0) && (cs_dist[i][j] < cmid))
				cmid = cs_dist[i][j];
			cmax = max(cs_dist[i][j],cmax);
		}
	fprintf(stderr,"cmid: %g, cmax: %g, max_size: %d\n",cmid,cmax,max_size);
	cmid = 1;
	fp = ffopen(xpm,"w");
	write_xpm3(fp,0,"Cluster size distribution","# clusters",timebuf,"Size",
			n_x,max_size,t_x,t_y,cs_dist,0,cmid,cmax,
			rlo,rmid,rhi,&nlevels);
	fclose(fp);
	cmax = 0.0;
	for(i=0; (i<n_x); i++)
		for(j=0; (j<max_size); j++) {
			cs_dist[i][j] *= (j+1);
			cmax = max(cs_dist[i][j],cmax);
		}
	fprintf(stderr,"cmid: %g, cmax: %g, max_size: %d\n",cmid,cmax,max_size);
	fp = ffopen(xpmw,"w");
	write_xpm3(fp,0,"Weighted cluster size distribution","Fraction",timebuf,"Size",
			n_x,max_size,t_x,t_y,cs_dist,0,cmid,cmax,
			rlo,rmid,rhi,&nlevels);
	fclose(fp);

	fp = xvgropen(histo,"Cluster size distribution","Cluster size","()");
	nhisto = 0;
	fprintf(fp,"%5d  %8.3f\n",0,0.0);
	for(j=0; (j<max_size); j++) {
		real nelem = 0;
		for(i=0; (i<n_x); i++)
			nelem += cs_dist[i][j];
		fprintf(fp,"%5d  %8.3f\n",j+1,nelem/n_x);
		nhisto += (int)((j+1)*nelem/n_x);
	}
	fprintf(fp,"%5d  %8.3f\n",j+1,0.0);
	fclose(fp);

	fprintf(stderr,"Total number of atoms in clusters =  %d\n",nhisto);

	sfree(clust_index);
	sfree(clust_size);
	sfree(index);
}

int gmx_clustsize(int argc,char *argv[])
{
	static char *desc[] = {
		"This program computes the size distributions of molecular/atomic clusters in",
		"the gas phase. The output is given in the form of a XPM file.",
		"The total number of clusters is written to a XVG file.[PAR]",
		"When the [TT]-mol[tt] option is given clusters will be made out of",
		"molecules rather than atoms, which allows clustering of large molecules.",
		"In this case an index file would still contain atom numbers",
		"or your calculcation will die with a SEGV.[PAR]",
		"When velocities are present in your trajectory, the temperature of",
		"the largest cluster will be printed in a separate xvg file assuming",
		"that the particles are free to move. If you are using constraints,",
		"please correct the temperature. For instance water simulated with SHAKE",
		"or SETTLE will yield a temperature that is 1.5 times too low. You can",
		"compensate for this with the -ndf option. Remember to take the removal",
		"of center of mass motion into account.[PAR]",
		"The [TT]-mc[tt] option will produce an index file containing the",
		"atom numbers of the largest cluster."
	};

	static real cutoff   = 0.35;
	static int  nskip    = 0;
	static int  nlevels  = 20;
	static int  ndf      = -1;
	static bool bMol     = FALSE;
	static bool bPBC     = TRUE;
	static rvec rlo      = { 1.0, 1.0, 0.0 };
	static rvec rhi      = { 0.0, 0.0, 1.0 };
	t_pargs pa[] = {
		{ "-cut",      FALSE, etREAL, {&cutoff},
		"Largest distance (nm) to be considered in a cluster" },
		{ "-mol",      FALSE, etBOOL, {&bMol},
		"Cluster molecules rather than atoms (needs tpr file)" },
		{ "-pbc",      FALSE, etBOOL, {&bPBC},
		"Use periodic boundary conditions" },
		{ "-nskip",    FALSE, etINT,  {&nskip},
		"Number of frames to skip between writing" },
		{ "-nlevels",  FALSE, etINT,  {&nlevels},
		"Number of levels of grey in xpm output" },
		{ "-ndf",      FALSE, etINT,  {&ndf},
		"Number of degrees of freedom of the entire system for temperature calculation. If not set the number of atoms times three is used." },
		{ "-rgblo",    FALSE, etRVEC, {rlo},
		"RGB values for the color of the lowest occupied cluster size" },
		{ "-rgbhi",    FALSE, etRVEC, {rhi},
		"RGB values for the color of the highest occupied cluster size" }
	};
#define NPA asize(pa)
	char       *fnNDX,*fnTPR;
	bool       bSQ,bRDF;
	t_rgb      rgblo,rgbhi;

	t_filenm   fnm[] = {
			{ efTRX, "-f",  NULL,         ffREAD  },
			{ efTPR, NULL,  NULL,         ffOPTRD },
			{ efNDX, NULL,  NULL,         ffOPTRD },
			{ efXPM, "-o", "csize",       ffWRITE },
			{ efXPM, "-ow","csizew",      ffWRITE },
			{ efXVG, "-nc","nclust",      ffWRITE },
			{ efXVG, "-mc","maxclust",    ffWRITE },
			{ efXVG, "-ac","avclust",     ffWRITE },
			{ efXVG, "-hc","histo-clust", ffWRITE },
			{ efXVG, "-temp","temp",     ffOPTWR },
			{ efNDX, "-mcn", "maxclust", ffOPTWR }
	};
#define NFILE asize(fnm)

	CopyRight(stderr,argv[0]);
	parse_common_args(&argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME | PCA_TIME_UNIT | PCA_BE_NICE,
			NFILE,fnm,NPA,pa,asize(desc),desc,0,NULL);

	fnNDX = ftp2fn_null(efNDX,NFILE,fnm);
	rgblo.r = rlo[XX],rgblo.g = rlo[YY],rgblo.b = rlo[ZZ];
	rgbhi.r = rhi[XX],rgbhi.g = rhi[YY],rgbhi.b = rhi[ZZ];

	fnTPR = ftp2fn_null(efTPR,NFILE,fnm);
	if (bMol && !fnTPR)
		gmx_fatal(FARGS,"You need a tpr file for the -mol option");

	clust_size(fnNDX, ftp2fn(efTRX,NFILE,fnm), opt2fn("-o",NFILE,fnm),
			opt2fn("-ow",NFILE,fnm),
			opt2fn("-nc",NFILE,fnm),opt2fn("-ac",NFILE,fnm),
			opt2fn("-mc",NFILE,fnm),opt2fn("-hc",NFILE,fnm),
			opt2fn("-temp",NFILE,fnm),opt2fn("-mcn",NFILE,fnm),
			bMol, bPBC, fnTPR,
			cutoff, nskip, nlevels, rgblo, rgbhi, ndf);

	thanx(stderr);

	return 0;
}

int
main(int argc, char *argv[])
{
  gmx_clustsize(argc,argv);
  return 0;
}
