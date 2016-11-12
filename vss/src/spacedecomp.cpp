#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "vss.h"

#define COLOUMB 332.0636


void patch_init(int natoms, double Roff, double const *cell, struct PSE *sel0, struct PSE *sel1, struct PSE *sel2, struct PAT *patch_data) {
  int i;
 
  // determine the number of bins in each dimension
  patch_data->flag_spacedecomp=1;
  for (i=0;i<3;i++) {
     patch_data->nbin[i]=cell[i*4+3]/Roff;
     if (patch_data->nbin[i]==1) patch_data->flag_spacedecomp=0;
     if (patch_data->nbin[i]==0) {printf("Error in patch_init: cell length = 0.\n"); printf("Do not use space decomposition.\n");  patch_data->flag_spacedecomp=0;} 
  }
  patch_data->nxy=patch_data->nbin[0]*patch_data->nbin[1];  
  patch_data->npatch=patch_data->nxy*patch_data->nbin[2];
  
  // allocate arrays
  patch_data->workcount = (int*) malloc(patch_data->npatch*sizeof(int));
  patch_data->aid2pid = (int*) malloc(natoms*sizeof(int));
  patch_data->tempxyzq = (double*) malloc(natoms*4*sizeof(double));  
  sel0->patchnatoms = (int*) malloc(patch_data->npatch*sizeof(int));
  sel1->patchnatoms = (int*) malloc(patch_data->npatch*sizeof(int));
  sel2->patchnatoms = (int*) malloc(patch_data->npatch*sizeof(int));
  sel0->patchdisp = (int*) malloc(patch_data->npatch*sizeof(int));
  sel1->patchdisp = (int*) malloc(patch_data->npatch*sizeof(int));
  sel2->patchdisp = (int*) malloc(patch_data->npatch*sizeof(int));
  sel0->patch2aid = (int*) malloc(sel0->natoms*sizeof(int));
  sel1->patch2aid = (int*) malloc(sel1->natoms*sizeof(int));
  sel2->patch2aid = (int*) malloc(sel2->natoms*sizeof(int));
  sel0->patchxyzq = (double*) malloc(sel0->natoms*4*sizeof(double));
  sel1->patchxyzq = (double*) malloc(sel1->natoms*4*sizeof(double));
  sel2->patchxyzq = (double*) malloc(sel2->natoms*4*sizeof(double));
  sel0->patchvdwp = (double*) malloc(sel0->natoms*4*sizeof(double));
  sel1->patchvdwp = (double*) malloc(sel1->natoms*4*sizeof(double));
  sel2->patchvdwp = (double*) malloc(sel2->natoms*4*sizeof(double));

}


void box_update(double const *cell, struct PAT *patch_data) {
  int i;

  // determine the patch size and box left end coordinates
  for (i=0;i<3;i++) {
     patch_data->patchsize[i]=cell[i*4+3]/patch_data->nbin[i];  // updated in each frame
     patch_data->boxlend[i]=cell[i]-0.5*cell[i*4+3];            // updated in each frame
     patch_data->boxrend[i]=cell[i]+0.5*cell[i*4+3];            // updated in each frame
  }

}


void binatoms(double const *cell, double *xyzq, struct VDW *vdw_data, struct PAT *patch_data, struct PSE *sel) {
// bin the atoms of a given selection, note that the xyz should be strictly confined in the box.
  int isel, iseldisp, pid, pdisp, quot[3];
  double length[3];

  for (int i=0;i<patch_data->npatch;i++) {
     sel->patchnatoms[i] = 0;
     sel->patchdisp[i] = 0;
     patch_data->workcount[i] = 0;
  }

  // put all atoms into the box
  for (int j=0;j<3;j++) {
     for (int i=0;i<sel->natoms;i++) {
        isel=sel->aid[i];
        iseldisp = 4*isel;
        length[j] = xyzq[iseldisp+j]-patch_data->boxlend[j];
        patch_data->tempxyzq[iseldisp+j]=xyzq[iseldisp+j];
        if (length[j]<0) {patch_data->tempxyzq[iseldisp+j] += cell[j*4+3]; }
        if (length[j]-cell[j*4+3] > 0) {patch_data->tempxyzq[iseldisp+j] -= cell[j*4+3]; }
     }
  }

  // bin the atoms, first build the patchnatoms
  for (int i=0;i<sel->natoms;i++) {
     isel = sel->aid[i];
     iseldisp = 4*isel;

     for (int j=0;j<3;j++) {
        length[j] = patch_data->tempxyzq[iseldisp+j]-patch_data->boxlend[j];
        quot[j] = length[j] / patch_data->patchsize[j];
        if (quot[j]<0) {quot[j] = 0;};   // For corrupted dcd, particles can be 1 box away!
        if (quot[j]>=patch_data->nbin[j]) {quot[j] = patch_data->nbin[j]-1;};  
     }

     pid = quot[0] + quot[1]*patch_data->nbin[0] + quot[2]*patch_data->nxy;

     patch_data->aid2pid[isel] = pid;
     patch_data->tempxyzq[iseldisp+3] = xyzq[iseldisp+3];

     sel->patchnatoms[pid]++;
  }

  // build patchdisp
  for (int i=1;i<patch_data->npatch;i++) {
     sel->patchdisp[i]=sel->patchdisp[i-1]+sel->patchnatoms[i-1];
  }

  // build patch2aid, patchxyzq, and patchvdwp
  for (int i=0;i<sel->natoms;i++) {
     isel = sel->aid[i];
     iseldisp = 4*isel;
     pid = patch_data->aid2pid[isel];
     pdisp = sel->patchdisp[pid] + patch_data->workcount[pid];
     sel->patch2aid[pdisp]=isel;
     for (int j=0;j<4;j++) {
        sel->patchxyzq[4*pdisp+j] = patch_data->tempxyzq[iseldisp+j];
     }
     for (int j=0;j<4;j++) {
        sel->patchvdwp[4*pdisp+j] = vdw_data->para[iseldisp+j];
     }
     patch_data->workcount[pid]++;
  }

}


void getneighbor(int pid, const double *cell, struct PAT *patch_data) {
// construct the neighbor patch list (neighbor,neighbor_nr) for a given pid
  int i, j, k, iwork, count, pindex[3], nindex[3], ndisp[3];

  iwork = pid%patch_data->nxy;
  pindex[2] = pid/patch_data->nxy;
  pindex[1] = iwork/patch_data->nbin[0];
  pindex[0] = iwork%patch_data->nbin[0];
 
  count=0;
  for (i=pindex[2]-1;i<pindex[2]+2;i++) {
      nindex[2]=i;
      ndisp[2]=0;
      if (i<0) {nindex[2]=i+patch_data->nbin[2];ndisp[2]=-1;}
      if (i==patch_data->nbin[2]) {nindex[2]=i-patch_data->nbin[2];ndisp[2]=1;}

     for (j=pindex[1]-1;j<pindex[1]+2;j++) {
        nindex[1]=j;
        ndisp[1]=0;
        if (j<0) {nindex[1]=j+patch_data->nbin[1];ndisp[1]=-1;}
        if (j==patch_data->nbin[1]) {nindex[1]=j-patch_data->nbin[1];ndisp[1]=1;}

        for (k=pindex[0]-1;k<pindex[0]+2;k++) {
            nindex[0]=k;
            ndisp[0]=0;
            if (k<0) {nindex[0]=k+patch_data->nbin[0];ndisp[0]=-1;}
            if (k==patch_data->nbin[0]) {nindex[0]=k-patch_data->nbin[0];ndisp[0]=1;}
            
            patch_data->neighbor[count] = nindex[0]+nindex[1]*patch_data->nbin[0]+nindex[2]*patch_data->nxy;
            patch_data->neighbor_nr[3*count]   = ndisp[0]*cell[3];
            patch_data->neighbor_nr[3*count+1] = ndisp[1]*cell[7];
            patch_data->neighbor_nr[3*count+2] = ndisp[2]*cell[11];
            count += 1;
        }
     }
  }

}


void compute_patch_vdwpmer(const double *cell, double ewald_factor, struct VDW *vdw_data, struct PAT *patch_data, struct PSE *sel, struct PSE *sel0, double *E_vdw, double *E_pmer) {
  double VDWE=0.0, PMERE=0.0;

  for (int I=0;I<patch_data->npatch;I++) {
     if (sel->patchnatoms[I]==0) continue;
     getneighbor(I,cell,patch_data);

     for (int J=0;J<27;J++) {
        int pid=patch_data->neighbor[J];
        int Jdisp=3*J;
   
        for (int i=0;i<sel->patchnatoms[I];i++) {
           int idisp=4*(sel->patchdisp[I]+i);
   
           for (int j=0;j<sel0->patchnatoms[pid];j++) {
              int jdisp=4*(sel0->patchdisp[pid]+j);
              double A, B, Ron2, Roff2, swdenom, swfactor, itar;
              double r, r2, r6, r12;
              Ron2 = vdw_data->Ron2;
              Roff2 = vdw_data->Roff2;
              swdenom = vdw_data->swdenom;
              patch_vdwAB_init(sel->patchvdwp[idisp],sel0->patchvdwp[jdisp],sel->patchvdwp[idisp+1],sel0->patchvdwp[jdisp+1],&A,&B);
              r2 = patch_getr2(&sel->patchxyzq[idisp],&sel0->patchxyzq[jdisp],&patch_data->neighbor_nr[Jdisp]);
              if (r2 > Roff2) continue;  // This will introduce a negligible error on electrostatic potential, but speed up the calculation. 
              r6 = r2*r2*r2;
              r12 = r6*r6;
              r = sqrt(r2);
              if (r2 > Ron2) {swfactor=(Roff2-r2)*(Roff2-r2)*(Roff2+2*r2-3*Ron2)/swdenom;}
              else {swfactor=1.0;}
              VDWE += (A/r12 - B/r6)*swfactor;
              itar = r * ewald_factor;
              PMERE += erfc(itar)*COLOUMB*sel->patchxyzq[idisp+3]*sel0->patchxyzq[jdisp+3]/r;
           }
        }
     }
  }

  *E_vdw=VDWE;
  *E_pmer=PMERE;

}


void patch_vdwAB_init(double epsilon_i, double epsilon_j, double R_i, double R_j, double *A, double *B) {
  double epsilon_ij, Rmin_ij, Rmin2, Rmin6, Rmin12;

  epsilon_ij = sqrt(epsilon_i*epsilon_j);
  Rmin_ij = R_i+R_j;
  Rmin2 = Rmin_ij * Rmin_ij;
  Rmin6 = Rmin2 * Rmin2 * Rmin2;
  Rmin12 = Rmin6 * Rmin6;

  *A = epsilon_ij * Rmin12;
  *B = 2.0 * epsilon_ij * Rmin6;
}


double patch_getr2(double *ixyz, double *jxyz, double *nr) {
// Note that j is shifted. (j+nr) is the correct position.
  double x, y, z, r2;

  x=ixyz[0]-jxyz[0]-nr[0];
  y=ixyz[1]-jxyz[1]-nr[1];
  z=ixyz[2]-jxyz[2]-nr[2];
  r2=x*x+y*y+z*z;

  return r2;
}


void patch_free(struct PAT *patch_data, struct PSE *sel0, struct PSE *sel1, struct PSE *sel2) {

  free(patch_data->workcount);
  free(patch_data->aid2pid);
  free(patch_data->tempxyzq);  
  free(sel0->patchnatoms);
  free(sel1->patchnatoms);
  free(sel2->patchnatoms);
  free(sel0->patchdisp);
  free(sel1->patchdisp);
  free(sel2->patchdisp);
  free(sel0->patch2aid);
  free(sel1->patch2aid);
  free(sel2->patch2aid);
  free(sel0->patchxyzq);
  free(sel1->patchxyzq);
  free(sel2->patchxyzq);
  free(sel0->patchvdwp);
  free(sel1->patchvdwp);
  free(sel2->patchvdwp);

}
