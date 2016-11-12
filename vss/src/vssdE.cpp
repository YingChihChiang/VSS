/*
 * Top level VSS potential routines
 *
 * $Id: vss.c,v 1.4 2005/07/20 15:37:39 johns Exp $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include "vss.h"
#include "pub3dfft.h"

#define M_PI 3.14159265358979323846
#define SQRT_PI 1.7724538509055160273  
#define COLOUMB 332.0636

void general_init(int natoms, double *xyz2, double *atoms, struct GEN *gen_data, struct SAM *SamDihed, struct SAM *SamAngle) {
  // Initialize data, note that when this routine is called, values of xyz in atoms are not yet assigned.
  int nrotate, jobcase;

  if (SamDihed->nentry!=0) {
     if (SamAngle->nentry!=0) {
        nrotate=sumint(SamDihed->nsam,SamDihed->nentry)*sumint(SamAngle->nsam,SamAngle->nentry);
        jobcase=1;
     } else {
        nrotate=sumint(SamDihed->nsam,SamDihed->nentry);
        jobcase=2;
     }
  } else {
     if (SamAngle->nentry!=0) {
        nrotate=sumint(SamAngle->nsam,SamAngle->nentry);
        jobcase=3;
     } else {
        nrotate=1;
        jobcase=4;
     }
  }

  gen_data->nrotate=nrotate;  
  gen_data->jobcase=jobcase;
  gen_data->sel2pdb = xyz2;
  gen_data->Eb = (double*) malloc(nrotate*sizeof(double));
  gen_data->DeltaE = (double*) malloc(nrotate*sizeof(double));

}


void general_free(struct GEN *gen_data) {

  free(gen_data->charge);
  free(gen_data->Eb);
  free(gen_data->DeltaE);

}


void vdw_data_init(double ewald_factor, struct VDW *vdw_data) {
  double Ron, Roff, Ron2, Roff2;
  double const pi=3.14159265359;

  Ron   = vdw_data->Ron;
  Roff  = vdw_data->Roff;
  Ron2  = Ron*Ron;
  Roff2 = Roff*Roff;  

  vdw_data->Ron2 = Ron2;
  vdw_data->Roff2 = Roff2;
  vdw_data->swdenom = (Roff2-Ron2)*(Roff2-Ron2)*(Roff2-Ron2);
  vdw_data->mycoeff = 2.0 * ewald_factor / sqrt(pi);
}


void vdwAB_init(int i, int j, struct VDW *vdw_data, double *A, double *B) {
  double epsilon_ij, Rmin_ij, Rmin2, Rmin6, Rmin12;

  epsilon_ij = sqrt(vdw_data->para[4*i]*vdw_data->para[4*j]);
  Rmin_ij = vdw_data->para[4*i+1]+vdw_data->para[4*j+1];
  Rmin2 = Rmin_ij * Rmin_ij;
  Rmin6 = Rmin2 * Rmin2 * Rmin2;
  Rmin12 = Rmin6 * Rmin6;

  *A = epsilon_ij * Rmin12;
  *B = 2.0 * epsilon_ij * Rmin6;
}


void sel_vdwAB_init(int i, int j, struct VDW *vdw_data, struct EXC *exc_data, double *A, double *B, double *excpair, double *excscale) {
  int k;
  double epsilon_ij, Rmin_ij, Rmin2, Rmin6, Rmin12;

  for (k=exc_data->excdisp3[i];k<exc_data->excdisp3[i+1];k++) {
     if (j==exc_data->exclist3[k]) {
        *A = 0.0;
        *B = 0.0;
        *excpair=1.0;
        *excscale=1.0;
        return;
     }
  }

  for (k=exc_data->excdisp4[i];k<exc_data->excdisp4[i+1];k++) {
     if (j==exc_data->exclist4[k]) {
        epsilon_ij = sqrt(vdw_data->para[4*i+2]*vdw_data->para[4*j+2]);
        Rmin_ij = vdw_data->para[4*i+3]+vdw_data->para[4*j+3];
        Rmin2 = Rmin_ij * Rmin_ij;
        Rmin6 = Rmin2 * Rmin2 * Rmin2;
        Rmin12 = Rmin6 * Rmin6;
        *A = epsilon_ij * Rmin12;
        *B = 2.0 * epsilon_ij * Rmin6;
        *excpair=0.0;
        *excscale=1.0;
        return;
     }
  }

  epsilon_ij = sqrt(vdw_data->para[4*i]*vdw_data->para[4*j]);
  Rmin_ij = vdw_data->para[4*i+1]+vdw_data->para[4*j+1];
  Rmin2 = Rmin_ij * Rmin_ij;
  Rmin6 = Rmin2 * Rmin2 * Rmin2;
  Rmin12 = Rmin6 * Rmin6;
  *A = epsilon_ij * Rmin12;
  *B = 2.0 * epsilon_ij * Rmin6;
  *excpair=0.0;
  *excscale=1.0;

}



void compute_intersel_vdwpmer(pmepot_data *data, const double *cell, int natoms, double *atoms, double *box_disp, struct VDW *vdw_data, struct EXC *exc_data,
                              struct PSE *sel, struct PSE *sel0, double *E_vdw, double *E_pme_r) {
// Evaluate VDW Energy (E_vdw) and PME Real Space Energy (E_pme_r) for pairs from different selections, eg. sel0-sel2.
  double VDWE=0.0, PMERE=0.0;
  double patchsize = 3.5*vdw_data->Roff2;  //3.5

  // Evaluate energy from the other boxes
  for (int i=0; i<sel->natoms; i++) {
     int isel=sel->aid[i];
     for (int j=0; j<sel0->natoms; j++) {
        double e_real;
        double A, B, Ron2, Roff2, swdenom;
        double r, r2, r6, r12, swfactor, itar;
        int jsel=sel0->aid[j];
        vdwAB_init(isel,jsel,vdw_data,&A,&B);
        Ron2 = vdw_data->Ron2;
        Roff2 = vdw_data->Roff2;
        swdenom = vdw_data->swdenom;
        for (int k=0; k<27; k++) {
           r2 = get_r2(atoms,isel,jsel,box_disp[3*k+0],box_disp[3*k+1],box_disp[3*k+2]);
           if (r2 > patchsize) continue;
           r6 = r2*r2*r2;
           r12 = r6*r6;
           r = sqrt(r2);
           if (r2 > Roff2) swfactor=0.0;
           else if (r2 <= Ron2) swfactor=1.0;
           else swfactor = (Roff2-r2)*(Roff2-r2)*(Roff2+2*r2-3*Ron2)/swdenom;
           VDWE += (A/r12 - B/r6)*swfactor; 
           itar = r * data->ewald_factor;
           e_real = erfc(itar);      
           e_real *= COLOUMB*atoms[4*isel+3]*atoms[4*jsel+3]/r;
           PMERE += e_real;
        }
     }
  }

  *E_vdw=VDWE;
  *E_pme_r=PMERE;

}


void compute_sel_vdwpmer(pmepot_data *data, const double *cell, int natoms, double *atoms, struct VDW *vdw_data, struct EXC *exc_data, struct PSE *sel, double *E_vdw, double *E_pme_r) {
// Evaluate VDW Energy (E_vdw) and PME Real Space Energy (E_pme_r) for a selected ligand, eg. sel1 or sel2.
  double VDWE=0.0, PMERE=0.0;

  for (int i=0; i<sel->natoms; i++) {
     int isel=sel->aid[i];
     for (int j=0; j<i; j++) {
        double e_real, excpair, excscale;
        double A, B, Ron2, Roff2, swdenom;
        double r, r2, r6, r12, swfactor, itar;
        int jsel=sel->aid[j];
        sel_vdwAB_init(isel,jsel,vdw_data,exc_data,&A,&B,&excpair,&excscale);
        Ron2 = vdw_data->Ron2;
        Roff2 = vdw_data->Roff2;
        swdenom = vdw_data->swdenom;
        r2 = get_r2(atoms,isel,jsel,0,0,0);
        r6 = r2*r2*r2;
        r12 = r6*r6;
        r = sqrt(r2);
        if (r2 > Roff2) swfactor=0.0;
        else if (r2 <= Ron2) swfactor=1.0;
        else swfactor = (Roff2-r2)*(Roff2-r2)*(Roff2+2*r2-3*Ron2)/swdenom;
        VDWE += (A/r12 - B/r6)*swfactor;              // A=B=0 if the pair is in the 1-3 excluded list, and A14 B14 is used if the pair is in 1-4 scaling list..
        itar = r * data->ewald_factor;
        e_real = erfc(itar) - excpair;        // Substract the contribtion to ELECT if the pair is in 1-3 excluded list.
        e_real *= COLOUMB*atoms[4*isel+3]*atoms[4*jsel+3]/r;
        PMERE += e_real;
     }
  }
  
  *E_vdw=VDWE;
  *E_pme_r=PMERE;

}


void compute_sel_nonb(pmepot_data *data, const double *cell, int natoms, double *atoms, struct VDW *vdw_data, struct EXC *exc_data, struct PSE *sel, double *E_vdw, double *E_elec) {
// Evaluate the nonbonded interactions (E_vdw, E_elec) for a selected ligand, e.g. sel1 or sel2. 
  *E_vdw=0.0;
  *E_elec=0.0;

  for (int i=0; i<sel->natoms; i++) {
     int isel=sel->aid[i];
     for (int j=0; j<i; j++) {
        double e_real, excpair, excscale;
        double A, B, Ron2, Roff2, swdenom;
        double r, r2, r6, r12, swfactor, itar;
        int jsel=sel->aid[j];
        sel_vdwAB_init(isel,jsel,vdw_data,exc_data,&A,&B,&excpair,&excscale);
        Ron2 = vdw_data->Ron2;
        Roff2 = vdw_data->Roff2;
        swdenom = vdw_data->swdenom;
        r2 = get_r2(atoms,isel,jsel,0,0,0);
        r6 = r2*r2*r2;
        r12 = r6*r6;
        r = sqrt(r2);
        if (r2 > Roff2) swfactor=0.0;
        else if (r2 <= Ron2) swfactor=1.0;
        else swfactor = (Roff2-r2)*(Roff2-r2)*(Roff2+2*r2-3*Ron2)/swdenom;
        *E_vdw  += (A/r12 - B/r6)*swfactor;              // A=B=0 if the pair is in the 1-3 excluded list, and A14 B14 is used if the pair is in 1-4 scaling list.
        e_real = 1.0 - excpair;                         // Substract the contribtion to ELECT if the pair is in 1-3 excluded list.
        e_real *= COLOUMB*atoms[4*isel+3]*atoms[4*jsel+3]/r;
        *E_elec += e_real;
     }
  }
}


double compute_pmek(pmepot_data *data, const double *cell, int natoms, double *atoms, double *charge,struct PSE *othersel) {
// Evaluate E_pmek for a selected system, e.g. sel1 + sel0.
  double rcell[12];
  double RecipE=0.0;

  // Sel (Set Other to 0)
  for (int i=0;i<natoms;i++) {
     atoms[4*i+3]=charge[i];
     for (int j=0;j<othersel->natoms;j++) {
        if (i==othersel->aid[j]) atoms[4*i+3]=0.0;
     }
  }

  fill_charges(data->dims,cell,natoms,atoms,data->q_arr,rcell,data->oddd);

  pubdz3d(1, data->dims[2], data->dims[1], data->dims[0],
    data->q_arr, data->dims[4], data->dims[3],
    data->fft_table, data->fft_ntable, data->fft_work);

  RecipE=compute_energy(data->q_arr, cell, rcell, data->dims, data->ewald_factor);
  RecipE *= COLOUMB;

  // Restore Charge 
  for (int i=0;i<natoms;i++) {
     atoms[4*i+3]=charge[i];
  }

  return RecipE;
}


double compute_pmeselfdiff(pmepot_data *data, double *atoms, struct PSE *sel1, struct PSE *sel2) {
// Evaluate K-space Ewald Self Energy difference between sel2 and sel1
  double SelfE=0.0;
  int i, aid;

  // sum over sel2
  SelfE=0.0;
  for (i=0;i<sel2->natoms;i++) {
     aid=sel2->aid[i];
     SelfE += atoms[4*aid+3]*atoms[4*aid+3];
  }
  // substrate over sel1
  for (i=0;i<sel1->natoms;i++) {
     aid=sel1->aid[i];
     SelfE -= atoms[4*aid+3]*atoms[4*aid+3];
  }
  SelfE *= -(data->ewald_factor * COLOUMB / SQRT_PI);

  return SelfE;
}


double get_r2(double const *xyzq, int i, int j, double box_disp_x, double box_disp_y, double box_disp_z) {
  double x,y,z;
  double r2;
  
  x = xyzq[4*i+0] - xyzq[4*j+0] - box_disp_x;
  y = xyzq[4*i+1] - xyzq[4*j+1] - box_disp_y;
  z = xyzq[4*i+2] - xyzq[4*j+2] - box_disp_z;
 
  r2=x*x+y*y+z*z;
  return r2;
}


void get_box_disp(double const *cell, double *box_disp) {
  int m1, m2, m3, mt,count;

  // main box's displacement is zero.
  box_disp[0]= 0.0;
  box_disp[1]= 0.0;
  box_disp[2]= 0.0;

  count=1;
  for (m1=-1;m1<2;m1++) {
     for (m2=-1;m2<2;m2++) {
        for (m3=-1;m3<2;m3++) {
           mt=m1*m1+m2*m2+m3*m3;
           if(mt==0) continue;
           box_disp[3*count+0]=m1*cell[3] + m2*cell[6] + m3*cell[9];
           box_disp[3*count+1]=m1*cell[4] + m2*cell[7] + m3*cell[10];
           box_disp[3*count+2]=m1*cell[5] + m2*cell[8] + m3*cell[11];
           count += 1;
        }
     }
  }

  return;
}

