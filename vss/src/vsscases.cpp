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

#define M_PI 3.14159265358979323846
#define SQRT_PI 1.7724538509055160273  
#define COLOUMB 332.0636

void vsscase1(const double *cell, double *atoms, struct GEN *gen_data, pmepot_data *data, struct VDW *vdw_data, 
           struct PSE *sel0, struct PSE *sel1, struct PSE *sel2, struct EXC *exc_data, struct MUT *mut_data, 
           struct PAT *patch_data, struct BOND *Dihed, struct BOND *Angle, struct SAM *SamDihed, struct SAM *SamAngle) 
{
  // Evaluate DeltaE for each frame, case 1, by Letitia: 06/08/2014, 12/11/2015, 23/01/2016
  //int i, j, k, myindex, nrotate, natoms;
  int j, k, nrotate, natoms;
  int dihed_id, angle_id;
  double E1_PME, E1_VDW, E2_PME, E2_VDW, SelfE_diff;
  double Evdw, Eelec;
  double *charge, *sel2pdb, *Eb, *DeltaE;
  double Ebmin=0.0;
  double E_vdw_sel, E_pme_sel, E_vdw_pair, E_pme_pair, E_work, box_disp[81];

  charge=gen_data->charge;
  sel2pdb=gen_data->sel2pdb;
  Eb=gen_data->Eb;
  DeltaE=gen_data->DeltaE;
  natoms=gen_data->natoms;

  if (patch_data->flag_spacedecomp) {
     box_update(cell,patch_data);
     binatoms(cell,atoms,vdw_data,patch_data,sel0);
     binatoms(cell,atoms,vdw_data,patch_data,sel1);
  } else {
     get_box_disp(cell,box_disp);
  }

  // Compute energy of sel1 + environment (sel0)
  if (sel0->natoms==0) {
     compute_sel_nonb(data,cell,natoms,atoms,vdw_data,exc_data,sel1,&Evdw,&Eelec);
     E1_PME = Eelec;
     E1_VDW = Evdw;
     SelfE_diff = 0.0;
  } else { 
     compute_sel_vdwpmer(data,cell,natoms,atoms,vdw_data,exc_data,sel1,&E_vdw_sel,&E_pme_sel);
     if (patch_data->flag_spacedecomp) compute_patch_vdwpmer(cell,data->ewald_factor,vdw_data,patch_data,sel1,sel0,&E_vdw_pair,&E_pme_pair);
     else compute_intersel_vdwpmer(data,cell,natoms,atoms,box_disp,vdw_data,exc_data,sel1,sel0,&E_vdw_pair,&E_pme_pair);
     E_work = compute_pmek(data,cell,natoms,atoms,charge,sel2);
     E1_PME = E_pme_sel+E_pme_pair+E_work;
     E1_VDW = E_vdw_sel+E_vdw_pair;
     // Compute the PME selfenergy difference between sel2 and sel1. (Constant for all geometry.)
     SelfE_diff = compute_pmeselfdiff(data,atoms,sel1,sel2);
  }

#ifdef Ham
  // Ghost ligand correction (Using a full molecule for ghost ligand)
  compute_sel_nonb(data,cell,natoms,atoms,vdw_data,exc_data,sel1,&Evdw,&Eelec);
  E1_PME -= Eelec;
  E1_VDW -= Evdw;
#endif

  nrotate=0;
  // SamDihed->nentry!=0, SamAngle->nentry!=0
  for (dihed_id=0;dihed_id<SamDihed->nentry;dihed_id++) {

     for (j=0;j<SamDihed->nsam[dihed_id];j++) {

        gen_sel2_geo(sel1->natoms,sel2->natoms,atoms,sel2pdb,mut_data);
        rot_sel2_dihed(sel2->natoms,dihed_id,SamDihed->samangle[j],&atoms[4*sel1->natoms],SamDihed);    

        for (angle_id=0;angle_id<SamAngle->nentry;angle_id++) {
           // Reset sel2 xyz for sampling of other angles (angle_id!=0)
           if (angle_id > 0) { 
              gen_sel2_geo(sel1->natoms,sel2->natoms,atoms,sel2pdb,mut_data); 
              rot_sel2_dihed(sel2->natoms,dihed_id,SamDihed->samangle[j],&atoms[4*sel1->natoms],SamDihed);
           }
           for (k=0;k<SamAngle->nsam[angle_id];k++) {
              rot_sel2_angle(sel2->natoms,angle_id,SamAngle->samangle[k],&atoms[4*sel1->natoms],SamAngle);
     
              Eb[nrotate] = calc_Eb(&atoms[4*sel1->natoms],Dihed,Angle);
              if (Eb[nrotate]-Ebmin < 0) {Ebmin=Eb[nrotate];}  // For numerical stability
        
              // For valid rotation: Update the patch_data for sel2
              if (patch_data->flag_spacedecomp) {binatoms(cell,atoms,vdw_data,patch_data,sel2);}
        
              // For valid rotation: Compute energy of sel2 + environment (sel0)
              if (sel0->natoms==0) {
                 compute_sel_nonb(data,cell,natoms,atoms,vdw_data,exc_data,sel2,&Evdw,&Eelec);
                 E2_PME = Eelec;
                 E2_VDW = Evdw;
              } else { 
                 compute_sel_vdwpmer(data,cell,natoms,atoms,vdw_data,exc_data,sel2,&E_vdw_sel,&E_pme_sel);
                 if (patch_data->flag_spacedecomp) compute_patch_vdwpmer(cell,data->ewald_factor,vdw_data,patch_data,sel2,sel0,&E_vdw_pair,&E_pme_pair);
                 else compute_intersel_vdwpmer(data,cell,natoms,atoms,box_disp,vdw_data,exc_data,sel2,sel0,&E_vdw_pair,&E_pme_pair);
                 E_work = compute_pmek(data,cell,natoms,atoms,charge,sel1);
                 E2_PME = E_pme_sel+E_pme_pair+E_work;
                 E2_VDW = E_vdw_sel+E_vdw_pair;
              }

#ifdef Ham        
              // Ghost ligand correction (Using a full molecule for ghost ligand)
              compute_sel_nonb(data,cell,natoms,atoms,vdw_data,exc_data,sel2,&Evdw,&Eelec);
              Eb[nrotate] += Evdw+Eelec;
              E2_PME -= Eelec;
              E2_VDW -= Evdw;
#endif

              // Assign DeltaE (rotated structure) and update nrotate
              DeltaE[nrotate] = (E2_PME + E2_VDW)-(E1_PME + E1_VDW)+SelfE_diff;
              nrotate++;
           }
        }
     }
  }

  // Return DeltaE, Eb, Ebmin via gen_data
  gen_data->Ebmin=Ebmin;
}


void vsscase2(const double *cell, double *atoms, struct GEN *gen_data, pmepot_data *data, struct VDW *vdw_data, 
           struct PSE *sel0, struct PSE *sel1, struct PSE *sel2, struct EXC *exc_data, struct MUT *mut_data, 
           struct PAT *patch_data, struct BOND *Dihed, struct BOND *Angle, struct SAM *SamDihed, struct SAM *SamAngle) 
{
  // Evaluate DeltaE for each frame, case 2, by Letitia: 06/08/2014, 12/11/2015, 23/01/2016
  //int i, j, k, myindex, nrotate, natoms;
  int j, nrotate, natoms;
  int dihed_id; //angle_id;
  double E1_PME, E1_VDW, E2_PME, E2_VDW, SelfE_diff;
  double Evdw, Eelec;
  double *charge, *sel2pdb, *Eb, *DeltaE;
  double Ebmin=0.0;
  double E_vdw_sel, E_pme_sel, E_vdw_pair, E_pme_pair, E_work, box_disp[81];

  charge=gen_data->charge;
  sel2pdb=gen_data->sel2pdb;
  Eb=gen_data->Eb;
  DeltaE=gen_data->DeltaE;
  natoms=gen_data->natoms;

  if (patch_data->flag_spacedecomp) {
     box_update(cell,patch_data);
     binatoms(cell,atoms,vdw_data,patch_data,sel0);
     binatoms(cell,atoms,vdw_data,patch_data,sel1);
  } else {
     get_box_disp(cell,box_disp);
  }

  // Compute energy of sel1 + environment (sel0)
  if (sel0->natoms==0) {
     compute_sel_nonb(data,cell,natoms,atoms,vdw_data,exc_data,sel1,&Evdw,&Eelec);
     E1_PME = Eelec;
     E1_VDW = Evdw;
     SelfE_diff = 0.0;
  } else { 
     compute_sel_vdwpmer(data,cell,natoms,atoms,vdw_data,exc_data,sel1,&E_vdw_sel,&E_pme_sel);
     if (patch_data->flag_spacedecomp) compute_patch_vdwpmer(cell,data->ewald_factor,vdw_data,patch_data,sel1,sel0,&E_vdw_pair,&E_pme_pair);
     else compute_intersel_vdwpmer(data,cell,natoms,atoms,box_disp,vdw_data,exc_data,sel1,sel0,&E_vdw_pair,&E_pme_pair);
     E_work = compute_pmek(data,cell,natoms,atoms,charge,sel2);
     E1_PME = E_pme_sel+E_pme_pair+E_work;
     E1_VDW = E_vdw_sel+E_vdw_pair;
     // Compute the PME selfenergy difference between sel2 and sel1. (Constant for all geometry.)
     SelfE_diff = compute_pmeselfdiff(data,atoms,sel1,sel2);
  }

/*
  // Here Letitia!!
  int shift;
  double E_vdw_sel0=0.0, E_pme_sel0=0.0, Eb_sel0;
  double patchsize = 3.5*vdw_data->Roff2;  //3.5
  shift = sel1->natoms+sel2->natoms;
  get_box_disp(cell,box_disp);
  for (int i=shift; i<natoms; i++) {
     //for (int j=shift; j<natoms; j++) {
     for (int j=shift; j<i; j++) {
        double e_real, excpair, excscale;
        double A, B, Ron2, Roff2, swdenom;
        double r, r2, r6, r12, swfactor, itar;
        excscale = 1.0;
        excpair = 0.0;
        // if (j==i) {excscale = 0.0; excpair =1.0;}
        if (j==i+1) {excscale = 0.0; excpair =1.0;}
        if (j==i+2) {excscale = 0.0; excpair =1.0;}
        vdwAB_init(i,j,vdw_data,&A,&B);
        Ron2 = vdw_data->Ron2;
        Roff2 = vdw_data->Roff2;
        swdenom = vdw_data->swdenom;
        for (int k=0; k<27; k++) {
           r2 = get_r2(atoms,i,j,box_disp[3*k+0],box_disp[3*k+1],box_disp[3*k+2]);
           if (r2 > patchsize) continue;
           r6 = r2*r2*r2;
           r12 = r6*r6;
           r = sqrt(r2);
           if (r2 > Roff2) swfactor=0.0;
           else if (r2 <= Ron2) swfactor=1.0;
           else swfactor = (Roff2-r2)*(Roff2-r2)*(Roff2+2*r2-3*Ron2)/swdenom;
           E_vdw_sel0 += (A/r12 - B/r6)*swfactor*excscale;              // A=B=0 if the pair is in the 1-3 excluded list, and A14 B14 is used if the pair is in 1-4 scaling list..
           itar = r * data->ewald_factor;
           e_real = erfc(itar) - excpair;        // Substract the contribtion to ELECT if the pair is in 1-3 excluded list.
           e_real *= COLOUMB*atoms[4*i+3]*atoms[4*j+3]/r;
           E_pme_sel0 += e_real;
        }
     }
  }
  // subtract self energy
  double SelfE=0.0;
  for (int i=shift; i<natoms; i++) {
     SelfE += atoms[4*i+3]*atoms[4*i+3];
  }
  for (int i=0; i<sel1->natoms; i++) {
     SelfE += atoms[4*i+3]*atoms[4*i+3];
  }
  SelfE *= -(data->ewald_factor * COLOUMB / SQRT_PI);

  printf("PME and VDW: %f %f \n",E1_PME + E_pme_sel0 - SelfE, E1_VDW + E_vdw_sel0);
  fflush(stdout);
  // Letitia, Here! 
*/
 
#ifdef Ham
  // Ghost ligand correction (Using a full molecule for ghost ligand)
  compute_sel_nonb(data,cell,natoms,atoms,vdw_data,exc_data,sel1,&Evdw,&Eelec);
  E1_PME -= Eelec;
  E1_VDW -= Evdw;
#endif

  nrotate=0;
  // SamDihed->nentry!=0, SamAngle->nentry=0
  for (dihed_id=0;dihed_id<SamDihed->nentry;dihed_id++) {

     for (j=0;j<SamDihed->nsam[dihed_id];j++) {

        gen_sel2_geo(sel1->natoms,sel2->natoms,atoms,sel2pdb,mut_data);
        rot_sel2_dihed(sel2->natoms,dihed_id,SamDihed->samangle[j],&atoms[4*sel1->natoms],SamDihed);

        Eb[nrotate] = calc_Eb(&atoms[4*sel1->natoms],Dihed,Angle);
        if (Eb[nrotate]-Ebmin < 0) {Ebmin=Eb[nrotate];}  // For numerical stability

        // For valid rotation: Update the patch_data for sel2
        if (patch_data->flag_spacedecomp) {binatoms(cell,atoms,vdw_data,patch_data,sel2);}

        // For valid rotation: Compute energy of sel2 + environment (sel0)
        if (sel0->natoms==0) {
           compute_sel_nonb(data,cell,natoms,atoms,vdw_data,exc_data,sel2,&Evdw,&Eelec);
           E2_PME = Eelec;
           E2_VDW = Evdw; 
        } else { 
           compute_sel_vdwpmer(data,cell,natoms,atoms,vdw_data,exc_data,sel2,&E_vdw_sel,&E_pme_sel);
           if (patch_data->flag_spacedecomp) compute_patch_vdwpmer(cell,data->ewald_factor,vdw_data,patch_data,sel2,sel0,&E_vdw_pair,&E_pme_pair);
           else compute_intersel_vdwpmer(data,cell,natoms,atoms,box_disp,vdw_data,exc_data,sel2,sel0,&E_vdw_pair,&E_pme_pair);
           E_work = compute_pmek(data,cell,natoms,atoms,charge,sel1);
           E2_PME = E_pme_sel+E_pme_pair+E_work;
           E2_VDW = E_vdw_sel+E_vdw_pair;
        }
#ifdef Ham
        // Ghost ligand correction (Using a full molecule for ghost ligand)
        compute_sel_nonb(data,cell,natoms,atoms,vdw_data,exc_data,sel2,&Evdw,&Eelec);
        Eb[nrotate] += Evdw+Eelec;
        E2_PME -= Eelec;
        E2_VDW -= Evdw;
#endif

        // Assign DeltaE (rotated structure) and update nrotate
        DeltaE[nrotate] = (E2_PME + E2_VDW)-(E1_PME + E1_VDW)+SelfE_diff;
        nrotate++;
       
     }
  }

  // Return DeltaE, Eb, Ebmin via gen_data
  gen_data->Ebmin=Ebmin;
}


void vsscase3(const double *cell, double *atoms, struct GEN *gen_data, pmepot_data *data, struct VDW *vdw_data, 
           struct PSE *sel0, struct PSE *sel1, struct PSE *sel2, struct EXC *exc_data, struct MUT *mut_data, 
           struct PAT *patch_data, struct BOND *Dihed, struct BOND *Angle, struct SAM *SamDihed, struct SAM *SamAngle) 
{
  // Evaluate DeltaE for each frame, case 3, by Letitia: 06/08/2014, 12/11/2015, 23/01/2016
  //int i, j, k, myindex, nrotate, natoms;
  int k, nrotate, natoms;
  //int dihed_id, angle_id;
  int angle_id;
  double E1_PME, E1_VDW, E2_PME, E2_VDW, SelfE_diff;
  double Evdw, Eelec;
  double *charge, *sel2pdb, *Eb, *DeltaE;
  double Ebmin=0.0;
  double E_vdw_sel, E_pme_sel, E_vdw_pair, E_pme_pair, E_work, box_disp[81];

  charge=gen_data->charge;
  sel2pdb=gen_data->sel2pdb;
  Eb=gen_data->Eb;
  DeltaE=gen_data->DeltaE;
  natoms=gen_data->natoms;

  if (patch_data->flag_spacedecomp) {
     box_update(cell,patch_data);
     binatoms(cell,atoms,vdw_data,patch_data,sel0);
     binatoms(cell,atoms,vdw_data,patch_data,sel1);
  } else {
     get_box_disp(cell,box_disp);
  }

  // Compute energy of sel1 + environment (sel0)
  if (sel0->natoms==0) {
     compute_sel_nonb(data,cell,natoms,atoms,vdw_data,exc_data,sel1,&Evdw,&Eelec);
     E1_PME = Eelec;
     E1_VDW = Evdw;
     SelfE_diff = 0.0;
  } else { 
     compute_sel_vdwpmer(data,cell,natoms,atoms,vdw_data,exc_data,sel1,&E_vdw_sel,&E_pme_sel);
     if (patch_data->flag_spacedecomp) compute_patch_vdwpmer(cell,data->ewald_factor,vdw_data,patch_data,sel1,sel0,&E_vdw_pair,&E_pme_pair);
     else compute_intersel_vdwpmer(data,cell,natoms,atoms,box_disp,vdw_data,exc_data,sel1,sel0,&E_vdw_pair,&E_pme_pair);
     E_work = compute_pmek(data,cell,natoms,atoms,charge,sel2);
     E1_PME = E_pme_sel+E_pme_pair+E_work;
     E1_VDW = E_vdw_sel+E_vdw_pair;
     // Compute the PME selfenergy difference between sel2 and sel1. (Constant for all geometry.)
     SelfE_diff = compute_pmeselfdiff(data,atoms,sel1,sel2);
  }

#ifdef Ham
  // Ghost ligand correction (Using a full molecule for ghost ligand)
  compute_sel_nonb(data,cell,natoms,atoms,vdw_data,exc_data,sel1,&Evdw,&Eelec);
  E1_PME -= Eelec;
  E1_VDW -= Evdw; 
#endif

  nrotate=0;
  // SamDihed->nentry=0, SamAngle->nentry!=0
  // generate sel2 xyz (Only Angle sampling)
  gen_sel2_geo(sel1->natoms,sel2->natoms,atoms,sel2pdb,mut_data);

  for (angle_id=0;angle_id<SamAngle->nentry;angle_id++) {
     // Reset sel2 xyz for sampling of other angles (angle_id!=0)
     if (angle_id > 0) {
        gen_sel2_geo(sel1->natoms,sel2->natoms,atoms,sel2pdb,mut_data);
     }
     for (k=0;k<SamAngle->nsam[angle_id];k++) {
        rot_sel2_angle(sel2->natoms,angle_id,SamAngle->samangle[k],&atoms[4*sel1->natoms],SamAngle);
 
        Eb[nrotate] = calc_Eb(&atoms[4*sel1->natoms],Dihed,Angle);
        if (Eb[nrotate]-Ebmin < 0) {Ebmin=Eb[nrotate];}  // For numerical stability
 
        // For valid rotation: Update the patch_data for sel2
        if (patch_data->flag_spacedecomp) { binatoms(cell,atoms,vdw_data,patch_data,sel2);}
 
        // For valid rotation: Compute energy of sel2 + environment (sel0)
        if (sel0->natoms==0) {
           compute_sel_nonb(data,cell,natoms,atoms,vdw_data,exc_data,sel2,&Evdw,&Eelec);
           E2_PME = Eelec;
           E2_VDW = Evdw;
        } else {
           compute_sel_vdwpmer(data,cell,natoms,atoms,vdw_data,exc_data,sel2,&E_vdw_sel,&E_pme_sel);
           if (patch_data->flag_spacedecomp) compute_patch_vdwpmer(cell,data->ewald_factor,vdw_data,patch_data,sel2,sel0,&E_vdw_pair,&E_pme_pair);
           else compute_intersel_vdwpmer(data,cell,natoms,atoms,box_disp,vdw_data,exc_data,sel2,sel0,&E_vdw_pair,&E_pme_pair);
           E_work = compute_pmek(data,cell,natoms,atoms,charge,sel1);
           E2_PME = E_pme_sel+E_pme_pair+E_work;
           E2_VDW = E_vdw_sel+E_vdw_pair;
        }

#ifdef Ham
        // Ghost ligand correction (Using a full molecule for ghost ligand)
        compute_sel_nonb(data,cell,natoms,atoms,vdw_data,exc_data,sel2,&Evdw,&Eelec);
        Eb[nrotate] += Evdw+Eelec;
        E2_PME -= Eelec;
        E2_VDW -= Evdw;
#endif 

        // Assign DeltaE (rotated structure) and update nrotate
        DeltaE[nrotate] = (E2_PME + E2_VDW)-(E1_PME + E1_VDW)+SelfE_diff;
        nrotate++;
     }
  }
 
  // Return DeltaE, Eb, Ebmin via gen_data
  gen_data->Ebmin=Ebmin;
}


void vsscase4(const double *cell, double *atoms, struct GEN *gen_data, pmepot_data *data, struct VDW *vdw_data, 
           struct PSE *sel0, struct PSE *sel1, struct PSE *sel2, struct EXC *exc_data, struct MUT *mut_data, 
           struct PAT *patch_data, struct BOND *Dihed, struct BOND *Angle, struct SAM *SamDihed, struct SAM *SamAngle) 
{
  // Evaluate DeltaE for each frame, case 4, by Letitia: 06/08/2014, 12/11/2015, 23/01/2016
  //int i, j, k, myindex, 
  int nrotate, natoms;
  //int dihed_id, angle_id;
  double E1_PME, E1_VDW, E2_PME, E2_VDW, SelfE_diff;
  double Evdw, Eelec;
  double *charge, *sel2pdb, *Eb, *DeltaE;
  double Ebmin=0.0;
  double E_vdw_sel, E_pme_sel, E_vdw_pair, E_pme_pair, E_work, box_disp[81];

  charge=gen_data->charge;
  sel2pdb=gen_data->sel2pdb;
  Eb=gen_data->Eb;
  DeltaE=gen_data->DeltaE;
  natoms=gen_data->natoms;

  if (patch_data->flag_spacedecomp) {
     box_update(cell,patch_data);
     binatoms(cell,atoms,vdw_data,patch_data,sel0);
     binatoms(cell,atoms,vdw_data,patch_data,sel1);
  } else {
     get_box_disp(cell,box_disp);
  }

  // Compute energy of sel1 + environment (sel0)
  if (sel0->natoms==0) {
     compute_sel_nonb(data,cell,natoms,atoms,vdw_data,exc_data,sel1,&Evdw,&Eelec);
     E1_PME = Eelec;
     E1_VDW = Evdw;
     SelfE_diff = 0.0;
  } else {   
     compute_sel_vdwpmer(data,cell,natoms,atoms,vdw_data,exc_data,sel1,&E_vdw_sel,&E_pme_sel);
     if (patch_data->flag_spacedecomp) compute_patch_vdwpmer(cell,data->ewald_factor,vdw_data,patch_data,sel1,sel0,&E_vdw_pair,&E_pme_pair);
     else compute_intersel_vdwpmer(data,cell,natoms,atoms,box_disp,vdw_data,exc_data,sel1,sel0,&E_vdw_pair,&E_pme_pair);
     E_work = compute_pmek(data,cell,natoms,atoms,charge,sel2);
     E1_PME = E_pme_sel+E_pme_pair+E_work;
     E1_VDW = E_vdw_sel+E_vdw_pair;
     // Compute the PME selfenergy difference between sel2 and sel1. (Constant for all geometry.)
     SelfE_diff = compute_pmeselfdiff(data,atoms,sel1,sel2);
  }

#ifdef Ham
  // Ghost ligand correction (Using a full molecule for ghost ligand)
  compute_sel_nonb(data,cell,natoms,atoms,vdw_data,exc_data,sel1,&Evdw,&Eelec);
  E1_PME -= Eelec;
  E1_VDW -= Evdw;
#endif

  nrotate=0;
  // SamDihed->nentry=0, SamAngle->nentry=0
  // generate sel2 xyz (No Dihedral nor Angle sampling)
  gen_sel2_geo(sel1->natoms,sel2->natoms,atoms,sel2pdb,mut_data);

  Eb[nrotate] = calc_Eb(&atoms[4*sel1->natoms],Dihed,Angle);
  if (Eb[nrotate]-Ebmin < 0) {Ebmin=Eb[nrotate];}  // For numerical stability
 
  // For valid rotation: Update the patch_data for sel2
  if (patch_data->flag_spacedecomp) {
     binatoms(cell,atoms,vdw_data,patch_data,sel2);
  }
 
  // For valid rotation: Compute energy of sel2 + environment (sel0)
  if (sel0->natoms==0) {
     compute_sel_nonb(data,cell,natoms,atoms,vdw_data,exc_data,sel2,&Evdw,&Eelec);
     E2_PME = Eelec;
     E2_VDW = Evdw;
  } else {  
     compute_sel_vdwpmer(data,cell,natoms,atoms,vdw_data,exc_data,sel2,&E_vdw_sel,&E_pme_sel);
     if (patch_data->flag_spacedecomp) compute_patch_vdwpmer(cell,data->ewald_factor,vdw_data,patch_data,sel2,sel0,&E_vdw_pair,&E_pme_pair);
     else compute_intersel_vdwpmer(data,cell,natoms,atoms,box_disp,vdw_data,exc_data,sel2,sel0,&E_vdw_pair,&E_pme_pair);
     E_work = compute_pmek(data,cell,natoms,atoms,charge,sel1);
     E2_PME = E_pme_sel+E_pme_pair+E_work;
     E2_VDW = E_vdw_sel+E_vdw_pair;
  }

#ifdef Ham
  // Ghost ligand correction (Using a full molecule for ghost ligand)
  compute_sel_nonb(data,cell,natoms,atoms,vdw_data,exc_data,sel2,&Evdw,&Eelec);
  Eb[nrotate] += Evdw+Eelec;
  E2_PME -= Eelec;
  E2_VDW -= Evdw;
#endif

  // Assign DeltaE (rotated structure) and update nrotate
  DeltaE[nrotate] = (E2_PME + E2_VDW)-(E1_PME + E1_VDW)+SelfE_diff;
  nrotate++;

  // Return DeltaE, Eb, Ebmin via gen_data
  gen_data->Ebmin=Ebmin;
}


