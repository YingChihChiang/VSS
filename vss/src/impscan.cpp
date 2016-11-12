#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "vss.h"

#define M_PI 3.14159265358979323846

void scandihed(int dihed_id, double RT, double *atoms, struct PSE *sel1, struct PSE *sel2, struct MUT *mut_data, struct BOND *Dihed, struct BOND *Angle, struct SAM *SamDihed) 
{
  // Perform a dihedral scan (every 1 degree) to determine the relevant dihedral angles for sel2.
  // by Letitia: 18/11/2015
  int i, j, disp, nsam, nfold=360, Ebindex;
  int *pickid;
  double samangle, *Eb, Ebmin, *sel1pdb, *sel2pdb;

  nsam=SamDihed->nsam[dihed_id];
  Eb= (double*) malloc(nfold*sizeof(double));
  pickid= (int*) malloc(nfold*sizeof(int));
  for (i=0;i<nfold;i++) {pickid[i]=i;}
  sel1pdb= (double*) malloc(3*sel1->natoms*sizeof(double));
  sel2pdb= (double*) malloc(3*sel2->natoms*sizeof(double));

  // create sel1pdb as a backup
  for (i=0;i<sel1->natoms;i++) {
     sel1pdb[3*i]=atoms[4*i];
     sel1pdb[3*i+1]=atoms[4*i+1];
     sel1pdb[3*i+2]=atoms[4*i+2];
  }

  // create sel2pdb
  for (i=0;i<sel2->natoms;i++) {
     disp=sel1->natoms+i;
     sel2pdb[3*i]=atoms[4*disp];
     sel2pdb[3*i+1]=atoms[4*disp+1];
     sel2pdb[3*i+2]=atoms[4*disp+2];
  }

  // replace sel1 xyz by sel2 xyz
  for (i=0;i<sel1->natoms;i++) {
     disp=sel1->natoms+i;
     atoms[4*i]=atoms[4*disp];
     atoms[4*i+1]=atoms[4*disp+1];
     atoms[4*i+2]=atoms[4*disp+2];
  }

  // Rotate the geometry. Construct the inner layer sampling for sel2.
  for (i=0;i<nfold;i++) {
     samangle = i*2*M_PI/nfold;
     gen_sel2_geo(sel1->natoms,sel2->natoms,atoms,sel2pdb,mut_data);
     rot_sel2_dihed(sel2->natoms,dihed_id,samangle,&atoms[4*sel1->natoms],SamDihed);
     Eb[i] = calc_Eb(&atoms[4*sel1->natoms],Dihed,Angle);
  }
 
  // pick the smallest energy and report
  for (i=0;i<nsam;i++) {
     Ebmin=Eb[i];
     Ebindex=i;
     for (j=i+1;j<nfold;j++) {
        if (Eb[j]<Ebmin) {
           Ebmin=Eb[j];
           Ebindex=j;
        }
     }
     Eb[Ebindex]=Eb[i];
     Eb[i]=Ebmin;
     disp=pickid[Ebindex];
     pickid[Ebindex]=pickid[i];
     pickid[i]=disp;;
  }
 
  // initialize the sampling angles for dihed_id in SamDihed.
  disp=sumint(SamDihed->nsam,dihed_id-1);
  printf("Important sampling for Dihedral %d, used angle: ",dihed_id);
  for (i=0;i<nsam;i++) {
     SamDihed->samangle[disp+i] = pickid[i]*2*M_PI/nfold;
     printf("%.2f  ",SamDihed->samangle[disp+i]);
  }
  printf("\n");

  // change the sel1 xyz back
  for (i=0;i<sel1->natoms;i++) {
     atoms[4*i]=sel1pdb[3*i];
     atoms[4*i+1]=sel1pdb[3*i+1];
     atoms[4*i+2]=sel1pdb[3*i+2];
  }

  // change the sel2 xyz back
  for (i=0;i<sel2->natoms;i++) {
     disp=sel1->natoms+i;
     atoms[4*disp]=sel2pdb[3*i];
     atoms[4*disp+1]=sel2pdb[3*i+1];
     atoms[4*disp+2]=sel2pdb[3*i+2];
  }


  free(Eb);
  free(pickid);
  free(sel1pdb);
  free(sel2pdb);
  return;

}


void scanangle(int angle_id, double RT, double *atoms, struct PSE *sel1, struct PSE *sel2, struct MUT *mut_data, struct BOND *Dihed, struct BOND *Angle, struct SAM *SamAngle) 
{
  // Perform an angle scan (every 1 degree) to determine the relevant dihedral angles for sel2.
  // by Letitia: 07/12/2015
  int i, j, disp, nsam, nfold=180, Ebindex;
  int *pickid;
  double samangle, *Eb, Ebmin, *sel1pdb, *sel2pdb;

  nsam=SamAngle->nsam[angle_id];
  Eb= (double*) malloc(nfold*sizeof(double));
  pickid= (int*) malloc(nfold*sizeof(int));
  for (i=0;i<nfold;i++) {pickid[i]=i;}
  sel1pdb= (double*) malloc(3*sel1->natoms*sizeof(double));
  sel2pdb= (double*) malloc(3*sel2->natoms*sizeof(double));

  // create sel1pdb as a backup
  for (i=0;i<sel1->natoms;i++) {
     sel1pdb[3*i]=atoms[4*i];
     sel1pdb[3*i+1]=atoms[4*i+1];
     sel1pdb[3*i+2]=atoms[4*i+2];
  }

  // create sel2pdb
  for (i=0;i<sel2->natoms;i++) {
     disp=sel1->natoms+i;
     sel2pdb[3*i]=atoms[4*disp];
     sel2pdb[3*i+1]=atoms[4*disp+1];
     sel2pdb[3*i+2]=atoms[4*disp+2];
  }

  // replace sel1 xyz by sel2 xyz
  for (i=0;i<sel1->natoms;i++) {
     disp=sel1->natoms+i;
     atoms[4*i]=atoms[4*disp];
     atoms[4*i+1]=atoms[4*disp+1];
     atoms[4*i+2]=atoms[4*disp+2];
  }

  // Rotate the geometry. Construct the inner layer sampling for sel2.
  for (i=0;i<nfold;i++) {
     samangle = i*M_PI/nfold;
     gen_sel2_geo(sel1->natoms,sel2->natoms,atoms,sel2pdb,mut_data);
     rot_sel2_angle(sel2->natoms,angle_id,samangle,&atoms[4*sel1->natoms],SamAngle);
     Eb[i] = calc_Eb(&atoms[4*sel1->natoms],Dihed,Angle);
  }
 
  // pick the smallest energy and report
  for (i=0;i<nsam;i++) {
     Ebmin=Eb[i];
     Ebindex=i;
     for (j=i+1;j<nfold;j++) {
        if (Eb[j]<Ebmin) {
           Ebmin=Eb[j];
           Ebindex=j;
        }
     }
     Eb[Ebindex]=Eb[i];
     Eb[i]=Ebmin;
     disp=pickid[Ebindex];
     pickid[Ebindex]=pickid[i];
     pickid[i]=disp;;
  }
 
  // initialize the sampling angles for angle_id in SamAngle.
  disp=sumint(SamAngle->nsam,angle_id-1);
  printf("Important sampling for Angle %d, used angle: ",angle_id);
  for (i=0;i<nsam;i++) {
     SamAngle->samangle[disp+i] = pickid[i]*M_PI/nfold;
     printf("%.2f  ",SamAngle->samangle[disp+i]);
  }
  printf("\n");
  

  // change the sel1 xyz back
  for (i=0;i<sel1->natoms;i++) {
     atoms[4*i]=sel1pdb[3*i];
     atoms[4*i+1]=sel1pdb[3*i+1];
     atoms[4*i+2]=sel1pdb[3*i+2];
  }

  // change the sel2 xyz back
  for (i=0;i<sel2->natoms;i++) {
     disp=sel1->natoms+i;
     atoms[4*disp]=sel2pdb[3*i];
     atoms[4*disp+1]=sel2pdb[3*i+1];
     atoms[4*disp+2]=sel2pdb[3*i+2];
  }


  free(Eb);
  free(pickid);
  free(sel1pdb);
  free(sel2pdb);
  return;

}


