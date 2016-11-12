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
#include <omp.h>
#include "vss.h"

#define M_PI 3.14159265358979323846
#define SQRT_PI 1.7724538509055160273  
#define COLOUMB 332.0636
#define BOLTZMAN 0.001987191

#define BUFSIZE 256

int vssmain(struct GEN *gen_data, struct DCDS *dcd_data, struct VDW *vdw_data, struct PSE *sel0, struct PSE *sel1, struct PSE *sel2, 
            struct EXC *exc_data, struct MUT *mut_data, struct BOND *Dihed, struct BOND *Angle, struct SAM *SamDihed, struct SAM *SamAngle) 
{
  // VSS by Letitia: 06/08/2014, 12/11/2015, 22/01/201, 22/01/2016

  // Reassign the general variables
  int natoms = gen_data->natoms;
  double RT = gen_data->RT;
  double ewaldfactor = gen_data->ewaldfactor; 
  double gridspacing = gen_data->gridspacing;

  // Reassign the dcd variables
  int *index1, *index3, fint, fstep, fend, tmax;
  char *dcdfile, *xstfile;
  double *dcd_xyz, *xyz1, *xyz2, *xyz3;
  void *dcd;
  FILE *xst;
  fint=dcd_data->fint; fstep=dcd_data->fstep; fend=dcd_data->fend; 
  dcdfile=dcd_data->dcdfile; xstfile=dcd_data->xstfile;
  index1=dcd_data->index1; index3=dcd_data->index3;
  xyz1 = dcd_data->xyz1; xyz2 = dcd_data->xyz2; xyz3 = dcd_data->xyz3;

  // Local variables
  int rc;
  int *txst, nframe=1;
  int *grid;
  double *cell, *atoms;
  double acc_exp_dE=0.0, dG[2];
  char namebuf[BUFSIZE];
  FILE *fd, *fd2, *fd3;
  pmepot_data **data;
  struct GEN *mygen_data;  
  struct PSE *mysel0, *mysel1, *mysel2;
  struct PAT *mypatch_data;

  // OpenMP initialization
  int nthreads=omp_get_num_procs();
  printf("nthreads = %d \n",nthreads);
  txst = (int *) malloc(nthreads*sizeof(int));
  grid = (int *) malloc(nthreads*3*sizeof(int));
  cell = (double *) malloc(nthreads*12*sizeof(double));
  atoms = (double *) malloc(nthreads*natoms*4*sizeof(double));
  data = (pmepot_data **) malloc(nthreads*sizeof(pmepot_data *));
  mygen_data = (struct GEN *) malloc(nthreads*sizeof(struct GEN));
  mysel0 = (struct PSE *) malloc(nthreads*sizeof(struct PSE));
  mysel1 = (struct PSE *) malloc(nthreads*sizeof(struct PSE));
  mysel2 = (struct PSE *) malloc(nthreads*sizeof(struct PSE));
  mypatch_data = (struct PAT *) malloc(nthreads*sizeof(struct PAT));
  for (int j=0;j<nthreads;j++) {
     mygen_data[j].charge = (double *) malloc(natoms*sizeof(double));
  }

  // Prepare for reading dcd and xst
  dcd = initDcdhandle(dcdfile, fint);
  if(!dcd){ printf(" \n"); printf("Cannot open dcdfile! \n"); printf(" \n"); return 1; }
  xst = fopen(xstfile, "r");
  if(!xst){ printf(" \n"); printf("Cannot open xstfile! \n"); printf(" \n"); return 1; }
  for(int i=0; i<3+fint; i++) fgets(namebuf, BUFSIZE, xst);  // Letitia: xst output frame 0 is NOT included in dcd file.
  rc = sscanf(namebuf, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", txst, cell+3, cell+4, cell+5, cell+6, cell+7, cell+8,cell+9, cell+10, cell+11, cell  , cell+1, cell+2);
  if(rc!=13){ printf(" \n"); printf("Incorrect format in xstfile! \n"); printf(" \n"); return 1; }

  // Initialize atoms, pmepot data, VDW data, Patch data (require cell and vdw_data_init), and general_data
  vdw_data_init(ewaldfactor,vdw_data);
  general_init(natoms,xyz2,atoms,gen_data,SamDihed,SamAngle); // evaluate gen_data->jobcase
  make_grid_from_cell(grid,cell,gridspacing);
  for (int j=0;j<nthreads;j++) {
     // mygen_data 
     for (int i=0;i<natoms;i++) { mygen_data[j].charge[i]=gen_data->charge[i];}
     general_init(natoms,xyz2,atoms,&mygen_data[j],SamDihed,SamAngle);
     mygen_data[j].natoms = natoms; 
     mygen_data[j].RT = RT; 
     mygen_data[j].ewaldfactor = ewaldfactor; 
     mygen_data[j].gridspacing = gridspacing; 
     // pmepot_data
     data[j] = pmepot_create(grid,ewaldfactor);
     // mysel0,1,2
     mysel0[j].natoms = sel0->natoms;
     mysel1[j].natoms = sel1->natoms;
     mysel2[j].natoms = sel2->natoms;
     mysel0[j].aid = sel0->aid;  // Not allocated. Don't free.   
     mysel1[j].aid = sel1->aid;  // Not allocated. Don't free. 
     mysel2[j].aid = sel2->aid;  // Not allocated. Don't free.
     patch_init(natoms,vdw_data->Roff,cell,&mysel0[j],&mysel1[j],&mysel2[j],&mypatch_data[j]); 
     // atoms(charge)
     for (int i=0;i<natoms;i++) { atoms[j*4*natoms+4*i+3]=gen_data->charge[i];}
  }

  // Preparing for processing time frames
  dcd_xyz = (double*)malloc(sizeof(double) * getNatoms(dcd) * 3);
  tmax = getNsets(dcd);
  if (fend == 0 || fend > tmax) {fend = tmax;}

  fd=fopen("vss.log","w");
  fd2=fopen("DeltaE.log","w");
  fd3=fopen("Eb.log","w");
  fprintf(fd,"#    Findex      Frame           dG    |  DeltaE \n"); 

#ifdef Ham
  printf("Ham. II is used! (ghost ligand with both bonded and nonbonded interactions) \n");
#endif

  // Determine which sampling is on. Run vss according to its case.
  if (gen_data->jobcase==1) {
     // Case 1: Both Dihedral and Angle Sampling
     printf("Case 1. \n");
     // Run VSS for each time frame
     for (int t=fint; t<fend; t+=fstep*nthreads){

        for (int j=0; j<nthreads; j++) {
           if ((t+j*fstep) >= fend) continue;
           double *mycell = &cell[j*12];
           double *myatoms = &atoms[j*natoms*4];
           // cell[] preparation
           fgets(namebuf, BUFSIZE, xst);
           rc = sscanf(namebuf, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", txst+j, mycell+3, mycell+4, mycell+5, mycell+6, mycell+7, mycell+8, mycell+9, mycell+10, mycell+11, mycell  , mycell+1, mycell+2);
           for (int i=0; i<fstep-1; i++) fgets(namebuf, BUFSIZE, xst);  // Letitia: dcd can jump directly to a given index; xst cannot. Read through xst.
           // atoms[] preparation
           getNextXYZ(dcd_xyz, dcd, fstep);
           cutXYZ(xyz1, dcd_xyz, index1, sel1->natoms);
           updateXYZofAtom(myatoms                              , xyz1, sel1->natoms);
           updateXYZofAtom(myatoms+4* sel1->natoms              , xyz2, sel2->natoms);
           if (sel0->natoms!=0) {
              cutXYZ(xyz3, dcd_xyz, index3, sel0->natoms);
              updateXYZofAtom(myatoms+4*(sel1->natoms+sel2->natoms), xyz3, sel0->natoms);
           }
        }

        // Dihedral & Angle scan to determine which angles should be used in sampling
        if (t==fint) {
           if (SamDihed->scanflag) for (int i=0;i<SamDihed->nentry;i++) { scandihed(i,RT,atoms,sel1,sel2,mut_data,Dihed,Angle,SamDihed);}
           if (SamAngle->scanflag) for (int i=0;i<SamAngle->nentry;i++) { scanangle(i,RT,atoms,sel1,sel2,mut_data,Dihed,Angle,SamAngle);}
        }

        # pragma omp parallel for num_threads(nthreads) schedule(static,1)
        for (int j=0; j<nthreads; j++) {
           if ((t+j*fstep) >= fend) continue;
           double *mycell = &cell[j*12];
           double *myatoms = &atoms[j*natoms*4];
           vsscase1(mycell,myatoms,&mygen_data[j],data[j],vdw_data,&mysel0[j],&mysel1[j],&mysel2[j],exc_data,mut_data,&mypatch_data[j],Dihed,Angle,SamDihed,SamAngle);
        }

        for (int j=0; j<nthreads; j++) {
           if ((t+j*fstep) >= fend) continue;
           for (int i=0;i<gen_data->nrotate;i++) {
              fprintf(fd2,"%15.7E ",mygen_data[j].DeltaE[i]);
              fprintf(fd3,"%15.7E ",mygen_data[j].Eb[i]);
           }
           fprintf(fd2,"\n");
           fprintf(fd3,"\n");
           calcdG(nframe,&acc_exp_dE,dG,&mygen_data[j]);
           fprintf(fd,"%10d %10d %15.7E %15.7E\n",t+j*fstep,txst[j],dG[0],dG[1]);
           //fprintf(fd,"%10d %15.7E %15.7E\n",txst[j],dG[0],dG[1]);
           nframe++;
        }
     }
  }

  if (gen_data->jobcase==2) {
     // Case 2: Only Dihedral Sampling
     printf("Case 2. \n");
     // Run VSS for each time frame
     for (int t=fint; t<fend; t+=fstep*nthreads){

        for (int j=0; j<nthreads; j++) {
           if ((t+j*fstep) >= fend) continue;
           double *mycell = &cell[j*12];
           double *myatoms = &atoms[j*natoms*4];
           // cell[] preparation
           fgets(namebuf, BUFSIZE, xst);
           rc = sscanf(namebuf, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", txst+j, mycell+3, mycell+4, mycell+5, mycell+6, mycell+7, mycell+8, mycell+9, mycell+10, mycell+11, mycell  , mycell+1, mycell+2);
           for (int i=0; i<fstep-1; i++) fgets(namebuf, BUFSIZE, xst);  // Letitia: dcd can jump directly to a given index; xst cannot. Read through xst.
           // atoms[] preparation
           getNextXYZ(dcd_xyz, dcd, fstep);
           cutXYZ(xyz1, dcd_xyz, index1, sel1->natoms);
           updateXYZofAtom(myatoms                              , xyz1, sel1->natoms);
           updateXYZofAtom(myatoms+4* sel1->natoms              , xyz2, sel2->natoms);
           if (sel0->natoms!=0) {
              cutXYZ(xyz3, dcd_xyz, index3, sel0->natoms);
              updateXYZofAtom(myatoms+4*(sel1->natoms+sel2->natoms), xyz3, sel0->natoms);
           }
        }

        // Dihedral & Angle scan to determine which angles should be used in sampling
        if (t==fint) {
           if (SamDihed->scanflag) for (int i=0;i<SamDihed->nentry;i++) { scandihed(i,RT,atoms,sel1,sel2,mut_data,Dihed,Angle,SamDihed);}
           if (SamAngle->scanflag) for (int i=0;i<SamAngle->nentry;i++) { scanangle(i,RT,atoms,sel1,sel2,mut_data,Dihed,Angle,SamAngle);}
        }

        # pragma omp parallel for num_threads(nthreads) schedule(static,1)
        for (int j=0; j<nthreads; j++) {
           if ((t+j*fstep) >= fend) continue;
           double *mycell = &cell[j*12];
           double *myatoms = &atoms[j*natoms*4];
           vsscase2(mycell,myatoms,&mygen_data[j],data[j],vdw_data,&mysel0[j],&mysel1[j],&mysel2[j],exc_data,mut_data,&mypatch_data[j],Dihed,Angle,SamDihed,SamAngle);
        }

        for (int j=0; j<nthreads; j++) { 
           if ((t+j*fstep) >= fend) continue;
           for (int i=0;i<gen_data->nrotate;i++) {
              fprintf(fd2,"%15.7E ",mygen_data[j].DeltaE[i]);
              fprintf(fd3,"%15.7E ",mygen_data[j].Eb[i]);
           }           
           fprintf(fd2,"\n");
           fprintf(fd3,"\n");
           calcdG(nframe,&acc_exp_dE,dG,&mygen_data[j]);
           fprintf(fd,"%10d %10d %15.7E %15.7E\n",t+j*fstep,txst[j],dG[0],dG[1]);
           //fprintf(fd,"%10d %15.7E %15.7E\n",txst[j],dG[0],dG[1]);
           nframe++;
        }
     }
  }

  if (gen_data->jobcase==3) {
     // Case 3: Only Angle Sampling
     printf("Case 3. \n");
     // Run VSS for each time frame
     for (int t=fint; t<fend; t+=fstep*nthreads){

        for (int j=0; j<nthreads; j++) {
           if ((t+j*fstep) >= fend) continue;
           double *mycell = &cell[j*12];
           double *myatoms = &atoms[j*natoms*4];
           // cell[] preparation
           fgets(namebuf, BUFSIZE, xst);
           rc = sscanf(namebuf, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", txst+j, mycell+3, mycell+4, mycell+5, mycell+6, mycell+7, mycell+8, mycell+9, mycell+10, mycell+11, mycell  , mycell+1, mycell+2);
           for (int i=0; i<fstep-1; i++) fgets(namebuf, BUFSIZE, xst);  // Letitia: dcd can jump directly to a given index; xst cannot. Read through xst.
           // atoms[] preparation
           getNextXYZ(dcd_xyz, dcd, fstep);
           cutXYZ(xyz1, dcd_xyz, index1, sel1->natoms);
           updateXYZofAtom(myatoms                              , xyz1, sel1->natoms);
           updateXYZofAtom(myatoms+4* sel1->natoms              , xyz2, sel2->natoms);
           if (sel0->natoms!=0) {
              cutXYZ(xyz3, dcd_xyz, index3, sel0->natoms);
              updateXYZofAtom(myatoms+4*(sel1->natoms+sel2->natoms), xyz3, sel0->natoms);
           }
        }

        // Dihedral & Angle scan to determine which angles should be used in sampling
        if (t==fint) {
           if (SamDihed->scanflag) for (int i=0;i<SamDihed->nentry;i++) { scandihed(i,RT,atoms,sel1,sel2,mut_data,Dihed,Angle,SamDihed);}
           if (SamAngle->scanflag) for (int i=0;i<SamAngle->nentry;i++) { scanangle(i,RT,atoms,sel1,sel2,mut_data,Dihed,Angle,SamAngle);}
        }

        # pragma omp parallel for num_threads(nthreads) schedule(static,1)
        for (int j=0; j<nthreads; j++) {
           if ((t+j*fstep) >= fend) continue;
           double *mycell = &cell[j*12];
           double *myatoms = &atoms[j*natoms*4];
           vsscase3(mycell,myatoms,&mygen_data[j],data[j],vdw_data,&mysel0[j],&mysel1[j],&mysel2[j],exc_data,mut_data,&mypatch_data[j],Dihed,Angle,SamDihed,SamAngle);
        }

        for (int j=0; j<nthreads; j++) {
           if ((t+j*fstep) >= fend) continue;
           for (int i=0;i<gen_data->nrotate;i++) {
              fprintf(fd2,"%15.7E ",mygen_data[j].DeltaE[i]);
              fprintf(fd3,"%15.7E ",mygen_data[j].Eb[i]);
           }
           fprintf(fd2,"\n");
           fprintf(fd3,"\n");
           calcdG(nframe,&acc_exp_dE,dG,&mygen_data[j]);
           fprintf(fd,"%10d %10d %15.7E %15.7E\n",t+j*fstep,txst[j],dG[0],dG[1]);
           //fprintf(fd,"%10d %15.7E %15.7E\n",txst[j],dG[0],dG[1]);
           nframe++;
        }
     }
  }

  if (gen_data->jobcase==4) {
     // Case 4: No improved Sampling
     printf("Case 4. \n");
     // Run VSS for each time frame
     for (int t=fint; t<fend; t+=fstep*nthreads){

        for (int j=0; j<nthreads; j++) {
           if ((t+j*fstep) >= fend) continue;
           double *mycell = &cell[j*12];
           double *myatoms = &atoms[j*natoms*4];
           // cell[] preparation
           fgets(namebuf, BUFSIZE, xst);
           rc = sscanf(namebuf, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", txst+j, mycell+3, mycell+4, mycell+5, mycell+6, mycell+7, mycell+8, mycell+9, mycell+10, mycell+11, mycell  , mycell+1, mycell+2);
           for (int i=0; i<fstep-1; i++) fgets(namebuf, BUFSIZE, xst);  // Letitia: dcd can jump directly to a given index; xst cannot. Read through xst.
           // atoms[] preparation
           getNextXYZ(dcd_xyz, dcd, fstep);
           cutXYZ(xyz1, dcd_xyz, index1, sel1->natoms);
           updateXYZofAtom(myatoms                              , xyz1, sel1->natoms);
           updateXYZofAtom(myatoms+4* sel1->natoms              , xyz2, sel2->natoms);
           if (sel0->natoms!=0) {
              cutXYZ(xyz3, dcd_xyz, index3, sel0->natoms);
              updateXYZofAtom(myatoms+4*(sel1->natoms+sel2->natoms), xyz3, sel0->natoms);
           }
        }

        // Dihedral & Angle scan to determine which angles should be used in sampling
        if (t==fint) {
           if (SamDihed->scanflag) for (int i=0;i<SamDihed->nentry;i++) { scandihed(i,RT,atoms,sel1,sel2,mut_data,Dihed,Angle,SamDihed);}
           if (SamAngle->scanflag) for (int i=0;i<SamAngle->nentry;i++) { scanangle(i,RT,atoms,sel1,sel2,mut_data,Dihed,Angle,SamAngle);}
        }

        # pragma omp parallel for num_threads(nthreads) schedule(static,1)
        for (int j=0; j<nthreads; j++) {
           if ((t+j*fstep) >= fend) continue;
           double *mycell = &cell[j*12];
           double *myatoms = &atoms[j*natoms*4];
           vsscase4(mycell,myatoms,&mygen_data[j],data[j],vdw_data,&mysel0[j],&mysel1[j],&mysel2[j],exc_data,mut_data,&mypatch_data[j],Dihed,Angle,SamDihed,SamAngle);
        }

        for (int j=0; j<nthreads; j++) {
           if ((t+j*fstep) >= fend) continue;
           for (int i=0;i<gen_data->nrotate;i++) {
              fprintf(fd2,"%15.7E ",mygen_data[j].DeltaE[i]);
              fprintf(fd3,"%15.7E ",mygen_data[j].Eb[i]);
           }
           fprintf(fd2,"\n");
           fprintf(fd3,"\n");
           calcdG(nframe,&acc_exp_dE,dG,&mygen_data[j]);
           fprintf(fd,"%10d %10d %15.7E %15.7E\n",t+j*fstep,txst[j],dG[0],dG[1]);
           //fprintf(fd,"%10d %15.7E %15.7E\n",txst[j],dG[0],dG[1]);
           nframe++;
        }
     }
  }

  closeDcdhandle(dcd);
  fclose(xst);
  fclose(fd);
  fclose(fd2);
  fclose(fd3);

  free(atoms);
  free(dcd_xyz);
  for (int j=0;j<nthreads;j++) {
     pmepot_destroy(data[j]);
     general_free(&mygen_data[j]);
     if (mypatch_data[j].flag_spacedecomp) patch_free(&mypatch_data[j],&mysel0[j],&mysel1[j],&mysel2[j]);
  }

  return 0;

}


void calcdG (int nframe, double *acc_exp_dE, double *dG, struct GEN *gen_data) {
// Evaluate weights and weighted dG for all structures in this step
  int i;  
  double norm=0.0;

  for (i=0;i<gen_data->nrotate;i++) {gen_data->Eb[i] = exp(-(gen_data->Eb[i]-gen_data->Ebmin)/gen_data->RT);}  // For numerical stability
  for (i=0;i<gen_data->nrotate;i++) {norm += gen_data->Eb[i];}
  for (i=0;i<gen_data->nrotate;i++) {gen_data->Eb[i] = gen_data->Eb[i]/norm;}
        
  dG[0] = 0.0;
  for (i=0;i<gen_data->nrotate;i++) {dG[0] += exp(-gen_data->DeltaE[i]/gen_data->RT)*gen_data->Eb[i];}
  dG[1] = -gen_data->RT*log(dG[0]);
  if (isinf(dG[1])) {dG[1]=1000000;}
  
  *acc_exp_dE += dG[0];
  dG[0] = -gen_data->RT*log(*acc_exp_dE/nframe);

}
