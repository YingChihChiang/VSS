#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "vss.h"

#if !defined(M_PI)
#define M_PI 3.14159265358979323846
#endif

void gen_sel2_geo(int nsel1, int nsel2, double *xyzq, double *sel2pdb, struct MUT *mut_data) {
// Generate sel2 geometry by placing functional groups (mutated atoms) on sel1.
// Letitia 2015.11.12
  int i, j, indexa1, indexb1, indexa2, indexb2, myindex, nsel;
  double theta1, theta2, phi1, phi2, length;
  double sintheta1, costheta1, sintheta2, costheta2, sinphi1, cosphi1, sinphi2, cosphi2;
  double vec1[3], vec2[3];

  // Assign sel1 xyz as default structure (Note that sel1 and sel2 should have a common base structure!)
  if (nsel1 > nsel2) {nsel=nsel2;} else {nsel=nsel1;}
  for (i=0;i<nsel;i++) {
     myindex = i+nsel1;
     xyzq[4*myindex]=xyzq[4*i];
     xyzq[4*myindex+1]=xyzq[4*i+1];
     xyzq[4*myindex+2]=xyzq[4*i+2];
  }

  // Assign mutated atoms' xyz
  myindex=0;
  for (i=0;i<mut_data->vecnum;i++) {
     // get vec1 vec2
     indexa1 = mut_data->veclist[2*i];
     indexb1 = mut_data->veclist[2*i+1];
     indexa2 = mut_data->veclist[2*mut_data->vecnum+2*i]-nsel1;
     indexb2 = mut_data->veclist[2*mut_data->vecnum+2*i+1]-nsel1;
     for (j=0;j<3;j++) {
        vec1[j] = xyzq[4*indexb1+j]-xyzq[4*indexa1+j];
        vec2[j] = sel2pdb[3*indexb2+j]-sel2pdb[3*indexa2+j];
     }

     // get theta and phi for vec1 and vec2
     length = sqrt(vec2[0]*vec2[0]+vec2[1]*vec2[1]);
     phi2   = atan2(vec2[1],vec2[0]);
     theta2 = atan(length/vec2[2]);
     length = sqrt(vec1[0]*vec1[0]+vec1[1]*vec1[1]);
     phi1   = atan2(vec1[1],vec1[0]);
     theta1 = atan(length/vec1[2]);
     if (theta1 < 0) {theta1 += M_PI;} // range {-PI/2, PI/2} --> {0, PI}
     if (theta2 < 0) {theta2 += M_PI;} // range {-PI/2, PI/2} --> {0, PI}
     sintheta1 = sin(theta1);
     costheta1 = cos(theta1);
     sinphi1 = sin(phi1);
     cosphi1 = cos(phi1);
     sintheta2 = sin(theta2);
     costheta2 = cos(theta2);
     sinphi2 = sin(phi2);
     cosphi2 = cos(phi2);

     // rotate vectors belong to vec2-group
     for (j=0;j<mut_data->atomnum[i];j++) {

        indexb2=mut_data->atomlist[myindex]-nsel1;

        vec2[0]=sel2pdb[3*indexb2+0]-sel2pdb[3*indexa2+0];
        vec2[1]=sel2pdb[3*indexb2+1]-sel2pdb[3*indexa2+1];
        vec2[2]=sel2pdb[3*indexb2+2]-sel2pdb[3*indexa2+2];

        vec2zaxis(vec2,sintheta2,costheta2,sinphi2,cosphi2);
        zaxis2vec(vec2,sintheta1,costheta1,sinphi1,cosphi1);

        // storage the mutated atoms' xyz
        indexb2 += nsel1;
        xyzq[4*indexb2+0]=vec2[0]+xyzq[4*indexa1+0];
        xyzq[4*indexb2+1]=vec2[1]+xyzq[4*indexa1+1];
        xyzq[4*indexb2+2]=vec2[2]+xyzq[4*indexa1+2];
        myindex++;
     }
  }

  // For debugging
  //for (i=0;i<nsel2;i++) {
  //   myindex=i+nsel1;
  //   printf("xyzq2 %d %f %f %f %f \n",myindex,xyzq[4*myindex+0],xyzq[4*myindex+1],xyzq[4*myindex+2],xyzq[4*myindex+3]);
  //}

}


void rot_sel2_dihed(int nsel2, int dihed_id, double sam_angle, double *xyzq, struct SAM *SamDihed) {
// Rotate one given dihedral angle of sel2 (xyzq start from sel2)
// Letitia 2015.11.12
  int i, j, indexa, indexb, disp;
  double theta, phi, length, implicit=0.0;
  double sintheta, costheta, sinphi, cosphi, sinangle, cosangle;
  double vec[3];

  // get vec
  indexa = SamDihed->list[4*dihed_id+1];
  indexb = SamDihed->list[4*dihed_id+2];
  vec[0] = xyzq[4*indexb]-xyzq[4*indexa];
  vec[1] = xyzq[4*indexb+1]-xyzq[4*indexa+1];
  vec[2] = xyzq[4*indexb+2]-xyzq[4*indexa+2];
     
  // get theta and phi
  length = sqrt(vec[0]*vec[0]+vec[1]*vec[1]);
  phi    = atan2(vec[1],vec[0]);
  theta  = atan(length/vec[2]);
  if (theta < 0) {theta += M_PI;} // range {-PI/2, PI/2} --> {0, PI}
  sintheta = sin(theta);
  costheta = cos(theta);
  sinphi = sin(phi);
  cosphi = cos(phi);

  // get the desired sampling angle (turn on the implicit if you want it to be corrected.)
  implicit=calc_dihed_angle(xyzq,SamDihed->list[4*dihed_id],SamDihed->list[4*dihed_id+1],SamDihed->list[4*dihed_id+2],SamDihed->list[4*dihed_id+3]);
  sinangle=sin(sam_angle-implicit);
  cosangle=cos(sam_angle-implicit);

  // rotate all atoms that is affected by this dihedral
  disp = sumint(SamDihed->naff,dihed_id-1);
  for (j=0;j<SamDihed->naff[dihed_id];j++) {
     indexb = SamDihed->afflist[j+disp];
     vec[0] = xyzq[4*indexb]-xyzq[4*indexa];
     vec[1] = xyzq[4*indexb+1]-xyzq[4*indexa+1];
     vec[2] = xyzq[4*indexb+2]-xyzq[4*indexa+2];

     vec2zaxis(vec,sintheta,costheta,sinphi,cosphi);
     rotate_around_zaxis(vec,sinangle,cosangle);
     zaxis2vec(vec,sintheta,costheta,sinphi,cosphi);

     // storage the new xyz value accordingly
     xyzq[4*indexb]   = vec[0]+xyzq[4*indexa];
     xyzq[4*indexb+1] = vec[1]+xyzq[4*indexa+1];
     xyzq[4*indexb+2] = vec[2]+xyzq[4*indexa+2];
  }

  // For debugging
  //for (i=0;i<nsel2;i++) {
  //   printf("xyzq2 %d %f %f %f %f \n",i,xyzq[4*i],xyzq[4*i+1],xyzq[4*i+2],xyzq[4*i+3]);
  //}
  
}


void rot_sel2_angle(int nsel2, int angle_id, double sam_angle, double *xyzq, struct SAM *SamAngle) {
// Rotate one given angle of sel2 (xyzq start from sel2)
// Letitia 2015.11.12
  int i, j, id1, id2, id3, disp;
  double angle, theta, phi, length;
  double sintheta, costheta, sinphi, cosphi, sinangle, cosangle;
  double vec1[3], vec2[3], vec3[3];

  // get vec1 vec2 vec3
  id1 = SamAngle->list[3*angle_id];
  id2 = SamAngle->list[3*angle_id+1];
  id3 = SamAngle->list[3*angle_id+2];
  for (j=0;j<3;j++) {
     vec1[j] = xyzq[4*id1+j]-xyzq[4*id2+j];
     vec2[j] = xyzq[4*id3+j]-xyzq[4*id2+j];
  }
  cross(vec1,vec2,vec3);

  // preparation for aligning vec3 to x-axis
  length = sqrt(vec3[0]*vec3[0]+vec3[1]*vec3[1]);
  phi    = atan2(vec3[1],vec3[0]);
  theta  = atan(length/vec3[2]);
  if (isnan(theta)) {theta=M_PI;}  // If norm(vec3)=0, theta=PI. (atan(0/0)=-nan)
  if (theta < 0) {theta += M_PI;} // range {-PI/2, PI/2} --> {0, PI}
  sintheta = sin(theta);
  costheta = cos(theta);
  sinphi = sin(phi);
  cosphi = cos(phi);

  // get original angle (under transformed coordinate, signed)
  vec2xaxis(vec1,costheta,sintheta,cosphi,sinphi);  
  vec2xaxis(vec2,costheta,sintheta,cosphi,sinphi);  
  angle=acos(innerproduct(vec1,vec2,3)/(norm(vec1,3)*norm(vec2,3)));
  if (isnan(angle)) {angle=M_PI;}  // Somehow acos sometimes fails to report acos(-1.00)=PI, but reports -nan.

  // get the desired sampling angle, rotating respect to the original angle.
  sinangle = sin(sam_angle-angle);
  cosangle = cos(sam_angle-angle);

  // rotate all atoms that is affected by this dihedral
  disp = sumint(SamAngle->naff,angle_id-1);
  for (j=0;j<SamAngle->naff[angle_id];j++) {
     id3 = SamAngle->afflist[j+disp];
     vec2[0] = xyzq[4*id3]-xyzq[4*id2];
     vec2[1] = xyzq[4*id3+1]-xyzq[4*id2+1];
     vec2[2] = xyzq[4*id3+2]-xyzq[4*id2+2];

     vec2xaxis(vec2,costheta,sintheta,cosphi,sinphi);
     rotate_around_xaxis(vec2,sinangle,cosangle);
     xaxis2vec(vec2,costheta,sintheta,cosphi,sinphi);

     // storage the new xyz value accordingly
     xyzq[4*id3]   = vec2[0]+xyzq[4*id2];
     xyzq[4*id3+1] = vec2[1]+xyzq[4*id2+1];
     xyzq[4*id3+2] = vec2[2]+xyzq[4*id2+2];
  }

  // For debugging
  //for (i=0;i<nsel2;i++) {
  //   printf("xyzq2 %d %f %f %f %f \n",i,xyzq[4*i],xyzq[4*i+1],xyzq[4*i+2],xyzq[4*i+3]);
  //}
}


