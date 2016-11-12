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
#include "vss.h"

#if !defined(M_PI)
#define M_PI 3.14159265358979323846
#endif

double calc_Eb(double *xyzq, struct BOND *Dihed, struct BOND *Angle) 
{
  // Calculate the bonded energy (only used for sel2)
  int i, j, k, l, m, n;
  double Eb, phi, delta_phi, phi_ub, delta_phi_ub, Eangle; 
  double rij[3], rjk[3], rkl[3], vecA[3], vecB[3];
  
  Eb=0.0;
  
  // Calculate the dihedral energy
  for (m = 0; m < Dihed->nentry; m++) {

     phi=calc_dihed_angle(xyzq,Dihed->list[4*m],Dihed->list[4*m+1],Dihed->list[4*m+2],Dihed->list[4*m+3]);

     if (Dihed->para[3*m+1] > 0) {
        delta_phi = Dihed->para[3*m+1]*phi-Dihed->para[3*m+2]*M_PI/180.0;
        Eb += Dihed->para[3*m]*(1+cos(delta_phi)); 
     } else {
        delta_phi = phi-Dihed->para[3*m+2]*M_PI/180.0;
        Eb += Dihed->para[3*m]*delta_phi*delta_phi;
     }

  }
 
  // Calculate the angle energy
  for (m = 0; m < Angle->nentry; m++) {
     // angle unit in rad
     phi = calc_angle(xyzq,Angle->list[3*m],Angle->list[3*m+1],Angle->list[3*m+2]);
     delta_phi = phi - Angle->para[4*m+1]*M_PI/180.0;
     phi_ub = calc_UB(xyzq,Angle->list[3*m],Angle->list[3*m+2]);     
     delta_phi_ub = phi_ub - Angle->para[4*m+3];
     Eb += Angle->para[4*m]*delta_phi*delta_phi + Angle->para[4*m+2]*delta_phi_ub*delta_phi_ub;
  }
  
  return Eb;

}


double calc_dihed_angle(double *xyzq, int i, int j, int k, int l) {
// Calculate the dihedral angle for a given index.
  int n;
  double xval, yval, phi;
  double rij[3], rjk[3], rkl[3], vecA[3], vecB[3];

  for (n=0;n<3;n++) {
      rij[n]=xyzq[4*j+n]-xyzq[4*i+n];
      rjk[n]=xyzq[4*k+n]-xyzq[4*j+n];
      rkl[n]=xyzq[4*l+n]-xyzq[4*k+n];
  }
  cross(rij,rjk,vecA);
  cross(rjk,rkl,vecB);
  // xval = cos_phi*|A|*|B|; yval = sin_phi*|A|*|B|
  xval = innerproduct(vecA,vecB,3);
  yval = innerproduct(vecA,rkl,3)*sqrt(innerproduct(rjk,rjk,3));
  phi=atan2(yval,xval);

  return phi;
}


double calc_angle(double *xyzq, int i, int j, int k) {
// Calculate the unsigned angle for a given index. (Valid for energy calculation only!)
  int n;
  double phi;
  double vec1[3], vec2[3];

  // get vec1 vec2 vec3
  for (n=0;n<3;n++) {
     vec1[n] = xyzq[4*i+n]-xyzq[4*j+n];
     vec2[n] = xyzq[4*k+n]-xyzq[4*j+n];
  }

  phi=acos(innerproduct(vec1,vec2,3)/(norm(vec1,3)*norm(vec2,3)));
  if (isnan(phi)) {phi=M_PI;}  // Somehow acos sometimes fails to report acos(-1.00)=PI, but reports -nan.

  return phi;
}


double calc_UB(double *xyzq, int i, int k) {
// Calculate the UB angle term for a given index
  int n;
  double phi;
  double vec[3];

  // get vec
  for (n=0;n<3;n++) {
     vec[n] = xyzq[4*i+n]-xyzq[4*k+n];
  }
  phi = sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);

  return phi;
}
