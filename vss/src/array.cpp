#include <stdio.h>
#include <math.h>
#include "vss.h"

int sumint(int *arr, int ndim) {
   int i, val=0;
 
   for (i=0;i<ndim;i++) {
      val += arr[i];
   }
   return val;
}


double sumdouble(double *arr, int ndim) {
   int i;
   double val=0.0;

   for (i=0;i<ndim;i++) {
      val += arr[i];
   }
   return val;
}


double norm(double *arr, int ndim) {
   int i;
   double val=0.0;
 
   for (i=0;i<ndim;i++) {
      val += arr[i]*arr[i];
   }
   val = sqrt(val);

   return val;
}

void cross(double *vecA, double *vecB, double *vecC) {
// Return C = A x B

  vecC[0] = vecA[1]*vecB[2] - vecA[2]*vecB[1];
  vecC[1] = vecA[2]*vecB[0] - vecA[0]*vecB[2];
  vecC[2] = vecA[0]*vecB[1] - vecA[1]*vecB[0];

  return;
}


double innerproduct(double *vecA, double *vecB, int ndim) {
// Return val = (A|B)
  int i;
  double val=0.0;

  for (i=0;i<ndim;i++) {
     val += vecA[i]*vecB[i];
  }

  return val;
}


void vec2zaxis(double *vec, double sintheta, double costheta, double sinphi, double cosphi) {
// From vec to z-axis.
  double x, y, z;

  x=vec[0];
  y=vec[1];
  z=vec[2];
  vec[0] =  costheta*cosphi*x + costheta*sinphi*y - sintheta*z;
  vec[1] = -sinphi*x + cosphi*y;
  vec[2] =  sintheta*cosphi*x + sintheta*sinphi*y + costheta*z;
}


void zaxis2vec(double *vec, double sintheta, double costheta, double sinphi, double cosphi) {
// From z-axis to vec, transpose of matrix from vec to z-axis.
  double x, y, z;

  x=vec[0];
  y=vec[1];
  z=vec[2];
  vec[0] =  costheta*cosphi*x - sinphi*y + sintheta*cosphi*z;
  vec[1] =  costheta*sinphi*x + cosphi*y + sintheta*sinphi*z;
  vec[2] = -sintheta*x + costheta*z;       
}


void rotate_around_zaxis(double *vec, double sinangle, double cosangle) {
  double x, y;

  x=vec[0];
  y=vec[1];
  vec[0] = cosangle*x - sinangle*y;
  vec[1] = sinangle*x + cosangle*y;
}


void rotate_around_xaxis(double *vec, double sinangle, double cosangle) {
  double y, z;

  y=vec[1];
  z=vec[2];
  vec[1] = cosangle*y - sinangle*z;
  vec[2] = sinangle*y + cosangle*z;
}


void vec2xaxis(double *vec, double costheta, double sintheta, double cosphi, double sinphi){
// From vec to x-axis.
  double x, y, z;

  x=vec[0];
  y=vec[1];
  z=vec[2];
  vec[0] =  sintheta*cosphi*x + sintheta*sinphi*y + costheta*z;
  vec[1] = -sinphi*x + cosphi*y;
  vec[2] = -costheta*cosphi*x - costheta*sinphi*y + sintheta*z;
}


void xaxis2vec(double *vec, double costheta, double sintheta, double cosphi, double sinphi){
// From x-axis to vec, transpose of matrix from vec to x-axis.
  double x,y,z;

  x=vec[0];
  y=vec[1];
  z=vec[2];
  vec[0] = sintheta*cosphi*x - sinphi*y - costheta*cosphi*z;
  vec[1] = sintheta*sinphi*x + cosphi*y - costheta*sinphi*z;
  vec[2] = costheta*x + sintheta*z;
}


