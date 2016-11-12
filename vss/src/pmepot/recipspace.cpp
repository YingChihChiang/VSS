/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/* $Id */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vss.h"

#if !defined(M_PI)
#define M_PI 3.14159265358979323846
#endif

static void dftmod(double *bsp_mod, double *bsp_arr, int nfft) {
  int j, k;
  double twopi, arg, sum1, sum2;
  double infft = 1.0/nfft;
/* Computes the modulus of the discrete fourier transform of bsp_arr, */
/*  storing it into bsp_mod */
  twopi =  2.0 * M_PI;

  for (k = 0; k <nfft; ++k) {
    sum1 = 0.;
    sum2 = 0.;
    for (j = 0; j < nfft; ++j) {
      arg = twopi * k * j * infft;
      sum1 += bsp_arr[j] * cos(arg);
      sum2 += bsp_arr[j] * sin(arg);
    }
    bsp_mod[k] = sum1*sum1 + sum2*sum2;
  }
}

void compute_b_spline(double *frac, double *M);

static void compute_b_moduli(double *bm, int K) {
  int i;
  double fr[3];

  double M[12];
  double *scratch = (double*) malloc(sizeof(double)*K);

  fr[0]=fr[1]=fr[2]=0.0;
  compute_b_spline(fr,M);  
  for (i=0; i<4; i++) bm[i] = M[i];
  for (i=4; i<K; i++) bm[i] = 0.0;
  dftmod(scratch, bm, K);
  for (i=0; i<K; i++) bm[i] = 1.0/scratch[i];

  free(scratch);
}

static void init_exp(double *xp, int K, double recip, double piob) {
  int i;
  double fac;
  fac = -piob*recip*recip;
  for (i=0; i<= K/2; i++)
    xp[i] = exp(fac*i*i);
} 

static double dot_product(const double *v1, const double *v2) {
  return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

static void cross_product(double *c, const double *v1, const double *v2) {
  c[0] = v1[1]*v2[2]-v2[1]*v1[2];
  c[1] = v2[0]*v1[2]-v1[0]*v2[2];
  c[2] = v1[0]*v2[1]-v2[0]*v1[1];
}

static double cell_volume(const double *cell) {
  double v, t[3];
  cross_product(t,&cell[3],&cell[6]);
  v = dot_product(t,&cell[9]);
  return abs(v);
}

static int is_parallel(const double *v1, const double *v2) {
  double t[3];
  cross_product(t,v1,v2);
  return ( t[0] == 0. && t[1] == 0. && t[2] == 0. );
}


double compute_energy(double *q_arr, const double *cell, const double *rcell,
			const int *dims, double ewald) {

  double energy = 0.0;
  double piob, i_pi_volume;
  double *bm1, *bm2, *bm3, *exp1, *exp2, *exp3;
  int pad2, pad3;
  int k1, k2, k3, ind;
  int K1, K2, K3, dim2, dim3;

  K1=dims[0]; K2=dims[1]; K3=dims[2]; dim2= dims[3]; dim3=dims[4];
  pad2 = (dim2-K2)*dim3;
  pad3 = dim3-K3-(K3 & 1 ? 1 : 2);

  bm1 = (double*) malloc(sizeof(double)*K1);
  bm2 = (double*) malloc(sizeof(double)*K2);
  bm3 = (double*) malloc(sizeof(double)*K3);

  exp1 = (double*) malloc(sizeof(double)*(K1/2 + 1));
  exp2 = (double*) malloc(sizeof(double)*(K2/2 + 1));
  exp3 = (double*) malloc(sizeof(double)*(K3/2 + 1));

  compute_b_moduli(bm1, K1);
  compute_b_moduli(bm2, K2);
  compute_b_moduli(bm3, K3);

  i_pi_volume = 1.0/(M_PI * cell_volume(cell));
  piob = M_PI/ewald;
  piob *= piob;

  /* for (n=0; n<6; virial[n++] = 0.0); */

  if ( cell[4] == 0. && cell[5] == 0. &&
       cell[6] == 0. && cell[8] == 0. &&
       cell[9] == 0. && cell[10] == 0. ) {

    double recipx = rcell[3];
    double recipy = rcell[7];
    double recipz = rcell[11];
    init_exp(exp1, K1, recipx, piob);
    init_exp(exp2, K2, recipy, piob);
    init_exp(exp3, K3, recipz, piob);

    for ( k1=0; k1<K1; ++k1 ) {
      double m1, m11, b1, xp1;
      int k1_s = k1<=K1/2 ? k1 : k1-K1;
      b1 = bm1[k1];
      m1 = k1_s*recipx;
      m11 = m1*m1;
      xp1 = i_pi_volume*exp1[abs(k1_s)];
      for ( k2=0; k2<K2; ++k2 ) {
        double m2, m22, b1b2, xp2;
        int k2_s = k2<=K2/2 ? k2 : k2-K2;
        b1b2 = b1*bm2[k2];
        m2 = k2_s*recipy;
        m22 = m2*m2;
        xp2 = exp2[abs(k2_s)]*xp1;
        if ( k1==0 && k2==0 ) {
          q_arr[0] = 0.0;
          q_arr[1] = 0.0;
          k3 = 1;
        } else {
          k3 = 0;
        }
        for ( ; k3<=K3/2; ++k3 ) {
          int mysize = (K3/2)+1;
          ind = 2*(k1*K2*mysize+k2*mysize+k3)+pad3*(k1*K2+k2)+pad2*k1;  // Letitia: remove index dependency;
          double m3, m33, xp3, msq, imsq, vir, fac;
          double theta3, theta, q2, qr, qc, C;
          theta3 = bm3[k3] *b1b2;
          m3 = k3*recipz;
          m33 = m3*m3;
          xp3 = exp3[k3];
          qr = q_arr[ind]; qc=q_arr[ind+1];
          q2 = 2*(qr*qr + qc*qc)*theta3;
          if ( (k3 == 0) || ( k3 == K3/2 && ! (K3 & 1) ) ) q2 *= 0.5;
          msq = m11 + m22 + m33;
          imsq = 1.0/msq;
          C = xp2*xp3*imsq;
          theta = theta3*C;
          q_arr[ind] *= theta;
          q_arr[ind+1] *= theta;
          vir = -2*(piob+imsq);
          fac = q2*C;
          energy += fac;
          /*
          virial[0] += fac*(1.0+vir*m11);
          virial[1] += fac*vir*m1*m2;
          virial[2] += fac*vir*m1*m3;
          virial[3] += fac*(1.0+vir*m22);
          virial[4] += fac*vir*m2*m3;
          virial[5] += fac*(1.0+vir*m33);
          */
        }
      }
    }

  } else if ( is_parallel(&cell[9],&rcell[9]) ) {

    const double *recip1, *recip2, *recip3;
    double recip3_x, recip3_y, recip3_z;
    recip1 = &rcell[3];
    recip2 = &rcell[6];
    recip3 = &rcell[9];
    recip3_x = recip3[0];
    recip3_y = recip3[1];
    recip3_z = recip3[2];
    init_exp(exp3, K3, sqrt(dot_product(recip3,recip3)), piob);

    for ( k1=0; k1<K1; ++k1 ) {
      double b1, m1[3];
      int k1_s = k1<=K1/2 ? k1 : k1-K1;
      b1 = bm1[k1];
      m1[0] = k1_s*recip1[0];
      m1[1] = k1_s*recip1[1];
      m1[2] = k1_s*recip1[2];
      /* xp1 = i_pi_volume*exp1[abs(k1_s)]; */
      for ( k2=0; k2<K2; ++k2 ) {
        double xp2, b1b2, m2_x, m2_y, m2_z;
        int k2_s = k2<=K2/2 ? k2 : k2-K2;
        b1b2 = b1*bm2[k2];
        m2_x = m1[0] + k2_s*recip2[0];
        m2_y = m1[1] + k2_s*recip2[1];
        m2_z = m1[2] + k2_s*recip2[2];
        /* xp2 = exp2[abs(k2_s)]*xp1; */
        xp2 = i_pi_volume*exp(-piob*(m2_x*m2_x+m2_y*m2_y+m2_z*m2_z));
        if ( k1==0 && k2==0 ) {
          q_arr[0] = 0.0;
          q_arr[1] = 0.0;
          k3 = 1;
        } else {
          k3 = 0;
        }
        for ( ; k3<=K3/2; ++k3 ) {
          int mysize=(K3/2)+1;
          ind = 2*(k1*K2*mysize+k2*mysize+k3)+pad3*(k1*K2+k2)+pad2*k1;  // Letitia: remove index dependency;
          double xp3, msq, imsq, vir, fac;
          double theta3, theta, q2, qr, qc, C;
          double m_x, m_y, m_z;
          theta3 = bm3[k3] *b1b2;
          m_x = m2_x + k3*recip3_x;
          m_y = m2_y + k3*recip3_y;
          m_z = m2_z + k3*recip3_z;
          msq = m_x*m_x + m_y*m_y + m_z*m_z;
          xp3 = exp3[k3];
          qr = q_arr[ind]; qc=q_arr[ind+1];
          q2 = 2*(qr*qr + qc*qc)*theta3;
          if ( (k3 == 0) || ( k3 == K3/2 && ! (K3 & 1) ) ) q2 *= 0.5;
          imsq = 1.0/msq;
          C = xp2*xp3*imsq;
          theta = theta3*C;
          q_arr[ind] *= theta;
          q_arr[ind+1] *= theta;
          vir = -2*(piob+imsq);
          fac = q2*C;
          energy += fac;
          /*
          virial[0] += fac*(1.0+vir*m_x*m_x);
          virial[1] += fac*vir*m_x*m_y;
          virial[2] += fac*vir*m_x*m_z;
          virial[3] += fac*(1.0+vir*m_y*m_y);
          virial[4] += fac*vir*m_y*m_z;
          virial[5] += fac*(1.0+vir*m_z*m_z);
          */
        }
      }
    }

  } else {

    const double *recip1, *recip2, *recip3;
    double recip3_x, recip3_y, recip3_z;
    recip1 = &rcell[3];
    recip2 = &rcell[6];
    recip3 = &rcell[9];
    recip3_x = recip3[0];
    recip3_y = recip3[1];
    recip3_z = recip3[2];

    for ( k1=0; k1<K1; ++k1 ) {
      double b1, m1[3];
      int k1_s = k1<=K1/2 ? k1 : k1-K1;
      b1 = bm1[k1];
      m1[0] = k1_s*recip1[0];
      m1[1] = k1_s*recip1[1];
      m1[2] = k1_s*recip1[2];
      /* xp1 = i_pi_volume*exp1[abs(k1_s)]; */
      for ( k2=0; k2<K2; ++k2 ) {
        double b1b2, m2_x, m2_y, m2_z;
        int k2_s = k2<=K2/2 ? k2 : k2-K2;
        b1b2 = b1*bm2[k2];
        m2_x = m1[0] + k2_s*recip2[0];
        m2_y = m1[1] + k2_s*recip2[1];
        m2_z = m1[2] + k2_s*recip2[2];
        /* xp2 = exp2[abs(k2_s)]*xp1; */
        if ( k1==0 && k2==0 ) {
          q_arr[0] = 0.0;
          q_arr[1] = 0.0;
          k3 = 1;
        } else {
          k3 = 0;
        }
        for ( ; k3<=K3/2; ++k3 ) {
          int mysize=(K3/2)+1;
          ind = 2*(k1*K2*mysize+k2*mysize+k3)+pad3*(k1*K2+k2)+pad2*k1;  // Letitia: remove index dependency;
          double xp3, msq, imsq, vir, fac;
          double theta3, theta, q2, qr, qc, C;
          double m_x, m_y, m_z;
          theta3 = bm3[k3] *b1b2;
          m_x = m2_x + k3*recip3_x;
          m_y = m2_y + k3*recip3_y;
          m_z = m2_z + k3*recip3_z;
          msq = m_x*m_x + m_y*m_y + m_z*m_z;
          /* xp3 = exp3[k3]; */
          xp3 = i_pi_volume*exp(-piob*msq);
          qr = q_arr[ind]; qc=q_arr[ind+1];
          q2 = 2*(qr*qr + qc*qc)*theta3;
          if ( (k3 == 0) || ( k3 == K3/2 && ! (K3 & 1) ) ) q2 *= 0.5;
          imsq = 1.0/msq;
          C = xp3*imsq;
          theta = theta3*C;
          q_arr[ind] *= theta;
          q_arr[ind+1] *= theta;
          vir = -2*(piob+imsq);
          fac = q2*C;
          energy += fac;
          /*
          virial[0] += fac*(1.0+vir*m_x*m_x);
          virial[1] += fac*vir*m_x*m_y;
          virial[2] += fac*vir*m_x*m_z;
          virial[3] += fac*(1.0+vir*m_y*m_y);
          virial[4] += fac*vir*m_y*m_z;
          virial[5] += fac*(1.0+vir*m_z*m_z);
          */
        }
      }
    }

  }

  free(bm1);
  free(bm2);
  free(bm3);
  
  free(exp1);
  free(exp2);
  free(exp3);

  /* for (n=0; n<6; ++n) virial[n] *= 0.5; */
  return 0.5*energy;
}

