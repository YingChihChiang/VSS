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
#include <time.h>
#include <math.h>
#include "vss.h"
#include "pub3dfft.h"

pmepot_data* pmepot_create(int *dims, double ewald_factor) {
  pmepot_data *data;
  int grid_size, max_dim;

  if ( dims[0] < 8 ) return 0;
  if ( dims[1] < 8 ) return 0;
  if ( dims[2] < 8 ) return 0;
  if ( dims[2] % 2 ) return 0;
  if ( ewald_factor <= 0. ) return 0;

  data = (typeof(data)) malloc(sizeof(pmepot_data));
  if ( ! data ) return 0;

  data->avg_count = 0;
  data->ewald_factor = ewald_factor;
  data->dims[0] = dims[0];
  data->dims[1] = dims[1];
  data->dims[2] = dims[2];
  data->dims[3] = dims[1];
  data->dims[4] = dims[2] + 2;
  grid_size = data->dims[0] * data->dims[3] * data->dims[4];
  data->grid_size = grid_size;
  max_dim = dims[0] > dims[1] ? dims[0] : dims[1];
  max_dim = max_dim > dims[2] ? max_dim : dims[2];
  data->max_dim = max_dim;
  data->fft_ntable = 4*max_dim+15;

  data->avg_potential = (double*) malloc(grid_size * sizeof(double));
  data->fft_table = (double*) malloc(3 * data->fft_ntable * sizeof(double));
  data->fft_work = (double*) malloc(2 * max_dim * sizeof(double));
  data->q_arr = (double*) malloc(data->grid_size * sizeof(double));
  if ( ! data->avg_potential || ! data->fft_table || ! data->fft_work ) {
    if ( data->avg_potential) free(data->avg_potential);
    if ( data->fft_table) free(data->fft_table);
    if ( data->fft_work) free(data->fft_work);
    free(data);
    return 0;
  }

  pubd3di(dims[2], dims[1], dims[0], data->fft_table, data->fft_ntable);

  return data;
}

void pmepot_destroy(pmepot_data *data) {
  free(data->avg_potential);
  free(data->fft_table);
  free(data->fft_work);
  free(data->q_arr);
  free(data);
}

double veclength(double vec[]){
    int i;
    double length=0.;
    for(i=0; i<3; i++)
	length += vec[i] * vec[i];
    return sqrt(length);
}

void make_grid_from_cell(int nr[], double cell[], double resolution){
    double a = veclength(cell+3) / resolution;
    double b = veclength(cell+6) / resolution;
    double c = veclength(cell+9) / resolution;

    nr[0] = good_fft_dim(a);
    nr[1] = good_fft_dim(b);
    nr[2] = good_fft_dim(c);
}

int good_fft_dim(double r){
    int max=512;
    int nl[25]={8,10,12,16,20,24,30,32,36,40,48,50,56,60,64,72,80,
	84,88,96,100,108,112,120,128};
    int ml[8]={1,2,3,4,5,6,8,10};
    int goodn = 8;
    int mi, ni, n;

    for(mi=0; mi<8; mi++){
	for(ni=0; ni<25; ni++){
	    n = ml[mi] * nl[ni];
	    if(n>max) return goodn;
	    if(n > goodn) goodn = n;
	    if(goodn >= r) return goodn;
	}
    }
    return goodn;
}
