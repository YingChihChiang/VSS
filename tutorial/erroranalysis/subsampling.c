#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

void chkfile(char *namefile, int *hdim, int *ddim);
void rdfile(char *namefile, int hdim, int ddim, int step, int *xdata, double *ydata);
double sum_array(double *ydata, int ddim);
void linear_regression(double *xdata, double *ydata, int ndim, double *alpha, double *beta);


void main() {
  int i, j, k, n, myindex, count, step, nstep, order;
  int ddim, ndim, mdim, blockdim, hdim[2], orgdim, repdim=100, *xdata;
  float temperature, rwork;
  double *ydata, *dGjdata, *dG_xdata, *dG_ydata, *var, *blockdata;
  double RT, tempa, alpha, beta, max_val_err, val_err;
  double const Boltzmann=0.001987191;
  FILE *fptr;
  char namefile[80], mystring[256], newstring[80]; 

  printf("Name of the file (frame  DeltaE)? \n");
  scanf("%s",namefile);
  printf("Temperature (K)? \n");
  scanf("%g",&temperature);
  printf("Take data every ? step \n");
  scanf("%d",&step);

  // Initialize the variables and read in vss.log
  RT = Boltzmann*temperature;
  chkfile(namefile,hdim,&orgdim);
  ddim = (orgdim/step);
  if(orgdim%step != 0) {ddim++;}
  xdata = malloc(ddim*sizeof(int));
  ydata = malloc(ddim*sizeof(double));
  blockdata = malloc(ddim*sizeof(double));
  dGjdata = malloc(repdim*ddim*sizeof(double)); 
  dG_xdata = malloc(ddim*sizeof(double));
  dG_ydata = malloc(ddim*sizeof(double));
  var = malloc(ddim*sizeof(double));

  rdfile(namefile,hdim[0],orgdim,step,xdata,ydata);

  // Turn Delta_E into Exp_(-beta*Delta_E)
  for (i=0;i<ddim;i++) {
    ydata[i] = exp(-ydata[i]/RT);
  }

  fptr = fopen("subsampling.log","w");
  fprintf(fptr,"# n*frame_step, dG(n), sigma^2(n), Uncertainty(n) \n");
  srand(time(NULL));
  count = 0;
  max_val_err=0;

  // Perform subsampling method (only for every possible value, note that the last one's error is poorly estimated.) 
  nstep=1;
  for (n=1;n<ddim+1;n+=nstep) {

     mdim = round(ddim/n);
     blockdim = repdim*mdim;

     for (i=0;i<ddim;i++) {blockdata[i]=ydata[i];} 

     // Construct the block data
     for (j=0;j<blockdim;j++) {

        // Sample ydata[0:n] without replacement --> shuffle the original data
        for (i=0;i<n;i++) {
           myindex = rand()%(ddim-i);
           tempa = blockdata[i];
           blockdata[i] = blockdata[i+myindex];
           blockdata[i+myindex] = tempa;
        }
           
        // Evaluate dGj (Fj) for the current block, average dGj gives dG (F) of the current n
        dGjdata[j] = sum_array(blockdata,n)/n;
        dGjdata[j] = -RT*log(dGjdata[j]);
     }
     
     dG_xdata[count] = sqrt(1.0/n); 
     dG_ydata[count] = sum_array(dGjdata,blockdim)/blockdim; 
     if (isinf(dG_ydata[count])) { continue; }

     // Evaluate the variance
     var[count] = 0;
     for (j=0;j<blockdim;j++) {
        var[count] += (dGjdata[j]-dG_ydata[count])*(dGjdata[j]-dG_ydata[count]);
     }
     var[count] = var[count]/blockdim;

     // Print dG(n), var(n), Error(n)
     val_err = 2*sqrt(n*var[count]/ddim);
     if (val_err - max_val_err > 0.0) { max_val_err = val_err;}

     fprintf(fptr,"%d %E %E %E \n",xdata[n-1],dG_ydata[count],var[count],val_err);

     if (xdata[n-1] >= 1000000) {
        nstep=25000; } 
     else if (xdata[n-1] >= 100000) {
        nstep=5000; }
     else if (xdata[n-1] >= 10000) {
        nstep=500; }
     else if (xdata[n-1] >= 1000) {
        nstep=50; }
     else if (xdata[n-1] >= 100) {
        nstep=5; }
     else {
        nstep=1;
     }

     count +=1;
  }
  fprintf(fptr,"%d %E %E %E \n",xdata[ddim-1],-RT*log(sum_array(ydata,ddim)/ddim),0.0,0.0);

  linear_regression(dG_xdata,dG_ydata,count,&alpha,&beta);
  fprintf(fptr,"# alpha and beta::  %E  %E \n",alpha,beta);
  fprintf(fptr,"# estimated dG upper bond:: %E \n",alpha);
  fprintf(fptr,"# Last error is highly underestimated because there is only one block! \n");
  fprintf(fptr,"# max. error::  %E \n",max_val_err);
  fclose(fptr);

  free(xdata);
  free(ydata);
  free(blockdata);
  free(dGjdata);
  free(dG_xdata);
  free(dG_ydata);
  free(var);

}


void chkfile(char *namefile, int *hdim, int *ddim) {
  int i, val;
  char mystring[256];
  FILE *fptr;

  fptr = fopen(namefile,"r");
  if (fptr == NULL) {printf("Error!! Cannot open the file. \n"); return;}

  hdim[0]=0;
  hdim[1]=0;
  *ddim=0;
  while(1) {
    fgets(mystring,256,fptr);
    if (feof(fptr)) break;
    if (strncmp(mystring,"#",1) == 0) {
      if (*ddim == 0) { 
        hdim[0] += 1;
      } else {
        hdim[1] += 1;
      }
    } else {
      *ddim += 1;
    }
  }

  fclose(fptr);
  return;
}


void rdfile(char *namefile, int hdim, int ddim, int step, int *xdata, double *ydata) {
  int i, j, ival;
  float dval;
  char mystring[256];
  FILE *fptr;

  fptr = fopen(namefile,"r");
  for (i=0;i<hdim;i++) {
    fgets(mystring,256,fptr);
  }

  for (i=0;i<ddim;i++) {
    fscanf(fptr,"%d  %E",&ival,&dval);
    if (i%step == 0) {
      xdata[(i/step)]=ival;
      ydata[(i/step)]=dval;
    }
  }
 
  fclose(fptr);
  return;
}


double sum_array(double *ydata, int ddim) {
  int i;
  double sum_val=0.0;
  
  for (i=0;i<ddim;i++) {
    sum_val += ydata[i];
  }
  return sum_val;
}


void linear_regression(double *xdata,double *ydata,int ndim,double *alpha,double *beta) {
  int i;
  double xbar, ybar, cov, var;

  xbar=sum_array(xdata,ndim)/ndim;
  ybar=sum_array(ydata,ndim)/ndim;

  cov=0.0;
  var=0.0;
  for (i=0;i<ndim;i++) {
     cov += (xdata[i]-xbar)*(ydata[i]-ybar);
     var += (xdata[i]-xbar)*(xdata[i]-xbar);
  }
  *beta  = cov/var;
  *alpha = ybar - *beta*xbar;

  return;
}

