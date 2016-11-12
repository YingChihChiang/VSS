#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

void chkfile(char *namefile, int *hdim, int *ddim);
void rdfile(char *namefile, int hdim, int ddim, double *xdata, double *ydata);
double sum_array(double *mydata, int ddim);
double eval_cov(double *data1, double *data2, int ddim);
void linear_regression(double *xdata,double *ydata,int ndim,double *alpha,double *beta); 
void minmax(double *mydata, int ddim, double *minval, double *maxval);

void main() {
  int i, j, k, ncut;
  int ddim, hdim[2], xcolumn, ycolumn;
  double *xdata, *ydata, *corr, cutoff, temperature;
  double avg, cov_val, var_val, alpha, beta, RT, minval, maxval, integral;
  double const Boltzmann=0.001987191;
  FILE *fptr;
  char namefile[80], mystring[256], newstring[80]; 

  printf("Name of the file (frame DeltaE)? \n");
  scanf("%s",namefile);
  printf("Temperature (K)? (0: use DeltaE instead of Exp(-beta*DeltaE)) \n");
  scanf("%lf",&temperature);
  printf("Cut off (0: take default value 0.05)? \n");
  scanf("%lf",&cutoff);

  // Initialize the variables and read in file
  RT=temperature*Boltzmann;
  chkfile(namefile,hdim,&ddim);
  xdata=malloc(ddim*sizeof(double));
  ydata=malloc(ddim*sizeof(double));
  corr=malloc(ddim*sizeof(double));
  rdfile(namefile,hdim[0],ddim,xdata,ydata);
  if (cutoff==0.0) {cutoff=0.05;}
 
  // Turn Energy into Exp(-beta*E)
  if (temperature != 0.0) { 
     minmax(ydata,ddim,&minval,&maxval);
     for (i=0;i<ddim;i++) {
        ydata[i] -= minval;    // Rescale the data for 'cov/var' stability. (Works only for Exp case!)
        ydata[i] = exp(-ydata[i]/RT);
     }
  }

  // Evaluate Average Y and Delta Y 
  avg=sum_array(ydata,ddim);
  avg=avg/ddim;
  for (i=0;i<ddim;i++) {ydata[i]-=avg;}

  // Evaluate Auto Correlation Function
  var_val=eval_cov(ydata,ydata,ddim);
  fptr=fopen("corr.log","w");
  for (i=0;i<ddim;i++) {
     corr[i]=eval_cov(ydata,&ydata[i],ddim-i);
     corr[i]=corr[i]/var_val;
     if (corr[i] > cutoff) {
        fprintf(fptr,"%g %g \n",xdata[i],corr[i]);
     } else {
        break;
     }
  }

  // Evaluate Relaxation Time via integration
  integral=0;
  for (j=0;j<i;j++) {
     integral += corr[j]; 
  }
  integral -= 0.5*(corr[0]+corr[i-1]);
  integral = integral*(xdata[1]-xdata[0]);
  fprintf(fptr,"# correlation time via integration: %g \n",integral);
  fprintf(fptr,"# statistical inefficiency: %g \n",1+2*integral);
  printf("# correlation time via integration: %g \n",integral);
  printf("# statistical inefficiency: %g \n",1+2*integral);
  
  fclose(fptr);

  free(xdata);
  free(ydata);
  free(corr);
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


void rdfile(char *namefile, int hdim, int ddim, double *xdata, double *ydata) {
  int i, j;
  float val[2];
  char mystring[256];
  FILE *fptr;

  fptr = fopen(namefile,"r");
  for (i=0;i<hdim;i++) {
    fgets(mystring,256,fptr);
  }

  for (i=0;i<ddim;i++) {
    fscanf(fptr,"%E  %E",&val[0],&val[1]);
    xdata[i]=val[0];
    ydata[i]=val[1];
  }
 
  fclose(fptr);
  return;
}


double sum_array(double *mydata, int ddim) {
  int i;
  double sum_val=0.0;
  
  for (i=0;i<ddim;i++) {
    sum_val += mydata[i];
  }
  return sum_val;
}


double eval_cov(double *data1, double *data2, int ddim) {
  int i;
  double cov=0.0;
  
  for (i=0;i<ddim;i++) {
     cov += data1[i]*data2[i];
  }

  return cov;
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


void minmax(double *mydata, int ddim, double *minval, double *maxval) {
  int i;

  *minval=mydata[0];
  *maxval=mydata[0];
  for (i=1;i<ddim;i++) {
    if(mydata[i]-*minval < 0.0) {*minval=mydata[i];}
    if(mydata[i]-*maxval > 0.0) {*maxval=mydata[i];}
  }
  return;
}

