#include "dcdplugin.c"
#define CHECK_ERR(RC, VERBAL) \
    if (RC < 0) {\
	if(VERBAL)\
	    fprintf(stderr, "skip/read_dcdstep failed: %d\n", RC);\
	return ERRCODE;\
    }

#define ERRCODE 0
void *initDcdhandle(char dcdfile[],  int firststep){
    int i, rc, natoms;
    void *v;
    dcdhandle *dcd;

    v = open_dcd_read(dcdfile, "dcd", &natoms);
    if (!v) {
	fprintf(stderr, "open_dcd_read failed for file %s\n", dcdfile);
	return ERRCODE;
    }
    dcd = (dcdhandle*) v;

    for(i=0; i<firststep; i++){
	rc = skip_dcdstep(dcd->fd, natoms, dcd->nfixed, dcd->charmm);
	dcd->first = 0;
	CHECK_ERR(rc,1);
    }
    return v;
}

#undef ERRCODE
#define ERRCODE 1

int getNextXYZ(double xyz[], void *v, int fstep){
    float unitcell[6] = {1.0f, 1.0f, 1.0f, 90.0f, 90.0f, 90.0f};
    int i, rc;
    dcdhandle *dcd = (dcdhandle*) v;

    rc = read_dcdstep(dcd->fd, dcd->natoms, dcd->x, dcd->y, dcd->z, unitcell,
	    dcd->nfixed, dcd->first, dcd->freeind, dcd->fixedcoords,
	    dcd->reverse, dcd->charmm);
    dcd->first = 0;
    CHECK_ERR(rc,1);

    //write data
    for(i=0; i<dcd->natoms; i++){
	xyz[3*i  ]=dcd->x[i];
	xyz[3*i+1]=dcd->y[i];
	xyz[3*i+2]=dcd->z[i];
    }

    //skip next (fstep - 1) frames
    for(i=1; i<fstep; i++){
	rc = skip_dcdstep(dcd->fd, dcd->natoms, dcd->nfixed, dcd->charmm);
	dcd->first = 0;
	CHECK_ERR(rc,0);
    }

    return 0;
}

void cutXYZ(double dest[], double src[], int index[], int n){
    int i;
    double *head = &(src[3 * index[0]]);
    int length=1;

    for(i=1; i<n; i++){
	if(index[i-1] + 1 == index[i]){
	    length++;
	}
	else{
	    memcpy(dest, head, sizeof(double) * length * 3);
	    dest += 3*length;
	    head = &(src[3 * index[i]]);
	    length = 1;
	}
    }
    memcpy(dest, head, sizeof(double) * length * 3);
}

void cutXYZQ(double xyzq[], double xyz[], double q[], int index[], int n){
    int i, j;

    for(i=0; i<n; i++){
	j = index[i];
	memcpy(&xyzq[4*i], &xyz[3*j], sizeof(double) * 3);
	xyzq[4*i+3] = q[j];
    }
}

void mergeXYZQ(double xyzq[], double xyz[], double q[], int n){
    int i;
    for(i=0; i<n; i++){
	memcpy(&xyzq[4*i], &xyz[3*i], sizeof(double) * 3);
	xyzq[4*i+3] = q[i];
    }
}

void updateXYZofAtom(double xyzq[], double xyz[], int n){
    int i;
    for(i=0; i<n; i++){
	memcpy(&xyzq[4*i], &xyz[3*i], sizeof(double) * 3);
    }
}

void closeDcdhandle(void *v){
    close_file_read((dcdhandle*)v);
}

int getNsets(void *v){
    return ((dcdhandle*)v)->nsets;
}

int getNatoms(void *v){
    return ((dcdhandle*)v)->natoms;
}

#ifdef TEST_READDCD
#include <stdio.h>
#include <stdlib.h>
int main(){
    double *dcd_xyz, tar_xyz[6*3];
    void *dcd;
    int n=6, tmax;
    int index[6]={6,7,8,10,11,13};
    int start=1, freq=3;

    int i, j;

    dcd = initDcdhandle("/home/andrew/vmd/VSS/data_ben2clb_trial/ben2chlben.dcd", start);
    if(!dcd){
	printf("initError!\n");
	return 1;
    }
    dcd_xyz = (double*)malloc(sizeof(double) * getNatoms(dcd) * 3);
    tmax = getNsets(dcd);
    printf("%d\n", tmax);

    for (i=start; i<tmax && i<10; i+=freq) {
	if(getNextXYZ(dcd_xyz, dcd, freq))
	    break;
	cutXYZ(tar_xyz, dcd_xyz, index, n);
	for(j=0; j<n; j++){
	    printf("%d %d %f %f %f\n", i, index[j], tar_xyz[3*j], tar_xyz[3*j+1], tar_xyz[3*j+2]);
	}
    }
    closeDcdhandle(dcd);

    return 0;
}
#endif
