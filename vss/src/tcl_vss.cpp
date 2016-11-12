/*
 * Tcl bindings for the Virtual Substitution Scan plugin
 *
 * $Id: tcl_vss.c,v 1.4 2005/07/20 15:37:40 johns Exp $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include "vss.h"

#include <tcl.h>

#define BUFSIZE 256

int tcl_vssrun(ClientData nodata, Tcl_Interp *interp,
			int objc, Tcl_Obj *const objv[]) {

  int i, j, itemp, idisp, natom, mycount, errcode=1;
  int data_count, sub_count, sub2_count;
  Tcl_Obj **data_list, **sub_list, **sub2_list;
  double *atoms;
  double d, RT, ewaldfactor, gridspacing;
  struct VDW vdw_data;
  struct PSE sel0, sel1, sel2;
  struct EXC exc_data;
  struct MUT mut_data;
  struct BOND Dihed, Angle;
  struct SAM SamDihed, SamAngle;
  struct GEN gen_data;

  //for dcd reading
  Tcl_Obj **dcd_list;
  //int rc, 
  int dcd_count, dcdfile_count, xstfile_count;
  struct DCDS dcd_data;
  

  if ( objc != 7 ) { Tcl_SetResult(interp,"args: gdata idata rdata samdata xyzq2 dcd",TCL_VOLATILE); return TCL_ERROR; }

  // GData  (Ron,Roff,RT,Ewaldfactor,gridspacing,nsel1,nsel2,nsel0)
  if ( Tcl_ListObjGetElements(interp,objv[1],&data_count,&data_list) != TCL_OK ) return TCL_ERROR;
  if ( Tcl_GetDoubleFromObj(interp,data_list[0],&vdw_data.Ron) != TCL_OK ) return TCL_ERROR;
  if ( Tcl_GetDoubleFromObj(interp,data_list[1],&vdw_data.Roff) != TCL_OK ) return TCL_ERROR;
  if ( Tcl_GetDoubleFromObj(interp,data_list[2],&RT) != TCL_OK ) return TCL_ERROR;
  if ( Tcl_GetDoubleFromObj(interp,data_list[3],&ewaldfactor) != TCL_OK ) return TCL_ERROR;
  if ( Tcl_GetDoubleFromObj(interp,data_list[4],&gridspacing) != TCL_OK ) return TCL_ERROR;
  if ( Tcl_GetIntFromObj(interp,data_list[5],&sel1.natoms) != TCL_OK ) return TCL_ERROR;
  if ( Tcl_GetIntFromObj(interp,data_list[6],&sel2.natoms) != TCL_OK ) return TCL_ERROR;
  if ( Tcl_GetIntFromObj(interp,data_list[7],&sel0.natoms) != TCL_OK ) return TCL_ERROR;

  natom=sel1.natoms+sel2.natoms+sel0.natoms;
  atoms = (double*) malloc(natom*4*sizeof(double));
  dcd_data.index1 = (int*)malloc(sel1.natoms*sizeof(int));
  dcd_data.index3 = (int*) malloc(sel0.natoms*sizeof(int));
  dcd_data.xyz1 = (double*) malloc(sel1.natoms*3*sizeof(double));
  dcd_data.xyz2 = (double*) malloc(sel2.natoms*3*sizeof(double));
  dcd_data.xyz3 = (double*) malloc(sel0.natoms*3*sizeof(double));
  gen_data.charge = (double*) malloc(natom*sizeof(double));

  // IData  ({excdisp3,exclist3}, {excdisp4,exclist4}, {sel1list,sel2list,sel0list}, {mutvec list}, {mutatom list}, {dihed list}, {angle list}, {sel1 index}, {sel0 index})
  if ( Tcl_ListObjGetElements(interp,objv[2],&data_count,&data_list) != TCL_OK ) return TCL_ERROR;
  if ( data_count != 9) {
    Tcl_SetResult(interp,"IData format: {excdisp3,exclist3}, {excdisp4,exclist4}, {sel1list,sel2list,sel0list}, {mutvec list}, {mutatom list}, {dihed list},  {angle list}, {sel1 index}, {sel0 index}",TCL_VOLATILE);
    free(atoms);
    return TCL_ERROR;
  }

  // IData: excluded list 3 (excdisp: position where excluded atoms' indices storaged in exclist.)
  if ( Tcl_ListObjGetElements(interp,data_list[0],&sub_count,&sub_list) != TCL_OK ) return TCL_ERROR;
  int seldisp=sel1.natoms+sel2.natoms+1;
  exc_data.excdisp3 = (int*) malloc(seldisp*sizeof(int));
  for ( i=0; i<seldisp; ++i) {
    if ( Tcl_GetIntFromObj(interp,sub_list[i],&exc_data.excdisp3[i]) != TCL_OK ) return TCL_ERROR;
  }
  exc_data.excnum3=exc_data.excdisp3[seldisp-1];
  exc_data.exclist3= (int*) malloc(exc_data.excnum3*sizeof(int));
  idisp = seldisp;
  for (i=0; i<exc_data.excnum3; ++i) {
    if ( Tcl_GetIntFromObj(interp,sub_list[idisp+i],&exc_data.exclist3[i]) != TCL_OK ) return TCL_ERROR;
  }

  // IData: excluded list 4 (excdisp: position where excluded atoms' indices storaged in exclist.)
  if ( Tcl_ListObjGetElements(interp,data_list[1],&sub_count,&sub_list) != TCL_OK ) return TCL_ERROR;
  exc_data.excdisp4 = (int*) malloc(seldisp*sizeof(int));
  for ( i=0; i<seldisp; ++i) {
    if ( Tcl_GetIntFromObj(interp,sub_list[i],&exc_data.excdisp4[i]) != TCL_OK ) return TCL_ERROR;
  }
  exc_data.excnum4=exc_data.excdisp4[seldisp-1];
  exc_data.exclist4= (int*) malloc(exc_data.excnum4*sizeof(int));
  idisp = seldisp;
  for (i=0; i<exc_data.excnum4; ++i) {
    if ( Tcl_GetIntFromObj(interp,sub_list[idisp+i],&exc_data.exclist4[i]) != TCL_OK ) return TCL_ERROR;
  } 

  // IData: sel1 list
  if ( Tcl_ListObjGetElements(interp,data_list[2],&sub_count,&sub_list) != TCL_OK ) return TCL_ERROR;
  sel1.aid = (int*) malloc(sel1.natoms*sizeof(int));
  for (i=0; i<sel1.natoms; ++i) {
    if ( Tcl_GetIntFromObj(interp,sub_list[i],&sel1.aid[i]) != TCL_OK ) return TCL_ERROR;
  }

  // IData: sel2 list
  sel2.aid = (int*) malloc(sel2.natoms*sizeof(int));
  idisp = sel1.natoms;
  for (i=0; i<sel2.natoms; ++i) {
    if ( Tcl_GetIntFromObj(interp,sub_list[idisp+i],&sel2.aid[i]) != TCL_OK ) return TCL_ERROR;
  }
  
  // IData: sel0 list
  sel0.aid = (int*) malloc(sel0.natoms*sizeof(int));
  idisp = sel1.natoms+sel2.natoms;
  for (i=0; i<sel0.natoms; ++i) {
    if ( Tcl_GetIntFromObj(interp,sub_list[idisp+i],&sel0.aid[i]) != TCL_OK ) return TCL_ERROR;
  }
 
  // IData: mutvec list
  if ( Tcl_ListObjGetElements(interp,data_list[3],&sub_count,&sub_list) != TCL_OK ) return TCL_ERROR;
  mut_data.vecnum = (sub_count/4);
  mut_data.veclist = (int*) malloc(sub_count*sizeof(int));
  for (i=0; i<sub_count;++i) {
     if ( Tcl_GetIntFromObj(interp,sub_list[i],&mut_data.veclist[i]) != TCL_OK ) return TCL_ERROR;
  }

  // IData: mutatom list
  if ( Tcl_ListObjGetElements(interp,data_list[4],&sub_count,&sub_list) != TCL_OK ) return TCL_ERROR;
  itemp = sub_count-mut_data.vecnum;
  idisp = mut_data.vecnum;
  mut_data.atomnum = (int*) malloc(mut_data.vecnum*sizeof(int));
  mut_data.atomlist = (int*) malloc(itemp*sizeof(int));
  for (i=0; i<mut_data.vecnum;++i) {
     if ( Tcl_GetIntFromObj(interp,sub_list[i],&mut_data.atomnum[i]) != TCL_OK ) return TCL_ERROR;
  }
  for (i=0; i<itemp;++i) {
     if ( Tcl_GetIntFromObj(interp,sub_list[idisp+i],&mut_data.atomlist[i]) != TCL_OK ) return TCL_ERROR;
  }

  // IData: dihedral list
  if ( Tcl_ListObjGetElements(interp,data_list[5],&sub_count,&sub_list) != TCL_OK ) return TCL_ERROR;
  Dihed.nentry = sub_count;
  Dihed.list = (int*) malloc(4*sub_count*sizeof(int));
  Dihed.para = (double*) malloc(3*sub_count*sizeof(double));
  for (i=0; i<sub_count;++i) {
     if ( Tcl_ListObjGetElements(interp,sub_list[i],&sub2_count,&sub2_list) != TCL_OK ) return TCL_ERROR;
     for ( j=0; j<4; ++j ) {
        if ( Tcl_GetIntFromObj(interp,sub2_list[j],&Dihed.list[4*i+j]) != TCL_OK )  return TCL_ERROR;
     }
  }

  // IData: angle list
  if ( Tcl_ListObjGetElements(interp,data_list[6],&sub_count,&sub_list) != TCL_OK ) return TCL_ERROR;
  Angle.nentry = sub_count;
  Angle.list = (int*) malloc(3*sub_count*sizeof(int));
  Angle.para = (double*) malloc(4*sub_count*sizeof(double));
  for (i=0; i<sub_count;++i) {
     if ( Tcl_ListObjGetElements(interp,sub_list[i],&sub2_count,&sub2_list) != TCL_OK ) return TCL_ERROR;
     for ( j=0; j<3; ++j ) {
        if ( Tcl_GetIntFromObj(interp,sub2_list[j],&Angle.list[3*i+j]) != TCL_OK )  return TCL_ERROR;
     }
  }

  // IData: sel1 index and sel0 index in original dcd file.
  if ( Tcl_ListObjGetElements(interp,data_list[7],&sub_count,&sub_list) != TCL_OK ) return TCL_ERROR;
  if (sub_count != sel1.natoms) { Tcl_SetResult(interp,"Number of sel1 atoms and number of sel1 index do not match.",TCL_VOLATILE); return TCL_ERROR; }
  for (i=0; i<sub_count;++i) {
     if ( Tcl_GetIntFromObj(interp,sub_list[i],&dcd_data.index1[i]) != TCL_OK )  return TCL_ERROR;
  }
  if ( Tcl_ListObjGetElements(interp,data_list[8],&sub_count,&sub_list) != TCL_OK ) return TCL_ERROR;
  if (sub_count != sel0.natoms) { Tcl_SetResult(interp,"Number of sel0 atoms and number of sel0 index do not match.",TCL_VOLATILE); return TCL_ERROR; }
  for (i=0; i<sub_count;++i) {
     if ( Tcl_GetIntFromObj(interp,sub_list[i],&dcd_data.index3[i]) != TCL_OK )  return TCL_ERROR;
  }


  // RData ({eps,Rmin/2,epsi_1-4,Rmin/2_1-4}, {k,n,delta}, {k,theda0,k_ub,theta_ub},{sel1 q},{sel0 q})
  if ( Tcl_ListObjGetElements(interp,objv[3],&data_count,&data_list) != TCL_OK ) return TCL_ERROR;
  vdw_data.para = (double*) malloc(natom*4*sizeof(double));
  for (i=0; i<natom; ++i) {
    if ( Tcl_ListObjGetElements(interp,data_list[i],&sub_count,&sub_list) != TCL_OK ) return TCL_ERROR;
    for ( j=0; j<4; ++j ) {
      if ( Tcl_GetDoubleFromObj(interp,sub_list[j],&vdw_data.para[4*i+j]) != TCL_OK ) return TCL_ERROR;
    }
  }
  if ( Tcl_ListObjGetElements(interp,data_list[natom],&sub_count,&sub_list) != TCL_OK ) return TCL_ERROR;
  if (sub_count != Dihed.nentry) { Tcl_SetResult(interp,"Number of dihedral list and Number of dihedral parameter does not match.",TCL_VOLATILE); return TCL_ERROR; }
  for (i=0; i<Dihed.nentry;++i) {
     if ( Tcl_ListObjGetElements(interp,sub_list[i],&sub2_count,&sub2_list) != TCL_OK ) return TCL_ERROR;
     for ( j=0; j<3; ++j ) {
        if ( Tcl_GetDoubleFromObj(interp,sub2_list[j],&Dihed.para[3*i+j]) != TCL_OK )  return TCL_ERROR;
     }
  }
  if ( Tcl_ListObjGetElements(interp,data_list[natom+1],&sub_count,&sub_list) != TCL_OK ) return TCL_ERROR;
  if (sub_count != Angle.nentry) { Tcl_SetResult(interp,"Number of angle list and Number of angle parameter does not match.",TCL_VOLATILE); return TCL_ERROR; }
  for (i=0; i<Angle.nentry;++i) {
     if ( Tcl_ListObjGetElements(interp,sub_list[i],&sub2_count,&sub2_list) != TCL_OK ) return TCL_ERROR;
     for ( j=0; j<4; ++j ) {
        if ( Tcl_GetDoubleFromObj(interp,sub2_list[j],&Angle.para[4*i+j]) != TCL_OK )  return TCL_ERROR;
     }
  }
  if ( Tcl_ListObjGetElements(interp,data_list[natom+2],&sub_count,&sub_list) != TCL_OK ) return TCL_ERROR;
  if (sub_count != sel1.natoms) { Tcl_SetResult(interp,"Number of sel1 atoms and number of sel1 charges do not match.",TCL_VOLATILE); return TCL_ERROR; }
  for (i=0; i<sub_count;++i) {
     if ( Tcl_GetDoubleFromObj(interp,sub_list[i],&gen_data.charge[i]) != TCL_OK )  return TCL_ERROR;
  }
  if ( Tcl_ListObjGetElements(interp,data_list[natom+3],&sub_count,&sub_list) != TCL_OK ) return TCL_ERROR;
  if (sub_count != sel0.natoms) { Tcl_SetResult(interp,"Number of sel0 atoms and number of sel0 charges do not match.",TCL_VOLATILE); return TCL_ERROR; }
  for (i=0; i<sub_count;++i) {
     idisp=i+sel1.natoms+sel2.natoms; 
     if ( Tcl_GetDoubleFromObj(interp,sub_list[i],&gen_data.charge[idisp]) != TCL_OK )  return TCL_ERROR;
  }

  // SamData ({SamDihed},{SamAngle})
  if ( Tcl_ListObjGetElements(interp,objv[4],&data_count,&data_list) != TCL_OK ) return TCL_ERROR;
  // Reading SamDihed ({nentry, scanflag, list, nsam, samangle, naff, afflist})
  if ( Tcl_GetIntFromObj(interp,data_list[0],&SamDihed.nentry) != TCL_OK )  return TCL_ERROR;
  if ( Tcl_GetIntFromObj(interp,data_list[1],&SamDihed.scanflag) != TCL_OK )  return TCL_ERROR;
  SamDihed.list = (int*) malloc(4*SamDihed.nentry*sizeof(int)); 
  SamDihed.nsam = (int*) malloc(SamDihed.nentry*sizeof(int)); 
  SamDihed.naff = (int*) malloc(SamDihed.nentry*sizeof(int)); 
  if ( Tcl_ListObjGetElements(interp,data_list[2],&sub_count,&sub_list) != TCL_OK ) return TCL_ERROR;
  for (i=0;i<SamDihed.nentry;i++) {
     for ( j=0; j<4; ++j ) {
        if ( Tcl_GetIntFromObj(interp,sub_list[i*4+j],&SamDihed.list[i*4+j]) != TCL_OK )  return TCL_ERROR;
     }
  }
  if ( Tcl_ListObjGetElements(interp,data_list[3],&sub_count,&sub_list) != TCL_OK ) return TCL_ERROR;
  for (i=0;i<SamDihed.nentry;i++) {
     if ( Tcl_GetIntFromObj(interp,sub_list[i],&SamDihed.nsam[i]) != TCL_OK )  return TCL_ERROR;
  }
  SamDihed.samangle = (double*) malloc(sumint(SamDihed.nsam,SamDihed.nentry)*sizeof(double));
  if ( Tcl_ListObjGetElements(interp,data_list[4],&sub_count,&sub_list) != TCL_OK ) return TCL_ERROR;
  itemp=sumint(SamDihed.nsam,SamDihed.nentry);
  for (i=0;i<itemp;i++) {
     if ( Tcl_GetDoubleFromObj(interp,sub_list[i],&SamDihed.samangle[i]) != TCL_OK )  return TCL_ERROR;
  }
  if ( Tcl_ListObjGetElements(interp,data_list[5],&sub_count,&sub_list) != TCL_OK ) return TCL_ERROR;
  for (i=0;i<SamDihed.nentry;i++) {
     if ( Tcl_GetIntFromObj(interp,sub_list[i],&SamDihed.naff[i]) != TCL_OK )  return TCL_ERROR;
  }
  SamDihed.afflist = (int*) malloc(sumint(SamDihed.naff,SamDihed.nentry)*sizeof(int));
  if ( Tcl_ListObjGetElements(interp,data_list[6],&sub_count,&sub_list) != TCL_OK ) return TCL_ERROR;
  itemp=sumint(SamDihed.naff,SamDihed.nentry);
  for (i=0;i<itemp;i++) {
     if ( Tcl_GetIntFromObj(interp,sub_list[i],&SamDihed.afflist[i]) != TCL_OK )  return TCL_ERROR;
  }
  // Reading SamAngle ({nentry, scanflag, list, nsam, samangle, naff, afflist})
  if ( Tcl_GetIntFromObj(interp,data_list[7],&SamAngle.nentry) != TCL_OK )  return TCL_ERROR;
  if ( Tcl_GetIntFromObj(interp,data_list[8],&SamAngle.scanflag) != TCL_OK )  return TCL_ERROR;
  SamAngle.list = (int*) malloc(3*SamAngle.nentry*sizeof(int));
  SamAngle.nsam = (int*) malloc(SamAngle.nentry*sizeof(int));
  SamAngle.naff = (int*) malloc(SamAngle.nentry*sizeof(int));
  if ( Tcl_ListObjGetElements(interp,data_list[9],&sub_count,&sub_list) != TCL_OK ) return TCL_ERROR;
  for (i=0;i<SamAngle.nentry;i++) {
     for ( j=0; j<3; ++j ) {
        if ( Tcl_GetIntFromObj(interp,sub_list[i*3+j],&SamAngle.list[i*3+j]) != TCL_OK )  return TCL_ERROR;
     }
  }
  if ( Tcl_ListObjGetElements(interp,data_list[10],&sub_count,&sub_list) != TCL_OK ) return TCL_ERROR;
  for (i=0;i<SamAngle.nentry;i++) {
     if ( Tcl_GetIntFromObj(interp,sub_list[i],&SamAngle.nsam[i]) != TCL_OK )  return TCL_ERROR;
  }
  SamAngle.samangle = (double*) malloc(sumint(SamAngle.nsam,SamAngle.nentry)*sizeof(double));
  if ( Tcl_ListObjGetElements(interp,data_list[11],&sub_count,&sub_list) != TCL_OK ) return TCL_ERROR;
  itemp=sumint(SamAngle.nsam,SamAngle.nentry);
  for (i=0;i<itemp;i++) {
     if ( Tcl_GetDoubleFromObj(interp,sub_list[i],&SamAngle.samangle[i]) != TCL_OK )  return TCL_ERROR;
  }
  if ( Tcl_ListObjGetElements(interp,data_list[12],&sub_count,&sub_list) != TCL_OK ) return TCL_ERROR;
  for (i=0;i<SamAngle.nentry;i++) {
     if ( Tcl_GetIntFromObj(interp,sub_list[i],&SamAngle.naff[i]) != TCL_OK )  return TCL_ERROR;
  }
  SamAngle.afflist = (int*) malloc(sumint(SamAngle.naff,SamAngle.nentry)*sizeof(int));
  if ( Tcl_ListObjGetElements(interp,data_list[13],&sub_count,&sub_list) != TCL_OK ) return TCL_ERROR;
  itemp=sumint(SamAngle.naff,SamAngle.nentry);
  for (i=0;i<itemp;i++) {
     if ( Tcl_GetIntFromObj(interp,sub_list[i],&SamAngle.afflist[i]) != TCL_OK )  return TCL_ERROR;
  }


  // xyzq2
  if ( Tcl_ListObjGetElements(interp,objv[5],&sub_count,&sub_list) != TCL_OK ) return TCL_ERROR;
  if (sub_count != sel2.natoms) { Tcl_SetResult(interp,"Number of sel2 atoms and number of sel2 xyzq do not match.",TCL_VOLATILE); return TCL_ERROR; }
  mycount = 4*sel1.natoms;
  for ( i=0; i<sub_count; ++i ) {
      if ( Tcl_ListObjGetElements(interp,sub_list[i],&sub2_count,&sub2_list) != TCL_OK ) { free(atoms); return TCL_ERROR; }
      if ( sub2_count != 4 ) { Tcl_SetResult(interp,"atom2 format: {{x y z q}...}",TCL_VOLATILE); free(atoms); return TCL_ERROR;}
      for ( j=0; j<3; ++j ) {
          if ( Tcl_GetDoubleFromObj(interp,sub2_list[j],&d) != TCL_OK ) { free(atoms); return TCL_ERROR; }
          dcd_data.xyz2[3*i+j] = d;
          atoms[mycount+4*i+j] = d;
      }
      if ( Tcl_GetDoubleFromObj(interp,sub2_list[3],&gen_data.charge[sel1.natoms+i]) != TCL_OK ) return TCL_ERROR;
  }


  // dcd
  if ( Tcl_ListObjGetElements(interp,objv[6],&dcd_count,&dcd_list) != TCL_OK ) return TCL_ERROR;
  if ( dcd_count != 5 ) { Tcl_SetResult(interp,"dcd format: {dcdfile xstfile fint fstep fend}",TCL_VOLATILE); return TCL_ERROR; }
  i=0;
  if ( !(dcd_data.dcdfile = Tcl_GetStringFromObj(dcd_list[i++], &dcdfile_count)) ) return TCL_ERROR;
  if ( !(dcd_data.xstfile = Tcl_GetStringFromObj(dcd_list[i++], &xstfile_count)) ) return TCL_ERROR;
  if ( Tcl_GetIntFromObj(interp,dcd_list[i++],&dcd_data.fint ) != TCL_OK ) return TCL_ERROR;
  if ( Tcl_GetIntFromObj(interp,dcd_list[i++],&dcd_data.fstep) != TCL_OK ) return TCL_ERROR;
  if ( Tcl_GetIntFromObj(interp,dcd_list[i++],&dcd_data.fend) != TCL_OK ) return TCL_ERROR;

  gen_data.natoms = natom;
  gen_data.RT = RT;
  gen_data.ewaldfactor = ewaldfactor;
  gen_data.gridspacing = gridspacing;
 
  // Run VSS
  errcode=vssmain(&gen_data,&dcd_data,&vdw_data,&sel0,&sel1,&sel2,&exc_data,&mut_data,&Dihed,&Angle,&SamDihed,&SamAngle);

  if ( errcode ) {
    Tcl_SetResult(interp,"Error: vss failed.",TCL_VOLATILE);
    return TCL_ERROR;
  }

  general_free(&gen_data);
  free(dcd_data.index1);
  free(dcd_data.index3);
  free(dcd_data.xyz1);
  free(dcd_data.xyz2);
  free(dcd_data.xyz3);
  free(exc_data.excdisp3);
  free(exc_data.excdisp4);
  free(exc_data.exclist3);
  free(exc_data.exclist4);
  free(sel1.aid);
  free(sel2.aid);
  free(sel0.aid);
  free(vdw_data.para);
  free(mut_data.veclist);
  free(mut_data.atomnum);
  free(mut_data.atomlist);
  free(Dihed.list);
  free(Dihed.para);
  free(Angle.list);
  free(Angle.para);
  free(SamDihed.list);
  free(SamDihed.nsam);
  free(SamDihed.samangle);
  free(SamDihed.naff);
  free(SamDihed.afflist);
  free(SamAngle.list);
  free(SamAngle.nsam);
  free(SamAngle.samangle);
  free(SamAngle.naff);
  free(SamAngle.afflist);
  return TCL_OK;

}


#if defined(VSSTCLDLL_EXPORTS) && defined(_WIN32)
#  undef TCL_STORAGE_CLASS
#  define TCL_STORAGE_CLASS DLLEXPORT

#define WIN32_LEAN_AND_MEAN // Exclude rarely-used stuff from Window s headers
#include <windows.h>

BOOL APIENTRY DllMain( HANDLE hModule, 
                       DWORD  ul_reason_for_call, 
                       LPVOID lpReserved
                                         )
{
    return TRUE;
}

EXTERN int Vss_Init(Tcl_Interp *interp) {

#else

extern "C" { 
  int Vss_Init(Tcl_Interp *interp) {

#endif

    Tcl_CreateObjCommand(interp,"vssrun",tcl_vssrun,
	(ClientData)NULL, (Tcl_CmdDeleteProc*)NULL);

 
    Tcl_PkgProvide(interp, "vss_core", "1.0.0");

    return TCL_OK;
  }
}
