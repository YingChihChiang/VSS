//# define nthreads 4

// pmepot_data_struct, from pmepot 
struct pmepot_data_struct {
  int dims[5];
  int grid_size;
  int max_dim;
  int fft_ntable;
  double ewald_factor;
  double oddd[12];
  int avg_count;
  double *avg_potential;
  double *fft_table;
  double *fft_work;
  double *q_arr;
};
typedef struct pmepot_data_struct pmepot_data;


// VSS data structures
struct VDW {
  int fftype;
  double A, B, Ron, Roff, Ron2, Roff2, swdenom, mycoeff;
  double *para;
};

struct EXC {
  int excnum3, excnum4;
  int *excdisp3, *exclist3, *excdisp4, *exclist4;
};

struct MUT {
  int vecnum;
  int *atomnum, *veclist, *atomlist;
};

struct BOND {
  int nentry;
  int *list;
  double *para;
};

struct SAM {
  int nentry, scanflag;
  int *nsam, *naff, *list, *afflist;
  double *samangle;
};

struct PAT {
  int flag_spacedecomp, npatch, nxy, nbin[3], neighbor[27];
  int *workcount, *aid2pid;
  double patchsize[3], boxlend[3], boxrend[3], neighbor_nr[81];
  double *tempxyzq;
};

struct PSE {
  int natoms, *aid;
  int *patchnatoms, *patchdisp, *patch2aid;
  double *patchxyzq, *patchvdwp;
};

struct GEN {
  int natoms, nrotate, jobcase;
  double RT, ewaldfactor, gridspacing, Ebmin;
  double *charge, *sel2pdb, *Eb, *DeltaE;
}; 

struct DCDS {
  int *index1, *index3, tmax, fint, fstep, fend;
  char *dcdfile, *xstfile;
  double *dcd_xyz, *xyz1, *xyz2, *xyz3;
  void *dcd;
  FILE *xst;
};

// pmepot functions header
pmepot_data* pmepot_create(int *dims, double ewald_factor);
void pmepot_destroy(pmepot_data *data);
int fill_charges(const int *dims, const double *cell, int natoms,const double *xyzq, double *q_arr, double *rcell, double *oddd);
double compute_energy(double *q_arr, const double *cell, const double *rcell,const int *dims, double ewald);
double veclength(double vec[]);
void make_grid_from_cell(int nr[], double cell[], double resolution);
int good_fft_dim(double r);


// VSS functions header
int vssmain(struct GEN *gen_data, struct DCDS *dcd_data, struct VDW *vdw_data, struct PSE *sel0, struct PSE *sel1, struct PSE *sel2, 
            struct EXC *exc_data, struct MUT *mut_data, struct BOND *Dihed, struct BOND *Angle, struct SAM *SamDihed, struct SAM *SamAngle);
void vsscase1(const double *cell, double *atoms, struct GEN *gen_data, pmepot_data *data, struct VDW *vdw_data, 
           struct PSE *sel0, struct PSE *sel1, struct PSE *sel2, struct EXC *exc_data, struct MUT *mut_data,
           struct PAT *patch_data, struct BOND *Dihed, struct BOND *Angle, struct SAM *SamDihed, struct SAM *SamAngle);
void vsscase2(const double *cell, double *atoms, struct GEN *gen_data, pmepot_data *data, struct VDW *vdw_data,
           struct PSE *sel0, struct PSE *sel1, struct PSE *sel2, struct EXC *exc_data, struct MUT *mut_data,
           struct PAT *patch_data, struct BOND *Dihed, struct BOND *Angle, struct SAM *SamDihed, struct SAM *SamAngle);
void vsscase3(const double *cell, double *atoms, struct GEN *gen_data, pmepot_data *data, struct VDW *vdw_data,
           struct PSE *sel0, struct PSE *sel1, struct PSE *sel2, struct EXC *exc_data, struct MUT *mut_data,
           struct PAT *patch_data, struct BOND *Dihed, struct BOND *Angle, struct SAM *SamDihed, struct SAM *SamAngle);
void vsscase4(const double *cell, double *atoms, struct GEN *gen_data, pmepot_data *data, struct VDW *vdw_data,
           struct PSE *sel0, struct PSE *sel1, struct PSE *sel2, struct EXC *exc_data, struct MUT *mut_data,
           struct PAT *patch_data, struct BOND *Dihed, struct BOND *Angle, struct SAM *SamDihed, struct SAM *SamAngle);
void general_init(int natoms, double *xyz2, double *atoms, struct GEN *gen_data, struct SAM *SamDihed, struct SAM *SamAngle);
void general_free(struct GEN *gen_data);
void calcdG (int nframe, double *acc_exp_dE, double *dG, struct GEN *gen_data);

// Van der Waals
void vdw_data_init(double ewald_factor, struct VDW *vdw_data);
void vdwAB_init(int i, int j, struct VDW *vdw_data, double *A, double *B);
void sel_vdwAB_init(int i, int j, struct VDW *vdw_data, struct EXC *exc_data, double *A, double *B, double *excpair, double *excscale);

// Energy
void compute_intersel_vdwpmer(pmepot_data *data, const double *cell, int natoms, double *atoms, double *box_disp, struct VDW *vdw_data, struct EXC *exc_data,
                              struct PSE *sel, struct PSE *sel0, double *E_vdw, double *E_pme_r);
double compute_pmek(pmepot_data *data, const double *cell,int natoms, double *atoms, double *charge, struct PSE *othersel);
void compute_sel_vdwpmer(pmepot_data *data, const double *cell, int natoms, double *atoms, struct VDW *vdw_data, struct EXC *exc_data, struct PSE *sel, double *E_vdw, double *E_pme_r);
void compute_sel_nonb(pmepot_data *data, const double *cell, int natoms, double *atoms, struct VDW *vdw_data, struct EXC *exc_data, struct PSE *sel, double *E_vdw, double *E_elec);
double compute_pmeselfdiff(pmepot_data *data, double *atoms, struct PSE *sel1, struct PSE *sel2);
double get_r2(double const *xyzq, int i, int j, double box_disp_x, double box_disp_y, double box_disp_z);
void get_box_disp(double const *cell, double *box_disp);

// Mutation 
void gen_sel2_geo(int nsel1, int nsel2, double *xyzq, double *sel2pdb, struct MUT *mut_data);
void rot_sel2_dihed(int nsel2, int dihed_id, double sam_angle, double *xyzq, struct SAM *SamDihed);
void rot_sel2_angle(int nsel2, int angle_id, double sam_angle, double *xyzq, struct SAM *SamAngle);

// Bonded Energy, Dihedral, Angle
double calc_Eb(double *xyzq, struct BOND *Dihed, struct BOND *Angle);
double calc_dihed_angle(double *xyzq, int i, int j, int k, int l);
double calc_angle(double *xyzq, int i, int j, int k);
double calc_UB(double *xyzq, int i, int k);
void scandihed(int dihed_id, double RT, double *atoms, struct PSE *sel1, struct PSE *sel2, struct MUT *mut_data, struct BOND *Dihed, struct BOND *Angle, struct SAM *SamDihed);
void scanangle(int angle_id, double RT, double *atoms, struct PSE *sel1, struct PSE *sel2, struct MUT *mut_data, struct BOND *Dihed, struct BOND *Angle, struct SAM *SamAngle);

// Patch
void patch_init(int natoms, double Roff, double const *cell, struct PSE *sel0, struct PSE *sel1, struct PSE *sel2, struct PAT *patch_data);
void box_update(double const *cell, struct PAT *patch_data);
void binatoms(double const *cell, double *xyzq, struct VDW *vdw_data, struct PAT *patch_data, struct PSE *sel);
void getneighbor(int pid, const double *cell, struct PAT *patch_data);
void compute_patch_vdwpmer(const double *cell, double ewald_factor, struct VDW *vdw_data, struct PAT *patch_data, struct PSE *sel, struct PSE *sel0, double *E_vdw, double *E_pmer);
void patch_vdwAB_init(double epsilon_i, double epsilon_j, double R_i, double R_j, double *A, double *B, int fftype);
double patch_getr2(double *ixyz, double *jxyz, double *nr);
void patch_free(struct PAT *patch_data, struct PSE *sel0, struct PSE *sel1, struct PSE *sel2);

// Array operation
int sumint(int *arr, int ndim);
double sumdouble(double *arr, int ndim);
double norm(double *arr, int ndim);
void cross(double *vecA, double *vecB, double *vecC);
double innerproduct(double *vecA, double *vecB, int ndim);
void vec2zaxis(double *vec, double sintheta, double costheta, double sinphi, double cosphi);
void zaxis2vec(double *vec, double sintheta, double costheta, double sinphi, double cosphi);
void rotate_around_zaxis(double *vec, double sinangle, double cosangle);
void rotate_around_xaxis(double *vec, double sinangle, double cosangle);
void vec2xaxis(double *vec, double costheta, double sintheta, double cosphi, double sinphi);
void xaxis2vec(double *vec, double costheta, double sintheta, double cosphi, double sinphi);

// ReadDcd
void* initDcdhandle(char dcdfile[],  int firststep);
int getNextXYZ(double xyz[], void *v, int fstep);
void cutXYZ(double dest[], double src[], int index[], int n);
void cutXYZQ(double xyzq[], double xyz[], double q[], int index[], int n);
void mergeXYZQ(double xyzq[], double xyz[], double q[], int n);
void updateXYZofAtom(double xyzq[], double xyz[], int n);
void closeDcdhandle(void *v);
int getNsets(void *v);
int getNatoms(void *v);

