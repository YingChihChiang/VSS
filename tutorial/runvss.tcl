# Tcl script for running vss, Letitia 22/09/2014, 25/01/2016
##############################################################################
# VSS usage:
# 4 possible cases: 
# Both Dihedral and Angle sampling: vss -sel1 $sel1 -file1 $file1 -file2 $file2 -parfile $parfile -mutvec $mutvec -samdihedlist $samdihedlist -samdihedval $samdihedval -samanglelist $samanglelist -samangleval $samangleval
# Only Dihedral sampling:           vss -sel1 $sel1 -file1 $file1 -file2 $file2 -parfile $parfile -mutvec $mutvec -samdihedlist $samdihedlist -samdihedval $samdihedval
# Only Angle sampling:              vss -sel1 $sel1 -file1 $file1 -file2 $file2 -parfile $parfile -mutvec $mutvec -samanglelist $samanglelist -samangleval $samangleval
# No Dihedral nor Angle sampling:   vss -sel1 $sel1 -file1 $file1 -file2 $file2 -parfile $parfile -mutvec $mutvec
#
# Only sel1, file1, sel2, file2, parfile, mutvec, are mandatory input arguments.
# Additionally assign initial findex (fint), step of findex (fstep), end of findex (fend), temperature (temp)  #E.g. -fint 0 -fstep 20 -fend 100 -temp 300.0
##############################################################################

# Set the parameters
set file1 "dcd/wb_benzene";     # lead compound's dcd, xst, psf (no extension)
set file2 "min_toluene";        # derivatives' pdb, psf (no extension)
set parfile "parfile.prm";      # used parameters (Charmm, CGenFF only)
set sel1 "resname BEN";         # ligand in dcd trajectory
set mutvec "0 1 0 1";           # site of mutation
set samdihedlist "2 0 1 12 ;";  # index of the dihedral angle to be sampled
set samdihedval "0 30 60 90 ;"; # value of the dihedral angle to be sampled
set samanglelist "0 1 12 ;";    # index of the angle to be sampled
set samangleval  "108.57 ;";    # value of the angle to be sampled

# Load the VSS plugin
source ./vss.tcl

# Run your VSS job
vss -fstep 1 -temp 300.0 -sel1 $sel1 -file1 $file1 -file2 $file2 -parfile $parfile -mutvec $mutvec -samdihedlist $samdihedlist -samdihedval $samdihedval

exit
