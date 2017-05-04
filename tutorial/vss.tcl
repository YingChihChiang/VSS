#
# Tcl/Tk front-end for Virtual Scan Substitution
#
# by YingChih Chiang (Letitia) 2014.02.27;
#
package provide vss 1.0.0

namespace eval ::VSS {

  package require vss_core 1.0

  namespace export vss

  variable default_grid_resolution 1.0
  variable default_grid_max_dim 512
  variable default_cell_padding 10.0
  variable default_ewald_factor 0.25795183763637
  variable default_Ron 10.0
  variable default_Roff 12.0
  variable Boltzmann 0.001987191
  variable M_PI 3.14159265358979323846 
  variable RData
}

proc vss { args } { return [eval ::VSS::vss $args] }

proc ::VSS::vss { args } {
  set nargs [llength $args]
  if {$nargs == 0 || $nargs % 2} {
    puts "usage: vss ?-arg val?..."
    puts "  -file1 <filename> file name of the sampled system (file1.dcd, file1.xst, file1.psf)"
    puts "  -file2 <filename> file name of the mutated molecule (file2.pdb, file2.psf)"
    puts "  -parfile <filename> parameter file contains both systems' VDW parameters"
    puts "  -sel1 <selection> molecule which contains beta=-1"
    puts "  -mutvec <selection> indices of atoms on the site of mutation, from sel1 and sel2"
    puts "  -fint <fint> evaluate FEP from frame fint"
    puts "  -fstep <fstep> evaluate FEP every fstep frame"
    puts "  -fend <fend> evaluate FEP till frame fend"
    puts "  -temp <temperature>"
    puts "  -grid <spacing> grid spacing for PME in Angstroms"
    puts "  -ewaldfactor <factor> specify ewald factor"
    puts "  -cutoff <{Ron Roff}> Giving switch on and cutt off"
    puts "  -impscandihed <integer> sample the most relevant n dihedral angles"
    puts "  -samdihedlist <index1 index2 index3 index4; ...> atoms that form the dihedral angle to be sampled"
    puts "  -samdihedval <float1 float2 ...; ...> values of the dihedral angle to be sampled (0 to 360)"
    puts "  -impscanangle <integer> sample the most relevant n angles"
    puts "  -samanglelist <index1 index2 index3;...> atoms that form the angle to be sampled"
    puts "  -samangleval <float1 float2 ...; ...> values of the angle to be sampled (0 to 180)"
    error "error: empty argument list or odd number of arguments: $args"
  }

  foreach {name val} $args {
    switch -- $name {
      -file1 {set arg(file1) $val}
      -file2 {set arg(file2) $val}
      -parfile {set arg(parfile) $val}
      -sel1 {set arg(sel1) $val}
      -mutvec {set arg(mutvec) $val}
      -fint {set arg(fint) $val}
      -fstep {set arg(fstep) $val}
      -fend {set arg(fend) $val}
      -temp { set arg(temp) $val }
      -grid { set arg(grid) $val }
      -ewaldfactor { set arg(ewaldfactor) $val }
      -cutoff { set arg(cutoff) $val }
      -impscandihed { set arg(impscandihed) $val }
      -impscanangle { set arg(impscanangle) $val }
      -samdihedlist { set arg(samdihedlist) $val }
      -samdihedval {set arg(samdihedval) $val}
      -samanglelist { set arg(samanglelist) $val }
      -samangleval {set arg(samangleval) $val}
      default { error "unknown argument: $name $val" }
    }
  }


  # Declare global variables and vss-global variables
  global vmd_frame
  variable RData
  variable M_PI

  # Check the needed inputs are read
  if [info exists arg(file1)] {
     append psffile $arg(file1) ".psf"
     append dcdfile $arg(file1) ".dcd"
     append xstfile $arg(file1) ".xst"
  } else {
     error "Variable file1 is not specified."
  }

  if [info exists arg(file2)] {
     append psffile2 $arg(file2) ".psf"
     append pdbfile2 $arg(file2) ".pdb"
  } else {
     error "Variable file2 is not specified."
  }

  if [info exists arg(parfile)] {
     set parfile $arg(parfile)
  } else {
     error "Variable parfile is not specified."
  }

  if [expr ![info exists arg(sel1)]] {error "Variable sel1 (base structure) is not specified."}
  if [expr ![info exists arg(mutvec)]] {error "Variable mutvec (mutation vectors) is not specified."}
  if [expr ![info exists arg(temp)]] {error "Variable temp (temparature) is not specified."}

  # Check how many mutation sites are in the mutation
  if {[llength [lindex $arg(mutvec) 0]]==0} {error "No mutation! Delta G = 0"}
  set mylength [llength $arg(mutvec)]
  if {[expr $mylength % 4] != 0} {
     error "Size of mutvec = 4 * mutation sites"
  } else { 
     set num_mutsite [expr $mylength/4]
  }


  # Load files of mutated molecule 
  mol load psf $psffile2 pdb $pdbfile2
  set sel2 [atomselect top "all"]
  set mid2 [$sel2 molid] 


  # Load files of the sampled system (only the first frame of dcd)
  mol new $psffile
  mol addfile $dcdfile waitfor 0
  set sel [atomselect top "all"]
  set sel1 [atomselect top $arg(sel1)]
  set mid [$sel1 molid]
  set nameofsel1 $arg(sel1)
  set sel0 [atomselect $mid "all and not (same fragment as $nameofsel1)"]
  set dual [atomselect $mid "(same fragment as $nameofsel1) and not index [$sel1 get index]"]

  # Get total number of frames and time step
  if [info exists arg(fstep)] {
     set fstep $arg(fstep)
  } else {
     set fstep 1
  }
  if [info exists arg(fint)] {
     set fint $arg(fint)
  } else {
     set fint 0
  }
  if {[info exists arg(fend)]} {
    set fend $arg(fend)
  } else {
    set fend 0
  }


  # Prepare xst file for readin cell (The first three lines are not needed.)
  set xstfd [open $xstfile r]
  set line1 [gets $xstfd]
  set line2 [gets $xstfd]
  set line3 [gets $xstfd]


  # Get temparature from tempfile
  set temperature $arg(temp)


  # Get grid dimensions
  if [info exists arg(grid)] {
    set gridspacing $arg(grid)
  } else {
    variable default_grid_resolution
    set gridspacing $default_grid_resolution 
  }


  # Get ewald factor
  if [info exists arg(ewaldfactor)] {
    set ewaldfactor $arg(ewaldfactor)
  } else {
    variable default_ewald_factor
    set ewaldfactor $default_ewald_factor
  }

  # Get switch Ron, cutoff Roff, and RT
  if [info exists arg(cutoff)] {
    set Ron [lindex $arg(cutoff) 0]
    set Roff [lindex $arg(cutoff) 1]
    if { $Ron-$Roff > 0.0 } {error "Switch > Cutoff. Input Error!"}
  } else {
    variable default_Ron
    variable default_Roff
    set Ron $default_Ron
    set Roff $default_Roff
  }
  variable Boltzmann
  set RT [expr $temperature*$Boltzmann]
  

  # Set the GData (General data)
  set GData [list $Ron $Roff $RT $ewaldfactor $gridspacing [$sel1 num] [$sel2 num] [$sel0 num]]


  # Get VDW parameters (Rsults append in RData directly)
  if { [catch {package require readcharmmpar 1.2} errmsg] } then {
     error "Could not load the readcharmmpar package needed to read VDW parameters"
  } else {
     set parfile $parfile
     set Allpara [::VSS::read_charmm_parameters $parfile]
     set VDWlist [$sel1 get type] 
     getVDWdata $VDWlist $Allpara
     set VDWlist [$sel2 get type]
     getVDWdata $VDWlist $Allpara
     set VDWlist [$sel0 get type]
     getVDWdata $VDWlist $Allpara
     unset Allpara VDWlist
  }

  # Get dihedral and improper angle list (only needed for sel2)
  set dihedral [get_dihedral_and_improper $sel2 $mid2 $parfile]
  lappend RData  [lindex $dihedral 1]

  # Get angle list
  set angle [get_angle $sel2 $mid2 $parfile]
  lappend RData [lindex $angle 1]
 
  # Append sel1 and sel0 charge to RData
  lappend RData [$sel1 get {charge}]
  lappend RData [$sel0 get {charge}]


  # Get excluded pair lists (Final data should be coherent.)
  set num_sel1 [$sel1 num]
  set num_sel2 [$sel2 num]
  set num_sel0 [$sel0 num]
  set num_dual [$dual num]
  set natom [$sel num]
  set sellist1 [$sel1 get index]; # incoherent
  set sellist2 [$sel2 get index]; # incoherent
  set sellist0 [$sel0 get index]; # incoherent
  if {$num_dual==0} {
     puts " "
     puts "No dual atoms. Pure MD."
     puts " "
     set num_dual 0
     set duallist {}
  } else {
     set num_dual [$dual num]
     set duallist [$dual get index]; # incoherent 
  }
  
  # Get excluded list 3, excluded list 4, and the corresponding displacement arrays.
  set exc_sel1 [getexclist $sel1 $mid 0 1 $duallist]
  set exc_sel2 [getexclist $sel2 $mid2 $num_sel1 0 $duallist]
  set exc_data3 [combineexc [lindex $exc_sel1 1] [lindex $exc_sel2 1] [lindex $exc_sel1 0] [lindex $exc_sel2 0]]
  set exc_data4 [combineexc [lindex $exc_sel1 3] [lindex $exc_sel2 3] [lindex $exc_sel1 2] [lindex $exc_sel2 2]]

  # Append package1 (displacement array and 1-3 excluded list) and package2 (displacement array and 1-4 scaled list) to IData
  lappend IData $exc_data3
  lappend IData $exc_data4
  unset exc_sel1 exc_sel2 exc_data3 exc_data4


  # Append sellist_lig (coherent sellist1 + coherent sellist2) to IData 
  for {set i 0} {$i < $num_sel1+$num_sel2+$num_sel0} {incr i} {
      lappend sellist_lig $i
  }
  lappend IData $sellist_lig


  # Get mutated vectors and append them to IData (no need for coherent)
  for {set i 0} {$i < $num_mutsite} {incr i} {
     lappend temp_mutvec1 [lindex $arg(mutvec) [expr 4*$i+0]]
     lappend temp_mutvec1 [lindex $arg(mutvec) [expr 4*$i+1]]
     lappend mutvec2 [lindex $arg(mutvec) [expr 4*$i+2]]
     lappend mutvec2 [lindex $arg(mutvec) [expr 4*$i+3]]
  }
  foreach {item} $temp_mutvec1 {
    for {set i 0} {$i < $num_sel1} {incr i} {
       if {$item==[lindex $sellist1 $i]} {
          lappend mutvec1 $i
          break
       }
    }
  }
  set mutvec2 [shift $mutvec2 [lindex $sellist2 0] $num_sel1]
  lappend IData [combine2 $mutvec1 $mutvec2]


  # Append mutated atom index (made coherent) in sel2 to IData
  for {set i 1} {$i <= $num_mutsite} {incr i} {
     set mysel [atomselect $mid2 "beta=$i"]
     if {[llength [$mysel get index]] == 0} {
        error "The mutated atoms are not labeled in pdb file."
     }
     set mymut [shift [$mysel get index] [lindex $sellist2 0] $num_sel1]
     set mysize [llength $mymut]
     for {set j 0} {$j < $mysize} {incr j} {
        lappend mutatom [lindex $mymut $j]
     }
     lappend mutsize $mysize
  }
  lappend IData [combine2 $mutsize $mutatom]  


  # Append dihedral list to IData
  lappend IData [lindex $dihedral 0]

  # Append angle list to IData
  lappend IData [lindex $angle 0]

  # Append sel1 index and sel0 index (in original dcd) to IData
  lappend IData [$sel1 get {index}]
  lappend IData [$sel0 get {index}]


  # Get SamDihed 
  if [info exists arg(samdihedlist)] {
     set samdihedlist_flag 1
  } else {
     set samdihedlist_flag 0
     set arg(samdihedlist) {}
  }

  if [info exists arg(impscandihed)] {
     set impscandihed_flag 1
  } else {
     set impscandihed_flag 0
     set arg(impscandihed) {}
  }

  if [info exists arg(samdihedval)] {
     set samdihedval_flag 1
  } else {
     set samdihedval_flag 0
     set arg(samdihedval) {}

  }

  set SamDihed [::VSS::init_SamData $mid2 4 $samdihedlist_flag $impscandihed_flag $samdihedval_flag $arg(samdihedlist) $arg(impscandihed) $arg(samdihedval)]

  # Get SamAngle
  if [info exists arg(samanglelist)] {
     set samanglelist_flag 1
  } else {
     set samanglelist_flag 0
     set arg(samanglelist) {}
  }
                            
  if [info exists arg(impscanangle)] {
     set impscanangle_flag 1
  } else {
     set impscanangle_flag 0
     set arg(impscanangle) {}
  }
                            
  if [info exists arg(samangleval)] {
     set samangleval_flag 1
  } else {
     set samangleval_flag 0
     set arg(samangleval) {}
  }

  set SamAngle [::VSS::init_SamData $mid2 3 $samanglelist_flag $impscanangle_flag $samangleval_flag $arg(samanglelist) $arg(impscanangle) $arg(samangleval)]

  # Get SamData
  set SamData [::VSS::combine2 $SamDihed $SamAngle]


  # Letitia: Run Job Here!
  set mytime [clock seconds]
  puts "[clock format $mytime -format "Job starts at %Y-%m-%d %H:%M:%S"]"
  set xyzq2 [$sel2 get {x y z charge}] 
  set dcd [list $dcdfile $xstfile $fint $fstep $fend]

  set vssdata [vssrun $GData $IData $RData $SamData $xyzq2 $dcd] 
 
  set mytime [clock seconds]
  puts "[clock format $mytime -format "Job ends at %Y-%m-%d %H:%M:%S"] "

  close $xstfd

  unset RData 
  return 
}


proc ::VSS::make_grid_from_cell { cell resolution {max 0} } {
  foreach { o a b c } $cell { break }
  set a [expr 1.0 * [veclength $a] / $resolution]
  set b [expr 1.0 * [veclength $b] / $resolution]
  set c [expr 1.0 * [veclength $c] / $resolution]
  set na [good_fft_dim $a $max]
  set nb [good_fft_dim $b $max]
  set nc [good_fft_dim $c $max]
  return [list $na $nb $nc]
}


proc ::VSS::good_fft_dim { r {max 0} } {
  variable default_grid_max_dim
  if { $max == 0 } { set max $default_grid_max_dim }
  if { $max < 8 } { error "max dimension of $max is too small" }
  set nl [list 8 10 12 16 20 24 30 32 36 40 48 50 56 60 64 72 80 \
		84 88 96 100 108 112 120 128]
  set goodn 8
  foreach mi {1 2 3 4 5 6 8 10} {
    foreach ni $nl {
      set n [expr $mi * $ni]
      if {$n > $max} { return $goodn }
      if {$n > $goodn} { set goodn $n }
      if {$goodn >= $r} { return $goodn }
    }
  }
  return $goodn
}


proc ::VSS::getVDWdata {VDWlist Allpara} {
# Read VDW parameter from parameter file and append it to RData
  variable RData
  foreach {item} $VDWlist {
     set temppara [::VSS::getvdwparam $Allpara  $item]
     lappend temppara2 [lindex $temppara 0]
     lappend temppara2 [lindex $temppara 1]
     for {set i 2} {$i < 4} {incr i} {
        if { [llength [lindex $temppara $i]] == 1 } {
           # scaled VDW parameters are given
           lappend temppara2 [lindex $temppara $i]
        } else {
           # scaled VDW parameters are not given
           lappend temppara2 [lindex $temppara [expr $i-2]]
        }
     }
     lappend RData $temppara2
     unset temppara2;
  }
}


proc ::VSS::combine2 {lista listb} {
# Combine 2 lists into one list, according to input order
  foreach {item} $lista {
    lappend totallist $item
  }
  foreach {item} $listb {
    lappend totallist $item
  }
  return $totallist
}


proc ::VSS::combine3 {lista listb listc} {
# Combine 3 lists into one list, according to input order
  foreach {item} $lista {
    lappend totallist $item
  }
  foreach {item} $listb {
    lappend totallist $item
  }
  foreach {item} $listc {
    lappend totallist $item
  }
  return $totallist
}


proc ::VSS::combineexc {disp1 disp2 list1 list2} {
# Combine 2 selections' data into one, according to the input order 
# Append one additional number (total num) at the end of disp
# Letitia 2015.07.18
  set num1 [llength $list1]
  set num2 [llength $list2]

  foreach {item} $disp1 {
    lappend totaldata $item
  }
  foreach {item} $disp2 {
    lappend totaldata [expr $item+$num1]
  }
  lappend totaldata [expr $num1+$num2]

  foreach {item} $list1 {
    lappend totaldata $item
  }
  foreach {item} $list2 {
    lappend totaldata $item
  }

  return $totaldata
}


proc ::VSS::shift {mylist argminus argplus} {
# Shift numbers in a list according to given arguments
  set shiftval [expr $argplus-$argminus]
  foreach {item} $mylist {
     lappend result [expr $item+$shiftval]
  }
  return $result
}


proc ::VSS::getexclist {sel mid seldisp flagsel1 duallist} {
# Return excluded list exclist3 and exclist4. 
# The list data are coherent and has no duallist's contribution. 
# Return also the displacement for exclist3 and exclist4. 
# The displacement is not yet coherent with the full sel1+sel2+sel0 system.
# Letitia 17.07.2015
  set natom [$sel num]
  set bondlist [$sel getbonds]
  set anglestring [molinfo $mid get angles]
  set dihedralstring [molinfo $mid get dihedrals]
  set sellist [$sel get index]
  set firstaid [lindex $sellist 0]
  set lastaid [lindex $sellist end]
  set ndual [llength $duallist]

  set bnum [llength [lindex $bondlist 0]]
  set anum [string length $anglestring]
  set dnum [string length $dihedralstring]

  # Get all excluded pairs from bonds (directly change the atom index to coherent aid)
  if {$flagsel1 && $ndual!=0} {
     for {set i 0} {$i < $natom} {incr i} {
        set item [lindex $bondlist $i]
        set nitem [llength $item]
        for {set j 0} {$j < $nitem} {incr j} {
           # Skip duallist's contribution and the other selections
           set flagsel 1
           set nid [lindex $item $j]
           for {set k 0} {$k<$ndual} {incr k} {
              set dualid [lindex $duallist $k]
              if {$nid == $dualid} {set flagsel 0}
              if {$nid == $dualid} {set flagsel 0}
           }
           if {$nid < $firstaid} {set flagsel 0}
           if {$nid > $lastaid} {set flagsel 0}
           if {$flagsel} {
              lappend exclist3 "[expr $i+$seldisp] [expr $nid+$seldisp-$firstaid]"
           } 
        }
     }
  } else {
     for {set i 0} {$i < $natom} {incr i} {
        set item [lindex $bondlist $i]
        set nitem [llength $item]
        for {set j 0} {$j < $nitem} {incr j} {
           # Skip duallist's contribution and the other selections
           set flagsel 1
           set nid [lindex $item $j]
           if {$nid < $firstaid} {set flagsel 0}
           if {$nid > $lastaid} {set flagsel 0}
           if {$flagsel} {
              lappend exclist3 "[expr $i+$seldisp] [expr $nid+$seldisp-$firstaid]"
           }
        }
     }
  }
  unset item nitem

  # Keep only non-repeated pairs
  if {[info exist exclist3] !=0} {
     set exclist3 [lsort -unique $exclist3]
  }

  # Get all excluded pairs from angles
  if {$flagsel1 && $ndual!=0} {
     if {$anum > 2} {
        set rep [list "\{" "" "\}" "" "\\" ""]
        set templist [split $anglestring]
        set tempstring [string map $rep $templist]
        set templist [split $tempstring]
        foreach {item1 item2 item3 item4} $templist {
           # Skip duallist's contribution and the other selections
           set flagsel 1
           for {set i 0} {$i<$ndual} {incr i} {
              set dualid [lindex $duallist $i]
              if {$item2 == $dualid} {set flagsel 0}
              if {$item4 == $dualid} {set flagsel 0}
           }
           if {$item2 < $firstaid} {set flagsel 0}
           if {$item2 > $lastaid} {set flagsel 0}
           if {$item4 < $firstaid} {set flagsel 0}
           if {$item4 > $lastaid} {set flagsel 0}
           if {$flagsel} {
              lappend exclist3 "[expr $item2+$seldisp-$firstaid] [expr $item4+$seldisp-$firstaid]"
              lappend exclist3 "[expr $item4+$seldisp-$firstaid] [expr $item2+$seldisp-$firstaid]"
           }
        }
        unset templist tempstring item1 item2 item3 item4
     }
  } else {
     if {$anum > 2} {
        set rep [list "\{" "" "\}" "" "\\" ""]
        set templist [split $anglestring]
        set tempstring [string map $rep $templist]
        set templist [split $tempstring]
        foreach {item1 item2 item3 item4} $templist {
           # Skip the other selections
           set flagsel 1
           if {$item2 < $firstaid} {set flagsel 0}
           if {$item2 > $lastaid} {set flagsel 0}
           if {$item4 < $firstaid} {set flagsel 0}
           if {$item4 > $lastaid} {set flagsel 0}
           if {$flagsel} {
              lappend exclist3 "[expr $item2+$seldisp-$firstaid] [expr $item4+$seldisp-$firstaid]"
              lappend exclist3 "[expr $item4+$seldisp-$firstaid] [expr $item2+$seldisp-$firstaid]"
           }
        }
        unset templist tempstring item1 item2 item3 item4
     }
  }
 
  # Keep only non-repeated pairs
  if {[info exist exclist3] !=0} {
     set exclist3 [lsort -unique $exclist3]
  }

  # Get all excluded pairs from dihedral angles
  if {$flagsel1 && $ndual!=0} {
     if {$dnum > 2} {
        set rep [list "\{" "" "\}" "" "\\" ""]
        set templist [split $dihedralstring]
        set tempstring [string map $rep $templist]
        set templist [split $tempstring]
        foreach {item1 item2 item3 item4 item5} $templist {
           # Skip duallist's contribution
           set flagsel 1
           for {set i 0} {$i<$ndual} {incr i} {
              set dualid [lindex $duallist $i]
              if {$item2 == $dualid} {set flagsel 0}
              if {$item5 == $dualid} {set flagsel 0}
           }
           if {$item2 < $firstaid} {set flagsel 0}
           if {$item2 > $lastaid} {set flagsel 0}
           if {$item5 < $firstaid} {set flagsel 0}
           if {$item5 > $lastaid} {set flagsel 0}
           if {$flagsel} {
              lappend exclist4 "[expr $item2+$seldisp-$firstaid]  [expr $item5+$seldisp-$firstaid]"
              lappend exclist4 "[expr $item5+$seldisp-$firstaid]  [expr $item2+$seldisp-$firstaid]"
           }
        }
        unset templist tempstring item1 item2 item3 item4 item5
     }
  } else {
     if {$dnum > 2} {
        set rep [list "\{" "" "\}" "" "\\" ""]
        set templist [split $dihedralstring]
        set tempstring [string map $rep $templist]
        set templist [split $tempstring]
        foreach {item1 item2 item3 item4 item5} $templist {
           # Skip the other selections
           set flagsel 1
           if {$item2 < $firstaid} {set flagsel 0}
           if {$item2 > $lastaid} {set flagsel 0}
           if {$item5 < $firstaid} {set flagsel 0}
           if {$item5 > $lastaid} {set flagsel 0}
           if {$flagsel} {
              lappend exclist4 "[expr $item2+$seldisp-$firstaid]  [expr $item5+$seldisp-$firstaid]"
              lappend exclist4 "[expr $item5+$seldisp-$firstaid]  [expr $item2+$seldisp-$firstaid]"
           }
        }
        unset templist tempstring item1 item2 item3 item4 item5
     }
  }

  # Keep only non-repeated pairs
  if {[info exist exclist4] !=0} {
     set exclist4 [lsort -unique $exclist4]
  }
 
  # Sort the data for exclist3 and excdisp3
  if {[info exist exclist3] !=0} {
     set mycount $seldisp
     foreach {item} $exclist3 {
        set item1 [lindex $item 0]
        set item2 [lindex $item 1]
        lappend Arr3($item1) $item2
     }
     set Arrint $seldisp
     set Arrend [expr $seldisp+$natom]
     set dispcount 0
     for {set i $Arrint} {$i < $Arrend} {incr i} {
        if {[info exist Arr3($i)] !=0} {
           set templist [lsort -integer $Arr3($i)]
           foreach {item} $templist {
              lappend trimexclist3 $item
           }
           lappend trimexcdisp3 $dispcount
           incr dispcount [llength $Arr3($i)]
        } else {
           lappend trimexcdisp3 $dispcount
        }
     }
  } else {
     for {set i 0} {$i < $natom} {incr i} {
        lappend trimexcdisp3 0
     }
  }

  # Sort the pairs for exclist4
  if {[info exist exclist4] !=0} {
     set mycount $seldisp
     foreach {item} $exclist4 {
        set item1 [lindex $item 0]
        set item2 [lindex $item 1]
        lappend Arr4($item1) $item2
     }
     set Arrint $seldisp
     set Arrend [expr $seldisp+$natom]
     set dispcount 0
     for {set i $Arrint} {$i < $Arrend} {incr i} {
        if {[info exist Arr4($i)] !=0} {
           set templist [lsort -integer $Arr4($i)]
           foreach {item} $templist {
              lappend trimexclist4 $item
           }
           lappend trimexcdisp4 $dispcount
           incr dispcount [llength $Arr4($i)]
        } else {
           lappend trimexcdisp4 $dispcount
        }
     }
  } else {
     for {set i 0} {$i < $natom} {incr i} {
        lappend trimexcdisp4 0
     }
  }

  # Set empty list if trimexclist3 and trimexclist4 do not exist
  if {[info exist trimexclist3]==0} { set trimexclist3 {} }
  if {[info exist trimexclist4]==0} { set trimexclist4 {} }

  return [list $trimexclist3 $trimexcdisp3 $trimexclist4 $trimexcdisp4]
}


proc ::VSS::get_dihedral_and_improper {sel mid parfile} {
# Return dihedral and improper list (with angle duplicated for periodicity)
# Return dihedral and improper parameters
# Letitia 13.04.2015

  # Get All parameters
  package require readcharmmpar 1.2
  set typelist [$sel get type]
  set Allpara [::VSS::read_charmm_parameters $parfile]
  set natom [$sel num]
  set rep [list  "{" "" "}" "" "unknown" ""]
  set outflag 1;   # Change to 2 when Improper is on!

  # Dihedral
  set mystring [molinfo $mid get dihedrals]
  if {![llength [lindex $mystring 0]]} {
     puts  "No dihedrals in molecule!"
     incr outflag -1
  } else {
     set templist [string map $rep $mystring]
     # Dihedral Parameters
     foreach {item1 item2 item3 item4} $templist {
        set mytype [list [lindex $typelist $item1] [lindex $typelist $item2] [lindex $typelist $item3] [lindex $typelist $item4]]
        set tempdata [::VSS::getdihedparam $Allpara $mytype]
        set ndihedral [expr [llength $tempdata]]
        for {set i 1} {$i < $ndihedral} {incr i} {
           lappend dihedlist [list $item1 $item2 $item3 $item4]
           lappend diheddata [lindex $tempdata $i]
        }
     }
  }
# # Improper 
# set mystring [molinfo $mid get impropers]
# if {![llength [lindex $mystring 0]]} { 
#    puts "No Impropers in molecule!"
#    incr outflag -1
# } else {
#    set templist [string map $rep $mystring]
#    # Improper Parameters (Append to dihedral)
#    foreach {item1 item2 item3 item4} $templist {
#       set mytype [list [lindex $typelist $item1] [lindex $typelist $item2] [lindex $typelist $item3] [lindex $typelist $item4]]
#       set tempdata [::VSS::getimproparam $Allpara $mytype]
#       set ndihedral [expr [llength $tempdata]]
#       for {set i 1} {$i < $ndihedral} {incr i} {
#          lappend dihedlist [list $item1 $item2 $item3 $item4]
#          lappend diheddata [lindex $tempdata $i]
#       }
#    }
# }

  # If there is no data, output an empty list
  if {$outflag} {
     return [list $dihedlist $diheddata]
  } else {
     return {}
  }
}


proc ::VSS::get_angle {sel mid parfile} {
# Return angle list (with angle duplicated for periodicity)
# Return angle parameters
# Letitia 10.11.2015

  # Get All parameters
  package require readcharmmpar 1.2
  set typelist [$sel get type]
  set Allpara [::VSS::read_charmm_parameters $parfile]
  set natom [$sel num]
  set rep [list  "{" "" "}" "" "unknown" ""]
  set outflag 1;

  # Angles
  set mystring [molinfo $mid get angles]
  if {![llength [lindex $mystring 0]]} {
     puts  "No angles in molecule!"
     incr outflag -1
  } else {
     set templist [string map $rep $mystring]
     # Angle Parameters
     foreach {item1 item2 item3} $templist {
        set mytype [list [lindex $typelist $item1] [lindex $typelist $item2] [lindex $typelist $item3]]
        set tempdata [::VSS::getangleparam $Allpara $mytype]
        lappend anglelist [list $item1 $item2 $item3]
        if {[llength [lindex $tempdata 2]]} {
           lappend angledata [list [lindex [lindex $tempdata 1] 0] [lindex [lindex $tempdata 1] 1] [lindex [lindex $tempdata 2] 0] [lindex [lindex $tempdata 2] 1]]
        } else {
           lappend angledata [list [lindex [lindex $tempdata 1] 0] [lindex [lindex $tempdata 1] 1] 0.0 0.0]
        }
     }
  }
 
  # If there is no data, output an empty list
  if {$outflag} {
     return [list $anglelist $angledata]
  } else {
     return {}
  }
}


proc ::VSS::init_SamData {mid2 nlist samdatalist_flag impscan_flag samdataval_flag samdatalist impscan samdataval} {
# Initiate SamData (Data for dihedral and/or angle sampling)
# Letitia 2015.12.06

  # Declare global variables and vss-global variables
  variable M_PI

  # SamData_list
  set SamData_list {}
  set SamData_nsam {}
  set SamData_samangle {}
  set SamData_naff {}
  set SamData_afflist {}
  if {$samdatalist_flag} {
     foreach {item} $samdatalist {
        if {$item==";"} { continue } else {lappend SamData_list $item}
     }
     set SamData_nentry [expr [llength $SamData_list]/$nlist]
  } else {
     set SamData_nentry 0
     set SamData_scanflag 0
     if {$nlist==4} {
           puts " "
           puts "No atom index of dihedral is given! No dihedral sampling is performed!"
        } else {
           puts "No atom index of angle is given! No angle sampling is performed!"
           puts " "
     }
     set SamData [list $SamData_nentry $SamData_scanflag $SamData_list $SamData_nsam $SamData_samangle $SamData_naff $SamData_afflist]
     return $SamData
  }

  # Assign number of sampling and sample angles for SamData
  if {$impscan_flag} { 
     if {$samdataval_flag} {
        set samdataval_flag 0
        puts " "
        if {$nlist==4} {
           puts "Both impscan and samsangle are set. Use impscan for dihedral sampling!"
        } else {
           puts "Both impscan and samsangle are set. Use impscan for angle sampling!"
        }
        puts " "
     }
     set SamData_scanflag 1
     for {set i 0} {$i<$SamData_nentry} {incr i} {lappend SamData_nsam $impscan}
     set totalnsam [expr $SamData_nentry*$impscan]
     for {set i 0} {$i<$totalnsam} {incr i} {lappend SamData_samangle 0}; # initialize the samangle for impscan, update in C.
  } else {
     if {!$samdataval_flag} {
        if {$nlist==4} {
           puts " "
           puts "No dihedral sampling detail is given. No dihedral sampling is performed!"
        } else {
           puts "No angle sampling detail is given. No angle sampling is performed!"
           puts " "
        }
        set SamData_nentry 0
        set SamData_scanflag 0
        set SamData [list $SamData_nentry $SamData_scanflag $SamData_list $SamData_nsam $SamData_samangle $SamData_naff $SamData_afflist]
        return $SamData
     }
  }

  if {$samdataval_flag} {
     set SamData_scanflag 0
     set nsam 0
     foreach {item} $samdataval {
        if {$item==";"} {
           lappend SamData_nsam $nsam
           set nsam 0
        } else {
           incr nsam
           lappend SamData_samangle [expr $item*$M_PI/180.0]
        }
     }
     if {[llength $SamData_nsam] != $SamData_nentry} {
        error "SamData: number of nsam != nentry. Note that samdataval input should end with ;."
     }
  }

  # Build the affected list for SamData
  if {$nlist==4} {
     if {$samdatalist_flag} {
        foreach {item1 item2 item3 item4} $SamData_list {
           set affdata [::VSS::getaff $mid2 $item2 $item3]
           set naff [lindex $affdata 0]
           set donelist [lindex $affdata 1]
           lappend SamData_naff $naff
           for {set i 0} {$i<$naff} {incr i} {lappend SamData_afflist [lindex $donelist $i]}
        }
     }
  } elseif {$nlist==3} {
     if {$samdatalist_flag} {
        foreach {item1 item2 item3} $SamData_list {
           set affdata [::VSS::getaff $mid2 $item2 $item3]
           set naff [lindex $affdata 0]
           set donelist [lindex $affdata 1]
           lappend SamData_naff $naff
           for {set i 0} {$i<$naff} {incr i} {lappend SamData_afflist [lindex $donelist $i]}
        }
     }
  } else {
    error "Bond list has $nlist atoms. Such bond type is not supported!"
  }

  set SamData [list $SamData_nentry $SamData_scanflag $SamData_list $SamData_nsam $SamData_samangle $SamData_naff $SamData_afflist]
  return $SamData
}


proc ::VSS::getaff {mid2 item2 item3} {
# Get affected atom list, called by init_SamData.
# Letitia 2015.12.06

  set donelist {}
  lappend donelist $item2 
  set todolist $item3
  while {[llength $todolist]!=0} {
     set affsel [atomselect $mid2 "index [lindex $todolist 0]"] 
     set afflist [lindex [$affsel getbonds] 0]
     set nafflist [llength $afflist]
     # update the todo list
     for {set i 0} {$i<$nafflist} {incr i} {
        set affid [lindex $afflist $i]
        set repeat [string first $affid $donelist]; # search the donelist to see if
        if {$repeat==-1} {lappend todolist $affid}; # not repeated, add to todolist
     }
     # move the first item of the todolist to the donelist
     lappend donelist [lindex $todolist 0]
     set todolist [lrange $todolist 1 end]
  }

  # The donelist is the afflist for this angle (except for the first element)
  set donelist [lrange $donelist 1 end]
  set naff [llength $donelist]

  return [list $naff $donelist]
}


# The following code reads xst file. They are taken from the
# orginal tcl script of pbctool

###########################################################
# Read unitcells from xst or xst file and return a list   #
# of vectors [list $cellorigin $cell1 $cell2 $cell3].     #
###########################################################

proc ::VSS::read_xst {file} {
   set cell {}
   set cell1 {}
   set cell2 {}
   set cell3 {}
   set cello {}

   set fd [open $file r]
      while {![eof $fd]} {
         set line [gets $fd]
         if {[string first \# $line]==-1 && [llength $line]>0} {
            # Get PBC vectors
            set cell1 [lrange $line 1 3]
            set cell2 [lrange $line 4 6]
            set cell3 [lrange $line 7 9]
            set cello [lrange $line 10 12]
            if {[lindex $line 0] > 0} {
               lappend cell [list $cello $cell1 $cell2 $cell3]
            }
         }
      }
      close $fd
   return $cell
}

# Below taken from readcharmmpar and modified for VSS
# Letitia 2017.05.02
#####################################################
### Get VDW parameters for atom with type $type.  ###
#####################################################

proc ::VSS::getvdwparam { params type } {
   set par {}

   # Check for each param set if it matches normal or reverse
   set found 0
   set foundparam {}
   foreach pset [lindex $params 4] {
      set ptype [lindex $pset 0]
      set found [string equal $ptype $type]
      if {$found} { set foundparam $pset; break }
   }

   if {!$found} {
      #puts "No parameters found for $type"
      return
   } else {
      set eps   [format "%6.4f" [lindex $foundparam 1 0]]
      set Rmin2 [format "%6.4f" [lindex $foundparam 1 1]]
      if {[llength [lindex $foundparam 2]]} {
         set eps_14   [format "%6.4f" [lindex $foundparam 2 0]]
         set Rmin2_14 [format "%6.4f" [lindex $foundparam 2 1]]
         set par [list $eps $Rmin2 $eps_14 $Rmin2_14]
      } else {
         set par [list $eps $Rmin2 {} {}]
      }
   }

   return $par
}

#########################################################
# Assign Parameters to an angle.                        #
# Returns something like                                #
# {{NR1 CPH1 HR3} {25.0 124.00} {20.00 2.14000}}        #
#########################################################

proc ::VSS::getangleparam { params types } {
   set par {}

   # Reversed typelist
   set revtypes [lrevert $types]

   # Check for each param set if it matches normal or reverse
   set found 0
   set foundparam {}
   # Modified by Letitia: there was a bug in this code
   foreach pset [lindex $params 1] {
      set ptypes [lindex $pset 0]
      set found [string equal $ptypes $types]
      if {$found} { set foundparam $pset; break }
      set found [string equal $ptypes $revtypes]
      if {$found} { 
         set foundparam $pset; 
         # 	 set angle [lrevert $angle]; 
         # 	 set names [lrevert $names]; 
         break
      }
   }

   if {!$found} {
      puts "No parameters found for $types / $revtypes"
      return {}
   } else {
      return [list [lindex $foundparam 0] [lindex $foundparam 1] [lindex $foundparam 2]]
   }
}

#########################################################
# Assign Parameters to a dihedral.                      #
# Returns something like                                #
# {{ON5 CN7 CN7 ON6B} {0.4 6 0.0} ...}                  #
#                                                       #
# Letitia: code fixed, runs and provides periodicity    #
# 13.04.2015                                            #
#########################################################

proc ::VSS::getdihedparam { params types } {  
   set par {}

   # Reversed typelist
   set revtypes [lrevert $types]
   # Get wildcard typelist
   set xtypes [list X [lindex $types 1] [lindex $types 2] X]
   # Reversed wildcard typelist
   set rxtypes [list X [lindex $revtypes 1] [lindex $revtypes 2] X]

   # Check for each param set if it matches normal, reverse or including wildcards
   set found 0
   set flagout 0
   set foundparam {}
   # Modified by Letitia: there was a bug in this code
   foreach pset [lindex $params 2] {
      set ptypes [lindex $pset 0]

      # Modified by Letitia: there could be periodicity!
      set found [string equal $ptypes $types]
      if {$found} { 
         lappend foundparam [lindex $pset 1]
         set typename [lindex $pset 0]
         set flagout 1
      }

      set found [string equal $ptypes $revtypes]
      if {$found} { 
         lappend foundparam [lindex $pset 1] 
         set typename [lindex $pset 0]
         set flagout 1
      }

      set found [string equal $ptypes $xtypes]
      if {$found} { 
         lappend foundparam [lindex $pset 1]
         set typename [lindex $pset 0]
         set flagout 1
      }

      set found [string equal $ptypes $rxtypes]
      if {$found} { 
         lappend foundparam $pset 
         set typename [lindex $pset 0]
         set flagout 1
      }
   }
  
   if {!$flagout} {
      puts "No parameters found for dihed $types"
      return {}
   } else {
      set foundparam [lsort -unique $foundparam]
      return [linsert $foundparam 0 $typename]
   }
}

########################################################################
### Read bond, angle, dihed, improper, nonbonded, nbifx, hbond       ###
### parameters from $parfile and returns them in a list.             ###
########################################################################

proc ::VSS::read_charmm_parameters { parfile } {
   set fd [open $parfile]

   # Read the parameter file:
   set bonds {}
   set angles {}
   set dihedrals {}
   set impropers {}
   set nonbonded {}
   set nbfix {};   # Explicit nonbonded pair interactions
   set hbonds {}
   set comment {}
   set section {}
   set remark {}

   set skip 0
   foreach line [split [read -nonewline $fd] \n] {
      set trimmed [string trim $line]
      if {[string index $trimmed 0]=="*"} {
         # Just a part of the header comment
         continue
      }
      if {[string index $trimmed 0]=="!"} {
         set comment $trimmed
         continue
      }
      if {[string length $trimmed]==0} {
         continue
      }
      set remidx [string first "!" $trimmed]
      if ($remidx>=0) {
         set remark [string range $trimmed $remidx end]
         set line [string range $trimmed 0 [expr $remidx-1]]
      }
      set keyword [lindex $line 0]

      if {$skip} { set skip 0; continue }

      # Look for beginning of next section
      if { [string equal BONDS     $keyword] } { set section bond; continue }
      if { [string equal ANGLES    $keyword] } { set section angl; continue }
      if { [string equal DIHEDRALS $keyword] } { set section dihe; continue }
      if { [string equal IMPROPERS  $keyword] } { set section impr; continue };   # Modified by Letitia: there was a bug    
      if { [string equal NONBONDED $keyword] } { 
         set section nonb; 
         # Look for continued line
         if {[lindex $line end] == "-"} {
            # Skip line:
            set skip 1
         }
         continue
      }
      if { [string equal NBFIX     $keyword] } { set section nbfi; continue }
      if { [string equal HBOND     $keyword] } { set section hbon; continue }

      if {$section=="nonb"} {
         if {[llength $nonbonded] && [llength $comment] && [string first ! $line]>1} {
            # append the second line comment if the ! is not at pos 0 or 1.
            lset nonbonded end end [join [list [lindex $nonbonded end end] $comment]]
         }
         if {[lindex $line 4] == 0.0} {
            #set rem [lrange $line 7 end]
            if {![llength $remark]} { set remark {{}} }
            lappend nonbonded [list [lindex $line 0] [lrange $line 2 3] [lrange $line 5 6] $remark]
         } else {
            #set rem [lrange $line 4 end]
            if {![llength $remark]} { set remark {{}} }
            lappend nonbonded [list [lindex $line 0] [lrange $line 2 3] {} $remark]
         }

      } elseif {$section=="bond"} {
         if {[llength $bonds] && [llength $comment]} {
            # append the second line comment
            lset bonds end end [join [list [lindex $bonds end end] $comment]]
         }
         lappend bonds     [list [lrange $line 0 1] [lrange $line 2 3] [lrange $line 4 end]]

      } elseif {$section=="angl"} {
         if {[llength $angles] && [llength $comment]} {
            # append the second line comment
            lset angles end end [join [list [lindex $angles end end] $comment]]
         }
         set ub [lrange $line 5 6]
         if {[string index $ub 0]=="!" || ![llength $ub]} {
            lappend angles    [list [lrange $line 0 2] [lrange $line 3 4] [list] [lrange $line 5 end]]
         } else  {
            lappend angles    [list [lrange $line 0 2] [lrange $line 3 4] $ub [lrange $line 7 end]]
         }

      } elseif {$section=="dihe"} {
         if {[llength $dihedrals] && [llength $comment]} {
            # append the second line comment
            lset dihedrals end 2 [join [list [lindex $dihedrals end 2] $comment]]
         }
         lappend dihedrals [list [lrange $line 0 3] [lrange $line 4 6] [lrange $line 7 end]]

      } elseif {$section=="impr"} {
         if {[llength $impropers] && [llength $comment]} {
            # append the second line comment
            lset impropers end 2 [join [list [lindex $impropers end 2] $comment]]
         }
         lappend impropers [list [lrange $line 0 3] [lrange $line 4 6] [lrange $line 7 end]]
      } elseif {$section=="nbfi"} {
         lappend nbfix     [list [lrange $line 0 1] [lrange $line 2 3] [lrange $line 4 end]]
      } elseif {$section=="hbon"} {
         lappend hbond     [list [lrange $line 0 1] [lrange $line 2 3] [lrange $line 4 end]]
      }
      set comment [list]
   }
   close $fd

   return [list $bonds $angles $dihedrals $impropers $nonbonded $nbfix $hbonds]
}

########################################
### Reverses the order of a list.    ###
########################################

proc ::VSS::lrevert { list } {
   set newlist {}
   for {set i [expr [llength $list]-1]} {$i>=0} {incr i -1} {
      lappend newlist [lindex $list $i]
   }
   return $newlist
}


