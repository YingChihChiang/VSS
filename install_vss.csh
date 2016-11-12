#!/bin/csh -e

# Set the original VMD plugin directory and VSS install directory
set VMDPLUGINDIR = ~/vmd/build/plugins
set TCLINCDIR = /usr/include/tcl8.5
set TCLLIBDIR = /usr/lib/x86_64-linux-gnu

# Check if VMDPLUGINDIR is valid
if ( ! -d "$VMDPLUGINDIR" || ! -d $VMDPLUGINDIR/LINUXAMD64 ) then
  echo "Your VMD plugin directory seems not right. Would you double check it?"
  exit 1
endif

# Compile and install vss and readcharmmpar
make clean
setenv TCLINC -I$TCLINCDIR
setenv TCLLIB -F$TCLLIBDIR
setenv PLUGINDIR $VMDPLUGINDIR
make LINUXAMD64 TCLINC=$TCLINC TCLLIB=$TCLLIB
make distrib


# Add vss to startup file
#touch ~/.vmdrc
#grep -q "vss" ~/.vmdrc || echo "source $VMDPLUGINDIR/LINUXAMD64/tcl/vss1.0/vss.tcl" >> ~/.vmdrc


# Final print
echo "    "
echo "    "
echo "VSS plugin has been installed to your VMD. Enjoy!"
echo "    "
echo "    "
