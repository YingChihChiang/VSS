# Makefile for VSS potential library

.SILENT:

.SUFFIXES: 

AR= ar
ARFLAGS = cr
RANLIB = ranlib

COMPILEDIR = ../compile
ARCHDIR=${COMPILEDIR}/lib_${ARCH}/tcl/vss1.0
INCDIR=-Isrc -Isrc/pmepot -Isrc/dcd
VPATH = src ${ARCHDIR}
MYCXXFLAGS = -m64 -O2 -fPIC -Wall -fopenmp
#MYCXXFLAGS = -m64 -O2 -fPIC -Wall -fopenmp -D Ham

##
## Only build vss if we have a Tcl library 
##
ifdef TCLLIB
ifdef TCLINC
ifdef TCLLDFLAGS
TARGETS = ${ARCHDIR} ${ARCHDIR}/libvss.so
endif
endif
endif


#
# Rules
#
bins:
win32bins:
dynlibs: ${TARGETS}
staticlibs: 
win32staticlibs: 

distrib:
	for localname in `find ../compile -name libvss.so -print` ; do \
		pluginname=`echo $$localname | sed s/..\\\/compile\\\/lib_// `; \
		dir=`dirname $(PLUGINDIR)/$$pluginname`; \
		mkdir -p $$dir; \
		cp $$localname $(PLUGINDIR)/$$pluginname; \
		cp src/pkgIndex.tcl $$dir ; \
		cp src/vss.tcl $$dir ; \
	done

${ARCHDIR}:
	mkdir -p ${ARCHDIR}

VSSOBJS = ${ARCHDIR}/vssmain.o \
        ${ARCHDIR}/vsscases.o \
        ${ARCHDIR}/vssdE.o \
	${ARCHDIR}/pmepot.o \
	${ARCHDIR}/gridcharges.o \
	${ARCHDIR}/recipspace.o \
	${ARCHDIR}/pub3dfft.o \
	${ARCHDIR}/impscan.o \
	${ARCHDIR}/calc_bondedE.o \
	${ARCHDIR}/mutation.o \
        ${ARCHDIR}/spacedecomp.o \
	${ARCHDIR}/array.o \
	${ARCHDIR}/readdcd.o \

${ARCHDIR}/libvss.so : ${VSSOBJS} ${ARCHDIR}/tcl_vss.o
	if [ -n "${TCLSHLD}" ]; \
	then ${TCLSHLD} $(LOPTO)$@ ${VSSOBJS} ${ARCHDIR}/tcl_vss.o ${TCLLIB} ${TCLLDFLAGS} ${LDFLAGS} -fopenmp ; \
	else ${SHLD} $(LOPTO)$@ ${VSSOBJS} ${ARCHDIR}/tcl_vss.o ${TCLLIB} ${TCLLDFLAGS} ${LDFLAGS} -fopenmp ; \
	fi

${ARCHDIR}/tcl_vss.o: tcl_vss.cpp vss.h
	${CXX} ${MYCXXFLAGS} ${TCLINC} ${INCDIR} -D_${ARCH} -DVSSTCLDLL_EXPORTS -c $< $(COPTO)$@

${ARCHDIR}/%.o: %.cpp
	${CXX} ${MYCXXFLAGS} ${INCDIR} -c $< $(COPTO)$@

${ARCHDIR}/%.o: pmepot/%.cpp
	${CXX} ${MYCXXFLAGS} ${INCDIR} -c $< $(COPTO)$@



