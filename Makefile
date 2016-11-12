##
## Plugin compilation and distribution rules
## Copyright 2002-2011, Board of Trustees of the University of Illinois
## Theoretical and Computational Biophysics Group
##
## Modifications Copyright 2015-2016, The Chinese University of Hong Kong (CUHK)
## Department of Physics, Biophysics Group
##
.SILENT:

default: make-arch-help

world:
	@echo "Building all supported targets..."
	csh -f build.csh

BUILDDIRS = \
  readcharmmpar \
  vss

distrib: 
	@echo "Populating distribution directory with compiled plugins"
	for dir in $(BUILDDIRS); do cd $$dir && $(MAKE) distrib && cd .. || exit 1 ; done

dynlibs:
	for dir in $(BUILDDIRS); do cd $$dir && $(MAKE) dynlibs && cd .. || exit 1 ; done

staticlibs:
	for dir in $(BUILDDIRS); do cd $$dir && $(MAKE) staticlibs && cd .. || exit 1 ; done

bins:
	for dir in $(BUILDDIRS); do cd $$dir && $(MAKE) bins && cd .. || exit 1 ; done


include Make-arch

clean:
	find compile \( -name \*.o -o -name \*.a -o -name \*.so -o -name \*.exp -o -name \*.lib -o -name \*.h \) -print | xargs rm -f
	find compile \( -name lib_\* \) -print | xargs rm -rf
	rm -f log.*

checkperms:
	@echo "Looking for bad file permissions..."
	find . ! -perm +004

