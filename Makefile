##################################################################################
# WHATEVER YOU DO, PLEASE DO IT TO Makefile.config AND LEAVE THIS MAKEFILE ALONE #
##################################################################################
include Makefile.config


#------------------------------------------------------------------#
# the menu                                                         #
#------------------------------------------------------------------#
MASTER_DEFINEFLAGS	=	$(DEFINEFLAGS)

#------------------------------------------------------------------#
# make settings available to all other Makefiles                   #
#------------------------------------------------------------------#
export CC
export FC
export OPTIMIZE
export CCFLAGS
export LNFLAGS
export MASTER_DEFINEFLAGS
export MAKE

# everything in src/
#===================
Cweb:	FORCE dirs
	cd src;\
	${MAKE} Cweb;\
	mv -f Cweb ../bin/Cweb

# everything in analysis/
#=======================
Cweb2ascii:	FORCE dirs
	cd analysis;\
	${MAKE} Cweb2ascii;\
	mv -f Cweb2ascii ../bin

Cwebhalos:	FORCE dirs
	cd analysis;\
	${MAKE} Cwebhalos;\
	mv -f Cwebhalos ../bin

#-------------------------------------------------------------------#
# "make clean" 
#-------------------------------------------------------------------#
clean:	FORCE
	@echo '';\
	echo '*=======================================================================*';\
	echo ' ${MAKE} clean: removing all unneeded files';\
	echo '*=======================================================================*';\
	echo '';\
	echo ''
	cd src; ${MAKE} clean

#-------------------------------------------------------------------#
# "make veryclean" 
#-------------------------------------------------------------------#
veryclean:	FORCE
	@echo '';\
	echo '*=======================================================================*';\
	echo ' ${MAKE} veryclean: removing all unneeded files (incl. bin/*)';\
	echo '*=======================================================================*';\
	echo '';\
	echo ''
	rm -f *~ *~.* bin/*;\
	cd src; ${MAKE} veryclean

FORCE:	

dirs:
	mkdir -p bin
