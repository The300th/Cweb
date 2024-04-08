#******************************************************************#
#                                                                  #
# Makefile for               Cweb                                  #
#                           ======                                 #
#                                                                  #
#******************************************************************#

.PHONY: clean veryclean

#*******************************************************************#
Cweb:			DEFINEFLAGS = ${MASTER_DEFINEFLAGS} -DAHF

export DEFINEFLAGS

# libraries in use
# ----------------
LIB_MATH 	=	m

# flags for compiling all individual routines
# -------------------------------------------
CFLAGS	=	${OPTIMIZE} ${CCFLAGS} ${DEFINEFLAGS}

# flags for linking all parts together
# ------------------------------------
LFLAGS	=	${LNFLAGS}

# all object files
# ----------------
OBJS	=	\
	comm.o \
	common.o\
	startrun.o
  
LIBS =    \
			libio/libio.a \
			libsfc/libsfc.a \
			libamr_serial/libamr_serial.a \
			libgravity/libgravity.a \
			libutility/libutility.a

# what's to make ?
# ----------------
Cweb:	main.o ${OBJS} ${LIBS}
	$(CC) ${PROFILE} -o Cweb main.o  ${OBJS} ${LIBS} $(LIBS) -l${LIB_MATH} ${LFLAGS}


clean:
	rm -f *.o Cweb* *.il
	cd libsfc ; ${MAKE} clean
	cd libio ; ${MAKE} clean
	cd libamr_serial ; ${MAKE} clean
	cd libgravity ; ${MAKE} clean
	cd libutility ; ${MAKE} clean

veryclean:
	rm -f *.o Cweb* *~ *~.* *.il
	cd libsfc ; ${MAKE} clean
	cd libio ; ${MAKE} clean
	cd libamr_serial ; ${MAKE} clean
	cd libgravity ; ${MAKE} clean
	cd libutility ; ${MAKE} clean


# dependencies for include files
# ------------------------------
common.h:	param.h tdef.h
tdef.h: 	param.h



# how to make the libraries?
# ------------------------------
libsfc/libsfc.a:
	cd libsfc ; ${MAKE} libsfc.a

libio/libio.a:
	cd libio ; ${MAKE} libio.a

libamr_serial/libamr_serial.a:
	cd libamr_serial ; ${MAKE} libamr_serial.a

libgravity/libgravity.a:
	cd libgravity ; ${MAKE} libgravity.a

libutility/libutility.a:
	cd libutility ; ${MAKE} libutility.a



#dependencies for .o files
# ------------------------
comm.o:	comm.c
	${CC} ${CFLAGS} -c comm.c

main.o:			main.c
	${CC} ${CFLAGS} -c main.c

startrun.o:		startrun.c
	${CC} ${CFLAGS} -c startrun.c
