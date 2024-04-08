.PHONY: clean veryclean

Cweb2ascii:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
Cwebhalos:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
Cwebfilaments:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}

CFLAGS  =	${OPTIMIZE} ${CCFLAGS} ${DEFINEFLAGS}
LFLAGS  =	${LNFLAGS}

export DEFINEFLAGS

LIBS = 	../src/libio/libio.a \
	../src/libutility/libutility.a \
	../src/libamr_serial/libamr_serial.a

# HOW TO MAKE THE ACTUAL PROGRAMMES
Cwebhalos:	Cwebhalos.c ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o Cwebhalos Cwebhalos.c ${LIBS} ${LIBS} -lm

Cweb2ascii:	Cweb2ascii.c ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o Cweb2ascii Cweb2ascii.c ${LIBS} ${LIBS} -lm

# HOW TO MAKE THE LIBRARIES
../src/libio/libio.a:
	cd ../src/libio ; ${MAKE} libio.a

../src/libamr_serial/libamr_serial.a:
	cd ../src/libamr_serial ; ${MAKE} libamr_serial.a

../src/libutility/libutility.a:
	cd ../src/libutility ; ${MAKE} libutility.a

clean:
	rm -f *.o

veryclean:
	rm -f *.o *~ *~.*
