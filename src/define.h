#ifndef DEFINE_INCLUDED
#define DEFINE_INCLUDED

/*=============================================================================
 * this is written into the logfile just for information
 *=============================================================================*/
#define VERSION 1.0
#define BUILD   0

/*=============================================================================
 * by default, Cweb only calculates the Vweb,
 * but you can also opt for Dweb and Pweb
 * Note, the Dweb and Pweb results are *not* written into the ASCII file
 *=============================================================================*/
//#define DWEB
//#define DWEB_MHD // Weiguang's way of calculating 2nd order (mixed) partial derivatives
//#define PWEB
//#define PWEB_MHD // Weiguang's way of calculating 2nd order (mixed) partial derivatives

/*=============================================================================
 * here we switch on/off various features of Cweb (i.e. DEFINEFLAGS)
 * (either switch them on here on using -DFEATURE in Makefile.config)
 *=============================================================================*/
//#define LAMBDA_TH       0.44  // this just dumps some info to stderr about web elements fractions at the end
//#define HUBBLE_Z        // use H(z) as the velocity shear tensor normalisation

//#define WRITE_ASCII     // writes ASCII output file (note, *only* Cweb will be written, DWEB/PWEB will only be written to binary outputs)

//#define READ_GRIDDATA    // instead of reading the simulation data, the code now simply reads the GRIDDATA file
//#define WRITE_GRIDDATA   // writes a file that contains the velocity and density field on the grid

//#define MULTIMASS            /* you MUST switch this on when the simulation features particles of different masses */

//#define GAS_PARTICLES        /* you MUST switch this on when the simulation contains gas and/or star particles     */
                               /* a few more words about this switch:
                                *   - for historical reasons it is called GAS_PARTICLES even though it actually
                                *     deals with baryons in the simulation
                                *   - if you do not switch */

//#define DMfocus           // ignore all gas and star particles, only keeping (all) dark matter particles (cf. main.c for the selection if() statement)
//#define USE_FFT3             /* if you have FFT3 installed, you might like to use that                             */
//#define WITH_MPI             /* switch on MPI domain decomposition                                                 */
//#define WITH_OPENMP          /* switch on OpenMP parallisation of for-loops                                        */

//#define TSC
//#define CIC                    /* assignment schemes */
//#define NGP

//#define BYTESWAP             /* forces a byteswap of the input file                                                */













/*---------------------------------------------------
 *                  STANDARD
 *         (best to not touch these!)
 *--------------------------------------------------*/
//#define TSC_SMOOTHING        /* smooth density and velocity fields on TSC nodes                                   */
//#define MHD_SMOOTHING        /* smooth density and velocity fields on MHD nodes                                   */
#define PERIODIC             /* use periodic boundary conditions                                                  */
#ifndef VERBOSE
#define VERBOSE              /* let the user know what's going on                                                 */
#endif
//#define VERBOSE2             /* dump as much runtime information as possible                                      */

//#define REFINE_BARYONIC_MASS         // use mass as refinement criterion for baryons (but number density for dark matter!)
//#define CHECK_RLIMIT_NOFILE          // uses system functions to increase file descriptor limitation if needed
#define FOPENCLOSE                   // open/close files, rather than opening multiple simulation files at the same time
#define BCASTHEADER                  // only one MPI task will read all the relevant header information and then broadcast
//#define NCPUREADING_EQ_NFILES        // this should speed up I/O of multiple snapshot files, but only works for this condition

#ifdef DWEB_MHD
#define DWEB
#endif

#ifdef PWEB_MHD
#define PWEB
#endif

/*--------------------------------------------------
 *           -DWITH_MPI or -DAHFrestart
 *--------------------------------------------------*/
#if (defined WITH_MPI || AHFrestart)
#  undef  REF_TEST
#  undef  VERBOSE
#endif


/*--------------------------------------------------
 *                    -DGADGET
 *--------------------------------------------------*/
#ifdef GADGET2
#define GADGET
#endif
#ifdef GADGET
 #ifndef MULTIMASS
  #define MULTIMASS	/* avoid unnecessary compiler warnings */
 #endif
#define GAS_PARTICLES
#endif

/*--------------------------------------------------
 *                    -DTIPSY
 *--------------------------------------------------*/
#ifdef TIPSY
#define  MULTIMASS
#define  GAS_PARTICLES
#endif

/*--------------------------------------------------
 *                    -DDEVA
 *--------------------------------------------------*/
#ifdef DEVA2
#define DEVA
//#define DEVA2_QHULL_FILE "../snapshots/cube00.qhull"
#endif

#ifdef DEVA
#define MULTIMASS
#define GAS_PARTICLES
#endif

/*--------------------------------------------------
 *                -DMARE_NOSTRUM
 *--------------------------------------------------*/
#ifdef MARE_NOSTRUM
#define MULTIMASS
#define GAS_PARTICLES
#endif

/*--------------------------------------------------
 *                -DMETALHACK
 *--------------------------------------------------*/
#ifdef METALHACK
#	ifndef MULTIMASS
#		define MULTIMASS
#	endif
#	ifndef GAS_PARTICLES
#		define GAS_PARTICLES;
#	endif
#	define METALDIE \
	fprintf(stderr,\
	        "There's metal in the air tonight, can you hear it call\n"\
	        "If you ain't got the balls, to take it you can\n"\
	        "Leave the hall\n");\
	exit(-666);
#endif


/*----------------------------------------------------------------------------
 *  misc definitions
 *----------------------------------------------------------------------------*/
#ifdef GAS_PARTICLES
#ifndef MULTIMASS
#define MULTIMASS
#endif
#endif

#ifdef NO_GAS
#undef GAS_PARTICLES
#endif

#ifdef VERBOSE2
#define VERBOSE
#endif

#ifdef VERBOSE
#define REF_TEST
#endif

#ifndef TSC
#ifndef CIC   /* forgotten to define mass assignemnt scheme ? => use TSC then... */
#ifndef NGP
#define TSC
#endif
#endif
#endif

/*--------------------------------------------
 * more transparent to read in source-code... 
 *--------------------------------------------*/
#ifndef CONTINUE
#define TERMINATE
#define TERMINATE2  /* used in leavers.c */
#define VERBOSELOG
#endif /* CONTINUE */

#ifdef PERIODIC
#define PERIODIC_X
#define PERIODIC_Y
#define PERIODIC_Z
#endif /* PERIODIC */

#endif

