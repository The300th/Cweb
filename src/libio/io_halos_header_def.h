#ifndef IO_HALOS_HEADER_DEF_H
#define IO_HALOS_HEADER_DEF_H

/**
 * \file io_halos_header_def.h
 *
 * Provides the structure definition for the AHF_HALOS header
 * structure. Including useful typedefs.
 */


/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include <inttypes.h>


/**********************************************************************\
 *    Global defines, structure definitions and typedefs              * 
\**********************************************************************/

/** Defines the length for the identifying header string */
#define HALOS_HEADER_HEADERSTRING 256

/**
 * The header structure itself
 */
struct io_halos_header_struct {
	char          header[HALOS_HEADER_HEADERSTRING];
	int32_t       multi_mass;
	long          no_part;
	long          no_species;
	double        total_mass;
	int32_t       no_timestep;
	double        boxsize;
	double        omega0;
	double        lambda0;
	double        pmass;
	double        a_initial;
	double        a_current;
};

/** Convenient typedef */
typedef struct io_halos_header_struct io_halos_header_struct_t;

/** Convenient typedef */
typedef io_halos_header_struct_t *io_halos_header_t;


#endif /* IO_HALOS_HEADER_DEF_H */
