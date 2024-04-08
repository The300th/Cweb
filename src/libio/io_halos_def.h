#ifndef IO_HALOS_DEF_H
#define IO_HALOS_DEF_H

/**
 * \file io_halos_def.h
 *
 * Provides the structure definition for the AHF_HALOS file structure.
 */


/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include <stdio.h>
#include <stdint.h>

#include "io_halos_header_def.h"
#include "io_file.h"


/**********************************************************************\
 *    Global defines, structure definitions and typedefs              * 
\**********************************************************************/

/**
 * The file structure itself
 */
struct io_halos_struct {
	/** Holds the filetype ID, *must* be first in the structure */
	io_file_type_t ftype;
#ifdef WITH_MPI
	/** The global rank of the process */
	int rank;
	/** The size of the global communicator */
	int size;
	/** Stores the communicator used for intra libio communication */
	MPI_Comm mycomm;
	/** The size of the intra-library communicator */
	int size_mycomm;
	/** The rank of the local process */
	int rank_mycomm;
#endif
	/** Holds the file stream */
	FILE *file;
	/** Holds the filename */
	char *fname;
	/** In which mode should this be opened */
	io_file_mode_t mode;
	/** Holds the header information */
	io_halos_header_t header;
	/** Stores the minimal weight of the particles */
	double minweight;
	/** Stores the maximal weight of the particles */
	double maxweight;
};

/** Convenient typedef */
typedef struct io_halos_struct io_halos_struct_t;

/** Convenient typedef */
typedef io_halos_struct_t *io_halos_t;


#endif /* IO_HALOS_DEF_H */
