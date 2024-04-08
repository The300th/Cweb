/**
 * \file io_halos_header.c
 *
 * Provides functions for reading and writing the header of AHF_HALOS
 * files.
 */


/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <stdint.h>

#include "io_halos_header.h"
#include "io_util.h"


/**********************************************************************\
 *    Local defines, structure definitions and typedefs               * 
\**********************************************************************/


/**********************************************************************\
 *    Prototypes of local functions                                   * 
\**********************************************************************/


/**********************************************************************\
 *    Implementation of global functions                              * 
\**********************************************************************/
extern io_halos_header_t
io_halos_header_get(io_logging_t log, io_halos_t f)
{
	io_halos_header_t dummy;
  FILE *fpin;
  char fileline[2048];
  
	/* Some sanity checks */
	if ((f == NULL) || (f->file == NULL))
		return NULL;

	/* Check if there already is a header, do nothing then */
	if (f->header != NULL)
		return f->header;

  /* Create the header structure */
	dummy = malloc(sizeof(io_halos_header_struct_t));
	if (dummy == NULL) {
		io_logging_memfatal(log, "AHF_HALOS header structure");
		return NULL;
	}
  
  /* Read all relevant parameters from separate file ahf_halos.info generated in io_parameters.c */
  fprintf(stderr,"o reading parameters from file ahf_halos.info\n");
  if((fpin = fopen("ahf_halos.info","r")) == NULL)
   {
    io_logging_fatal(log,"Could not open ahf_halos.info containing the following information:");
    io_logging_fatal(log,"52817495                HALOS_N");
    io_logging_fatal(log,"20.0                    HALOS_BOXSIZE [Mpc/h]");
    io_logging_fatal(log,"0.74                    HALOS_REDSHIFT");
    io_logging_fatal(log,"Please generate this file. Aborting now!");
    free(dummy);
    return NULL;
   }
  fgets(fileline, 2048, fpin);
  sscanf(fileline,"%ld",&(dummy->no_part));
  fgets(fileline, 2048, fpin);
  sscanf(fileline,"%lf",&(dummy->boxsize));
  fgets(fileline, 2048, fpin);
  sscanf(fileline,"%lf",&(dummy->a_current));
  dummy->a_current = 1./dummy->a_current - 1.; // conversion to redshift
  fclose(fpin);

	/* Go to the header */
	rewind(f->file);


	/* over-read the header line */
	fscanf(f->file, "%s", dummy->header);
  
  /* set the remaining header parameters */
  // TODO: do we actually need any of these for Cweb?????
  dummy->multi_mass  = 1;
  dummy->no_species  = dummy->no_part;
  dummy->total_mass  = 0.;
  dummy->no_timestep = 100;
  dummy->omega0      = 1.;
  dummy->lambda0     = 0.;
  dummy->pmass       = 1.e10;
  dummy->a_initial   = 1000.;
  
	/* Set the header */
	f->header = dummy;

	/* And return it */
	return dummy;
}

extern io_halos_header_t
io_halos_header_new(io_logging_t log)
{
	io_halos_header_t dummy;

	dummy = malloc(sizeof(io_halos_header_struct_t));
    if (dummy == NULL) {
		io_logging_memfatal(log, "AHF_HALOS header structure");
		return NULL;
	}

	return dummy;
}

extern void
io_halos_header_del(io_logging_t log, io_halos_header_t *header)
{
	if ( (header == NULL) || (*header == NULL) )
		return;

	free(*header);

	*header = NULL;

	return;
}

extern void
io_halos_header_write(io_logging_t log,
                      io_halos_header_t header,
                      io_halos_t f)
{
	/* NOT IMPLEMENTED */

	return;
}

extern void
io_halos_header_log(io_logging_t log, io_halos_header_t header)
{
	io_logging_msg(log, INT32_C(5),
	               "Headerobject information:");
	io_logging_msg(log, INT32_C(5),
	               "  Header:                      %s",
	               header->header);
	io_logging_msg(log, INT32_C(5),
	               "  Multimass:                   %" PRIi32,
	               header->multi_mass);
	io_logging_msg(log, INT32_C(5),
	               "  Number of particles:         %li",
	               header->no_part);
	io_logging_msg(log, INT32_C(5),
	               "  Number of mass species:      %li",
	               header->no_species);
	io_logging_msg(log, INT32_C(5),
	               "  Total mass:                  %e",
	               header->total_mass);
	io_logging_msg(log, INT32_C(5),
	               "  Number of timestep         : %" PRIi32,
	               header->no_timestep);
	io_logging_msg(log, INT32_C(5),
	               "  Boxsize                    : %e",
	               header->boxsize);
	io_logging_msg(log, INT32_C(5),
	               "  Omega0                     : %e",
	               header->omega0);
	io_logging_msg(log, INT32_C(5),
	               "  lambda0                    : %e",
	               header->lambda0);
	io_logging_msg(log, INT32_C(5),
	               "  pmass                      : %e",
	               header->pmass);
	io_logging_msg(log, INT32_C(5),
	               "  a_initial                  : %e",
	               header->a_initial);
	io_logging_msg(log, INT32_C(5),
	               "  a_current                  : %e",
	               header->a_current);

	return;
}


/**********************************************************************\
 *    Implementation of local functions                               * 
\**********************************************************************/
