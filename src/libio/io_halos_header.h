#ifndef IO_HALOS_HEADER_H
#define IO_HALOS_HEADER_H

/**
 * \file io_halos_header.h
 *
 * Provides functions for reading and writing the header of AHF_HALOS
 * files.
 */


/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include "io_halos_header_def.h"
#include "io_halos_def.h"
#include "io_logging.h"


/**********************************************************************\
 *    Global defines, structure definitions and typedefs              * 
\**********************************************************************/


/**********************************************************************\
 *    Prototypes of global functions                                  * 
\**********************************************************************/

/**
 * \brief Reads a AHF_HALOS header.
 *
 * The function will first rewind the file pointer to the start
 * of the header and then read in everything. Then it will
 * position the file pointer to the start of the particle date.
 *
 * \param log  A logging object.
 * \param f    An AHF_HALOS file object.
 *
 * \return A freshly filled header, or NULL, in case of memory
 *         problems.
 */
extern io_halos_header_t
io_halos_header_get(io_logging_t log, io_halos_t f);

/**
 * \brief Generates an empty header object.
 *
 * \param log  A logging object.
 *
 * \return A freshly allocated header, or NULL, in case of memory
 *         problems.
 */
extern io_halos_header_t
io_halos_header_new(io_logging_t log);

/**
 * \brief This will delete an halos_header object.
 *
 * \param log      The logging object.
 * \param *header  A pointer to the variable holding the header object.
 *                 This variable will be set to NULL.
 *
 * \return Nothing.
 */
extern void
io_halos_header_del(io_logging_t log, io_halos_header_t *header);

/**
 * \brief Writes the header to the file
 *
 * \param log     A logging object.
 * \param header  The header to write.
 * \param f       The file the header will be written to.
 *
 * \return Nothing.
 */
extern void
io_halos_header_write(io_logging_t log,
                      io_halos_header_t header,
                      io_halos_t f);

/**
 * \briefs Writes the header information to the logfile.
 *
 * \param log     The logging object.
 * \param header  The header object.
 *
 * \return Nothing.
 */
extern void
io_halos_header_log(io_logging_t log, io_halos_header_t header);


#endif /* IO_HALOS_HEADER_H */
