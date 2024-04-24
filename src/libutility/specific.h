#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "../common.h"
#include "../param.h"
#include "../tdef.h"

#ifndef SPECIFIC_INCLUDED
#define SPECIFIC_INCLUDED

/* AMIGA specific routines */
MINMAX  MinMax               (double,double,double);
MINMAX  MinMaxBound          (double,double,double,double);
void    read_amiga_header    (FILE *infile, info_io *io, int *SWAPBYTES);
void    sanity_check         ();
void    ic_unit_conversion   ();
double  init_header_masses   ();
void    binning_parameter    (HALO halo, int *nbins, double *dist_min, double *dist_max);
void    get_c2fslope         (double func[3][3][3], double slope[3]);
int     get_axes             (double itensor[3][3], double *axis1, double *axis2, double *axis3);
int     idx_inv              (long unsigned *idx, int numHalos, int j);
double  f1mod                (double x, double y);
double  Laplace_pot          (nptr tsc_nodes[3][3][3], double spacing);
double  Laplace_temp1        (nptr tsc_nodes[3][3][3], double spacing);
void    write_filename       (char *f_name, char *prefix, unsigned l1dim);
void    write_logfile        (double timecounter, double timestep, int no_timestep);
void    write_parameterfile  ();
double  calc_cNFW            (double V2_max, double V2_vir);
extern int
cmp_sfckey_part(const void *p1, const void *p2);
int     check_for_MPIboundary(double x, double y, double z);
int     check_for_zerodens   (nptr tsc_nodes[3][3][3], double dens_threshold);

#define WRITEAHFLOGO(fx)  \
fprintf(fx,"=======================================================================\n");\
fprintf(fx,"\t       CCCCCC                                   \n");\
fprintf(fx,"\t    CCCC                                         \n");\
fprintf(fx,"\t  CCC                                           \n");\
fprintf(fx,"\t CC                W           W  EEEEEE  BBBBB    \n");\
fprintf(fx,"\t CC                 W         W   E       B    B   \n");\
fprintf(fx,"\t  CCC          ====  W   W   W    EEEE    BBBBB    \n");\
fprintf(fx,"\t    CCCC              W W W W     E       B    B   \n");\
fprintf(fx,"\t       CCCCCC          W   W      EEEEEE  BBBBB    (v%3.1f/%03d)\n",VERSION,BUILD);\
fprintf(fx,"=======================================================================\n\n");


#endif
