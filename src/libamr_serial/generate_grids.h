#include "../param.h"
#include "../tdef.h"

#ifndef GENERATE_GRIDS_INCLUDED
#define GENERATE_GRIDS_INCLUDED

gridls  *gen_domgrids     (int *no_grids);
boolean  gen_refgrid      (gridls **grid_list, int *no_grids);
boolean  gen_AMRhierarchy (gridls **grid_list, int *no_grids);

#ifdef READ_GRIDDATA
void     read_griddata    (gridls **grid_list, int *no_grids, char *filename);
#else
void     write_griddata   (gridls **grid_list, int *no_grids, char *filename);
#endif

#endif
