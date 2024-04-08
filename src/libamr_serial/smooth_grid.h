#include "../param.h"
#include "../tdef.h"

#ifndef SMOOTH_GRID_INCLUDED
#define SMOOTH_GRID_INCLUDED

void     Gaussian_smooth_gridFFT   (gridls *cur_grid, double Rsmooth);
void     TSCMHDsmooth_grid         (gridls *cur_grid);

#endif

