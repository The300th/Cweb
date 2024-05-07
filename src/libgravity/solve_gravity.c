#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "gravity.h"
#include "../libutility/utility.h"
#include "../libamr_serial/amr_serial.h"

/*
 * there are two ways to stop the GS iteration procedure:
 *
 * 1. the residuals are smaller than the truncation error
 *
 * 2. the residuals are smaller than LIMIT given below
 *
 */
#define LIMIT     0.00001

#define ONE_THIRD 0.33333333333333333333

#include "gravity.h"

# ifdef PWEB

/*==============================================================================
 * solve on coarsest grid
 *==============================================================================*/
void solve_cg(gridls *cur_grid)
{
  flouble *dens_array;         /* density array pointer */
  pqptr    cur_pquad;          /* current pquad         */
  cqptr    cur_cquad;          /* current cquad         */
  nqptr    cur_nquad;          /* current nquad         */
  nptr     cur_node;           /* current node          */
  long     i, j, k, l1dim, FFTarray_length;
  double   FourPiGa;
  
  /* conversion factor to go from densito to source term */
  FourPiGa = simu.FourPiG*calc_super_a(cur_grid->timecounter);
  
  /* array dimension */
  l1dim           = cur_grid->l1dim;
  FFTarray_length = 2*l1dim*l1dim*l1dim;
  
  /* generate complex (!) density array for FFT */
  if((dens_array = (flouble *) calloc(FFTarray_length, sizeof(flouble))) == NULL)
   {
    fprintf(io.logfile,"solve_cg: could not allocate density array for FFT\n");
    fflush(io.logfile);
    fclose(io.logfile);
    exit(1);
   }
  
  /* fill density array for FFT ... no need for quad-ll's !! */
  cur_pquad = cur_grid->pquad;
  for(k = 0, cur_cquad = cur_pquad->loc; k < l1dim; k++, cur_cquad++)
    for(j = 0, cur_nquad = cur_cquad->loc; j < l1dim; j++, cur_nquad++)
      for(i = 0, cur_node = cur_nquad->loc; i < l1dim; i++, cur_node++)
       {
        dens_array[Re(i,j,k,l1dim)] = cur_node->dens * FourPiGa;  /* real part      */
        dens_array[Im(i,j,k,l1dim)] = 0.0;                        /* imaginary part */
       }
  
  /* solve by FFT */
  fft_potential(dens_array, l1dim);
  
  /* fill node potential values */
  for(k = 0, cur_cquad = cur_pquad->loc; k < cur_grid->l1dim; k++, cur_cquad++)
    for(j = 0, cur_nquad = cur_cquad->loc; j < cur_grid->l1dim; j++, cur_nquad++)
      for(i = 0, cur_node = cur_nquad->loc; i < cur_grid->l1dim; i++,cur_node++)
        cur_node->pot = dens_array[Re(i,j,k,l1dim)];
  
  /* destroy memory assigned to dens_array */
  free(dens_array);
  
#ifdef TEST_POT
   {
      nptr     tsc_nodes[3][3][3];
      double   Delta_pot, Source;
      FILE     *fp;
      
      fp = fopen("test_pot.dat","w");
      cur_pquad = cur_grid->pquad;
      for(k = 0, cur_cquad = cur_pquad->loc; k < l1dim; k++, cur_cquad++)
         for(j = 0, cur_nquad = cur_cquad->loc; j < l1dim; j++, cur_nquad++)
            for(i = 0, cur_node = cur_nquad->loc; i < l1dim; i++, cur_node++)
            {
               tsc_nodes[1][1][1] = cur_node;
               get_TSCnodes(cur_grid, cur_pquad, cur_cquad, cur_nquad, tsc_nodes, &k, &j, &i);
               
               Delta_pot = Laplace_pot(tsc_nodes, pow2(cur_grid->spacing));
               Source    = cur_node->dens * FourPiGa;

               fprintf(fp,"%d %d %d %g %g %g\n",i,j,k,Delta_pot,Source,Delta_pot-Source);
            }
      fclose(fp);
   }
#endif
  
}


/*============================================================================
 * test convergence of current grid
 *============================================================================*/
boolean converged(gridls *cur_grid)
{
  double trunc_error;
  double residual;
  
  residual     =             cur_grid->cur_resid;
  trunc_error  = ONE_THIRD * cur_grid->trunc_err;
  
  if     (residual <= LIMIT)
    return TRUE;
  else if(residual <= (trunc_error * CONVCRIT)) 
    return TRUE;
  else
    return FALSE;
}

/*=============================================================================
 * test for slow convergence
 *=============================================================================*/
boolean slow_conv(gridls *cur_grid)
{
#ifdef NO_MULTIGRID
  return FALSE;
#endif
  if(cur_grid->cur_resid > (ETA * cur_grid->old_resid))
    return TRUE;
  else
    return FALSE;
}

/*=============================================================================
 * solve for potential on domain grids
 *=============================================================================*/
void solve_dom_gravity(gridls *grid_list)
{
  int     grid_no;             /* the current grid number     */
  int     lstgrid_no;          /* no_grids minus one          */
  int     i;                   /* index for GS sweep loop     */
  gridls *cur_grid;            /* pointer to the current grid */
  boolean cg_conv;             /* has current grid converged? */
  int     gs_sweeps;
  
  gs_sweeps = DOMSWEEPS;
  
  solve_cg(global.dom_grid);
  return;
  
}


/*==============================================================================
 * control routine for multi grid solver
 *==============================================================================*/
void solve_gravity(gridls *grid_list, int curgrid_no)
{
  solve_dom_gravity(grid_list);
  
  // TODO: also solve on refinment grids...
}

#endif //PWEB
