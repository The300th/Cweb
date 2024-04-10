#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "amr_serial.h"
#include "../libutility/utility.h"

/*===============================================================================
 * gen_domgrids: generate the domain grids
 *===============================================================================*/
gridls *gen_domgrids(int *no_grids)
{
  gridls *grid_list;       /* pointer to array of grid entries          */
  gridls *cur_grid;        /* pointer to current grid                   */
  double  mass2dens;        /* conversion factor to get mass density     */
  double  mass2partdens;    /* conversion factor to get particle density */
  double  frac, dl1dim, dno_part;
  int     exp;
  int     l1dim;
  
  /* calculate the total number of domain grids */
  frac        = (double) simu.NGRID_DOM / (double) simu.NGRID_MIN;
  exp         = (float) (log(frac)/log(2.));
  *no_grids   = exp + 1;
  
  /*======================================
   * there are '*no_grids' to be created:
   *        a) simu.NGRID_MIN^3 grid
   *            ...
   *        x) simu.NGRID_DOM^3 grid
   *=======================================*/
  
  /* allocate memory for all grids */
  grid_list = (gridls *) calloc(*no_grids, sizeof(gridls));
  
  /*
   * 'grid_list' is the pointer to the first gridls-structure
   * out of '*no_grids' of such structures
   */
  
  /* generate the all grids and initialize the structure with the corresponding data */
  for(cur_grid = grid_list, l1dim = simu.NGRID_MIN; l1dim <= simu.NGRID_DOM; l1dim *= 2, cur_grid++)
  {
    /* mass2dens <=> 1/rho_bar [code units] */
    dl1dim        = (double)l1dim;
    dno_part      = (double)simu.no_part;
    mass2dens     = pow3(dl1dim)/simu.no_vpart;
    mass2partdens = pow3(dl1dim)/dno_part;
    
    /*===========================
     * fill in grid_list details
     *===========================*/
    cur_grid->timecounter    = global.super_t;
    cur_grid->l1dim          = l1dim;
    cur_grid->spacing        = (double)1.0 / (double) cur_grid->l1dim;
    cur_grid->spacing2       = pow2(cur_grid->spacing);
    
    cur_grid->masstodens     = mass2dens;
    cur_grid->masstopartdens = mass2partdens;
#ifdef MULTIMASS
    cur_grid->critdens       = simu.Nth_dom * mass2partdens;
#else
    cur_grid->critdens       = simu.Nth_dom * mass2dens;
#endif
    cur_grid->old_resid      = 0.0;
    cur_grid->cur_resid      = 0.0;
    cur_grid->no_sweeps      = 0;
    
    cur_grid->time.potential = 0;
    cur_grid->time.density   = 0;
    cur_grid->time.DK        = 0;
    cur_grid->time.grid      = 0;
    cur_grid->time.hydro     = 0;
    
    cur_grid->multistep             = 0;
    cur_grid->next                  = FALSE;
    
    /*========================================================
     * generate the actual pquad, cquad, and nquad structures
     *========================================================*/
    alloc_quads(cur_grid, l1dim);
  }
  
  /*===========================
   * keep track of domain grid
   *===========================*/
  global.dom_grid   = grid_list + (*no_grids - 1);
  global.domgrid_no = (*no_grids - 1);
  
  return(grid_list);
}


/*==============================================================================
 * gen_refgrid:
 *--------------------
 * 1. test current grid for refinement, and generate adaptive grid (if required)
 * 2. re-link particles
 * 3. ajdust densities as neccessary.
 *==============================================================================*/
boolean gen_refgrid(gridls **grid_list, int *no_grids)
{
  gridls *coa_grid, *fin_grid;     /* old/new grid pointer                    */
  boolean refined;                 /* refined YES/NO flag                     */
  boolean mg;                      /* does refinement hold enough particles   */
  int     idim;
  
#ifdef REF_TEST
  fprintf(stderr,"  - gen_refgrid:         current coarse grid         = %ld\n",
          (*grid_list + (*no_grids-1))->l1dim);
#endif
  
  /* set coarse grid pointer */
  coa_grid = *grid_list + *no_grids - 1;
  fin_grid = *grid_list + *no_grids;
  
  /* timing... */
  coa_grid->time.grid -= time(NULL);
  
  if(coa_grid->next == FALSE)  /* fin_grid does not exist ? */
  {
    /* reallocate the grid list array (including space for one additional grid) */
    (*grid_list) = (gridls *) realloc(*grid_list, (*no_grids + 1) * sizeof(gridls));
    
    if((*grid_list) == NULL)
    {
      fprintf(io.logfile,"gen_refgrid: error reallocating grid list\n");
      fflush(io.logfile);
      fclose(io.logfile);
      exit(1);
    }
    
    /* grid_list block might have moved in memory */
    global.dom_grid = *grid_list + global.domgrid_no;
    coa_grid        = *grid_list + *no_grids - 1;
    fin_grid        = *grid_list + *no_grids;
    
    /* fill in new grid's entries */
    fin_grid->masstodens     = coa_grid->masstodens     * CRITMULTI;
    fin_grid->masstopartdens = coa_grid->masstopartdens * CRITMULTI;
#ifdef MULTIMASS
    fin_grid->critdens       = simu.Nth_ref*fin_grid->masstopartdens;
#else
    fin_grid->critdens       = simu.Nth_ref*fin_grid->masstodens;
#endif
    fin_grid->timecounter    = coa_grid->timecounter;
    fin_grid->timestep       = coa_grid->timestep/2;
    fin_grid->l1dim          = coa_grid->l1dim * 2;
    fin_grid->spacing        = (double)1.0 / (double)fin_grid->l1dim;
    fin_grid->spacing2       = fin_grid->spacing * fin_grid->spacing;
    
    fin_grid->no_pquad       = 0;
    fin_grid->pquad_array    = NULL;
    
    /* we measure the time throughout the two DKD cycles of this fin_grid */
    fin_grid->time.potential = 0;
    fin_grid->time.density   = 0;
    fin_grid->time.DK        = 0;
    fin_grid->time.grid      = 0;
    fin_grid->time.hydro     = 0;
    
    fin_grid->no_sweeps      = 0;
    fin_grid->cur_resid      = 0.;
    
    /* remember that fin_grid exists from now on */
    fin_grid->next = FALSE;
    coa_grid->next = TRUE;
  }
  
  
  /* (re-)set values that are changing within the two DKD cycles the grid lives through */
  fin_grid->old_resid             = 0.0;
  fin_grid->cur_resid             = 0.0;
  
  fin_grid->multistep             = 0;
  
  
  /*------------------------------
   * now call refinement rountine
   *------------------------------*/
  refined = refine_grid(fin_grid, coa_grid);
  
  if(refined)
  {
#ifdef REF_TEST
    fprintf(stderr,"                         placed new refinement       = %ld\n", fin_grid->l1dim);
#endif
    
    /* update number of grids */
    *no_grids += 1;
    
    /* now relink particles to new grid */
    mg = relink(coa_grid, fin_grid);        /* performs: unassign_part */
    
    /* now assign mass to new grid: THIS IS NEEDED AS IT RETURNS "mg" ! */
    zero_dens        (fin_grid);
    mg = assign_npart(fin_grid);
    
#ifdef REF_TEST
    fprintf(stderr,"                      grid %6ld:  nnodes=%ld npart=%ld => %g (<! %g)\n",
            fin_grid->l1dim,fin_grid->size.no_nodes,fin_grid->size.no_part,
            (double)fin_grid->size.no_nodes/(double)fin_grid->size.no_part, CRITMULTI);
#endif
    
    /* only use grid if it holds enough particles */
    if(mg == FALSE)
    {
      refined = FALSE;
      relink_back(coa_grid, fin_grid);
      fin_grid->size.no_nodes         = 0;
      fin_grid->size.no_part          = 0;
      
      free_grid(fin_grid, no_grids);
      
#ifdef REF_TEST
      fprintf(stderr,"                        destroyed new refinement    = %ld (%g)\n",
              fin_grid->l1dim, fin_grid->critdens);
#endif
    }
  }
  
  else
  {
#ifdef REF_TEST
    fprintf(stderr,"                     NO new refinement grid !\n\n");
#endif
    /*
     * no need to call free_grid()
     * -> temporary pquads, cquads, nquads already destroyed within ref_grid()
     *
     * no need to realloc grid_list
     * -> grid will be kept in memory from now on...
     */
  }
  
  /* grid_list block might have changed position in memory */
  global.dom_grid = *grid_list + global.domgrid_no;
  
  /* ...timing */
  coa_grid->time.grid += time(NULL);
  
  return(refined);
}

/*================================================================================
 * recursively generate the AMR hierarchy
 *================================================================================*/
boolean gen_AMRhierarchy(gridls **grid_list, int *no_grids)
{
  gridls *cur_grid;                   /* current grid pointer                  */
  gridls *fin_grid;                   /* fine grid pointer                     */
  gridls *for_grid;                   /* sometimes needed for debugging...     */
  int     fingrid_no, curgrid_no;     /* index of finest and current grid      */
  int     grid_no;
  
  boolean runAHF;  /* only run AHF when the finest grid level has been reached */
  boolean refined; /* refinement flag                                          */
  
  /* up to now we plan to run AHF */
  runAHF = TRUE;
  
#ifdef READ_GRIDDATA
  runAHF = FALSE;    // we read everything from a file and hence no need to for those calls below to get density perfectly right on all levels!
#endif
  
  /*=========================================================================
   * try to refine finest grid:
   *
   * 1. get (number!) density field on currently finest grid right
   * 2. try to refine grid using 'number of particles per node' criterion
   *=========================================================================*/
  
  /* currently the finest grid... */
  fin_grid = *grid_list + *no_grids - 1;
  
  if(!(((global.fst_cycle == FALSE) && (fin_grid->l1dim == global.fin_l1dim)) ||
       ((fin_grid->l1dim) == simu.NGRID_MAX)))
  {
    timing.genrefgrids  -= time(NULL);
    
    /* get number density: "npart/node" */
    fin_grid->time.grid -= time(NULL);
    zero_dens   (fin_grid);
    assign_npart(fin_grid);
    fin_grid->time.grid += time(NULL);
    
    /* try to refine grid under "npart/node" criterion */
    refined = gen_refgrid(grid_list, no_grids);
    
    timing.genrefgrids  += time(NULL);
  }
  else
  {
    refined = FALSE;
  }
  
  /*--------------------------------------------------
   * when refined we make a recursive call to step()
   *--------------------------------------------------*/
  if(refined == TRUE)
  {
    /* recursive call to gen_AMRhierarchy() */
    runAHF = gen_AMRhierarchy(grid_list, no_grids);
    
    curgrid_no      = *no_grids - 1;
    fingrid_no      = *no_grids - 1;
    global.dom_grid = *grid_list + global.domgrid_no;
    cur_grid        = *grid_list + curgrid_no;
  }
  
  /*---------------------------------------------
   * otherwise stay with the actual finest grid
   *---------------------------------------------*/
  else /* refined == FALSE */
  {
    /* no new refinement => stay with current finest grid */
    curgrid_no      = *no_grids - 1;
    fingrid_no      = *no_grids - 1;
    global.dom_grid = *grid_list + global.domgrid_no;
    cur_grid        = *grid_list + curgrid_no;
    
    /* store force resolution in header info */
    io.header.cur_reflevel = (float)(curgrid_no-global.domgrid_no);
    io.header.cur_frcres   = (float)(simu.boxsize/(double)cur_grid->l1dim);
    io.header.cur_frcres  *= 3000.;
    
    
    energy.time -= time(NULL);
    energy.time += time(NULL);
    
    if(global.fst_cycle == TRUE)
    {
      global.fin_l1dim = cur_grid->l1dim;
      global.fst_cycle = FALSE;
    }
  } /* if(refined) */
  
  
  
  /*=========================================================================
   *
   *         WE NOW HAVE ACCESS TO THE FULL BLOWN GRID HIERACHY
   *
   *=========================================================================*/
  
  
  /*=========================================================================
   *          here are those parts dealing with AHFpotcentre,
   *              the actual AHF is happening in main.c
   *=========================================================================*/
  ahf.time -= time(NULL);
  
  /* we only run AHF when having reached the end of the grid hierarchy */
  if(runAHF)
  {
#ifdef VERBOSE
    fprintf(stderr,"  - Obtaining correct density and velocity fields on all grid levels ... ");
#endif
    
    /* get the density absolutely right on all amr levels */
    /*====================================================*/
    timing.densrecovery -= time(NULL);
    /* re-assign density right on all grids */
    for(for_grid = cur_grid; for_grid >= global.dom_grid; for_grid--)
      restore_dens(for_grid);
    if(curgrid_no != fingrid_no)
      refill_dens(cur_grid);
    if(cur_grid != global.dom_grid)
      stack_dens(cur_grid, global.dom_grid+1);
    timing.densrecovery += time(NULL);
    
#ifdef VERBOSE
    fprintf(stderr,"done\n");
#endif
  }
  ahf.time += time(NULL);
  
  /* return FALSE indicating that AHF has been run */
  return(FALSE);
}


/*================================================================================
 * read regular grid data into the finest domain grid
 *================================================================================*/
void read_griddata(gridls **grid_list, int *no_grids, char *filename)
{
  gridls       *dom_grid;                   /* domain grid pointer                     */
  FILE          *gridfile;
  partptr       cur_part;
  pqptr         cur_pquad;
  cqptr         cur_cquad, icur_cquad;
  nqptr         cur_nquad, icur_nquad;
  nptr          cur_node;
  long          ipart, ipart_node, x, y, z, nd[3];
  double        cur_shift;
#ifdef PWEB
  float         tempd[5];
#else
  float         tempd[4];
#endif
  char          line[MAXSTRING], infile[MAXSTRING];
  int           notboundary;
  
#ifdef TSC
  sprintf(infile,"%s.TSCmesh", filename);
#endif
#ifdef CIC
  sprintf(infile,"%s.TSCmesh", filename);
#endif
#ifdef NGP
  sprintf(infile,"%s.NGPmesh", filename);
#endif
  sprintf(infile,"%s-%05d",infile,simu.NGRID_DOM);
  
#ifdef WITH_MPI
  sprintf(infile,"%s-MPIrank%d",infile,global_mpi.rank);
#endif
  
  // the domain grid pointer
  dom_grid = *grid_list + *no_grids - 1;
  
  /* shift of cell centre as compared to edge of box [grid units] */
  cur_shift = 0.5/(double)dom_grid->l1dim;
  
  fprintf(stderr,"read_griddata(): reading mesh data from file %s\n",infile);
  
  if((gridfile = fopen(infile,"rb")) == NULL)
  {
    fprintf(stderr,"could not open %s\n", infile);
    exit(1);
  }
  
  fread(&nd, sizeof(float), 3, gridfile);
  if ((int) nd[0] != (int) simu.NGRID_DOM)
  {
    fprintf(stderr,"mesh dimension not fit!! %ld\t%d\n", nd[0],simu.NGRID_DOM);
    // exit(1);
  }
  
  fread(&global.z,     sizeof(double), 1, gridfile);
  fread(&simu.boxsize, sizeof(double), 1, gridfile);
  fread(&simu.omega0,  sizeof(double), 1, gridfile);
  fread(&simu.lambda0, sizeof(double), 1, gridfile);
  fread(&simu.pmass,   sizeof(double), 1, gridfile);
  fread(&simu.no_part, sizeof(long unsigned), 1, gridfile);  
  
  /* loop over all nodes */
  for(cur_pquad=dom_grid->pquad; cur_pquad!=NULL; cur_pquad=cur_pquad->next)
  {
    for(cur_cquad = cur_pquad->loc, z = cur_pquad->z;
        cur_cquad < cur_pquad->loc + cur_pquad->length;
        cur_cquad++, z++)
    {
      for(icur_cquad=cur_cquad; icur_cquad!=NULL; icur_cquad=icur_cquad->next)
      {
        for(cur_nquad = icur_cquad->loc, y = icur_cquad->y;
            cur_nquad < icur_cquad->loc + icur_cquad->length;
            cur_nquad++, y++)
        {
          for(icur_nquad=cur_nquad; icur_nquad!=NULL; icur_nquad=icur_nquad->next)
          {
            for(cur_node = icur_nquad->loc, x = icur_nquad->x;
                cur_node < icur_nquad->loc + icur_nquad->length;
                cur_node++, x++)
              
            {
#ifdef WITH_MPI
              notboundary = check_for_MPIboundary(((double)x+0.5)/(double)dom_grid->l1dim, ((double)y+0.5)/(double)dom_grid->l1dim, ((double)z+0.5)/(double)dom_grid->l1dim);
#else
              notboundary = TRUE;
#endif
              /* only consider this cell if not within the MPI boundary zone
               *    -> otherwise it will be read by some other MPI task already */
              if(notboundary == TRUE)
              {
#ifdef PWEB
                fread(&tempd, sizeof(float), 5, gridfile);
                cur_node->dens     = tempd[0];  // its density contrast in the data right now
                cur_node->densV[X] = tempd[1];
                cur_node->densV[Y] = tempd[2];
                cur_node->densV[Z] = tempd[3];
                cur_node->pot      = tempd[4];
#else
                fread(&tempd, sizeof(float), 4, gridfile);
                cur_node->dens     = tempd[0];  // its density contrast in the data right now
                cur_node->densV[X] = tempd[1];
                cur_node->densV[Y] = tempd[2];
                cur_node->densV[Z] = tempd[3];
#endif
              }
            }
          }
        }
      }
    }
  }
  fclose(gridfile);
}

/*================================================================================
 * write simulation grid data into the file
 *================================================================================*/
void write_griddata(gridls **grid_list, int *no_grids, char *filename)
{
  gridls        *dom_grid;                   /* domain grid pointer                     */
  FILE          *gridfile;
  partptr       cur_part;
  pqptr         cur_pquad;
  cqptr         cur_cquad, icur_cquad;
  nqptr         cur_nquad, icur_nquad;
  nptr          cur_node;
  long          ipart, ipart_node, x, y, z;
  double        cur_shift;
  char          outfile[MAXSTRING];
  float         tmp_float;
  int           notboundary;
  
  // the domain grid pointer
  dom_grid = *grid_list + *no_grids - 1;
  
  /* shift of cell centre as compared to edge of box [grid units] */
  cur_shift = 0.5/(double)dom_grid->l1dim;
  
#ifdef TSC
  sprintf(outfile,"%s.TSCmesh", filename);
#endif
#ifdef CIC
  sprintf(outfile,"%s.TSCmesh", filename);
#endif
#ifdef NGP
  sprintf(outfile,"%s.NGPmesh", filename);
#endif
  sprintf(outfile,"%s-%05d",outfile,simu.NGRID_DOM);
  
#ifdef WITH_MPI
  sprintf(outfile,"%s-MPIrank%d",outfile,global_mpi.rank);
#endif
  
  fprintf(stderr,"open %s to write mesh data. \n", outfile);
  if((gridfile = fopen(outfile,"wb")) == NULL)
  {
    fprintf(stderr,"could not open %s\n", outfile);
    //exit(1);
  }
  
  fwrite(&simu.NGRID_DOM, sizeof(float), 1, gridfile);
  fwrite(&simu.NGRID_DOM, sizeof(float), 1, gridfile);
  fwrite(&simu.NGRID_DOM, sizeof(float), 1, gridfile);
  
  fwrite(&global.z,     sizeof(double), 1, gridfile);
  fwrite(&simu.boxsize, sizeof(double), 1, gridfile);
  fwrite(&simu.omega0,  sizeof(double), 1, gridfile);
  fwrite(&simu.lambda0, sizeof(double), 1, gridfile);
  fwrite(&simu.pmass,   sizeof(double), 1, gridfile);
  fwrite(&simu.no_part, sizeof(long unsigned), 1, gridfile);
  
  /* loop over all nodes */
  for(cur_pquad=dom_grid->pquad; cur_pquad!=NULL; cur_pquad=cur_pquad->next)
  {
    for(cur_cquad = cur_pquad->loc, z = cur_pquad->z;
        cur_cquad < cur_pquad->loc + cur_pquad->length;
        cur_cquad++, z++)
    {
      for(icur_cquad=cur_cquad; icur_cquad!=NULL; icur_cquad=icur_cquad->next)
      {
        for(cur_nquad = icur_cquad->loc, y = icur_cquad->y;
            cur_nquad < icur_cquad->loc + icur_cquad->length;
            cur_nquad++, y++)
        {
          for(icur_nquad=cur_nquad; icur_nquad!=NULL; icur_nquad=icur_nquad->next)
          {
            for(cur_node = icur_nquad->loc, x = icur_nquad->x;
                cur_node < icur_nquad->loc + icur_nquad->length;
                cur_node++, x++)
              
            {
#ifdef WITH_MPI
              notboundary = check_for_MPIboundary(((double)x+0.5)/(double)dom_grid->l1dim, ((double)y+0.5)/(double)dom_grid->l1dim, ((double)z+0.5)/(double)dom_grid->l1dim);
#else
              notboundary = TRUE;
#endif
              /* only consider this cell if not within the MPI boundary zone
               *    -> otherwise it will be written by some other MPI task already */
              if(notboundary == TRUE)
              {
                tmp_float = cur_node->dens;
                fwrite(&tmp_float, sizeof(float), 1, gridfile);
                
                tmp_float = cur_node->densV[X];
                fwrite(&tmp_float, sizeof(float), 1, gridfile);
                
                tmp_float = cur_node->densV[Y];
                fwrite(&tmp_float, sizeof(float), 1, gridfile);
                
                tmp_float = cur_node->densV[Z];
                fwrite(&tmp_float, sizeof(float), 1, gridfile);
                
#ifdef PWEB
                tmp_float = cur_node->pot;
                fwrite(&tmp_float, sizeof(float), 1, gridfile);
#endif
              }
            }
          }
        }
      }
    }
  }
  fclose(gridfile);
}
