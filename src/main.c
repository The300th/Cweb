#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdint.h>

#include "common.h"
#ifdef WITH_MPI
#	include <mpi.h>
#	include "libutility/loadbalance.h"
#	include "comm.h"
#endif

#include "libsfc/sfc.h"
#include "startrun.h"
#include "libutility/utility.h"
#include "libamr_serial/amr_serial.h"
#include "libgravity/gravity.h"


#ifdef WITH_MPI
static void
local_communicate_all(void);
#endif

#ifdef AHFrfocus
#include <assert.h>
static void
local_focusSphere(void);
#ifdef AHFrfocusRescale
static void
local_focusSphereRescale(void);
#endif
#endif

/*==============================================================================
 * MAIN: where everything starts ....
 *==============================================================================*/
int main(int argc, char **argv)
{
  gridls  *grid_list;        /* pointer to list of grids            */
  
  int     no_grids;          /* total number of grids               */
  int     no_timestep;       /* number of coarse grid timesteps     */
  int     no_first_timestep; /* number of initial timestep          */
  
  double  timecounter;       /* time variable                       */
  double  timestep;          /* timestep size                       */
  double  timecounter_final; /* for all sorts of tests...           */
  double  Hz, Hz2;
  double  Rsmooth;
  double  Rsmooth_min, Rsmooth_max; // in case we want to implement a for-loop over various smoothing scales
  int     Nsmooth, ismooth;
  
  char     AMIGA_input[MAXSTRING];
  
  /* all the Vwen stuff */
  gridls*       cur_grid;
  pqptr         cur_pquad;
  cqptr         cur_cquad, icur_cquad;
  nqptr         cur_nquad, icur_nquad;
  nptr          cur_node;
  nptr          tsc_nodes[3][3][3];
  nptr          mhd_nodes[5][5][5];
  partptr       cur_part;
  long          ipart, x, y, z, npart;
  double        dx, dy, dz;
  double        twodx, twody, twodz;
  double        Vx, Vy, Vz;
  double        Vx_x2,Vx_x0, Vy_x2,Vy_x0, Vz_x2, Vz_x0;
  double        Vx_y2,Vx_y0, Vy_y2,Vy_y0, Vz_y2, Vz_y0;
  double        Vx_z2,Vx_z0, Vy_z2,Vy_z0, Vz_z2, Vz_z0;
  double        dVxdx, dVxdy, dVxdz;
  double        dVydx, dVydy, dVydz;
  double        dVzdx, dVzdy, dVzdz;
  double        meanV[NDIM],sumweight,weight;
  double        x_fac, v_fac, phi_fac, pot_fac, rho_fac;
  double        local_shear[3][3];
  double        vorticity[3];
  double        itensor[3][3];
    int         converged;
  double        lambda1, lambda2, lambda3;
  long          i, j, k;
  double        Wg, distance2;
  double        dFdx[3][3][3], dFdy[3][3][3], dFdz[3][3][3];
#ifdef DWEB
  double        dtensor[3][3], dambda1, dambda2, dambda3;
#endif
#ifdef PWEB
  double        ptensor[3][3], pambda1, pambda2, pambda3;
#endif
  char          outfile[MAXSTRING], infile[MAXSTRING];
  FILE         *fpout;
  float         tmp_float, tmp_x, tmp_y, tmp_z;
  uint64_t      tmp_uint64;
  int32_t       one=1;
  
  uint64_t      Nknots, Nfilaments, Nsheets, Nvoids, Nnodes_written;
  double        LAMBDA_THRESHOLD;
#ifdef WITH_MPI
  uint64_t      Nknots_all, Nfilaments_all, Nsheets_all, Nvoids_all, Nnodes_written_all;
#endif
  
#ifdef WITH_MPI
  uint64_t newparts;
#endif
  
  int haveneighbours, notboundary;
  
#ifdef DEBUG_CWEB
  double VxMean, VyMean, VzMean, MeanNorm;
  double DensMean;
  FILE *fdebug;
  fdebug = fopen("CwebDEBUG.dat","w");
#endif
  
  /*============================================================
   * we always read the relevant parameters from an input file!
   *===========================================================*/
  if(argc<2)
  {
    fprintf(stderr,"usage:    %s Cweb.input\n", argv[0]);
    fprintf(stderr,"       or %s --parameterfile\n", argv[0]);
    exit(1);
  }
  
  /*============================================================
   * maybe the user only wants the parameterfile?
   *===========================================================*/
  if(strcmp(argv[1],"--parameterfile") == 0)
  {
    global_io.params                 = (io_parameter_t) calloc(1,sizeof(io_parameter_struct_t));
    global_io.params->outfile_prefix = (char *) calloc(MAXSTRING,sizeof(char));
    global.a                         = 1;
    strcpy(global_io.params->outfile_prefix,"Cweb");
    write_parameterfile();
    exit(0);
  }
  else
  {
    strcpy(AMIGA_input, argv[1]);
  }
  
  
#if (!defined WITH_MPI)
  WRITEAHFLOGO(stderr);
#endif
  
  /* how much memory per node and particle for this particular run */
  global.bytes_node = sizeof(struct node);
  global.bytes_part = sizeof(struct particle);
  
#	ifdef WITH_MPI
  /* Initialize the MPI environment */
  common_initmpi(&argc, &argv);
#		ifdef MPI_TIMING
  global_mpi.start = MPI_Wtime();
#		endif
#	endif
  
#ifndef READ_GRIDDATA
  /*========================================================
   * startrun:    input the initial data from infile
   *========================================================*/
  timing.io       -= time(NULL);
  
  timing.startrun -= time(NULL);
  startrun((argc > 1) ? argv[1] : NULL, &timecounter, &timestep, &no_first_timestep);
  timing.startrun += time(NULL);
  
#ifdef VERBOSE
  fprintf(stderr,"Finished startrun()\n");
#endif
  
#ifdef DEBUG_STARTRUN
  /*===========================================================
   * DEBUG_STARTRUN:
   * we simply check if the particles have been read correctly
   *===========================================================*/
  {
    FILE *fpout;
    char outname[MAXSTRING];
    partptr cur_part;
    
#ifdef WITH_MPI
    sprintf(outname,"test-%d.ascii",global_mpi.rank);
#else
    sprintf(outname,"test.ascii");
#endif
    
    fpout = fopen(outname,"w");
    
    for(cur_part=global_info.fst_part; cur_part<(global_info.fst_part+global_info.no_part); cur_part++)
      fprintf(fpout,"%e %e %e\n",cur_part->pos[X]*simu.boxsize,cur_part->pos[Y]*simu.boxsize,cur_part->pos[Z]*simu.boxsize);
    
    fclose(fpout);
#ifdef WITH_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    //exit(0);
  }
#endif /* DEBUG_STARTRUN */
  
  /*==========================================================================================
   * DMfocus:
   *
   * keep dark matter particles for the analysis
   *
   *==========================================================================================*/
#if (defined DMfocus && defined MULTIMASS && defined GAS_PARTICLES)
  /* global_info.no_part
   * global_info.fst_part
   *                       => the no. of particles and relevant pointer for this CPU */
  timing.ptfocus -= time(NULL);
  {
    long unsigned no_part;
    partptr       fst_part, cur_part, new_part;
    int           ikeep;
    
    fprintf(stderr,"\n==================================================================\n");
    fprintf(stderr,"                          DMfocus\n");
    fprintf(stderr,"               ? ARE YOU SURE ABOUT THIS FLAG ?\n");
    fprintf(stderr,"==================================================================\n");
    fprintf(stderr,"AHF will now remove all particles whose type is 0 and 4, i.e. only keeping the dark matter particles\n");
    fprintf(stderr,"starting with %ld particles -> ",global_info.no_part);
    
    /* 1. count number of particles to keep */
    no_part  = 0;
    for(cur_part=global_info.fst_part; cur_part<(global_info.fst_part+global_info.no_part); cur_part++)
    {
      /* only keep the dm particles */
      if(cur_part->u < 0 && fabs(cur_part->u-PSTAR) > ZERO)
        no_part++;
    }
    
    /* allocate memory for new particles */
    fst_part = c_part(no_part);
    
    /* 2. remove all other particles */
    new_part = fst_part;
    for(cur_part=global_info.fst_part; cur_part<(global_info.fst_part+global_info.no_part); cur_part++)
    {
      ikeep = 0;
      
      /* only keep the dm particles */
      if(cur_part->u < 0 && fabs(cur_part->u-PSTAR) > ZERO)
        ikeep = 1;
      
      if(ikeep)
      {
        new_part->pos[X] = cur_part->pos[X];
        new_part->pos[Y] = cur_part->pos[Y];
        new_part->pos[Z] = cur_part->pos[Z];
        new_part->mom[X] = cur_part->mom[X];
        new_part->mom[Y] = cur_part->mom[Y];
        new_part->mom[Z] = cur_part->mom[Z];
        new_part->weight = cur_part->weight;
        new_part->u      = cur_part->u;
#if (!(defined AHF_NO_PARTICLES && defined AHFlean))
        new_part->id     = cur_part->id;
#endif
        new_part++;
      }
    }
    
    /* erase old particle list and store new one */
    free(global_info.fst_part);
    global_info.fst_part = fst_part;
    
    /* update global.no_part parameter */
    global_info.no_part  = no_part;
    fprintf(stderr,"ended with %ld particles\n\n",global_info.no_part);
  }
  timing.ptfocus += time(NULL);
#endif /* DMfocus */
  
  
  /*==========================================================================================
   * AHFptfocus:
   *
   * only use a certain type of particles ("pt") and focus ("focus") the AHF analysis on them
   *
   *==========================================================================================*/
#if (defined AHFptfocus && defined MULTIMASS && defined GAS_PARTICLES)
  /* global_info.no_part
   * global_info.fst_part
   *                       => the no. of particles and relevant pointer for this CPU */
  timing.ptfocus -= time(NULL);
  {
    long unsigned no_part;
    partptr       fst_part, cur_part, new_part;
    int           ikeep;
    
    fprintf(stderr,"\n==================================================================\n");
    fprintf(stderr,"                          AHFptfocus\n");
    fprintf(stderr,"               ? ARE YOU SURE ABOUT THIS FLAG ?\n");
    fprintf(stderr,"==================================================================\n");
    fprintf(stderr,"AHF will now remove all particles whose type is not %d\n",AHFptfocus);
    fprintf(stderr,"starting with %ld particles -> ",global_info.no_part);
    
    /* 1. count number of particles to keep */
    no_part  = 0;
    for(cur_part=global_info.fst_part; cur_part<(global_info.fst_part+global_info.no_part); cur_part++)
    {
      /* we only want to keep those particles with type AHFptfocus */
      if(AHFptfocus == 0)
      {
        if(cur_part->u >= AHFptfocus)
          no_part++;
      }
      else
      {
        if(fabs(cur_part->u+AHFptfocus) < ZERO)
          no_part++;
      }
      
      /* only keep the high-resolution particles */
      // if(cur_part->u >= 0 || cur_part->u == PDM || cur_part->u == PSTAR)
      //   no_part++;
      
    }
    
    /* allocate memory for new particles */
    fst_part = c_part(no_part);
    
    /* 2. remove all other particles */
    new_part = fst_part;
    for(cur_part=global_info.fst_part; cur_part<(global_info.fst_part+global_info.no_part); cur_part++)
    {
      ikeep = 0;
      
      /* we only want ot keep those particles with type AHFptfocus */
      if(AHFptfocus == 0)
      {
        if(cur_part->u >= AHFptfocus)
          ikeep = 1;
      }
      else
      {
        if(fabs(cur_part->u+AHFptfocus) < ZERO)
          ikeep = 1;
      }
      
      /* only keep the high-resolution particles */
      // if(cur_part->u >= 0 || cur_part->u == PDM || cur_part->u == PSTAR)
      //   ikeep = 1;
      
      if(ikeep)
      {
        new_part->pos[X] = cur_part->pos[X];
        new_part->pos[Y] = cur_part->pos[Y];
        new_part->pos[Z] = cur_part->pos[Z];
        new_part->mom[X] = cur_part->mom[X];
        new_part->mom[Y] = cur_part->mom[Y];
        new_part->mom[Z] = cur_part->mom[Z];
        new_part->weight = cur_part->weight;
        new_part->u      = cur_part->u;
#if (!(defined AHF_NO_PARTICLES && defined AHFlean))
        new_part->id     = cur_part->id;
#endif
        new_part++;
      }
    }
    
    /* erase old particle list and store new one */
    free(global_info.fst_part);
    global_info.fst_part = fst_part;
    
    /* update global.no_part parameter */
    global_info.no_part  = no_part;
    fprintf(stderr,"ended with %ld particles\n\n",global_info.no_part);
  }
  timing.ptfocus += time(NULL);
#endif /* AHFptfocus */
  
  
#ifdef AHFrfocus
  /*====================================================================
   * This is for focussing on a Sphere defined in param.h
   * This assumes that periodicity can be neglected for deciding
   * whether a particle is inside the selected sphere or not.
   *====================================================================*/
  timing.rfocus -= time(NULL);
  local_focusSphere();
#ifdef AHFrfocusRescale
  local_focusSphereRescale();
#endif
  timing.rfocus += time(NULL);
#endif
  
#		if (defined WITH_MPI && defined MPI_TIMING)
  global_mpi.stop = MPI_Wtime();
  io_logging_msg(global_io.log, INT32_C(1), "Startrun done in %fs", global_mpi.stop-global_mpi.start);
  global_mpi.start = global_mpi.stop;
#		endif
  
#ifdef VERBOSE
  fprintf(stderr,"Generating sfc keys ... ");
#endif
  
#ifdef WITH_MPI
  timing.loadbalance -= time(NULL);
  /* Sort the particles in a particle block structure */
  io_logging_section(global_io.log, "Initial Load-Balancing and Particle Distribution");
  io_logging_subsection(global_io.log, "Loadbalancing");
  loadbalance_update(global_io.log, global_info.loadbal, global_info.fst_part, global_info.no_part);
#			ifdef MPI_TIMING
  global_mpi.stop = MPI_Wtime();
  io_logging_msg(global_io.log, INT32_C(1), "Loadbalance done in %fs", global_mpi.stop-global_mpi.start);
  global_mpi.start = global_mpi.stop;
#			else
  io_logging_msg(global_io.log, INT32_C(1), "Loadbalance done.");
#			endif
  loadbalance_log(global_io.log, global_info.loadbal);
  timing.loadbalance += time(NULL);
#		else /* WITH _MPI */
  /* Generate the SFC keys for all particles */
  timing.sfckey -= time(NULL);
  //	for (uint64_t i=0; i<global_info.no_part; i++) {
  //		partptr part=global_info.fst_part+i;
  //		part->sfckey = sfc_curve_calcKey(global_info.ctype,
  //		                                 (double)(part->pos[0]),
  //		                                 (double)(part->pos[1]),
  //		                                 (double)(part->pos[2]),
  //		                                 BITS_PER_DIMENSION);
  //	}
  //	/* Sorting all particles to have fast access later on */
  //	qsort(global_info.fst_part,
  //	      global_info.no_part,
  //	      sizeof(part),
  //	      &cmp_sfckey_part);
  timing.sfckey += time(NULL);
#endif /* WITH_MPI*/
  
#ifdef WITH_MPI
  timing.distribution -= time(NULL);
  
  /* Do a first sort of the particles, required for distributing */
  io_logging_subsection(global_io.log, "Sorting particles");
  qsort(global_info.fst_part,
        global_info.no_part,
        sizeof(part),
        &cmp_sfckey_part);
#			ifdef MPI_TIMING
  global_mpi.stop = MPI_Wtime();
  io_logging_msg(global_io.log, INT32_C(1), "Sorting done in %fs", global_mpi.stop-global_mpi.start);
  global_mpi.start = global_mpi.stop;
#			else
  io_logging_msg(global_io.log, INT32_C(1), "Sorting done.");
#			endif
  
  /* Distribute the particles */
  io_logging_subsection(global_io.log, "Distributing particles");
  io_logging_msg(global_io.log, INT32_C(0), "Currently having %"PRIu64" particles.", global_info.no_part);
  comm_dist_part(global_io.log,
                 &(global_info.fst_part),
                 &(global_info.no_part),
                 global_info.loadbal);
#			ifdef MPI_TIMING
  global_mpi.stop = MPI_Wtime();
  io_logging_msg(global_io.log, INT32_C(1), "Distributing done in %fs", global_mpi.stop-global_mpi.start);
  global_mpi.start = global_mpi.stop;
#			else
  io_logging_msg(global_io.log, INT32_C(1), "Distributing done.");
#			endif
  io_logging_msg(global_io.log, INT32_C(0), "Having %"PRIu64" particles!", global_info.no_part);
  
  /* Do the AHF distribution*/
  io_logging_subsection(global_io.log, "AHF distribution (duplicating)");
  newparts = comm_dist_part_ahf(global_io.log,
                                &(global_info.fst_part),
                                &(global_info.no_part),
                                global_info.loadbal);
  io_logging_msg(global_io.log, INT32_C(0), "Received %"PRIu64" new particles.", newparts);
  /* We need to sort the particles again */
  qsort(global_info.fst_part, global_info.no_part, sizeof(part), &cmp_sfckey_part);
#				ifdef MPI_TIMING
  global_mpi.stop = MPI_Wtime();
  io_logging_msg(global_io.log, INT32_C(1), "AHF distribution done in %fs", global_mpi.stop-global_mpi.start);
  global_mpi.start = global_mpi.stop;
#				else
  io_logging_msg(global_io.log, INT32_C(1), "AHF distribution done.");
#				endif
  
  timing.distribution += time(NULL);
#endif /* WITH_MPI */
#ifdef VERBOSE
  fprintf(stderr,"done\n");
#endif
  
  
#if AHFstep_split_only
  /*====================================================================
   * we only split the data using the SFC and
   * dump the data into multilpe files
   *====================================================================*/
  {
    io_file_t dumpf;
    io_file_strg_struct_t strg;
    char *fname;
    
    /* Start tge section */
    io_logging_section(global_io.log, "Dumping AHF chunk to file");
    
    /* First generate the filename */
    fname = (char *)malloc( sizeof(char) *( strlen(global_io.params->outfile_prefix)+30));
    if (fname == NULL) {
      io_logging_memfatal(global_io.log, "filename string");
      common_terminate(EXIT_FAILURE);
    }
    sprintf(fname, "%s.chunk.%04i.dump", global_io.params->outfile_prefix, global_mpi.rank);
    io_logging_msg(global_io.log, UINT32_C(0), "Used filename: %s", fname);
    
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    
    /* Assign particles to structure */
    strg.posx.val = (void *)(global_info.fst_part->pos);
    strg.posx.stride =   (char *)((global_info.fst_part+1)->pos  ) - (char *)(global_info.fst_part->pos);
    strg.posy.val = (void *)(global_info.fst_part->pos+1);
    strg.posy.stride =   (char *)((global_info.fst_part+1)->pos+1) - (char *)(global_info.fst_part->pos+1);
    strg.posz.val = (void *)(global_info.fst_part->pos+2);
    strg.posz.stride =   (char *)((global_info.fst_part+1)->pos+2) - (char *)(global_info.fst_part->pos+2);
    strg.momx.val = (void *)(global_info.fst_part->mom);
    strg.momx.stride =   (char *)((global_info.fst_part+1)->mom  ) - (char *)(global_info.fst_part->mom);
    strg.momy.val = (void *)(global_info.fst_part->mom+1);
    strg.momy.stride =   (char *)((global_info.fst_part+1)->mom+1) - (char *)(global_info.fst_part->mom+1);
    strg.momz.val = (void *)(global_info.fst_part->mom+2);
    strg.momz.stride =   (char *)((global_info.fst_part+1)->mom+2) - (char *)(global_info.fst_part->mom+2);
#	ifdef MULTIMASS
    strg.weight.val = (void *)&(global_info.fst_part->weight);
    strg.weight.stride =   (char *)&((global_info.fst_part+1)->weight) - (char *)&(global_info.fst_part->weight);
#	else
    strg.weight.val = NULL;
    strg.weight.stride = (ptrdiff_t)0;
#	endif /* MULTIMASS */
#	ifdef GAS_PARTICLES
    strg.u.val = (void *)&(global_info.fst_part->u);
    strg.u.stride =   (char *)&((global_info.fst_part+1)->u) - (char *)&(global_info.fst_part->u);
#	else
    strg.u.val = NULL;
    strg.u.stride = (ptrdiff_t)0;
#	endif /* GAS_PARTICLE */
#	if (defined AHFlean && defined AHF_NO_PARTICLES)
    strg.id.val = NULL;
    strg.id.stride = (ptrdiff_t)0;
#	else
    strg.id.val = &(global_info.fst_part->id);
    strg.id.stride =   (char *)&((global_info.fst_part+1)->id) - (char *)&(global_info.fst_part->id);
#	endif
    strg.bytes_float = sizeof(global_info.fst_part->pos[0]);
#	if (defined AHFlean && defined AHF_NO_PARTICLES)
    strg.bytes_int = 0;
#	else
    strg.bytes_int = sizeof(global_info.fst_part->id);
#	endif
    
    /* Open the dump file now */
    dumpf = io_file_open(global_io.log, fname, IO_FILE_ARES, IO_FILE_UNKOWN_SWAPPING, IO_FILE_WRITE, 0);
    
    /* Write the particles */
    io_file_writepart(global_io.log, dumpf, 0, global_info.no_part, strg);
    
    /* Set the header values */
    ((io_ares_t)dumpf)->header->no_part = (uint64_t)simu.no_part;
    ((io_ares_t)dumpf)->header->no_species = UINT64_C(0);
    ((io_ares_t)dumpf)->header->no_vpart = simu.no_vpart;
    ((io_ares_t)dumpf)->header->boxsize = simu.boxsize;
    ((io_ares_t)dumpf)->header->omega0 = simu.omega0;
    ((io_ares_t)dumpf)->header->lambda0 = simu.lambda0;
    ((io_ares_t)dumpf)->header->pmass = simu.pmass;
    ((io_ares_t)dumpf)->header->minweight = simu.min_weight;
    ((io_ares_t)dumpf)->header->maxweight = simu.max_weight;
    ((io_ares_t)dumpf)->header->a_initial = simu.a_initial;
    ((io_ares_t)dumpf)->header->a_current = global.a;
    ((io_ares_t)dumpf)->header->timestep = timestep;
    ((io_ares_t)dumpf)->header->minkey = global_info.loadbal->fstkey[global_mpi.rank];
    ((io_ares_t)dumpf)->header->maxkey = global_info.loadbal->lstkey[global_mpi.rank];
    ((io_ares_t)dumpf)->header->lb_level = global_info.loadbal->level;
    ((io_ares_t)dumpf)->header->rank = global_mpi.rank;
    ((io_ares_t)dumpf)->header->size = global_mpi.size;
    
    /* Log the file */
    io_file_log(global_io.log, dumpf);
    
    /* Close the file and clean up*/
    io_file_close(global_io.log, &dumpf);
    free(fname);
  }
  common_terminate(EXIT_SUCCESS);
#endif /*  AHFstep_split_only */
  
  
#ifdef AHF_DUMP_AFTER_READ_TO_ASCII
  /*====================================================================
   * write an ASCII file of the data just read
   *====================================================================*/
  {
    FILE *dumpf;
    char *fname;
    
    /* First generate the filename */
    fname = (char *)malloc( sizeof(char) *( strlen(global_io.params->outfile_prefix)+35));
    if (fname == NULL) {
      io_logging_memfatal(global_io.log, "filename string");
      common_terminate(EXIT_FAILURE);
    }
#ifdef WITH_MPI
    sprintf(fname, "%s.chunk.%04i.ascii", global_io.params->outfile_prefix, global_mpi.rank);
#else
    sprintf(fname, "%s.DUMP.ascii", global_io.params->outfile_prefix);
#endif
    io_logging_msg(global_io.log, UINT32_C(0), "Used filename: %s", fname);
    fflush(NULL);
#ifdef WITH_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
    dumpf = fopen(fname, "w");
    fprintf(dumpf, "# x y z  vx vy vz  ID\n");
    for (uint64_t i=0L; i<global_info.no_part; i++) {
#ifdef GAS_PARTICLES
      fprintf(dumpf, "%15e %15e %15e   %15e %15e %15e  %15e  %16.8g %lu\n",
              global_info.fst_part[i].pos[0],
              global_info.fst_part[i].pos[1],
              global_info.fst_part[i].pos[2],
              global_info.fst_part[i].mom[0],
              global_info.fst_part[i].mom[1],
              global_info.fst_part[i].mom[2],
              global_info.fst_part[i].u,
              global_info.fst_part[i].weight,
              (unsigned long)global_info.fst_part[i].id);
#else
      fprintf(dumpf, "%15e %15e %15e   %15e %15e %15e  %lu\n",
              global_info.fst_part[i].pos[0],
              global_info.fst_part[i].pos[1],
              global_info.fst_part[i].pos[2],
              global_info.fst_part[i].mom[0],
              global_info.fst_part[i].mom[1],
              global_info.fst_part[i].mom[2],
              (unsigned long)global_info.fst_part[i].id);
#endif
    }
    fclose(dumpf);
    common_terminate(EXIT_SUCCESS);
  }
#endif /* AHF_DUMP_AFTER_READ_TO_ASCII */
  
  
#ifdef WITH_MPI
  loadbalance_minimalMemory(global_io.log, global_info.loadbal);
#endif
  io_logging_msg(global_io.log, INT32_C(5), "main:  running with %" PRIu64 " particles", global_info.no_part);
  io_logging_part(global_io.log, "Handing over logging to Cweb");
  
  timing.io       += time(NULL);
  
  
  
  /*=====================================================================
   * at this point we completely read in the data file
   * and are ready to proceed with generating the
   * grid hierarchy or what else we plan to do...
   *=====================================================================*/
  //fprintf(stderr,"\n");
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#else // READ_GRIDDATA
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // we need to set a lot of internally used values that are usually initialized in startrun()
  fprintf(stderr,"Not reading the simulation, but instead a previously calculated grid file\n");
  fprintf(stderr,"Initialiazing parameters from %s...",argv[1]);
  
  local_startrunParams(argv[1]);
  
  // these values are found in the Cweb.input file:
  simu.NGRID_DOM     = global_io.params->NGRID_DOM;
  simu.NGRID_MIN     = simu.NGRID_DOM;
  simu.Nth_dom       = global_io.params->Nth_dom;
  simu.Nth_ref       = global_io.params->Nth_ref;
  simu.lb_level      = global_io.params->lb_level;
  simu.MaxGatherRad  = global_io.params->MaxGatherRad;
  simu.UserDvir      = global_io.params->UserDvir;
  simu.UseRhoBack    = global_io.params->UseRhoBack;
  simu.NGRID_MAX     = global_io.params->NGRID_MAX;
  simu.AHF_MINPART   = global_io.params->AHF_MINPART;
  simu.RSMOOTH       = global_io.params->RSMOOTH;
  simu.GADGET_m2Msunh= global_io.params->GADGET_m2Msunh;
  simu.GADGET_l2Mpch = global_io.params->GADGET_l2Mpch;
  
#endif // READ_GRIDDATA
  
  
  /*====================================================================
   *  GENERATE THE FULL BLOWN AMR HIERARCHY AND ORGANIZE IT INTO A TREE
   *====================================================================*/
  
  
  /*=====================================================================
   * generate the domain grids: simu.NGRID_MIN^3, ...., simu.NGRID_DOM^3
   *=====================================================================*/
#ifdef VERBOSE
  fprintf(stderr,"Generating domain grid ... ");
#endif
  timing.gendomgrids -= time(NULL);
  grid_list = gen_domgrids(&no_grids);
  timing.gendomgrids += time(NULL);
#ifdef VERBOSE
  fprintf(stderr,"done\n");
#endif
  
#ifndef READ_GRIDDATA
  /*=====================================================================
   * build initial linked list
   *=====================================================================*/
#ifdef VERBOSE
  fprintf(stderr,"Building linked-list on domain grid ... ");
#endif
  timing.ll -= time(NULL);
  ll(global_info.no_part, global_info.fst_part, global.dom_grid);
  global.fst_part = global_info.fst_part;
  global.no_part  = global_info.no_part;
  timing.ll += time(NULL);
#ifdef VERBOSE
  fprintf(stderr,"done\n");
#endif
  
  /*================================================================
   * assign particles to the domain grid with simu.NGRID_DOM^3 nodes
   *================================================================*/
#ifdef VERBOSE
  fprintf(stderr,"Assigning particles to the domain grid ... ");
#endif
  zero_dens(global.dom_grid);
  assign_npart(global.dom_grid);
#ifdef VERBOSE
  fprintf(stderr,"done\n");
#endif
#endif // READ_GRIDDATA
  
  /*================================================================
   * initialize some counters
   *================================================================*/
  no_timestep         = no_first_timestep+1;  /* count total number of integration steps */
  global.total_time   = 0.;                   /* cumulative total time for simulation    */
  global.output_count = 0;                    /* count the number of outputs             */
  
  /* make *current* time step available to AHF/etc. routines */
  global.no_timestep = no_first_timestep;
  
#ifndef READ_GRIDDATA
  /*=========================================================================================
   * recursively call gen_AMRhierarchy() to generate the AMR hierarchy...
   *=========================================================================================*/
  // generate AMR hierarchy
  // CAREFUL: this sets various other parameters and values associated with grids that are required...
  // ...therefore, do not remove this call even when dealing with only one domain grid!
#ifdef VERBOSE
  fprintf(stderr,"Generate AMR grids (from LgridDomain=%d to LgridMax=%d):\n",simu.NGRID_DOM,simu.NGRID_MAX);
#endif
  global.fst_cycle = TRUE;
  
  gen_AMRhierarchy(&grid_list, &no_grids);
#endif // READ_GRIDDATA
  
  /*=========================================================================================
   * do we need the potential?
   *=========================================================================================*/
#ifdef PWEB
#ifndef READ_GRIDDATA
  /* solve for potential on domain grid(s) */
  fprintf(stderr,"PWEB:   solving for potential on domain grid %lu ... ",global.dom_grid->l1dim);
  solve_gravity(grid_list, 0);
  fprintf(stderr,"done\n");
  
  /* TODO: solve for potential on refinement grids */
  
#endif // READ_GRIDDATA
#endif /* PWEB */
  
  
  /*=========================================================================================
   * at this point we are able to read griddata into grid_list
   * (only works for the serial version as a single files is required)
   *=========================================================================================*/
#ifdef READ_GRIDDATA
  // now we can read the grid data
  read_griddata(&grid_list, &no_grids, global_io.params->icfile_name);
  //  write_griddata(&grid_list, &no_grids, global_io.params->icfile_name);
  //  exit(0);
  
  // we need to set some parameters that are usually initialized during the reading of the simulation snapshot
  global.a              = 1./(1.+global.z);
  global.super_t        = calc_super_t(global.a);
  global.t              = calc_t(global.a);
  simu.no_vpart         = simu.no_part;
  simu.no_species       = 1;
  simu.a_initial        = 0.001;
  simu.a_final          = 1.0;
  create_timeline(simu.a_initial, simu.a_final, &simu.timeline);
  simu.z_initial        = (double)1.0/simu.a_initial - (double)1.0;
  simu.z_final          = (double)1.0/simu.a_final   - (double)1.0;
  simu.multi_mass       = 0.0;
  simu.min_weight       = 1.0;
  simu.max_weight       = 1.0;
  simu.t_unit           = 1.0/H0;
  simu.np_limit         = TRUE;
  simu.mean_dens        = (double) 1.0;
  simu.double_precision = 0;
  simu.mmfocus          = 0;
  simu.hydro            = 0;
  simu.magneto          = 0;
  simu.SHIFT            = ((double)0.5000000/(double) simu.NGRID_DOM);
  simu.super_t_initial  = calc_super_t(simu.a_initial);
  simu.super_t_final    = calc_super_t(simu.a_final);
  simu.t_initial        = calc_t(simu.a_initial);
  simu.t_final          = calc_t(simu.a_final);
  simu.gamma            = 0.0;
  simu.omegab           = 0.0;
  simu.omegaDM          = simu.omega0;
  simu.FourPiG          = 1.5*simu.omega0;
  simu.f_b              = 0.0;
  simu.H_frac           = 0.0;
  simu.T_init           = 0.0;
  simu.B_init           = 0.0;
  simu.e_init           = 0.0;
  simu.no_halos         = 0;
  simu.med_weight       = simu.max_weight;
  simu.l_unit           = 0.0;
  simu.m_unit           = 0.0;
  simu.no_gas           = 0;
  simu.no_stars         = 0;

#endif
  
  
  /*=========================================================================================
   * or we might like to write the standard TSC fields to file (before smoothing)
   * (only works for the serial version as a single files is written)
   *=========================================================================================*/
#ifdef WRITE_GRIDDATA
  write_griddata(&grid_list, &no_grids, global_io.params->icfile_name);
#ifdef WRITE_GRIDDATA_TERMINATE
  fprintf(stderr,"mesh data written, aborting now\n");
  exit(0);
#endif
#endif
  
  
  
  
  
  
  
  
  
  
  
  
  
  /*=========================================================================================
   *=========================================================================================
   *=========================================================================================
   *=========================================================================================
   *               eventually perform Cweb analysis on AMR hierarchy
   *=========================================================================================*
   *=========================================================================================*
   *=========================================================================================*
   *=========================================================================================*/
  ahf.time -= time(NULL);
  
  x_fac   = simu.boxsize;
  v_fac   = simu.boxsize/simu.t_unit/global.a;
  phi_fac = pow2(H0);  //Grav * simu.pmass / (simu.boxsize * global.a);
  pot_fac = pow2(H0*simu.boxsize);
  rho_fac = simu.pmass / pow3(simu.boxsize);
#ifdef HUBBLE_Z
  Hz      = calc_Hubble(global.a);
#else
  Hz      = H0;
#endif
  Hz2     = pow2(Hz);
  Rsmooth = simu.RSMOOTH;
  
  fprintf(stderr,"Conversion factors to physical coordinates: x_fac = %lf,  v_fac = %lf\n",x_fac,v_fac);
  
  /* loop over all nodes on all grids */
  for(cur_grid=global.dom_grid;cur_grid<global.dom_grid+no_grids;cur_grid++)
  {
    
    // we are smoothing all densities with a Gaussian filter (which only makes sense when the filter is larger than the grid-spacing)
    if(Rsmooth > 2*simu.boxsize/(double)cur_grid->l1dim) {
#ifdef WITH_MPI
      sprintf(outfile,"%s.%05ld.Rs=%4.2lf.MPIrank%04d.Cweb",global_io.params->outfile_prefix,cur_grid->l1dim,Rsmooth,global_mpi.rank);
#else
      sprintf(outfile,"%s.%05ld.Rs=%4.2lf.Cweb",global_io.params->outfile_prefix,cur_grid->l1dim,Rsmooth);
#endif
      fprintf(stderr,"Gaussian smoothing %ld grid (Rsmooth=%f [Mpc/h])... ",cur_grid->l1dim,Rsmooth);
      Gaussian_smooth_gridFFT(cur_grid, Rsmooth);
      fprintf(stderr,"done\n");
    }
    else {
      fprintf(stderr,"  -> you are trying to smooth on a scale (Rsmooth=%lf [Mpc/h]) smaller than the 2x grid spacing (B/L=%lf [Mpc/h]): will not smooth!\n",Rsmooth,simu.boxsize/(double)cur_grid->l1dim);
#ifdef WITH_MPI
      sprintf(outfile,"%s.%05ld.Rs=0.%04d.Cweb",global_io.params->outfile_prefix,cur_grid->l1dim,global_mpi.rank);
#else
      sprintf(outfile,"%s.%05ld.Rs=0.Cweb",global_io.params->outfile_prefix,cur_grid->l1dim);
#endif
    }
    
    /* count number of nodes written to file */
    Nnodes_written   = 0l;
    
#ifdef WRITE_ASCII
    strcat(outfile,"-ascii");
    fpout=fopen(outfile,"w");
    fprintf(fpout,"#x(1) y(2) z(3) dens(4) Vx(5) Vy(6) Vz(7) Wx(8) Wy(9) Wz(10) lambda1(11) lambda2(12) lambda3(13) e1(14-16) e2(17-19) e3(20-22)\n");
#else
    fpout=fopen(outfile,"wb");
    fwrite(&one, sizeof(int32_t), 1, fpout);
    tmp_uint64  = (uint64_t)Nnodes_written;  // placeholder for the time being, will be overwritten at end
    FWRITE_TMP_UINT64;
    tmp_uint64  = (uint64_t)cur_grid->l1dim;
    FWRITE_TMP_UINT64;
    tmp_float   = simu.boxsize;
    FWRITE_TMP_FLOAT;
    // we could/should also add tags that indicate whether the file contains DWEB/PWEB or not...
#endif
    
    fprintf(stderr,"Calculating Vweb on %8ld grid (writing results to %s) ... ",cur_grid->l1dim,outfile);
#ifdef DWEB
    fprintf(stderr,"\n    -> getting Dweb at the same time ... ");
#endif
#ifdef PWEB
    fprintf(stderr,"\n    -> getting Pweb at the same time ... ");
#endif

    /* for centred finite differences (in physical units) */
    dx    = cur_grid->spacing*x_fac;
    dy    = cur_grid->spacing*x_fac;
    dz    = cur_grid->spacing*x_fac;
    twodx = 2.0*dx;
    twody = 2.0*dy;
    twodz = 2.0*dz;

    Nknots           = 0;
    Nfilaments       = 0;
    Nsheets          = 0;
    Nvoids           = 0;
#ifdef LAMBDA_TH
    LAMBDA_THRESHOLD = LAMBDA_TH;
#else
    LAMBDA_THRESHOLD = 0.0;
#endif
    
    for(cur_pquad=cur_grid->pquad; cur_pquad!=NULL; cur_pquad=cur_pquad->next)
    {
      for(cur_cquad = cur_pquad->loc, z = cur_pquad->z; cur_cquad < cur_pquad->loc + cur_pquad->length; cur_cquad++, z++)
      {
        for(icur_cquad=cur_cquad; icur_cquad!=NULL; icur_cquad=icur_cquad->next)
        {
          for(cur_nquad = icur_cquad->loc, y = icur_cquad->y; cur_nquad < icur_cquad->loc + icur_cquad->length; cur_nquad++, y++)
          {
            for(icur_nquad=cur_nquad; icur_nquad!=NULL; icur_nquad=icur_nquad->next)
            {
              for(cur_node = icur_nquad->loc, x = icur_nquad->x; cur_node < icur_nquad->loc + icur_nquad->length; cur_node++, x++)
              {
#ifdef WITH_MPI
                notboundary = check_for_MPIboundary(((double)x+0.5)/(double)cur_grid->l1dim, ((double)y+0.5)/(double)cur_grid->l1dim, ((double)z+0.5)/(double)cur_grid->l1dim);
#else
                notboundary = TRUE;
#endif
                /* 1. only consider this cell if not within the MPI boundary zone
                 *    -> otherwise it will be written by some other MPI task already */
                if(notboundary == TRUE)
                {
                   // simply reset tsc_nodes[][][] -> will be filled by get_TSCnodes() below....
                   for(i = 0; i <= 2; i++)
                      for(j = 0; j <= 2; j++)
                         for(k = 0; k <= 2; k++){
                            tsc_nodes[i][j][k]=NULL;
                         }
                   tsc_nodes[1][1][1] = cur_node;
                   get_TSCnodes(cur_grid, cur_pquad, icur_cquad, icur_nquad, tsc_nodes, &z, &y, &x);
                   
                   /* 2. only deal with current cell if the tensors can be properly calculated, i.e. the cell and its neighbours have (dens+mean_dens)>0 */
                   //fprintf(stderr,"\n dens=%f (ll=%ld)",cur_node->dens+simu.mean_dens,cur_node->ll);
                   if(check_for_zerodens(tsc_nodes, CWEB_DENSTHRESHOLD) == 1)
                   {
//                      /* 3. only use nodes that have neighbours for the finite difference calculation of dVi/drj
//                       *    -> otherwise ignore as it is a boundary cell on a refinement level */
//                      haveneighbours = test_tsc(tsc_nodes);
//                      if(haveneighbours == TRUE)            // this test only makes sense to AMD grids, which are currently not implemented...
                      {
                         
                         /* we need to divide the momentum density by the actual mass density to get the velocity */
                         Vx_x2 = tsc_nodes[1][1][2]->densV[X]/(tsc_nodes[1][1][2]->dens+simu.mean_dens);
                         Vx_x0 = tsc_nodes[1][1][0]->densV[X]/(tsc_nodes[1][1][0]->dens+simu.mean_dens);
                         Vx_y2 = tsc_nodes[1][2][1]->densV[X]/(tsc_nodes[1][2][1]->dens+simu.mean_dens);
                         Vx_y0 = tsc_nodes[1][0][1]->densV[X]/(tsc_nodes[1][0][1]->dens+simu.mean_dens);
                         Vx_z2 = tsc_nodes[2][1][1]->densV[X]/(tsc_nodes[2][1][1]->dens+simu.mean_dens);
                         Vx_z0 = tsc_nodes[0][1][1]->densV[X]/(tsc_nodes[0][1][1]->dens+simu.mean_dens);
                         
                         Vy_x2 = tsc_nodes[1][1][2]->densV[Y]/(tsc_nodes[1][1][2]->dens+simu.mean_dens);
                         Vy_x0 = tsc_nodes[1][1][0]->densV[Y]/(tsc_nodes[1][1][0]->dens+simu.mean_dens);
                         Vy_y2 = tsc_nodes[1][2][1]->densV[Y]/(tsc_nodes[1][2][1]->dens+simu.mean_dens);
                         Vy_y0 = tsc_nodes[1][0][1]->densV[Y]/(tsc_nodes[1][0][1]->dens+simu.mean_dens);
                         Vy_z2 = tsc_nodes[2][1][1]->densV[Y]/(tsc_nodes[2][1][1]->dens+simu.mean_dens);
                         Vy_z0 = tsc_nodes[0][1][1]->densV[Y]/(tsc_nodes[0][1][1]->dens+simu.mean_dens);
                         
                         Vz_x2 = tsc_nodes[1][1][2]->densV[Z]/(tsc_nodes[1][1][2]->dens+simu.mean_dens);
                         Vz_x0 = tsc_nodes[1][1][0]->densV[Z]/(tsc_nodes[1][1][0]->dens+simu.mean_dens);
                         Vz_y2 = tsc_nodes[1][2][1]->densV[Z]/(tsc_nodes[1][2][1]->dens+simu.mean_dens);
                         Vz_y0 = tsc_nodes[1][0][1]->densV[Z]/(tsc_nodes[1][0][1]->dens+simu.mean_dens);
                         Vz_z2 = tsc_nodes[2][1][1]->densV[Z]/(tsc_nodes[2][1][1]->dens+simu.mean_dens);
                         Vz_z0 = tsc_nodes[0][1][1]->densV[Z]/(tsc_nodes[0][1][1]->dens+simu.mean_dens);
                         
                         /* finite differences (in physical units) */
                         dVxdx = v_fac*(Vx_x2-Vx_x0)/twodx;
                         dVxdy = v_fac*(Vx_y2-Vx_y0)/twody;
                         dVxdz = v_fac*(Vx_z2-Vx_z0)/twodz;
                         
                         dVydx = v_fac*(Vy_x2-Vy_x0)/twodx;
                         dVydy = v_fac*(Vy_y2-Vy_y0)/twody;
                         dVydz = v_fac*(Vy_z2-Vy_z0)/twodz;
                         
                         dVzdx = v_fac*(Vz_x2-Vz_x0)/twodx;
                         dVzdy = v_fac*(Vz_y2-Vz_y0)/twody;
                         dVzdz = v_fac*(Vz_z2-Vz_z0)/twodz;
                         // TODO: not sure about the correct direction/sign here [2]-[0] or [0]-[2] !?
                         
                         /* vorticity */
                         vorticity[X] = dVzdy-dVydz;
                         vorticity[Y] = dVxdz-dVzdx;
                         vorticity[Z] = dVydx-dVxdy;
                         
                         /* local shear tensor (units according to Hoffman et al. 2012) */
                         local_shear[X][X] = -0.5/Hz*( (dVxdx) + (dVxdx) );  // dVx/dx + dVx/dx
                         local_shear[Y][Y] = -0.5/Hz*( (dVydy) + (dVydy) );  // dVy/dy + dVy/dy
                         local_shear[Z][Z] = -0.5/Hz*( (dVzdz) + (dVzdz) );  // dVz/dz + dVz/dz
                         
                         local_shear[X][Y] = -0.5/Hz*( (dVxdy) + (dVydx) );  // dVx/dy + dVy/dx
                         local_shear[X][Z] = -0.5/Hz*( (dVxdz) + (dVzdx) );  // dVx/dz + dVz/dx
                         local_shear[Y][Z] = -0.5/Hz*( (dVydz) + (dVzdy) );  // dVy/dz + dVz/dy
                         
                         local_shear[Y][X] = local_shear[X][Y];
                         local_shear[Z][X] = local_shear[X][Z]; // symmetric tensor
                         local_shear[Z][Y] = local_shear[Y][Z];
                         
                         /* get eigenvalues and eigenvectors (leave local_shear[][] unspoiled for debugging) */
                         itensor[X][X] = local_shear[X][X];
                         itensor[Y][Y] = local_shear[Y][Y];
                         itensor[Z][Z] = local_shear[Z][Z];
                         itensor[X][Y] = local_shear[X][Y];
                         itensor[X][Z] = local_shear[X][Z];
                         itensor[Y][Z] = local_shear[Y][Z];
                         itensor[Y][X] = local_shear[Y][X];
                         itensor[Z][X] = local_shear[Z][X];
                         itensor[Z][Y] = local_shear[Z][Y];
                         lambda1 = -1000.0;
                         lambda2 = -1000.0;
                         lambda3 = -1000.0;
                         converged = get_axes(itensor, &lambda1, &lambda2, &lambda3);
                         // NOTE: the eigenvalues and -vectors are ordered upon return axis1>axis2>axis3
                         
#ifdef IGNORE_JACOBI_NONCONVERGENCE
                         if(converged == 0){
                            fprintf(stderr," ---> (Vweb) at this grid position x=%ld y=%ld z=%ld dens=%g densVx=%g densVy=%g densVz=%g\n",
                                    x,y,z,cur_node->dens,cur_node->densV[X],cur_node->densV[Y],cur_node->densV[Z]);
                            lambda1 = lambda2 = lambda3 = -3.1415;
                            itensor[X][X] = itensor[X][Y] = itensor[X][Z] = -3.1415;
                            itensor[X][Y] = itensor[Y][Y] = itensor[Y][Z] = -3.1415;
                            itensor[X][Z] = itensor[Y][Z] = itensor[Z][Z] = -3.1415;
                         }
#endif
                         
                         /* count web elements */
                         if(lambda3 > LAMBDA_THRESHOLD)                               Nknots++;
                         if(lambda2 > LAMBDA_THRESHOLD && lambda3 < LAMBDA_THRESHOLD) Nfilaments++;
                         if(lambda1 > LAMBDA_THRESHOLD && lambda2 < LAMBDA_THRESHOLD) Nsheets++;
                         if(lambda1 < LAMBDA_THRESHOLD)                               Nvoids++;
                         
#if defined(DWEB) || defined(PWEB)
                         /* Now we use mhd nodes to make proper derivations */
                         for(i = 0; i <= 4; i++)
                            for(j = 0; j <= 4; j++)
                               for(k = 0; k <= 4; k++){
                                  mhd_nodes[i][j][k]=NULL;
                               }
                         get_MHDnodes(cur_grid, cur_pquad, z, y, x, mhd_nodes);
//                         haveneighbours = test_mhd(mhd_nodes);
                         haveneighbours = TRUE; // this test only makes sense for AMR grids, which are currently not implemented
#endif

#ifdef DWEB  // doese DWEB actually make sense at all!?!?!
                         for(i = 0; i <= 2; i++)
                            for(j = 0; j <= 2; j++)
                               dtensor[i][j] = 0;
                         dambda1 = -1000.0;
                         dambda2 = -1000.0;
                         dambda3 = -1000.0;

                         converged = 0;

#ifdef DWEB_AK // my own version --> not accurate enough!
                         /* Old version only using the direct neighbours only, no fancy derivatives... */
                         dtensor[0][0] = ((tsc_nodes[1][1][2]->dens - tsc_nodes[1][1][1]->dens) - (tsc_nodes[1][1][1]->dens - tsc_nodes[1][1][0]->dens))/dx/dx;
                         dtensor[1][1] = ((tsc_nodes[1][2][1]->dens - tsc_nodes[1][1][1]->dens) - (tsc_nodes[1][1][1]->dens - tsc_nodes[1][0][1]->dens))/dy/dy; // this is 2nd order accurate
                         dtensor[2][2] = ((tsc_nodes[2][1][1]->dens - tsc_nodes[1][1][1]->dens) - (tsc_nodes[1][1][1]->dens - tsc_nodes[0][1][1]->dens))/dz/dz;
                         
                         /* directly use the second order derivative to calculate d^2 rho/dx/dy. */
                         dtensor[0][1] = ((tsc_nodes[1][2][2]->dens - tsc_nodes[1][0][2]->dens)/twody - (tsc_nodes[1][2][0]->dens - tsc_nodes[1][0][0]->dens)/twody)/twodx;
                         dtensor[0][2] = ((tsc_nodes[2][1][2]->dens - tsc_nodes[0][1][2]->dens)/twodz - (tsc_nodes[2][1][0]->dens - tsc_nodes[0][1][0]->dens)/twody)/twodx; // this is 2nd order accurate (A.10 in https://onlinelibrary.wiley.com/doi/pdf/10.1002/9781119083405.app1)
                         dtensor[1][2] = ((tsc_nodes[2][2][1]->dens - tsc_nodes[0][2][1]->dens)/twodz - (tsc_nodes[2][0][1]->dens - tsc_nodes[0][0][1]->dens)/twody)/twody;
                         
                         /* using the mean of first order derivative to calculate d^2 rho/dx/dy. */
                         //                       dtensor[0][1] = ((tsc_nodes[1][2][1]->dens - tsc_nodes[0][2][1]->dens) - (tsc_nodes[1][0][1]->dens - tsc_nodes[0][0][1]->dens) +
                         //                                        (tsc_nodes[2][2][1]->dens - tsc_nodes[1][2][1]->dens) - (tsc_nodes[2][0][1]->dens - tsc_nodes[1][0][1]->dens))/dx/dy/4.0;
                         //                       dtensor[0][2] = ((tsc_nodes[1][1][2]->dens - tsc_nodes[0][1][2]->dens) - (tsc_nodes[1][1][0]->dens - tsc_nodes[0][1][0]->dens) +
                         //                                        (tsc_nodes[2][1][2]->dens - tsc_nodes[1][1][2]->dens) - (tsc_nodes[2][1][0]->dens - tsc_nodes[1][1][0]->dens))/dx/dz/4.0;
                         //                       dtensor[1][2] = ((tsc_nodes[1][1][2]->dens - tsc_nodes[1][0][2]->dens) - (tsc_nodes[1][1][0]->dens - tsc_nodes[1][0][0]->dens) +
                         //                                        (tsc_nodes[1][2][2]->dens - tsc_nodes[1][1][2]->dens) - (tsc_nodes[1][2][0]->dens - tsc_nodes[1][1][0]->dens))/dx/dy/4.0;
                         
                         dtensor[1][0] = dtensor[0][1];
                         dtensor[2][0] = dtensor[0][2]; // the Hessian matrix is symmetric
                         dtensor[2][1] = dtensor[1][2];
                         
                         // add the missing conversion to physical units (note, dx/dy/dz are already in physical units!)
                         for(i = 0; i <= 2; i++)
                            for(j = 0; j <= 2; j++)
                               dtensor[i][j] *= rho_fac;
                         
                         converged = get_axes(dtensor, &dambda1, &dambda2, &dambda3);
#else // DWEB_AK
                         Wg        = 0.0;
                         
                         if(haveneighbours == TRUE)  //Note here that this assumes tsc_nodes always have neighbours when mhd_nodes have neighbours.
                         {
                            for(i = 0; i <= 2; i++)
                               for(j = 0; j <= 2; j++)
                                  for(k = 0; k <= 2; k++){
                                     dFdx[i][j][k]=0;
                                     dFdy[i][j][k]=0;
                                     dFdz[i][j][k]=0;
                                  }
                            /* calculate the first order derivation*/
                            for(i = 1; i <= 3; i++)
                               for(j = 1; j <= 3; j++)
                                  for (k = 1; k<=3; k++){
                                     dFdx[i-1][j-1][k-1] = rho_fac*(mhd_nodes[i][j][k+1]->dens - mhd_nodes[i][j][k-1]->dens) / twodx;
                                     dFdy[i-1][j-1][k-1] = rho_fac*(mhd_nodes[i][j+1][k]->dens - mhd_nodes[i][j-1][k]->dens) / twody;
                                     dFdz[i-1][j-1][k-1] = rho_fac*(mhd_nodes[i+1][j][k]->dens - mhd_nodes[i-1][j][k]->dens) / twodz;
                                  }
                            /* calculate the second order derivation*/
                            for(i = 0; i <= 2; i++)
                               for(j = 0; j <= 2; j++){
                                  //                            distance2 = pow((i - 1) * (i - 1) + (j - 1) * (j - 1), 0.5);
                                  //                            distance2 = (1 + distance2) * (1 + distance2);  //weight = 1/(1+d)^2
                                  //                            Wg += 1./distance2;
                                  dtensor[0][0] += (dFdx[i][j][2] - dFdx[i][j][0]) / twodx;// / distance2 * x_fac;  // d^2F/dx^2
                                  dtensor[1][1] += (dFdy[i][2][j] - dFdy[i][0][j]) / twody;// / distance2 * x_fac;  // d^2F/dy^2
                                  dtensor[2][2] += (dFdz[2][i][j] - dFdz[0][i][j]) / twodz;// / distance2 * x_fac;  // d^2F/dz^2
                                  dtensor[0][1] += (dFdx[i][2][j] - dFdx[i][0][j]) / twody;// / distance2 * x_fac;  // d(dF/dx)/dy
                                  dtensor[0][2] += (dFdx[2][i][j] - dFdx[0][i][j]) / twodz;// / distance2 * x_fac;  // d(dF/dx)/dz
                                  dtensor[1][0] += (dFdy[i][j][2] - dFdy[i][j][0]) / twodx;// / distance2 * x_fac;  // d(dF/dy)/dx
                                  dtensor[1][2] += (dFdy[2][i][j] - dFdy[0][i][j]) / twodz;// / distance2 * x_fac;  // d(dF/dy)/dz
                                  dtensor[2][0] += (dFdz[i][j][2] - dFdz[i][j][0]) / twodx;// / distance2 * x_fac;  // d(dF/dz)/dx
                                  dtensor[2][1] += (dFdz[i][2][j] - dFdz[i][0][j]) / twody;// / distance2 * x_fac;  // d(dF/dz)/dy
                               }
                            //                        for(i = 0; i <= 2; i++)
                            //                          for(j = 0; j <= 2; j++)
                            //                            dtensor[i][j] /= Wg;
                            converged = get_axes(dtensor, &dambda1, &dambda2, &dambda3);
                            
                         }  //if have neighbours for mhd_nodes
                         
#endif // DWEB_AK
                         
#ifdef IGNORE_JACOBI_NONCONVERGENCE
                         if(converged == 0){
                            fprintf(stderr," ---> (Dweb) at this grid position x=%ld y=%ld z=%ld dens=%g densVx=%g densVy=%g densVz=%g\n",
                                    x,y,z,cur_node->dens,cur_node->densV[X],cur_node->densV[Y],cur_node->densV[Z]);
                            dambda1 = dambda2 = dambda3 = -3.1415;
                            dtensor[X][X] = dtensor[X][Y] = dtensor[X][Z] = -3.1415;
                            dtensor[X][Y] = dtensor[Y][Y] = dtensor[Y][Z] = -3.1415;
                            dtensor[X][Z] = dtensor[Y][Z] = dtensor[Z][Z] = -3.1415;
                         }
#endif
                         
#endif // DWEB
                         
                         
#ifdef PWEB
                         for(i = 0; i <= 2; i++)
                            for(j = 0; j <= 2; j++)
                               ptensor[i][j] = 0;
                         pambda1 = -1000.0;
                         pambda2 = -1000.0;
                         pambda3 = -1000.0;
                         
                         converged = 0;

#ifdef PWEB_AK // my own version --> not accurate enough!
                         
                         // Old version only using the direct neighbours, no fancy derivatives...(does not give reasonable results though!?)
                         ptensor[0][0] = ((tsc_nodes[1][1][2]->pot - tsc_nodes[1][1][1]->pot) - (tsc_nodes[1][1][1]->pot - tsc_nodes[1][1][0]->pot))/dx/dx;
                         ptensor[1][1] = ((tsc_nodes[1][2][1]->pot - tsc_nodes[1][1][1]->pot) - (tsc_nodes[1][1][1]->pot - tsc_nodes[1][0][1]->pot))/dy/dy; // this is 2nd order accurate
                         ptensor[2][2] = ((tsc_nodes[2][1][1]->pot - tsc_nodes[1][1][1]->pot) - (tsc_nodes[1][1][1]->pot - tsc_nodes[0][1][1]->pot))/dz/dz;
                         
                         /* directly use the second order derivative to calculate d^2 rho/dx/dy. */
                         ptensor[0][1] = ((tsc_nodes[1][2][2]->pot - tsc_nodes[1][0][2]->pot)/twody - (tsc_nodes[1][2][0]->pot - tsc_nodes[1][0][0]->pot)/twody)/twodx;
                         ptensor[0][2] = ((tsc_nodes[2][1][2]->pot - tsc_nodes[0][1][2]->pot)/twodz - (tsc_nodes[2][1][0]->pot - tsc_nodes[0][1][0]->pot)/twody)/twodx; // this is 2nd order accurate (A.10 in https://onlinelibrary.wiley.com/doi/pdf/10.1002/9781119083405.app1)
                         ptensor[1][2] = ((tsc_nodes[2][2][1]->pot - tsc_nodes[0][2][1]->pot)/twodz - (tsc_nodes[2][0][1]->pot - tsc_nodes[0][0][1]->pot)/twody)/twody;
                         
                         /* using the mean of first order derivative to calculate d^2 rho/dx/dy. */
                         //                       ptensor[0][1] = ((tsc_nodes[1][2][1]->pot - tsc_nodes[0][2][1]->pot) - (tsc_nodes[1][0][1]->pot - tsc_nodes[0][0][1]->pot) +
                         //                                        (tsc_nodes[2][2][1]->pot - tsc_nodes[1][2][1]->pot) - (tsc_nodes[2][0][1]->pot - tsc_nodes[1][0][1]->pot))/dx/dy/4.0;
                         //                       ptensor[0][2] = ((tsc_nodes[1][1][2]->pot - tsc_nodes[0][1][2]->pot) - (tsc_nodes[1][1][0]->pot - tsc_nodes[0][1][0]->pot) +
                         //                                        (tsc_nodes[2][1][2]->pot - tsc_nodes[1][1][2]->pot) - (tsc_nodes[2][1][0]->pot - tsc_nodes[1][1][0]->pot))/dx/dz/4.0;
                         //                       ptensor[1][2] = ((tsc_nodes[1][1][2]->pot - tsc_nodes[1][0][2]->pot) - (tsc_nodes[1][1][0]->pot - tsc_nodes[1][0][0]->pot) +
                         //                                        (tsc_nodes[1][2][2]->pot - tsc_nodes[1][1][2]->pot) - (tsc_nodes[1][2][0]->pot - tsc_nodes[1][1][0]->pot))/dx/dy/4.0;
                         
                         ptensor[1][0] = ptensor[0][1];
                         ptensor[2][0] = ptensor[0][2]; // the Hessian matrix is symmetric
                         ptensor[2][1] = ptensor[1][2];
                         
                         // add the missing conversion to physical units (note, dx/dy/dz are already in physical units!)
                         for(i = 0; i <= 2; i++)
                            for(j = 0; j <= 2; j++)
                               ptensor[i][j] *= pot_fac/Hz2;
                         
                         converged = get_axes(ptensor, &pambda1, &pambda2, &pambda3);
                         
#else // PWEB_AK
                         Wg = 0.0;
                         
                         if(haveneighbours == TRUE)  //Note here that this assumes tsc_nodes always have neighbours when mhd_nodes have neighbours.
                         {
                            for(i = 0; i <= 2; i++)
                               for(j = 0; j <= 2; j++)
                                  for(k = 0; k <= 2; k++){
                                     dFdx[i][j][k]=0;
                                     dFdy[i][j][k]=0;
                                     dFdz[i][j][k]=0;
                                  }
                            /* calculate the first order derivation*/
                            for(i = 1; i <= 3; i++)
                               for(j = 1; j <= 3; j++)
                                  for (k = 1; k<=3; k++){
                                     dFdx[i-1][j-1][k-1] = (mhd_nodes[i][j][k+1]->pot - mhd_nodes[i][j][k-1]->pot) / twodx; // pot in unit of grid which needs to * x_fac^2
                                     dFdy[i-1][j-1][k-1] = (mhd_nodes[i][j+1][k]->pot - mhd_nodes[i][j-1][k]->pot) / twody; // dx,dy,dz with length of 2 grids are already in physical units!
                                     dFdz[i-1][j-1][k-1] = (mhd_nodes[i+1][j][k]->pot - mhd_nodes[i-1][j][k]->pot) / twodz; //  [+1] - [-1] includes grids (2-0, 3-1, 4-2), [i,j,k] only selects central grids [1, 2, 3] for the second derivative
                                  }
                            /* calculate the second order derivation*/
                            for(i = 0; i <= 2; i++)
                               for(j = 0; j <= 2; j++){
                                  //                            distance2 = pow((i - 1) * (i - 1) + (j - 1) * (j - 1), 0.5);
                                  //                            distance2 = (1 + distance2) * (1 + distance2);  //weight = 1/(1+d)^2
                                  //                            Wg += 1./distance2;
                                  ptensor[0][0] += (dFdx[i][j][2] - dFdx[i][j][0]) / twodx;  // d^2F/dx^2
                                  ptensor[2][2] += (dFdz[2][i][j] - dFdz[0][i][j]) / twodz;  // d^2F/dz^2
                                  ptensor[0][1] += (dFdx[i][2][j] - dFdx[i][0][j]) / twody;  // d(dF/dx)/dy
                                  ptensor[0][2] += (dFdx[2][i][j] - dFdx[0][i][j]) / twodz;  // d(dF/dx)/dz
                                  ptensor[1][1] += (dFdy[i][2][j] - dFdy[i][0][j]) / twody;  // d^2F/dy^2
                                  ptensor[1][0] += (dFdy[i][j][2] - dFdy[i][j][0]) / twodx;  // d(dF/dy)/dx
                                  ptensor[1][2] += (dFdy[2][i][j] - dFdy[0][i][j]) / twodz;  // d(dF/dy)/dz
                                  ptensor[2][0] += (dFdz[i][j][2] - dFdz[i][j][0]) / twodx;  // d(dF/dz)/dx
                                  ptensor[2][1] += (dFdz[i][2][j] - dFdz[i][0][j]) / twody;  // d(dF/dz)/dy
                               }
                            //                        for(i = 0; i <= 2; i++)
                            //                          for(j = 0; j <= 2; j++)
                            //                            dtensor[i][j] /= Wg;
                            
                            // missing unit factor for potential, making tensor[][] dimensionless, too!?
                            for(i = 0; i <= 2; i++)
                               for(j = 0; j <= 2; j++)
                                  ptensor[i][j] *= pot_fac/Hz2;
                            
                            // the Hessian matrix should be symmetric:
                            //fprintf(stderr,"%lf =? %lf\n",ptensor[0][1],ptensor[1][0]); AK: these components are *not* identical!?
                            
                            // the trace of ptensor[][] should be the density (the source term, cf. solve_gravity.c)
                            //                        fprintf(stderr,"%g = %g ... ",
                            //                                cur_node->dens*simu.FourPiG*calc_super_a(cur_grid->timecounter),
                            //                                (ptensor[0][0]+ptensor[1][1]+ptensor[2][2])*Hz2/pot_fac*pow2(x_fac));
                            
                            converged = get_axes(ptensor, &pambda1, &pambda2, &pambda3);
                         }  //if have neighbours for mhd_nodes
                         
#endif // PWEB_AK
                         
#ifdef IGNORE_JACOBI_NONCONVERGENCE
                         if(converged == 0){
                            fprintf(stderr," ---> (Pweb) at this grid position x=%ld y=%ld z=%ld dens=%g densVx=%g densVy=%g densVz=%g\n",
                                    x,y,z,cur_node->dens,cur_node->densV[X],cur_node->densV[Y],cur_node->densV[Z]);
                            pambda1 = pambda2 = pambda3 = -3.1415;
                            ptensor[X][X] = ptensor[X][Y] = ptensor[X][Z] = -3.1415;
                            ptensor[X][Y] = ptensor[Y][Y] = ptensor[Y][Z] = -3.1415;
                            ptensor[X][Z] = ptensor[Y][Z] = ptensor[Z][Z] = -3.1415;
                         }
#endif
                         
#endif // PWEB
                         
                         /* write information to output file */
                         Nnodes_written++;
                         
                         tmp_x = ((float)x+0.5)/cur_grid->l1dim*x_fac;
                         tmp_y = ((float)y+0.5)/cur_grid->l1dim*x_fac;
                         tmp_z = ((float)z+0.5)/cur_grid->l1dim*x_fac;
                         
#ifdef AHFrfocusRescale
                         tmp_x /= simu.boxsize;
                         tmp_y /= simu.boxsize;
                         tmp_z /= simu.boxsize;
                         
                         tmp_x *= (2*AHFrfocusR)/simu.boxsize;
                         tmp_y *= (2*AHFrfocusR)/simu.boxsize;
                         tmp_z *= (2*AHFrfocusR)/simu.boxsize;
                         
                         tmp_x += (AHFrfocusX-AHFrfocusR)/simu.boxsize;
                         tmp_y += (AHFrfocusY-AHFrfocusR)/simu.boxsize;
                         tmp_z += (AHFrfocusZ-AHFrfocusR)/simu.boxsize;
                         
                         tmp_x *= simu.boxsize;
                         tmp_y *= simu.boxsize;
                         tmp_z *= simu.boxsize;
#endif
                         
#ifdef WRITE_ASCII
                         fprintf(fpout,"%16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8g %16.8g %16.8g %16.8lf %16.8lf %16.8lf %16.8lf %16.8lf %16.8lf %16.8lf %16.8lf %16.8lf\n",
                                 tmp_x,
                                 tmp_y,
                                 tmp_z,
                                 cur_node->dens+simu.mean_dens, // density in terms of background density
                                 cur_node->densV[X]/(cur_node->dens+simu.mean_dens)*v_fac,
                                 cur_node->densV[Y]/(cur_node->dens+simu.mean_dens)*v_fac, // velocity in cell
                                 cur_node->densV[Z]/(cur_node->dens+simu.mean_dens)*v_fac,
                                 vorticity[X],
                                 vorticity[Y], // already in the correct units
                                 vorticity[Z],
                                 lambda1,lambda2,lambda3,
                                 itensor[0][0],itensor[1][0],itensor[2][0],
                                 itensor[0][1],itensor[1][1],itensor[2][1],
                                 itensor[0][2],itensor[1][2],itensor[2][2]
                                 );
                         fflush(fpout);
#else
                         tmp_float = tmp_x;
                         FWRITE_TMP_FLOAT;
                         tmp_float = tmp_y;
                         FWRITE_TMP_FLOAT;
                         tmp_float = tmp_z;
                         FWRITE_TMP_FLOAT;
                         tmp_float = cur_node->dens+simu.mean_dens; // density in terms of background density
                         FWRITE_TMP_FLOAT;
                         tmp_float = cur_node->densV[X]/(cur_node->dens+simu.mean_dens)*v_fac;
                         FWRITE_TMP_FLOAT;
                         tmp_float = cur_node->densV[Y]/(cur_node->dens+simu.mean_dens)*v_fac; // velocity in cell
                         FWRITE_TMP_FLOAT;
                         tmp_float = cur_node->densV[Z]/(cur_node->dens+simu.mean_dens)*v_fac;
                         FWRITE_TMP_FLOAT;
                         tmp_float = vorticity[X];
                         FWRITE_TMP_FLOAT;
                         tmp_float = vorticity[Y];
                         FWRITE_TMP_FLOAT;
                         tmp_float = vorticity[Z];
                         FWRITE_TMP_FLOAT;
                         tmp_float = lambda1;
                         FWRITE_TMP_FLOAT;
                         tmp_float = lambda2;
                         FWRITE_TMP_FLOAT;
                         tmp_float = lambda3;
                         FWRITE_TMP_FLOAT;
                         tmp_float = itensor[0][0];
                         FWRITE_TMP_FLOAT;
                         tmp_float = itensor[1][0];
                         FWRITE_TMP_FLOAT;
                         tmp_float = itensor[2][0];
                         FWRITE_TMP_FLOAT;
                         tmp_float = itensor[0][1];
                         FWRITE_TMP_FLOAT;
                         tmp_float = itensor[1][1];
                         FWRITE_TMP_FLOAT;
                         tmp_float = itensor[2][1];
                         FWRITE_TMP_FLOAT;
                         tmp_float = itensor[0][2];
                         FWRITE_TMP_FLOAT;
                         tmp_float = itensor[1][2];
                         FWRITE_TMP_FLOAT;
                         tmp_float = itensor[2][2];
                         FWRITE_TMP_FLOAT;
#ifdef DWEB
                         tmp_float = dambda1;
                         FWRITE_TMP_FLOAT;
                         tmp_float = dambda2;
                         FWRITE_TMP_FLOAT;
                         tmp_float = dambda3;
                         FWRITE_TMP_FLOAT;
                         tmp_float = dtensor[0][0];
                         FWRITE_TMP_FLOAT;
                         tmp_float = dtensor[1][0];
                         FWRITE_TMP_FLOAT;
                         tmp_float = dtensor[2][0];
                         FWRITE_TMP_FLOAT;
                         tmp_float = dtensor[0][1];
                         FWRITE_TMP_FLOAT;
                         tmp_float = dtensor[1][1];
                         FWRITE_TMP_FLOAT;
                         tmp_float = dtensor[2][1];
                         FWRITE_TMP_FLOAT;
                         tmp_float = dtensor[0][2];
                         FWRITE_TMP_FLOAT;
                         tmp_float = dtensor[1][2];
                         FWRITE_TMP_FLOAT;
                         tmp_float = dtensor[2][2];
                         FWRITE_TMP_FLOAT;
#endif
#ifdef PWEB
                         tmp_float = cur_node->pot*pot_fac;
                         FWRITE_TMP_FLOAT;
                         tmp_float = pambda1;
                         FWRITE_TMP_FLOAT;
                         tmp_float = pambda2;
                         FWRITE_TMP_FLOAT;
                         tmp_float = pambda3;
                         FWRITE_TMP_FLOAT;
                         tmp_float = ptensor[0][0];
                         FWRITE_TMP_FLOAT;
                         tmp_float = ptensor[1][0];
                         FWRITE_TMP_FLOAT;
                         tmp_float = ptensor[2][0];
                         FWRITE_TMP_FLOAT;
                         tmp_float = ptensor[0][1];
                         FWRITE_TMP_FLOAT;
                         tmp_float = ptensor[1][1];
                         FWRITE_TMP_FLOAT;
                         tmp_float = ptensor[2][1];
                         FWRITE_TMP_FLOAT;
                         tmp_float = ptensor[0][2];
                         FWRITE_TMP_FLOAT;
                         tmp_float = ptensor[1][2];
                         FWRITE_TMP_FLOAT;
                         tmp_float = ptensor[2][2];
                         FWRITE_TMP_FLOAT;
#endif
#ifdef UonGrid
                         tmp_float = cur_node->u;
                         FWRITE_TMP_FLOAT;
#endif
#endif
                      } // if(haveneighbours)
                   } // if(check_for_zerodens())
                   else
                   {
                      /* count and write empty cells as void */
                      Nvoids++;
                      
                      /* write information to output file */
                      Nnodes_written++;
                      
                      tmp_x = ((float)x+0.5)/cur_grid->l1dim*x_fac;
                      tmp_y = ((float)y+0.5)/cur_grid->l1dim*x_fac;
                      tmp_z = ((float)z+0.5)/cur_grid->l1dim*x_fac;
                      
#ifdef AHFrfocusRescale
                      tmp_x /= simu.boxsize;
                      tmp_y /= simu.boxsize;
                      tmp_z /= simu.boxsize;
                      
                      tmp_x *= (2*AHFrfocusR)/simu.boxsize;
                      tmp_y *= (2*AHFrfocusR)/simu.boxsize;
                      tmp_z *= (2*AHFrfocusR)/simu.boxsize;
                      
                      tmp_x += (AHFrfocusX-AHFrfocusR)/simu.boxsize;
                      tmp_y += (AHFrfocusY-AHFrfocusR)/simu.boxsize;
                      tmp_z += (AHFrfocusZ-AHFrfocusR)/simu.boxsize;
                      
                      tmp_x *= simu.boxsize;
                      tmp_y *= simu.boxsize;
                      tmp_z *= simu.boxsize;
#endif
#ifdef WRITE_ASCII
                      fprintf(fpout,"%16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8g %16.8g %16.8g %16.8lf %16.8lf %16.8lf %16.8lf %16.8lf %16.8lf %16.8lf %16.8lf %16.8lf\n",
                              tmp_x,
                              tmp_y,
                              tmp_z,
                              cur_node->dens+simu.mean_dens, // density in terms of background density (should be < CWEB_DENSTHRESHOLD)
                              -1.,
                              -1.,
                              -1.,
                              -1.,
                              -1.,
                              -1.,
                              -1.,-1.,-1.,
                              -1.,-1.,-1.,
                              -1.,-1.,-1.,
                              -1.,-1.,-1.
                              );
                      fflush(fpout);
#else
                      tmp_float = tmp_x;
                      FWRITE_TMP_FLOAT;
                      tmp_float = tmp_y;
                      FWRITE_TMP_FLOAT;
                      tmp_float = tmp_z;
                      FWRITE_TMP_FLOAT;
                      tmp_float = cur_node->dens+simu.mean_dens; // density in terms of background density (should be < CWEB_DENSTHRESHOLD)
                      FWRITE_TMP_FLOAT;
                      tmp_float = -1.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -1.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -1.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -1.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -1.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -1.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -1.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -1.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -1.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -1.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -1.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -1.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -1.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -1.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -1.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -1.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -1.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -1.;
                      FWRITE_TMP_FLOAT;
#ifdef DWEB
                      tmp_float = -10.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -10.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -10.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -10.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -10.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -10.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -10.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -10.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -10.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -10.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -10.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -10.;
                      FWRITE_TMP_FLOAT;
#endif
#ifdef PWEB
                      tmp_float = -10.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -10.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -10.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -10.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -10.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -10.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -10.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -10.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -10.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -10.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -10.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -10.;
                      FWRITE_TMP_FLOAT;
                      tmp_float = -10.;
                      FWRITE_TMP_FLOAT;
#endif
#ifdef UonGrid
                      tmp_float = -1.;
                      FWRITE_TMP_FLOAT;
#endif
#endif
                   } // else(fabs(cur_node->dens+simu.mean_dens) > CWEB_DENSTHRESHOLD)
                   
                } // if(notboundary)
              } // for(cur_node)
            }
          }
        }
      }
    }
    fprintf(stderr,"done\n");
    
#ifndef WRITE_ASCII
    /* correct number of written nodes */
    fseek(fpout,(long)sizeof(int32_t),SEEK_SET); // skip the "one"
    tmp_uint64 = Nnodes_written;
    FWRITE_TMP_UINT64;
#endif
    
    fclose(fpout);
    
#ifdef LAMBDA_TH
#ifdef WITH_MPI
    // sum up all individual web counts
    //    MPI_Reduce(&Nknots,     &Nknots_all,     8, MPI_BYTE, MPI_SUM, 0, MPI_COMM_WORLD);
    //    MPI_Reduce(&Nfilaments, &Nfilaments_all, 8, MPI_BYTE, MPI_SUM, 0, MPI_COMM_WORLD);
    //    MPI_Reduce(&Nsheets,    &Nsheets_all,    8, MPI_BYTE, MPI_SUM, 0, MPI_COMM_WORLD);
    //    MPI_Reduce(&Nvoids,     &Nvoids_all,     8, MPI_BYTE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Nknots,     &Nknots_all,     1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Nfilaments, &Nfilaments_all, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Nsheets,    &Nsheets_all,    1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Nvoids,     &Nvoids_all,     1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
#endif

      /* dump filling fractions */
#ifdef WITH_MPI
      if(global_mpi.rank == 0)
#endif
      {
          fprintf(stderr,"\nVolume Filling Fractions: (grid=%ld, LambdaThreshold=%g)\n",cur_grid->l1dim,LAMBDA_THRESHOLD);
          fprintf(stderr,"=========================\n");
          fprintf(stderr,"knots     = %g\n",(double)Nknots/(double)pow3(cur_grid->l1dim));
          fprintf(stderr,"filaments = %g\n",(double)Nfilaments/(double)pow3(cur_grid->l1dim));
          fprintf(stderr,"sheets    = %g\n",(double)Nsheets/(double)pow3(cur_grid->l1dim));
          fprintf(stderr,"voids     = %g\n",(double)Nvoids/(double)pow3(cur_grid->l1dim));
      }
#endif //LAMBDA_TH

  } // for(cur_grid)
  ahf.time += time(NULL);
  
  /*=========================================================================================
   * update logfile and say bye-bye
   *=========================================================================================*/
#ifndef READ_GRIDDATA
  write_logfile(timecounter, timestep, no_timestep);
  
  /* free all allocated memory... */
  fprintf(stderr,"Clearing up...");
  free(grid_list);
  
  /*============================================================================
   *                                   BYE BYE!
   *============================================================================*/
  free(io.icfile_name);
  free(io.dumpfile_name);
  free(io.logfile_name);
  free(io.outfile_prefix);
  free(global.termfile_name);
  
  free(global.fst_part);
  if(global.fst_gas)
    free(global.fst_gas);
  if(global.fst_star)
    free(global.fst_star);
  
  
  fprintf(io.logfile, "==========================================================\n");
  fprintf(io.logfile, "                       FINISHED (v%3.1f/%03d)\n",VERSION,BUILD);
  fprintf(io.logfile, "==========================================================\n");
  fclose(io.logfile);
  
#	ifdef WITH_MPI
  /* Gracefully terminate MPI */
  MPI_Finalize();
#	endif
  
#ifdef DEBUG_CWEB
  fclose(fdebug);
#endif
#endif //READ_GRIDDATA
  
  fprintf(stderr,"done\n");
  
  return EXIT_SUCCESS;
}

#ifdef WITH_MPI
/*====================================================================
 *  MPI communications
 *====================================================================*/
static void
local_communicate_all(void)
{
  if (global_mpi.rank == 0) {
    /*********************\
     *      SENDING      *
     \*********************/
    /** This is where to send to */
    int target;
    
    for (target = 1; target < global_mpi.size; target++) {
      MPI_Send(&global, sizeof(struct info_global),
               MPI_BYTE,
               target, 0, MPI_COMM_WORLD);
    }
    /*********************\
     *    DONE SENDING   *
     \*********************/
  } else {
    /*********************\
     *     RECIEVING     *
     \*********************/
    /** The status */
    MPI_Status stat;
    
    MPI_Recv(&global, sizeof(struct info_global),
             MPI_BYTE,
             0, 0, MPI_COMM_WORLD, &stat);
    
    /*********************\
     *  DONE  RECIEVING  *
     \*********************/
  }
  
  return;
}
#endif

#ifdef AHFrfocus
/*====================================================================
 *  remove all particles outside a spherical region about centre[]
 *
 *  (though called AHFrfocus it is a multi-purpose routine)
 *====================================================================*/
static void
local_focusSphere(void)
{
  uint64_t oldNumPart = global_info.no_part;
  uint64_t newNumPart = UINT64_C(0);
  double center[3];
  double radSqr;
  int8_t   *tags;
  uint64_t i, j;
  partptr newParts;
  
  /* The sphere to focus on, give in Mpc/h (cf. param.h) */
  center[X] = AHFrfocusX;
  center[Y] = AHFrfocusY;
  center[Z] = AHFrfocusZ;
  radSqr    = AHFrfocusR;
  
  fprintf(stderr, "ATTENTION!  This is AHFrfocus calling:\n");
  fprintf(stderr, "This will remove all particles outside this sphere:\n");
  fprintf(stderr, "Center: (%lf %lf %lf) Mpc/h",
          center[X], center[Y], center[Z]);
  fprintf(stderr, "Radius: %lf Mpc/h\n", radSqr);
  fprintf(stderr,"starting with %ld particles -> ",global_info.no_part);
  
  /* Convert center and radius to AHF units */
  center[X] /= simu.boxsize;
  center[Y] /= simu.boxsize;
  center[Z] /= simu.boxsize;
  radSqr /= simu.boxsize;
  radSqr *= radSqr;
  
  /* Keeps a record whether the particle is kept or not. */
  tags = malloc(sizeof(int8_t) * oldNumPart);
  
  /* Check which particles to keep (ignoring periodicity) */
  for (i=UINT64_C(0); i<oldNumPart; i++) {
    double dpos[3];
    double rSqr;
    
    dpos[X] = global_info.fst_part[i].pos[X] - center[X];
    dpos[Y] = global_info.fst_part[i].pos[Y] - center[Y];
    dpos[Z] = global_info.fst_part[i].pos[Z] - center[Z];
    rSqr  = dpos[X] * dpos[X] + dpos[Y] * dpos[Y] + dpos[Z] * dpos[Z];
    if (rSqr <= radSqr)
    {
      tags[i] = 1;
      newNumPart++;
    } else {
      tags[i] = 0;
    }
  }
  
  /* Get new particles and copy the old ones over. */
  newParts = c_part(newNumPart);
  for (j = i = UINT64_C(0); i<oldNumPart; i++) {
    if (tags[i] == 1) {
      memcpy(newParts + j, global_info.fst_part + i,
             sizeof(struct particle));
      j++;
    }
  }
  assert(j == newNumPart);
  
  /* Clean */
  free(tags);
  free(global_info.fst_part);
  
  /* Activate the new particles. */
  global_info.fst_part = newParts;
  global_info.no_part  = newNumPart;
  fprintf(stderr,"ended with %ld particles\n\n",global_info.no_part);
}

#ifdef AHFrfocusRescale
/*====================================================================
 *  mimic a box of sidelength 2*AHFrfocusR
 *====================================================================*/
static void
local_focusSphereRescale(void)
{
  uint64_t i;
  double   xshift, yshift, zshift;
  
  xshift = (AHFrfocusX-AHFrfocusR)/simu.boxsize;
  yshift = (AHFrfocusY-AHFrfocusR)/simu.boxsize;
  zshift = (AHFrfocusZ-AHFrfocusR)/simu.boxsize;
  
  for(i=0; i<global_info.no_part; i++) {
    
    // shift new box to lower left corner
    global_info.fst_part[i].pos[X] -= xshift;
    global_info.fst_part[i].pos[Y] -= yshift;
    global_info.fst_part[i].pos[Z] -= zshift;
    
    // scale box
    global_info.fst_part[i].pos[X] = global_info.fst_part[i].pos[X]*simu.boxsize/(2*AHFrfocusR);
    global_info.fst_part[i].pos[Y] = global_info.fst_part[i].pos[Y]*simu.boxsize/(2*AHFrfocusR);
    global_info.fst_part[i].pos[Z] = global_info.fst_part[i].pos[Z]*simu.boxsize/(2*AHFrfocusR);
  }
}
#endif

#endif


