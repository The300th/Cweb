#include <stddef.h>
#include <math.h>
#ifdef WITH_OPENMP
#include <omp.h>
#endif

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "amr_serial.h"
#include "../libutility/utility.h"

/*=============================================================================
 * Gaussian smoothing
 *
 * done in Fourier space...
 * NOTE, we smooth all possible arrays in the same loop!
 *       This is highly memory consuming, but much faster
 *       If memoy becomes an issue, it should be easy to do the smoothing of the fields sequentially
 *=============================================================================*/
#ifdef PWEB
void smooth_arrays(flouble *parray, flouble *darray, flouble *xarray, flouble *yarray, flouble *zarray, long l1dim, double Rsmooth)
#else
void smooth_arrays(flouble *darray, flouble *xarray, flouble *yarray, flouble *zarray, long l1dim, double Rsmooth)
#endif
{
  double        Filter, sf, sf2;
  double        FFTnorm;
  double        kx, ky, kz, ksquare, kf;
  long          *k;
  long          k_dim, k_nyq;
  long          i, l, m, n;
  int           forward=1, backward=-1;
  unsigned long nn[NDIM];
  double        A, B; // short-cuts for constant factors
  
  // Gaussian variables
  double        sigma        = 1.0;                               // sigma of Gaussian filter in box units
  double        Rsmooth_box  = Rsmooth/simu.boxsize;              // smoothing scale in box units
  double        s            = 2.0 * pow2(sigma*Rsmooth_box);
  double        a            = 1/s;
  
  /* dimensions of grid */
  nn[0] = l1dim;
  nn[1] = l1dim;
  nn[2] = l1dim;
  k_dim = l1dim;
  k_nyq = l1dim/2;

  /* this will help extracting information from array^ */
  k = (long*) calloc(l1dim,sizeof(long));
  for(i=0; i<l1dim; i++)
  {
    if(i <= k_nyq)
      k[i] = i;
    else
      k[i] = -(l1dim-i);
  }
  
  /* fourn normalisation factor */
  FFTnorm = (double)pow3(l1dim);
  
  /* fundamental wave in box units (B=1): kf = 2*pi */
  kf  = TWOPI;
  
  /* short-cuts */
  sf  = PI/(double)l1dim/kf;
  sf2 = pow2(sf);
  A   = sqrt(PI/a);
  B   = pow2(PI)/(sf*a); // a is in inverse-k and also needs to be adjusted to internal units via 'sf'
  
  /* forward FFT of all arrays */
  fourn(darray-1, nn-1, NDIM, forward);
  fourn(xarray-1, nn-1, NDIM, forward);
  fourn(yarray-1, nn-1, NDIM, forward);
  fourn(zarray-1, nn-1, NDIM, forward);
#ifdef PWEB
  fourn(parray-1, nn-1, NDIM, forward);
#endif
  
  /* loop over all cells, applying the smoothing function */
  for(n = 0; n < k_dim; n++)
  {
    kz = kf * k[n];
    
    for(m = 0; m < k_dim; m++)
    {
      ky = kf * k[m];
      
      for(l = 0; l < k_dim; l++)
      {
        kx = kf * k[l];
        
        // k^2 in box units
        ksquare = sf2 * (kx * kx + ky * ky + kz * kz);
        
        // Gaussian filter (internal units):
        //    f(x) = exp(-a x^2)
        // => f^(k)= sqrt(pi/a) * exp(-pi^2 k^2 /a)
//        Filter = A * exp(-B*ksquare);   // the factor A cancels when considering the missing 1/sqrt(2pi sigma^2) normalisation factor in f(x)
        Filter = exp(-B*ksquare);

//        fprintf(stderr,"%g %g (%g)\n",sqrt(ksquare),Filter,TWOPI/Rsmooth_box*sf);
        
        // apply filter to all arrays
        darray[Re(l,m,n,l1dim)] *= Filter;
        darray[Im(l,m,n,l1dim)] *= Filter;
        xarray[Re(l,m,n,l1dim)] *= Filter;
        xarray[Im(l,m,n,l1dim)] *= Filter;
        yarray[Re(l,m,n,l1dim)] *= Filter;
        yarray[Im(l,m,n,l1dim)] *= Filter;
        zarray[Re(l,m,n,l1dim)] *= Filter;
        zarray[Im(l,m,n,l1dim)] *= Filter;
#ifdef PWEB
        parray[Re(l,m,n,l1dim)] *= Filter;
#endif
      }
    }
  }
  
  
  /* backward FFT */
  fourn(darray-1, nn-1, NDIM, backward);
  fourn(xarray-1, nn-1, NDIM, backward);
  fourn(yarray-1, nn-1, NDIM, backward);
  fourn(zarray-1, nn-1, NDIM, backward);
#ifdef PWEB
  fourn(parray-1, nn-1, NDIM, backward);
#endif

  /* normalize FFT output */
  for(i=0; i<2*l1dim*l1dim*l1dim; i++)
  {
    darray[i] /= FFTnorm;
    xarray[i] /= FFTnorm;
    yarray[i] /= FFTnorm;
    zarray[i] /= FFTnorm;
#ifdef PWEB
    parray[i] /= FFTnorm;
#endif
  }
  
  /* delete k-array again */
  free(k);
}


/*=============================================================================
 * Gaussian_smooth_gridFFT()
 *
 * wrapper for smooting all density fields
 * note, we are smooting all fields using the same for-loops,
 *       but this means to allocate as many temporary l1dim^3 grids as needed
 *       if memory becomes and issue, change the structure here...
 *=============================================================================*/
void Gaussian_smooth_gridFFT(gridls *cur_grid, double Rsmooth)
{
  flouble *dens_array, *densVx_array, *densVy_array, *densVz_array; /* density array pointers */
#ifdef PWEB
  flouble *pot_array;
#endif
  pqptr    cur_pquad;          /* current pquad         */
  cqptr    cur_cquad;          /* current cquad         */
  nqptr    cur_nquad;          /* current nquad         */
  nptr     cur_node;           /* current node          */
  long     i, j, k, l1dim, FFTarray_length;
  
  /* array dimension */
  l1dim           = cur_grid->l1dim;
  FFTarray_length = 2*l1dim*l1dim*l1dim;
  
  /* generate complex (!) density arrays for FFT */
  if((dens_array = (flouble *) calloc(FFTarray_length, sizeof(flouble))) == NULL)
  {
    fprintf(io.logfile,"Gaussian_smooth_gridFFT: could not allocate density array for FFT\n");
    fflush(io.logfile);
    fclose(io.logfile);
    exit(1);
  }
  /* generate complex (!) density arrays for FFT */
  if((densVx_array = (flouble *) calloc(FFTarray_length, sizeof(flouble))) == NULL)
  {
    fprintf(io.logfile,"Gaussian_smooth_gridFFT: could not allocate density array for FFT\n");
    fflush(io.logfile);
    fclose(io.logfile);
    exit(1);
  }
  if((densVy_array = (flouble *) calloc(FFTarray_length, sizeof(flouble))) == NULL)
  {
    fprintf(io.logfile,"Gaussian_smooth_gridFFT: could not allocate density array for FFT\n");
    fflush(io.logfile);
    fclose(io.logfile);
    exit(1);
  }
  if((densVz_array = (flouble *) calloc(FFTarray_length, sizeof(flouble))) == NULL)
  {
    fprintf(io.logfile,"Gaussian_smooth_gridFFT: could not allocate density array for FFT\n");
    fflush(io.logfile);
    fclose(io.logfile);
    exit(1);
  }
#ifdef PWEB
  if((pot_array = (flouble *) calloc(FFTarray_length, sizeof(flouble))) == NULL)
  {
    fprintf(io.logfile,"Gaussian_smooth_gridFFT: could not allocate density array for FFT\n");
    fflush(io.logfile);
    fclose(io.logfile);
    exit(1);
  }
#endif

  /* fill density array for FFT ... no need for quad-ll's !! */
  cur_pquad = cur_grid->pquad;
  for(k = 0, cur_cquad = cur_pquad->loc; k < l1dim; k++, cur_cquad++)
    for(j = 0, cur_nquad = cur_cquad->loc; j < l1dim; j++, cur_nquad++)
      for(i = 0, cur_node = cur_nquad->loc; i < l1dim; i++, cur_node++)
      {
        dens_array[Re(i,j,k,l1dim)]   = cur_node->dens;     /* real part      */
        dens_array[Im(i,j,k,l1dim)]   = 0.0;                /* imaginary part */
        densVx_array[Re(i,j,k,l1dim)] = cur_node->densV[X]; /* real part      */
        densVx_array[Im(i,j,k,l1dim)] = 0.0;                /* imaginary part */
        densVy_array[Re(i,j,k,l1dim)] = cur_node->densV[Y]; /* real part      */
        densVy_array[Im(i,j,k,l1dim)] = 0.0;                /* imaginary part */
        densVz_array[Re(i,j,k,l1dim)] = cur_node->densV[Z]; /* real part      */
        densVz_array[Im(i,j,k,l1dim)] = 0.0;                /* imaginary part */
#ifdef PWEB
        pot_array[Re(i,j,k,l1dim)]    = cur_node->pot;      /* real part      */
        pot_array[Im(i,j,k,l1dim)]    = 0.0;                /* imaginary part */
#endif
      }
  
  /* solve by FFT */
#ifdef PWEB
  smooth_arrays(pot_array, dens_array, densVx_array, densVy_array, densVz_array, l1dim, Rsmooth);
#else
  smooth_arrays(dens_array, densVx_array, densVy_array, densVz_array, l1dim, Rsmooth);
#endif
  
  /* fill node potential values */
  for(k = 0, cur_cquad = cur_pquad->loc; k < cur_grid->l1dim; k++, cur_cquad++)
    for(j = 0, cur_nquad = cur_cquad->loc; j < cur_grid->l1dim; j++, cur_nquad++)
      for(i = 0, cur_node = cur_nquad->loc; i < cur_grid->l1dim; i++,cur_node++) {
        cur_node->dens     = dens_array[Re(i,j,k,l1dim)];
        cur_node->densV[X] = densVx_array[Re(i,j,k,l1dim)];
        cur_node->densV[Y] = densVy_array[Re(i,j,k,l1dim)];
        cur_node->densV[Z] = densVz_array[Re(i,j,k,l1dim)];
#ifdef PWEB
        cur_node->pot      = pot_array[Re(i,j,k,l1dim)];
#endif
      }
  
  /* destroy memory assigned to densV_arrays */
  free(dens_array);
  free(densVx_array);
  free(densVy_array);
  free(densVz_array);
#ifdef PWEB
  free(pot_array);
#endif
}

/*=============================================================================
 * unsophisticated average smoothing
 *=============================================================================*/
void TSCMHDsmooth_grid(gridls *cur_grid)
{
#ifdef WITH_OPENMP
  long          ipquad;
#endif

  pqptr         cur_pquad;
  cqptr         cur_cquad, icur_cquad;
  nqptr         cur_nquad, icur_nquad;
  nptr          cur_node;
  long          cur_x, cur_y, cur_z;
  nptr          tsc_nodes[3][3][3];
  nptr          mhd_nodes[5][5][5];
  int           i, j, k, nsmooth;

  /*==============================================================
   * three loops over cur_grid:
   *  1. reset temporary values
   *  2. smooth storing results in temporary storage
   *  3. copy temporary storage over to actual variables
   *==============================================================*/




  /*----------------
   * 1. init-loop
   *----------------*/
#ifdef WITH_OPENMP
#pragma omp parallel private(ipquad, cur_pquad,  cur_cquad, icur_cquad, cur_nquad, icur_nquad, cur_node, cur_x, cur_y, cur_z)  shared(cur_grid)
#pragma omp for schedule(static)
  for(ipquad=0; ipquad<cur_grid->no_pquad; ipquad++) {
    cur_pquad=cur_grid->pquad_array[ipquad];
#else
  for(cur_pquad=cur_grid->pquad; cur_pquad != NULL; cur_pquad=cur_pquad->next) {
#endif
      cur_z = cur_pquad->z;

      for(cur_cquad = cur_pquad->loc;
          cur_cquad < cur_pquad->loc + cur_pquad->length;
          cur_cquad++, cur_z++)
       {
        for(icur_cquad  = cur_cquad;
            icur_cquad != NULL;
            icur_cquad  = icur_cquad->next)
         {
          cur_y = icur_cquad->y;

          for(cur_nquad = icur_cquad->loc;
              cur_nquad < icur_cquad->loc + icur_cquad->length;
              cur_nquad++, cur_y++)
           {
            for(icur_nquad  = cur_nquad;
                icur_nquad != NULL;
                icur_nquad  = icur_nquad->next)
             {
              cur_x = icur_nquad->x;

              for(cur_node = icur_nquad->loc;
                  cur_node < icur_nquad->loc + icur_nquad->length;
                  cur_node++, cur_x++)
               {
                cur_node->densVtmp[X]   = 0.0;
                cur_node->densVtmp[Y]   = 0.0;
                cur_node->densVtmp[Z]   = 0.0;
                cur_node->force.temp[0] = 0.0;
#ifdef UonGrid
                cur_node->force.temp[1] = 0.0;
#endif
#ifdef PWEB
                cur_node->force.temp[2] = 0.0;
#endif
               }
             }
           }
         }
       }
    }




    /*----------------
     * 2. smooth-loop
     *----------------*/
#ifdef WITH_OPENMP
#pragma omp parallel private(ipquad, cur_pquad,  cur_cquad, icur_cquad, cur_nquad, icur_nquad, cur_node, cur_x, cur_y, cur_z)  shared(cur_grid)
#pragma omp for schedule(static)
  for(ipquad=0; ipquad<cur_grid->no_pquad; ipquad++) {
    cur_pquad=cur_grid->pquad_array[ipquad];
#else
  for(cur_pquad=cur_grid->pquad; cur_pquad != NULL; cur_pquad=cur_pquad->next) {
#endif
      cur_z = cur_pquad->z;

      for(cur_cquad = cur_pquad->loc;
          cur_cquad < cur_pquad->loc + cur_pquad->length;
          cur_cquad++, cur_z++)
       {
         for(icur_cquad  = cur_cquad;
             icur_cquad != NULL;
             icur_cquad  = icur_cquad->next)
          {
           cur_y = icur_cquad->y;

           for(cur_nquad = icur_cquad->loc;
               cur_nquad < icur_cquad->loc + icur_cquad->length;
               cur_nquad++, cur_y++)
            {
              for(icur_nquad  = cur_nquad;
                  icur_nquad != NULL;
                  icur_nquad  = icur_nquad->next)
               {
                cur_x = icur_nquad->x;

                for(cur_node = icur_nquad->loc;
                    cur_node < icur_nquad->loc + icur_nquad->length;
                    cur_node++, cur_x++)
                 {
                  /* keep track of how many cells have been used during the smoothing */
                  nsmooth = 0;

                  /*  smooth over neighbours */
#ifdef MHD_SMOOTHING
                  get_MHDnodes(cur_grid, cur_pquad, cur_z, cur_y, cur_x, mhd_nodes);
                  for(k = 0; k < 5; k++){
                    for(j = 0; j < 5; j++){
                      for(i = 0; i < 5; i++){
                        if(mhd_nodes[k][j][i] != NULL) {
                          nsmooth++;
                          cur_node->densVtmp[X]   += mhd_nodes[k][j][i]->densV[X];
                          cur_node->densVtmp[Y]   += mhd_nodes[k][j][i]->densV[Y];
                          cur_node->densVtmp[Z]   += mhd_nodes[k][j][i]->densV[Z];
                          cur_node->force.temp[0] += mhd_nodes[k][j][i]->dens;
#ifdef UonGrid
                          cur_node->force.temp[1] += mhd_nodes[k][j][i]->u;
#endif
#ifdef PWEB
                          cur_node->force.temp[2] += mhd_nodes[k][j][i]->pot;
#endif
                        }
                      }
                    }
                  }
#else // MHD_SMOOTHING
                  tsc_nodes[1][1][1] = cur_node;
                  get_TSCnodes(cur_grid, cur_pquad, icur_cquad, icur_nquad, tsc_nodes, &cur_z, &cur_y, &cur_x);
                  for(k = 0; k < 3; k++){
                    for(j = 0; j < 3; j++){
                      for(i = 0; i < 3; i++){
                        if(tsc_nodes[k][j][i] != NULL) {
                          nsmooth++;
                          cur_node->densVtmp[X]   += tsc_nodes[k][j][i]->densV[X];
                          cur_node->densVtmp[Y]   += tsc_nodes[k][j][i]->densV[Y];
                          cur_node->densVtmp[Z]   += tsc_nodes[k][j][i]->densV[Z];
                          cur_node->force.temp[0] += tsc_nodes[k][j][i]->dens;
#ifdef UonGrid
                          cur_node->force.temp[1] += tsc_nodes[k][j][i]->u;
#endif
#ifdef PWEB
                          cur_node->force.temp[2] += tsc_nodes[k][j][i]->pot;
#endif
                        }
                      }
                    }
                  }
#endif // MHD_SMOOTHING

                  /* store smoothed momentum density in densVtmp[] and mass density in force.forces[X] */
                  cur_node->densVtmp[X]   /= (double)(nsmooth);
                  cur_node->densVtmp[Y]   /= (double)(nsmooth);
                  cur_node->densVtmp[Z]   /= (double)(nsmooth);
                  cur_node->force.temp[0] /= (double)(nsmooth);
#ifdef UonGrid
                  cur_node->force.temp[1] /= (double)(nsmooth);
#endif
#ifdef PWEB
                  cur_node->force.temp[2] /= (double)(nsmooth);
#endif
                 }
               }
            }
          }
       }
     }



    /*----------------
     * 3. copy-loop
     *----------------*/
#ifdef WITH_OPENMP
#pragma omp parallel private(ipquad, cur_pquad,  cur_cquad, icur_cquad, cur_nquad, icur_nquad, cur_node, cur_x, cur_y, cur_z)  shared(cur_grid)
#pragma omp for schedule(static)
    for(ipquad=0; ipquad<cur_grid->no_pquad; ipquad++) {
      cur_pquad=cur_grid->pquad_array[ipquad];
#else
    for(cur_pquad=cur_grid->pquad; cur_pquad != NULL; cur_pquad=cur_pquad->next) {
#endif
        cur_z = cur_pquad->z;

        for(cur_cquad = cur_pquad->loc;
            cur_cquad < cur_pquad->loc + cur_pquad->length;
            cur_cquad++, cur_z++)
         {
          for(icur_cquad  = cur_cquad;
              icur_cquad != NULL;
              icur_cquad  = icur_cquad->next)
           {
            cur_y = icur_cquad->y;

            for(cur_nquad = icur_cquad->loc;
                cur_nquad < icur_cquad->loc + icur_cquad->length;
                cur_nquad++, cur_y++)
             {
              for(icur_nquad  = cur_nquad;
                  icur_nquad != NULL;
                  icur_nquad  = icur_nquad->next)
               {
                cur_x = icur_nquad->x;

                for(cur_node = icur_nquad->loc;
                    cur_node < icur_nquad->loc + icur_nquad->length;
                    cur_node++, cur_x++)
                 {
                  cur_node->densV[X]  = cur_node->densVtmp[X];
                  cur_node->densV[Y]  = cur_node->densVtmp[Y];
                  cur_node->densV[Z]  = cur_node->densVtmp[Z];
                  cur_node->dens      = cur_node->force.temp[0];
#ifdef UonGrid
                  cur_node->u         = cur_node->force.temp[1];
#endif
#ifdef PWEB
                  cur_node->pot       = cur_node->force.temp[2];
#endif
                 }
               }

             }
           }
         }
      }


}
    
    
void createFilter(double gKernel[5][5][5])
{
  // set standard deviation to 1.0
  double sigma = 1.0;
  double r, s = 2.0 * sigma * sigma;
  
  // sum is for normalization
  double sum = 0.0;
  
  // generate 5x5 kernel
  for (int x = -2; x <= 2; x++) {
    for(int y = -2; y <= 2; y++) {
      for(int z = -2; y <= 2; y++) {
        r = sqrt(x*x + y*y + z*z);
        gKernel[x + 2][y + 2][z + 2] = (exp(-(r*r)/s))/(M_PI * s);
        sum += gKernel[x + 2][y + 2][z + 2];
      }
    }
  }
  
  // normalize the Kernel
  for(int i = 0; i < 5; ++i)
    for(int j = 0; j < 5; ++j)
      for(int k = 0; j < 5; ++j)
        gKernel[i][j][k] /= sum;
  
}
    



