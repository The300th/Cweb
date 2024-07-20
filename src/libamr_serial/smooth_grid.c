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
 * TopHat smoothing
 *
 * done in Fourier space...
 * NOTE, we smooth all possible arrays in the same loop!
 *       This is highly memory consuming, but much faster
 *       If memoy becomes an issue, it should be easy to do the smoothing of the fields sequentially
 *=============================================================================*/
#ifdef PWEB
void TopHat_smooth_arrays(flouble *parray, flouble *darray, flouble *vxarray, flouble *vyarray, flouble *vzarray, long l1dim, double Rsmooth)
#else
void TopHat_smooth_arrays(flouble *darray, flouble *vxarray, flouble *vyarray, flouble *vzarray, long l1dim, double Rsmooth)
#endif
{
   double        Filter;
   double        FFTnorm;
   double        kx, ky, kz, ksquare, k_box, kf, x;
   long          *k;
   long          k_dim, k_nyq;
   long          i, l, m, n;
   int           forward=1, backward=-1;
   unsigned long nn[NDIM];
   
   // TopHat variables
   double        Rsmooth_box = Rsmooth/simu.boxsize; // smoothing scale in box units
   
   /* dimensions of grid */
   nn[0] = l1dim;
   nn[1] = l1dim;
   nn[2] = l1dim;
   k_dim = l1dim;
   k_nyq = l1dim/2;
   
   /* k-modes in box units */
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
   kf = TWOPI;
      
   /* forward FFT of all arrays */
   fourn(darray-1, nn-1, NDIM, forward);
   fourn(vxarray-1, nn-1, NDIM, forward);
   fourn(vyarray-1, nn-1, NDIM, forward);
   fourn(vzarray-1, nn-1, NDIM, forward);
#ifdef PWEB
   fourn(parray-1, nn-1, NDIM, forward);
#endif
   
   /* loop over all cells, applying the smoothing function */
#ifdef WITH_OPENMP
#pragma omp parallel for private(n,m,l,kz,ky,kx,ksquare,k_box,Filter) shared(darray,vxarray,vyarray,vzarray,l1dim,k_dim,kf,k) schedule(static)
#endif
   for(n = 0; n < k_dim; n++)
   {
      kz = kf * k[n];
      
      for(m = 0; m < k_dim; m++)
      {
         ky = kf * k[m];
         
         for(l = 0; l < k_dim; l++)
         {
            kx = kf * k[l];
            
            // k in box units
            ksquare = (kx * kx + ky * ky + kz * kz);
            k_box   = sqrt(ksquare);
            
            // argument of window function
            x = k_box*Rsmooth_box;
            
            // TopHat filter:
            if(k_box > 0.0)
              Filter = 3.0*(sin(x)-x*cos(x))/pow3(x);
            else
              Filter = 1.0;
            
//            fprintf(stderr,"%g %g (%g %g %g)\n",k_box,Filter,kx,ky,kz);
                        
            // apply filter to all arrays
            darray[Re(l,m,n,l1dim)] *= Filter;
            darray[Im(l,m,n,l1dim)] *= Filter;
            
            vxarray[Re(l,m,n,l1dim)] *= Filter;
            vxarray[Im(l,m,n,l1dim)] *= Filter;
            
            vyarray[Re(l,m,n,l1dim)] *= Filter;
            vyarray[Im(l,m,n,l1dim)] *= Filter;
            
            vzarray[Re(l,m,n,l1dim)] *= Filter;
            vzarray[Im(l,m,n,l1dim)] *= Filter;
#ifdef PWEB
            parray[Re(l,m,n,l1dim)] *= Filter;
            parray[Im(l,m,n,l1dim)] *= Filter;
#endif
         }
      }
   }
   
   
   /* backward FFT */
   fourn(darray-1, nn-1, NDIM, backward);
   fourn(vxarray-1, nn-1, NDIM, backward);
   fourn(vyarray-1, nn-1, NDIM, backward);
   fourn(vzarray-1, nn-1, NDIM, backward);
#ifdef PWEB
   fourn(parray-1, nn-1, NDIM, backward);
#endif
   
   /* normalize FFT output */
   for(i=0; i<2*l1dim*l1dim*l1dim; i++)
   {
      darray[i] /= FFTnorm;
      vxarray[i] /= FFTnorm;
      vyarray[i] /= FFTnorm;
      vzarray[i] /= FFTnorm;
#ifdef PWEB
      parray[i] /= FFTnorm;
#endif
   }
   
   /* delete k-array again */
   free(k);
}


/*=============================================================================
 * Gaussian smoothing
 *
 * done in Fourier space...
 * NOTE, we smooth all possible arrays in the same loop!
 *       This is highly memory consuming, but much faster
 *       If memoy becomes an issue, it should be easy to do the smoothing of the fields sequentially
 *=============================================================================*/
#ifdef PWEB
void Gauss_smooth_arrays(flouble *parray, flouble *darray, flouble *vxarray, flouble *vyarray, flouble *vzarray, long l1dim, double Rsmooth)
#else
void Gauss_smooth_arrays(flouble *darray, flouble *vxarray, flouble *vyarray, flouble *vzarray, long l1dim, double Rsmooth)
#endif
{
   double        Filter;
   double        FFTnorm;
   double        kx, ky, kz, ksquare, kf;
   long          *k;
   long          k_dim, k_nyq;
   long          i, l, m, n;
   int           forward=1, backward=-1;
   unsigned long nn[NDIM];
   double        A, B; // short-cuts for constant factors
   
   // Gaussian variables
   double        sigma        = 1.0/(double)l1dim;                 // sigma of Gaussian filter in box units
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
   
   A   = sqrt(PI/a);
   B   = pow2(PI)/(a);
   
   /* forward FFT of all arrays */
   fourn(darray-1, nn-1, NDIM, forward);
   fourn(vxarray-1, nn-1, NDIM, forward);
   fourn(vyarray-1, nn-1, NDIM, forward);
   fourn(vzarray-1, nn-1, NDIM, forward);
#ifdef PWEB
   fourn(parray-1, nn-1, NDIM, forward);
#endif
   
   /* loop over all cells, applying the smoothing function */
#ifdef WITH_OPENMP
#pragma omp parallel for private(n,m,l,kz,ky,kx,ksquare,Filter) shared(darray,vxarray,vyarray,vzarray,l1dim,k_dim,kf,k,B) schedule(static)
#endif
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
            ksquare = (kx * kx + ky * ky + kz * kz);
            
            // Gaussian filter (internal units):
            //    f(x) = exp(-a x^2)
            // => f^(k)= sqrt(pi/a) * exp(-pi^2 k^2 /a)
            //        Filter = A * exp(-B*ksquare);   // the factor A cancels when considering the missing 1/sqrt(2pi sigma^2) normalisation factor in f(x)
            Filter = exp(-B*ksquare);
            
            //        fprintf(stderr,"%g %g (%g)\n",sqrt(ksquare),Filter,TWOPI/Rsmooth_box);
            
            // apply filter to all arrays
            darray[Re(l,m,n,l1dim)] *= Filter;
            darray[Im(l,m,n,l1dim)] *= Filter;
            
            vxarray[Re(l,m,n,l1dim)] *= Filter;
            vxarray[Im(l,m,n,l1dim)] *= Filter;
            
            vyarray[Re(l,m,n,l1dim)] *= Filter;
            vyarray[Im(l,m,n,l1dim)] *= Filter;
            
            vzarray[Re(l,m,n,l1dim)] *= Filter;
            vzarray[Im(l,m,n,l1dim)] *= Filter;
#ifdef PWEB
            parray[Re(l,m,n,l1dim)] *= Filter;
            parray[Im(l,m,n,l1dim)] *= Filter;
#endif
         }
      }
   }
   
   
   /* backward FFT */
   fourn(darray-1, nn-1, NDIM, backward);
   fourn(vxarray-1, nn-1, NDIM, backward);
   fourn(vyarray-1, nn-1, NDIM, backward);
   fourn(vzarray-1, nn-1, NDIM, backward);
#ifdef PWEB
   fourn(parray-1, nn-1, NDIM, backward);
#endif
   
   /* normalize FFT output */
   for(i=0; i<2*l1dim*l1dim*l1dim; i++)
   {
      darray[i] /= FFTnorm;
      vxarray[i] /= FFTnorm;
      vyarray[i] /= FFTnorm;
      vzarray[i] /= FFTnorm;
#ifdef PWEB
      parray[i] /= FFTnorm;
#endif
   }
   
   /* delete k-array again */
   free(k);
}


/*=============================================================================
 * Smooth_gridFFT()
 *
 * wrapper for smooting all density fields
 * note, we are smooting all fields using the same for-loops,
 *       but this means to allocate as many temporary l1dim^3 grids as needed
 *       if memory becomes and issue, change the structure here...
 *=============================================================================*/
void Smooth_gridFFT(gridls *cur_grid, double Rsmooth)
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
      fprintf(io.logfile,"Smooth_gridFFT: could not allocate density array for FFT\n");
      fflush(io.logfile);
      fclose(io.logfile);
      exit(1);
   }
   /* generate complex (!) density arrays for FFT */
   if((densVx_array = (flouble *) calloc(FFTarray_length, sizeof(flouble))) == NULL)
   {
      fprintf(io.logfile,"Smooth_gridFFT: could not allocate density array for FFT\n");
      fflush(io.logfile);
      fclose(io.logfile);
      exit(1);
   }
   if((densVy_array = (flouble *) calloc(FFTarray_length, sizeof(flouble))) == NULL)
   {
      fprintf(io.logfile,"Smooth_gridFFT: could not allocate density array for FFT\n");
      fflush(io.logfile);
      fclose(io.logfile);
      exit(1);
   }
   if((densVz_array = (flouble *) calloc(FFTarray_length, sizeof(flouble))) == NULL)
   {
      fprintf(io.logfile,"Smooth_gridFFT: could not allocate density array for FFT\n");
      fflush(io.logfile);
      fclose(io.logfile);
      exit(1);
   }
#ifdef PWEB
   if((pot_array = (flouble *) calloc(FFTarray_length, sizeof(flouble))) == NULL)
   {
      fprintf(io.logfile,"Smooth_gridFFT: could not allocate density array for FFT\n");
      fflush(io.logfile);
      fclose(io.logfile);
      exit(1);
   }
#endif
   
   /* fill density array for FFT ... no need for quad-ll's !! */
   cur_pquad = cur_grid->pquad;
#ifdef WITH_OPENMP
#ifdef PWEB
#pragma omp parallel private(k, j, i, cur_cquad, cur_nquad, cur_node) shared(cur_pquad, dens_array, densVx_array, densVy_array, densVz_array, l1dim, pot_array)
#else
#pragma omp parallel private(k, j, i, cur_cquad, cur_nquad, cur_node) shared(cur_pquad, dens_array, densVx_array, densVy_array, densVz_array, l1dim)
#endif
#pragma omp for schedule(static)
#endif
   for(k = 0; k < l1dim; k++) {
      cur_cquad = (cur_pquad->loc)+k;
      for(j = 0, cur_nquad = cur_cquad->loc; j < l1dim; j++, cur_nquad++) {
         for(i = 0, cur_node = cur_nquad->loc; i < l1dim; i++, cur_node++) {
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
      }
   }
   
   
   /* the actual smoothing is done on the temporary grids */
#ifdef GAUSSIAN_SMOOTHING
#ifdef PWEB
   Gauss_smooth_arrays(pot_array, dens_array, densVx_array, densVy_array, densVz_array, l1dim, Rsmooth);
#else
   Gauss_smooth_arrays(dens_array, densVx_array, densVy_array, densVz_array, l1dim, Rsmooth);
#endif
#else // GAUSSIAN_SMOOTHING
#ifdef PWEB
   TopHat_smooth_arrays(pot_array, dens_array, densVx_array, densVy_array, densVz_array, l1dim, Rsmooth);
#else
   TopHat_smooth_arrays(dens_array, densVx_array, densVy_array, densVz_array, l1dim, Rsmooth);
#endif
#endif // GAUSSIAN_SMOOTHING
   
   /* copy back to original grid */
#ifdef WITH_OPENMP
#ifdef PWEB
#pragma omp parallel private(k, j, i, cur_cquad, cur_nquad, cur_node) shared(cur_pquad, dens_array, densVx_array, densVy_array, densVz_array, l1dim, pot_array)
#else
#pragma omp parallel private(k, j, i, cur_cquad, cur_nquad, cur_node) shared(cur_pquad, dens_array, densVx_array, densVy_array, densVz_array, l1dim)
#endif
#pragma omp for schedule(static)
#endif
   for(k = 0; k < cur_grid->l1dim; k++) {
      cur_cquad = (cur_pquad->loc)+k;
      for(j = 0, cur_nquad = cur_cquad->loc; j < cur_grid->l1dim; j++, cur_nquad++) {
         for(i = 0, cur_node = cur_nquad->loc; i < cur_grid->l1dim; i++,cur_node++) {
            cur_node->dens     = dens_array[Re(i,j,k,l1dim)];
            cur_node->densV[X] = densVx_array[Re(i,j,k,l1dim)];
            cur_node->densV[Y] = densVy_array[Re(i,j,k,l1dim)];
            cur_node->densV[Z] = densVz_array[Re(i,j,k,l1dim)];
#ifdef PWEB
            cur_node->pot      = pot_array[Re(i,j,k,l1dim)];
#endif
         }
      }
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
 
 
// routine once used to perform a smoothing in real-space over 5x5x5 cells using this (Gaussian) kernel
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
    



