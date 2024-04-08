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

/* this constructs coa_dens values by simply restricting fin_dens */
//#define NEW_F2C_SOURCE

/*==============================================================================
*
* This file includes everything for the restriction operations
*
* visible functions:
*
* 1. f2c_dens       -> interpolate density         from fine to coarse grid
*
* 2. f2c_pot        -> interpolate potential       from fine to coarse grid
*                      (store it in force.temp[1] -> only for trunc_err())
*
* 3. f2c_volume     -> correct coarse grid volume for potential energy calculation:
*                      force.temp[0]:
*                       > 2.0 + do not use this node for volume integration
*                       < 2.0 = use this node for volume integration
*                       < 1.0 = use modified volume element
*                      force.temp[1]:     modified volume element
*
* 4. fc_overlap     -> temporarily stores the density contribution
*                      from particles spilling across the boundary
*                      from the coarse into the fine grid in force.temp[0]
*
*==============================================================================*/


/*=============================================================================
 * fine_to_coarse interpolation of density
 *=============================================================================*/
void f2c_dens(gridls *fin_grid)
{
   gridls       *coa_grid;
   
   long          ipquad;
   
   pqptr         coa_pquad, tcoa_pquad;
   cqptr         coa_cquad, tcoa_cquad;
   nqptr         coa_nquad, tcoa_nquad;
   nptr          coa_node;
   long          coa_x, coa_y, coa_z;
   nptr          tsc_coanodes[3][3][3];
   
   pqptr         fin_pquad;
   cqptr         fin_cquad, ifin_cquad;
   nqptr         fin_nquad, ifin_nquad;
   nptr          fin_node;
   long          fin_x, fin_y, fin_z;
  
   /* it's your responsibility to make sure that a coarse grid actually exists */
   coa_grid = fin_grid-1;
   
   /*==============================================================
    * loop over fine grid (and simultanesouly over coarse grid...)
    *==============================================================*/
#ifdef WITH_OPENMP
#pragma omp parallel private(ipquad, fin_pquad,  fin_cquad, ifin_cquad, fin_nquad, ifin_nquad, fin_node, fin_x, fin_y, fin_z, tcoa_pquad, coa_cquad, tcoa_cquad, coa_nquad, tcoa_nquad, coa_node, coa_x, coa_y, coa_z)  shared(fin_grid, coa_grid, simu)
#pragma omp for schedule(static)
   for(ipquad=0; ipquad<fin_grid->no_pquad; ipquad++)
     {
     fin_pquad=fin_grid->pquad_array[ipquad];
#else
   for(fin_pquad=fin_grid->pquad; fin_pquad != NULL; fin_pquad=fin_pquad->next)
     {
#endif
      fin_z = fin_pquad->z;
      coa_z = fin_z/2;
      
      /* find correct coa_pquad */
      for(tcoa_pquad = coa_grid->pquad; coa_z > tcoa_pquad->z + tcoa_pquad->length; tcoa_pquad = tcoa_pquad->next)
         ;
      
      /* jump to correct cquad */
      coa_cquad = tcoa_pquad->loc + (coa_z - tcoa_pquad->z);
      
      
      for(fin_cquad = fin_pquad->loc;
          fin_cquad < fin_pquad->loc + fin_pquad->length; 
          fin_cquad++, fin_z++)  
        {  
         for(ifin_cquad  = fin_cquad; 
             ifin_cquad != NULL; 
             ifin_cquad  = ifin_cquad->next)
           {
            fin_y = ifin_cquad->y;
            coa_y = fin_y/2;
            
            /* find correct coa_cquad */
            for(tcoa_cquad = coa_cquad; coa_y > tcoa_cquad->y + tcoa_cquad->length; tcoa_cquad = tcoa_cquad->next)
               ;
            
            /* jump to correct nquad */
            coa_nquad = tcoa_cquad->loc + (coa_y - tcoa_cquad->y);
            
            for(fin_nquad = ifin_cquad->loc;  
                fin_nquad < ifin_cquad->loc + ifin_cquad->length; 
                fin_nquad++, fin_y++) 
              { 
               for(ifin_nquad  = fin_nquad; 
                   ifin_nquad != NULL; 
                   ifin_nquad  = ifin_nquad->next)
                 {
                  fin_x = ifin_nquad->x;
                  coa_x = fin_x/2;
                  
                  /* find correct coarse nquad */
                  for(tcoa_nquad = coa_nquad; coa_x > tcoa_nquad->x + tcoa_nquad->length; tcoa_nquad = tcoa_nquad->next)
                     ;
                  
                  /* jump to correct (first!) node */
                  coa_node = tcoa_nquad->loc + (coa_x - tcoa_nquad->x);
                  
                  for(fin_node = ifin_nquad->loc; 
                      fin_node < ifin_nquad->loc + ifin_nquad->length; 
                      fin_node++, fin_x++)
                    {
                     
                     /* we now have access to both the fin_node and all surrounding mother coa_nodes */
                     
                     /* average the mass of fin_nodes */
                     coa_node->dens     += 0.125*(fin_node->dens+simu.mean_dens);
                     coa_node->densV[X] += 0.125*(fin_node->densV[X]);
                     coa_node->densV[Y] += 0.125*(fin_node->densV[Y]);
                     coa_node->densV[Z] += 0.125*(fin_node->densV[Z]);
#ifdef UonGrid
                     coa_node->u        += 0.125*(fin_node->u);
#endif
                     
                     /* move to next coa_node */
                     if(is_even(fin_x) == FALSE)
                       {
                        coa_node++;
                        coa_x++;
                       }
                     
                    }
                 }
               
               /* move to next coa_nquad */
               if(is_even(fin_y) == FALSE)
                 {
                  coa_nquad++;
                  coa_y++;
                 }
              }
           }
         
         /* move to next coa_cquad */
         if(is_even(fin_z) == FALSE)
           {
            coa_cquad++;
            coa_z++;
           }
        }
     }
}

/*=============================================================================
* fine_coarse density overlap 
* (this is actually a c2f interpolation but placed into fine_to_coarse.c)
*=============================================================================*/
void fc_overlap(gridls *fin_grid)
{
   gridls       *coa_grid;
   
   long          i, j, k, ipquad;
   
   pqptr         coa_pquad, tcoa_pquad;
   cqptr         coa_cquad, tcoa_cquad;
   nqptr         coa_nquad, tcoa_nquad;
   nptr          coa_node;
   long          coa_x, coa_y, coa_z;
   double        coa_dens;
   nptr          tsc_coanodes[3][3][3];
   
   pqptr         fin_pquad;
   cqptr         fin_cquad, ifin_cquad;
   nqptr         fin_nquad, ifin_nquad;
   nptr          fin_node;
   long          fin_x, fin_y, fin_z;
   nptr          tsc_finnodes[3][3][3];
   
   /* it's your responsibility to make sure that a coarse grid actually exists */
   coa_grid = fin_grid-1;
   
   
   /*==============================================================
    * loop over fine grid (and simultanesouly over coarse grid...)
    *==============================================================*/
#ifdef WITH_OPENMP
#pragma omp parallel private(i, j, k, ipquad, tcoa_pquad, coa_cquad, tcoa_cquad, coa_nquad, tcoa_nquad, coa_node, coa_x, coa_y, coa_z, tsc_coanodes, coa_dens, fin_pquad,  fin_cquad, ifin_cquad, fin_nquad, ifin_nquad, fin_node, fin_x, fin_y, fin_z, tsc_finnodes) shared(fin_grid, coa_grid)
#pragma omp for schedule(static)
   for(ipquad=0; ipquad<fin_grid->no_pquad; ipquad++)
     {
     fin_pquad=fin_grid->pquad_array[ipquad];
#else
   for(fin_pquad=fin_grid->pquad; fin_pquad != NULL; fin_pquad=fin_pquad->next)
     {
#endif
      fin_z = fin_pquad->z;
      coa_z = fin_z/2;
      
      /* find correct coa_pquad */
      for(tcoa_pquad = coa_grid->pquad; coa_z > tcoa_pquad->z + tcoa_pquad->length; tcoa_pquad = tcoa_pquad->next)
         ;
      
      /* jump to correct cquad */
      coa_cquad = tcoa_pquad->loc + (coa_z - tcoa_pquad->z);
      
      
      for(fin_cquad = fin_pquad->loc;
          fin_cquad < fin_pquad->loc + fin_pquad->length; 
          fin_cquad++, fin_z++)  
        {  
         for(ifin_cquad  = fin_cquad; 
             ifin_cquad != NULL; 
             ifin_cquad  = ifin_cquad->next)
           {
            fin_y = ifin_cquad->y;
            coa_y = fin_y/2;
            
            /* find correct coa_cquad */
            for(tcoa_cquad = coa_cquad; coa_y > tcoa_cquad->y + tcoa_cquad->length; tcoa_cquad = tcoa_cquad->next)
               ;
            
            /* jump to correct nquad */
            coa_nquad = tcoa_cquad->loc + (coa_y - tcoa_cquad->y);
            
            for(fin_nquad = ifin_cquad->loc;  
                fin_nquad < ifin_cquad->loc + ifin_cquad->length; 
                fin_nquad++, fin_y++) 
              { 
               for(ifin_nquad  = fin_nquad; 
                   ifin_nquad != NULL; 
                   ifin_nquad  = ifin_nquad->next)
                 {
                  fin_x = ifin_nquad->x;
                  coa_x = fin_x/2;
                  
                  /* find correct coarse nquad */
                  for(tcoa_nquad = coa_nquad; coa_x > tcoa_nquad->x + tcoa_nquad->length; tcoa_nquad = tcoa_nquad->next)
                     ;
                  
                  /* jump to correct (first!) node */
                  coa_node = tcoa_nquad->loc + (coa_x - tcoa_nquad->x);
                  
                  for(fin_node = ifin_nquad->loc; 
                      fin_node < ifin_nquad->loc + ifin_nquad->length; 
                      fin_node++, fin_x++)
                    {
                     /* get actual coa_dens */
                     coa_dens = coa_node->dens + simu.mean_dens;
                     
                     /* we only need to consider those coa_node that carry some remnant density */
                     if(fabs(coa_dens) > ZERO)
                       {
                        /*
                         * do not interpolate the coa_dens to fin_node 
                         * as we cannot calculate a slope for coa_dens!
                         *
                         * => some of the neighboring coa_nodes 
                         *    do not carry any density at all yet...
                         */
                        
                        /* copy remnant coa_dens density to fin_node */
                        fin_node->force.temp[0] = coa_dens;
                        fin_node->densVtmp[X]   = coa_node->densV[X];
                        fin_node->densVtmp[Y]   = coa_node->densV[Y];
                        fin_node->densVtmp[Z]   = coa_node->densV[Z];
#ifdef UonGrid
                        fin_node->force.temp[1] = coa_node->u;
#endif
                       }
                     else
                       {
                        fin_node->force.temp[0] = 0.0;
                        fin_node->densVtmp[X]   = 0.0;
                        fin_node->densVtmp[Y]   = 0.0;
                        fin_node->densVtmp[Z]   = 0.0;
#ifdef UonGrid
                        fin_node->force.temp[1] = 0.0;
#endif
                       }
                     
                     
                     /* move to next coa_node */
                     if(is_even(fin_x) == FALSE)
                       {
                        coa_node++;
                        coa_x++;
                       }
                     
                    }
                 }
               
               /* move to next coa_nquad */
               if(is_even(fin_y) == FALSE)
                 {
                  coa_nquad++;
                  coa_y++;
                 }
              }
           }
         
         /* move to next coa_cquad */
         if(is_even(fin_z) == FALSE)
           {
            coa_cquad++;
            coa_z++;
           }
        }
     }
}



