#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
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


/*===================================================================================
* coarse_to_fine interpolation of coa_node->dens: ONLY USED FOR DEBUGGING OCTSPLIT
*                                   (however, it can't harm to have the source...)
* BTW, fc_overlap() in fine_to_coarse.c uses the same principle as c2f_dens()!
*===================================================================================*/
void c2f_dens(gridls *coa_grid, gridls *fin_grid)
{
   long          i, j, k, idim, ipquad;
   
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
   
   double        fc_sep[NDIM], slope[NDIM], func[NDIM][NDIM][NDIM];
   
   /*==============================================================
    * loop over fine grid (and simultanesouly over coarse grid...)
    *==============================================================*/
#ifdef WITH_OPENMP
#pragma omp parallel private(ipquad, i, j, k, idim, coa_pquad, tcoa_pquad, coa_cquad, tcoa_cquad, coa_nquad, tcoa_nquad, coa_node, coa_x, coa_y, coa_z, tsc_coanodes, fin_pquad, fin_cquad, ifin_cquad, fin_nquad, ifin_nquad, fin_node, fin_x, fin_y, fin_z, fc_sep, slope, func) shared(fin_grid, coa_grid)
#pragma omp for schedule(static)
   for(ipquad=0; ipquad<fin_grid->no_pquad; ipquad++)
     {
       fin_pquad = fin_grid->pquad_array[ipquad];
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
                  
                  /* jump to correct node */
                  coa_node = tcoa_nquad->loc + (coa_x - tcoa_nquad->x);
                  
                  for(fin_node = ifin_nquad->loc; 
                      fin_node < ifin_nquad->loc + ifin_nquad->length; 
                      fin_node++, fin_x++)
                    {                     
                                          
                     /* distance of fin_node to mother coa_node (in coa_grid units!) */
                     fc_sep[X] = ((double)fin_x - (double)(2*coa_x) - 0.5) / 2.;
                     fc_sep[Y] = ((double)fin_y - (double)(2*coa_y) - 0.5) / 2.;
                     fc_sep[Z] = ((double)fin_z - (double)(2*coa_z) - 0.5) / 2.;
                     
                     /* find all 26 neighbouring coarse nodes */
                     tsc_coanodes[1][1][1] = coa_node;
                     get_TSCnodes(coa_grid, tcoa_pquad, tcoa_cquad, tcoa_nquad, tsc_coanodes, &coa_z, &coa_y, &coa_x);
                     
                     /* calculate slopes in all three dimensions (note: the spacing is 1 in coa_grid units!) */
                     for(k=0; k<NDIM; k++)
                        for(j=0; j<NDIM; j++)
                           for(i=0; i<NDIM; i++)
                              func[k][j][i] = tsc_coanodes[k][j][i]->dens;
                     get_c2fslope(func, slope);
                     
                     /* linear extrapolation to fin_node */
                     fin_node->dens = coa_node->dens + ( slope[X]*fc_sep[X]
                                                       + slope[Y]*fc_sep[Y]
                                                       + slope[Z]*fc_sep[Z]);
                                          
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


