#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <libgen.h>

//================================================================================
//
// structure of Cweb binary file:
//--------------------------------
//int32_t     dummy
//uint64_t    Nnodes (=L^3 for the serial version)
//uint64_t    L
//float       BoxSize in Mpc/h
//
//this repeats Nnodes times:
//22 floats    x,y,z,dens,Vx,Vy,Vz,Wx,Wy,Wz,lambda1,lambda2,lambda3,local_shear[]
//
//================================================================================


//=================================================================================
// DEFINEFLAGS
//=================================================================================
//#define WRITE_FULL_SET   // converts all properties found in the binary file to ASCII

//===================================================================
// the headers for the ASCII output files
//===================================================================
#ifdef WRITE_FULL_SET
#define HEADER_STRING    "#x(1) y(2) z(3) dens(4) Vx(5) Vy(6) Vz(7) Wx(8) Wy(9) Wz(10) lambda1(11) lambda2(12) lambda3(13) e1(14-16) e2(17-19) e3(20-22)\n"
#else
#define HEADER_STRING    "#only a selection of the properties in the binary file have been converted to ASCII, check Cweb2ascii.c for details\n"
#endif

//===================================================================
// global variables
//===================================================================
uint64_t Nknots;
uint64_t Nfilaments;
uint64_t Nsheets;
uint64_t Nvoids;
double   lambda_threshold;
uint64_t l1dim;

//===================================================================
// structures
//===================================================================
typedef struct {
  float x;
  float y;    // cell-centre [kpc/h]
  float z;
  float dens; // overdensity
  float Vx;
  float Vy;   // velocity [km/sec]
  float Vz;
  float Wx;
  float Wy;   // vorticity [km/sec]
  float Wz;
  float lambda1;
  float lambda2;
  float lambda3;
  float local_shear[3][3];
#ifdef PWEB
  float pot;
  float pambda1;
  float pambda2;
  float pambda3;
  float pocal_shear[3][3];
#endif
#ifdef UonGrid
  float u;
#endif
} Cweb_t;

//===================================================================
// PROTOTYPES
//===================================================================
void     convert                           (int , char **, char *);
void     convert_Cweb                      (FILE *, FILE *);
int      get_CwebID                        (Cweb_t);

//===================================================================
// routines copied from AHF to keep this code independent
//===================================================================
int     ReadUInt64         (FILE *, uint64_t      *, int);
int     ReadFloat          (FILE *, float         *, int);
#define TRUE  1
#define FALSE 0
#define pow3(x) ((x)*(x)*(x))

/*=================================================================================================
 * write_Cweb(): FEEL FREE TO FILTER UNWANTED PROPERTIES BY COMMENTING THEM OUT
 *=================================================================================================*/
void write_Cweb(Cweb_t Cweb, FILE *fpout)
{
  int CwebID;
  
  CwebID = get_CwebID(Cweb);
  
#ifdef WRITE_FULL_SET
  // the full set
  fprintf(fpout,"%f ",Cweb.x);
  fprintf(fpout,"%f ",Cweb.y);    // cell-centre [kpc/h]
  fprintf(fpout,"%f ",Cweb.z);
  fprintf(fpout,"%f ",Cweb.dens); // overdensity
  fprintf(fpout,"%f ",Cweb.Vx);
  fprintf(fpout,"%f ",Cweb.Vy);   // velocity [km/sec]
  fprintf(fpout,"%f ",Cweb.Vz);
  fprintf(fpout,"%f ",Cweb.Wx);
  fprintf(fpout,"%f ",Cweb.Wy);   // vorticity [km/sec]
  fprintf(fpout,"%f ",Cweb.Wz);
  fprintf(fpout,"%f ",Cweb.lambda1);
  fprintf(fpout,"%f ",Cweb.lambda2);
  fprintf(fpout,"%f ",Cweb.lambda3);
  fprintf(fpout,"%f ",Cweb.local_shear[0][0]);
  fprintf(fpout,"%f ",Cweb.local_shear[1][0]);
  fprintf(fpout,"%f ",Cweb.local_shear[2][0]);
  fprintf(fpout,"%f ",Cweb.local_shear[0][1]);
  fprintf(fpout,"%f ",Cweb.local_shear[1][1]);
  fprintf(fpout,"%f ",Cweb.local_shear[2][1]);
  fprintf(fpout,"%f ",Cweb.local_shear[0][2]);
  fprintf(fpout,"%f ",Cweb.local_shear[1][2]);
  fprintf(fpout,"%f ",Cweb.local_shear[2][2]);
#ifdef UonGrid
  fprintf(fpout,"%f ",Cweb.u);
#endif
#else
  // only write a subset
  fprintf(fpout,"%f ",Cweb.x);
  fprintf(fpout,"%f ",Cweb.y);
  fprintf(fpout,"%f ",Cweb.z);
  fprintf(fpout,"%f ",Cweb.dens);
  fprintf(fpout,"%f ",Cweb.lambda1);
  fprintf(fpout,"%f ",Cweb.lambda2);
  fprintf(fpout,"%f ",Cweb.lambda3);
#ifdef UonGrid
  fprintf(fpout,"%f ",Cweb.u);
#endif
  fprintf(fpout,"%d ",CwebID);
#endif
  
  fprintf(fpout,"\n");
}

/*=================================================================================================
 * main()
 *=================================================================================================*/
int main(argc,argv)
int argc;
char **argv;
{   
  char   outfile[2048], **infile;
  char   prefix[2048];

  int    i, slen, nfiles;

  fprintf(stderr,"======================================================================\n");
  fprintf(stderr,"              convert Cweb binary files to ASCII files\n");
  fprintf(stderr,"     (the code also concatenates all files when run in MPI mode\n");
  fprintf(stderr,"      and gives a quick statistic for volume filling fractions)\n");
  fprintf(stderr,"======================================================================\n");

  if(argc<4)
    {
      fprintf(stderr,"usage: %s prefix Nfiles lambda_threshold\n",*argv);
      exit(1);
    }
  
  //===================================================================
  // deal with command line
  //===================================================================
  strcpy(prefix,argv[1]);
  nfiles  = atoi(argv[2]);
  lambda_threshold = atof(argv[3]);
  
  // prepare array holding filenames
  infile = (char **) calloc(nfiles, sizeof(char *));
  for (i=0; i<nfiles; i++) {
    infile[i] = (char *) calloc(2048, sizeof(char));
  }

  // construct filenames
  if(nfiles == 1) {
    sprintf(infile[0],"%s.Cweb",prefix);
  }
  else {
    for(i=0; i<nfiles; i++) {
      sprintf(infile[i],"%s.%04d.Cweb",prefix,i);
    }
  }
  sprintf(outfile,"%s.Cweb-ascii",basename(prefix));
  
  // be verbose
  fprintf(stderr,"o reading from:\n");
  for(i=0; i<nfiles; i++) {
    fprintf(stderr,"   %s\n",infile[i]);
  }
  fprintf(stderr,"o writing to: %s\n",outfile);
  
  //===================================================================
  // wrapper for file-type conversion routines
  //===================================================================
  convert(nfiles, infile, outfile);
}
  
/*=================================================================================================
 * convert_halos()
 *=================================================================================================*/
void  convert(int nfiles, char **infile, char *outfile)
{
  int i;
  FILE *fpin, *fpout;
  
  // open output file
  fpout = fopen(outfile,"w");
  if(fpout == NULL) {
    fprintf(stderr,"FATAL: cannot open %s for writing\n",outfile);
    exit(0);
  }
  
  // write the header line
  fprintf(fpout,HEADER_STRING);
  
  // for volume fraction statistics
  Nknots           = 0;
  Nfilaments       = 0;
  Nsheets          = 0;
  Nvoids           = 0;
  l1dim            = 0;
  
  // loop over all input files
  for(i=0; i<nfiles; i++) {
    
    // open input file
    fpin = fopen(infile[i],"rb");
    if(fpin == NULL) {
      fprintf(stderr,"FATAL: cannot open %s for writing\n",infile[i]);
      exit(0);
    }
    
    convert_Cweb(fpin,fpout);
    fclose(fpin);
  }
  
  // close output file
  fclose(fpout);
  
  // dump filling fractions
  fprintf(stderr,"\nVolume Filling Fractions: (grid=%"PRIu64", LambdaThreshold=%g)\n",l1dim,lambda_threshold);
  fprintf(stderr,"=========================\n");
  fprintf(stderr,"knots     = %lf\n",(double)Nknots/(double)pow3(l1dim));
  fprintf(stderr,"filaments = %lf\n",(double)Nfilaments/(double)pow3(l1dim));
  fprintf(stderr,"sheets    = %lf\n",(double)Nsheets/(double)pow3(l1dim));
  fprintf(stderr,"voids     = %lf\n",(double)Nvoids/(double)pow3(l1dim));
  fprintf(stderr,"knots     = %"PRIu64"\n",Nknots);
  fprintf(stderr,"filaments = %"PRIu64"\n",Nfilaments);
  fprintf(stderr,"sheets    = %"PRIu64"\n",Nsheets);
  fprintf(stderr,"voids     = %"PRIu64"\n",Nvoids);
}

/*=================================================================================================
 * convert_Cweb()
 *=================================================================================================*/
void  convert_Cweb(FILE *fpin, FILE *fpout)
{
  uint64_t numHalos;
  uint32_t numColumns;
  uint64_t i;
  int32_t  one, Pweb;
  int      swap=0;
  Cweb_t   Cweb;
  uint64_t Nnodes;
  uint64_t L;
  float    BoxSize;

  // figure out swap status
  fread(&one, sizeof(int32_t), 1, fpin);
  if(one == 1)    swap = 0;
  else            swap = 1;

  fread(&Pweb, sizeof(int32_t), 1, fpin);
  ReadUInt64(fpin, &Nnodes,   swap);
  ReadUInt64(fpin, &L,        swap);
  ReadFloat (fpin, &BoxSize,  swap);
  
  // minimal consistency check
  if(l1dim != 0 && l1dim != L) {
    fprintf(stderr,"the files you are merging are not for the same grid level: L=%ld vs. l1dim=%"PRIu64"\nABORTING\n",L,l1dim);
    exit(0);
  }
  
  // we can safely use L
  l1dim = (uint64_t) L;
  
#ifdef VERBOSE
  fprintf(stderr,"o reading %ld cells from file (swap=%d,l1dim=%ld)\n",Nnodes,swap,l1dim);
#endif

  // read in Cweb properties
  for(i=0; i<Nnodes; i++) {
    ReadFloat(fpin, &Cweb.x,           swap);
    ReadFloat(fpin, &Cweb.y,           swap);
    ReadFloat(fpin, &Cweb.z,           swap);
    ReadFloat(fpin, &Cweb.dens,        swap);
    ReadFloat(fpin, &Cweb.Vx,          swap);
    ReadFloat(fpin, &Cweb.Vy,          swap);
    ReadFloat(fpin, &Cweb.Vz,          swap);
    ReadFloat(fpin, &Cweb.Wx,          swap);
    ReadFloat(fpin, &Cweb.Wy,          swap);
    ReadFloat(fpin, &Cweb.Wz,          swap);
    ReadFloat(fpin, &Cweb.lambda1,     swap);
    ReadFloat(fpin, &Cweb.lambda2,     swap);
    ReadFloat(fpin, &Cweb.lambda3,     swap);
    ReadFloat(fpin, &(Cweb.local_shear[0][0]),     swap);
    ReadFloat(fpin, &(Cweb.local_shear[1][0]),     swap);
    ReadFloat(fpin, &(Cweb.local_shear[2][0]),     swap);
    ReadFloat(fpin, &(Cweb.local_shear[0][1]),     swap);
    ReadFloat(fpin, &(Cweb.local_shear[1][1]),     swap);
    ReadFloat(fpin, &(Cweb.local_shear[2][1]),     swap);
    ReadFloat(fpin, &(Cweb.local_shear[0][2]),     swap);
    ReadFloat(fpin, &(Cweb.local_shear[1][2]),     swap);
    ReadFloat(fpin, &(Cweb.local_shear[2][2]),     swap);
#ifdef PWEB
    if(Pweb) {
      ReadFloat(fpin, &Cweb.pot,     swap);
      ReadFloat(fpin, &Cweb.pambda1,     swap);
      ReadFloat(fpin, &Cweb.pambda2,     swap);
      ReadFloat(fpin, &Cweb.pambda3,     swap);
      ReadFloat(fpin, &(Cweb.pocal_shear[0][0]),     swap);
      ReadFloat(fpin, &(Cweb.pocal_shear[1][0]),     swap);
      ReadFloat(fpin, &(Cweb.pocal_shear[2][0]),     swap);
      ReadFloat(fpin, &(Cweb.pocal_shear[0][1]),     swap);
      ReadFloat(fpin, &(Cweb.pocal_shear[1][1]),     swap);
      ReadFloat(fpin, &(Cweb.pocal_shear[2][1]),     swap);
      ReadFloat(fpin, &(Cweb.pocal_shear[0][2]),     swap);
      ReadFloat(fpin, &(Cweb.pocal_shear[1][2]),     swap);
      ReadFloat(fpin, &(Cweb.pocal_shear[2][2]),     swap);
    }
#endif
#ifdef UonGrid
    ReadFloat(fpin, &(Cweb.u),         swap);
#endif

    //=================================================================================
    // write Cweb to ASCII file
    //=================================================================================
    write_Cweb(Cweb, fpout);
    
    
    /* count web elements */
    if(Cweb.lambda3 > lambda_threshold)                                    Nknots++;
    if(Cweb.lambda2 > lambda_threshold && Cweb.lambda3 < lambda_threshold) Nfilaments++;
    if(Cweb.lambda1 > lambda_threshold && Cweb.lambda2 < lambda_threshold) Nsheets++;
    if(Cweb.lambda1 < lambda_threshold)                                    Nvoids++;
  
  } // for(Nnodes)
}


////////////////////////////////////////////////////////////////////////////////////////


/*
 Read a possibly byte swapped unsigned long integer
 */
int ReadUInt64(FILE *fptr, uint64_t *n,int swap)
{
  unsigned char *cptr,tmp;
  
  if(sizeof(uint64_t) == 4)
   {
    if (fread(n,4,1,fptr) != 1)
      return(FALSE);
    if (swap) {
      cptr = (unsigned char *)n;
      tmp     = cptr[0];
      cptr[0] = cptr[3];
      cptr[3] = tmp;
      tmp     = cptr[1];
      cptr[1] = cptr[2];
      cptr[2] = tmp;
    }
   }
  else if(sizeof(uint64_t) == 8)
   {
    if (fread(n,8,1,fptr) != 1)
      return(FALSE);
    if (swap) {
      cptr = (unsigned char *)n;
      tmp     = cptr[0];
      cptr[0] = cptr[7];
      cptr[7] = tmp;
      tmp     = cptr[1];
      cptr[1] = cptr[6];
      cptr[6] = tmp;
      tmp     = cptr[2];
      cptr[2] = cptr[5];
      cptr[5] = tmp;
      tmp     = cptr[3];
      cptr[3] = cptr[4];
      cptr[4] = tmp;
    }
   }
  else
   {
    fprintf(stderr,"ReadUInt64: something wrong...cannot read long\n");
    exit(0);
   }
  
  return(TRUE);
}

/*
 Read a possibly byte swapped floating point number
 Assume IEEE format
 */
int ReadFloat(FILE *fptr,float *n, int swap)
{
  unsigned char *cptr,tmp;
  
  if(sizeof(float) != 4)
   {
    fprintf(stderr,"ReadFloat: sizeof(float)=%ld and not 4\n",sizeof(float));
    exit(0);
   }
  
  if (fread(n,4,1,fptr) != 1)
    return(FALSE);
  if (swap)
   {
    cptr = (unsigned char *)n;
    tmp     = cptr[0];
    cptr[0] = cptr[3];
    cptr[3] = tmp;
    tmp     = cptr[1];
    cptr[1] = cptr[2];
    cptr[2] = tmp;
   }
  return(TRUE);
}

/*=================================================================================================
 * get_CwebID()
 *=================================================================================================*/
int get_CwebID(Cweb_t Cweb)
{
  if(Cweb.lambda3 > lambda_threshold)
    return(3);
  if(Cweb.lambda2 > lambda_threshold && Cweb.lambda3 < lambda_threshold)
    return(2);
  if(Cweb.lambda1 > lambda_threshold && Cweb.lambda2 < lambda_threshold)
    return(1);
  if(Cweb.lambda1 < lambda_threshold)
    return(0);

}

