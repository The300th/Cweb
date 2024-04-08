#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <libgen.h>
#include <math.h>

//=================================================================================
// DEFINEFLAGS
//=================================================================================
#define AHF_ASCII  // we still allow AHF_halos to be in ASCII format
//#define DEBUG
#define READ_FOFHALOS

//===================================================================
// the headers for the ASCII output files
//===================================================================
#define HEADER_STRING    "#ID(1) CwebType(2)"

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
#ifdef UonGrid
  float u;
#endif
  int   type; // this is being determined while reading the Cweb file
} Cweb_t;

typedef struct halo {
  uint32_t numColumns;
  uint64_t ID;
  uint64_t hostHalo;
  uint32_t numSubStruct;
  float    Mvir;
  uint32_t npart;
  float    Xc;
  float    Yc;
  float    Zc;
  float    VXc;
  float    VYc;
  float    VZc;
  float    Rvir;
  float    Rmax;
  float    r2;
  float    mbp_offset;
  float    com_offset;
  float    Vmax;
  float    v_esc;
  float    sigV;
  float    lambda;
  float    lambdaE;
  float    Lx;
  float    Ly;
  float    Lz;
  float    b;
  float    c;
  float    Eax;
  float    Eay;
  float    Eaz;
  float    Ebx;
  float    Eby;
  float    Ebz;
  float    Ecx;
  float    Ecy;
  float    Ecz;
  float    ovdens;
  uint32_t nbins;
  float    fMhires;
  float    Ekin;
  float    Epot;
  float    SurfP;
  float    Phi0;
  float    cNFW;
  uint32_t n_gas;
  float    M_gas;
  float    lambda_gas;
  float    lambdaE_gas;
  float    Lx_gas;
  float    Ly_gas;
  float    Lz_gas;
  float    b_gas;
  float    c_gas;
  float    Eax_gas;
  float    Eay_gas;
  float    Eaz_gas;
  float    Ebx_gas;
  float    Eby_gas;
  float    Ebz_gas;
  float    Ecx_gas;
  float    Ecy_gas;
  float    Ecz_gas;
  float    Ekin_gas;
  float    Epot_gas;
  uint32_t n_star;
  float    M_star;
  float    lambda_star;
  float    lambdaE_star;
  float    Lx_star;
  float    Ly_star;
  float    Lz_star;
  float    b_star;
  float    c_star;
  float    Eax_star;
  float    Eay_star;
  float    Eaz_star;
  float    Ebx_star;
  float    Eby_star;
  float    Ebz_star;
  float    Ecx_star;
  float    Ecy_star;
  float    Ecz_star;
  float    Ekin_star;
  float    Epot_star;
  float    mean_z_gas;
  float    mean_z_star;
} halo_t;

//===================================================================
// global variables
//===================================================================
uint64_t Nknots;
uint64_t Nfilaments;
uint64_t Nsheets;
uint64_t Nvoids;
double   lambda_threshold;
uint64_t l1dim;
float    BoxSize;
float    ahf2Cweb;

Cweb_t  *Cweb;
uint64_t CwebNodes;
halo_t  *halo;
uint64_t Nhalos;

//===================================================================
// prototypes and the likes
//===================================================================
#define CwebIndex(i,j,k)  ((i)+(j)*l1dim+(k)*l1dim*l1dim)

void     get_Cweb     (int , char **);
void     get_AHFhalos (int , char **);
void     read_Cweb    (FILE *);
void     read_halos   (FILE *);
void     read_FOFhalos(FILE *);
void     get_ijk      (float, float, float, uint64_t *, uint64_t *, uint64_t *);

//===================================================================
// routines copied from AHF to keep this code independent
//===================================================================
int     ReadUInt32         (FILE *, uint32_t *,      int);
int     ReadUInt64         (FILE *, uint64_t      *, int);
int     ReadFloat          (FILE *, float         *, int);
#define TRUE  1
#define FALSE 0
#define pow3(x) ((x)*(x)*(x))

/*=================================================================================================
 * main()
 *=================================================================================================*/
int main(argc,argv)
int argc;
char **argv;
{   
  char     outfile[2048], **infile;
  char     Cwebprefix[2048], AHFprefix[2048];
  FILE   *fpout;
  int      slen, Cwebnfiles, AHFnfiles;
  int      WebType;
  uint64_t i,j,k,ihalo;

  fprintf(stderr,"======================================================================\n");
  fprintf(stderr,"                    Cweb classification of AHF_halos\n");
  fprintf(stderr,"              (currently only works for regular Cweb grids!)\n");
  fprintf(stderr,"======================================================================\n");

  if(argc<7)
    {
      fprintf(stderr,"usage: %s Cwebprefix Cwebnfiles lambda_threshold AHFprefix AHFnfiles ahf2Cweb\n",*argv);
      exit(1);
    }
  
  //===================================================================
  // deal with command line
  //===================================================================
  strcpy(Cwebprefix,      argv[1]);
  Cwebnfiles       = atoi(argv[2]);
  lambda_threshold = atof(argv[3]);
  strcpy(AHFprefix,       argv[4]);
  AHFnfiles        = atoi(argv[5]);
  ahf2Cweb         = atof(argv[6]);
  
  
  //===================================================================
  // read all Cweb files into memory
  //===================================================================
  // prepare array holding filenames
  infile = (char **) calloc(Cwebnfiles, sizeof(char *));
  for (i=0; i<Cwebnfiles; i++) {
    infile[i] = (char *) calloc(2048, sizeof(char));
  }

  // construct filenames
  if(Cwebnfiles == 1) {
    sprintf(infile[0],"%s.Cweb",Cwebprefix);
  }
  else {
    for(i=0; i<Cwebnfiles; i++) {
      sprintf(infile[i],"%s.%04d.Cweb",Cwebprefix,(int)i);
    }
  }
  
  // be verbose
  fprintf(stderr,"o reading Cweb from:\n");
  for(i=0; i<Cwebnfiles; i++) {
    fprintf(stderr,"   %s\n",infile[i]);
  }
  
  // pass filenames to general purpose reading routine
  get_Cweb(Cwebnfiles, infile);

  // free infile[]
  for (i=0; i<Cwebnfiles; i++) {
    free(infile[i]);
  }
  free(infile);
  
  //===================================================================
  // read all AHF_halos files into memory
  //===================================================================
  // prepare array holding filenames
  infile = (char **) calloc(AHFnfiles, sizeof(char *));
  for (i=0; i<AHFnfiles; i++) {
    infile[i] = (char *) calloc(2048, sizeof(char));
  }
  
  // construct filenames
#ifdef READ_FOFHALOS
  if(AHFnfiles == 1) {
    sprintf(infile[0],"%s",AHFprefix);
  }
  else {
    for(i=0; i<AHFnfiles; i++) {
      sprintf(infile[i],"%s.%04d",AHFprefix,(int)i);
    }
  }
#else
  if(AHFnfiles == 1) {
    sprintf(infile[0],"%s.AHF_halos",AHFprefix);
  }
  else {
    for(i=0; i<AHFnfiles; i++) {
      sprintf(infile[i],"%s.%04d.AHF_halos",AHFprefix,(int)i);
    }
  }
#endif
  
  // be verbose
  fprintf(stderr,"o reading AHF_halos from:\n");
  for(i=0; i<AHFnfiles; i++) {
    fprintf(stderr,"   %s\n",infile[i]);
  }
  
  // pass filenames to general purpose reading routine
  get_AHFhalos(AHFnfiles, infile);
  
  // free infile[]
  for (i=0; i<AHFnfiles; i++) {
    free(infile[i]);
  }
  free(infile);
  
  //===================================================================
  // cross-correlate Cweb and AHF_halos
  //===================================================================
  // prepare filename for output file
  sprintf(outfile,"%s.AHF_Cweb%4.2f",basename(AHFprefix),lambda_threshold);
  fprintf(stderr,"\no writing %"PRIu64" halos to: %s\n",Nhalos,outfile);

  // open output file
  fpout = fopen(outfile,"w");
  if(fpout == NULL) {
    fprintf(stderr,"FATAL: cannot open %s for writing\n",outfile);
    exit(0);
  }

  // classify AHF_halos according to Cweb simultaneously writing output file
  fprintf(fpout,"%s\n",HEADER_STRING);
  for(ihalo=0; ihalo<Nhalos; ihalo++) {
    
#ifdef DEBUG
    fprintf(stderr,"ihalo=%"PRIu64" %f %f %f\n",ihalo,halo[ihalo].Xc,halo[ihalo].Yc,halo[ihalo].Zc);
#endif
    
    get_ijk(halo[ihalo].Xc,halo[ihalo].Yc,halo[ihalo].Zc,&i,&j,&k);
    fprintf(fpout,"%"PRIu64" %d\n",halo[ihalo].ID,Cweb[CwebIndex(i,j,k)].type);
  }
  
  // close output file
  fclose(fpout);
  
  //===================================================================
  // cleaning up
  //===================================================================
  free(halo);
  free(Cweb);

  fprintf(stderr,"END\n");
}
  
/*=================================================================================================
 * get_Cweb()
 *=================================================================================================*/
void  get_Cweb(int nfiles, char **infile)
{
  int i;
  FILE *fpin;
  
  // for volume fraction statistics
  Nknots           = 0;
  Nfilaments       = 0;
  Nsheets          = 0;
  Nvoids           = 0;
  l1dim            = 0;
  
  // loop over all input files
  Cweb    = NULL;
  BoxSize = -1.;
  for(i=0; i<nfiles; i++) {
    
    // open input file
    fpin = fopen(infile[i],"rb");
    if(fpin == NULL) {
      fprintf(stderr,"FATAL: cannot open %s for reading\n",infile[i]);
      exit(0);
    }
    
    // read cosmic web updating Cweb and CwebNodes
    read_Cweb(fpin);
    fclose(fpin);
  }
  
  // double check for inconsistencies
  if(l1dim*l1dim*l1dim != CwebNodes) {
    fprintf(stderr,"you read a different number of nodes than expected: l1dim*l1dim*l1dim=%"PRIu64" vs. CwebNodes=%"PRIu64"\nABORTING\n",l1dim*l1dim*l1dim,CwebNodes);
  }
  
  // be verbose
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
  fprintf(stderr,"BoxSize   = %f\n\n",BoxSize);
}

/*=================================================================================================
 * get_AHFhalos()
 *=================================================================================================*/
void  get_AHFhalos(int nfiles, char **infile)
{
  int i;
  uint64_t ihalo;
  FILE *fpin;
  
  // loop over all input files
  halo   = NULL;
  Nhalos = 0;
  for(i=0; i<nfiles; i++) {
    
    // open input file
#ifdef AHF_ASCII
    fpin = fopen(infile[i],"r");
#else
    fpin = fopen(infile[i],"rb");
#endif
    if(fpin == NULL) {
      fprintf(stderr,"FATAL: cannot open %s for reading\n",infile[i]);
      exit(0);
    }
    
    // read AHF_halos
#ifdef READ_FOFHALOS
    read_FOFhalos(fpin);
#else
    read_halos(fpin);
#endif
    fclose(fpin);
  }
  
#ifdef DEBUG
  for(ihalo=0; ihalo<Nhalos; ihalo++) {
    fprintf(stderr,"get_AHFhalos: %"PRIu64" %f %f %f\n",ihalo,halo[ihalo].Xc,halo[ihalo].Yc,halo[ihalo].Zc);
  }
  exit(0);
#endif
}

/*=================================================================================================
 * read_Cweb()
 *=================================================================================================*/
void read_Cweb(FILE *fpin)
{
  uint64_t numHalos;
  uint32_t numColumns;
  uint64_t i,j,k,n,idx;
  int32_t  one;
  int      swap=0;
  uint64_t Nnodes;
  uint64_t L;
  Cweb_t   Cweb_tmp;
  float    BoxSize_tmp;

  // figure out swap status
  fread(&one, sizeof(int32_t), 1, fpin);
  if(one == 1)    swap = 0;
  else            swap = 1;

  ReadUInt64(fpin, &Nnodes,       swap);
  ReadUInt64(fpin, &L,            swap);
  ReadFloat (fpin, &BoxSize_tmp,  swap);

  // minimal consistency check
  if(l1dim != 0 && l1dim != L) {
    fprintf(stderr,"the files you are merging are not for the same grid level: L=%"PRIu64" vs. l1dim=%"PRIu64"\nABORTING\n",L,l1dim);
    exit(0);
  }
  if(BoxSize > 0. && BoxSize != BoxSize_tmp) {
    fprintf(stderr,"the files you are merging are not for the same boxsize: BoxSize=%f vs. BoxSize_tmp=%f\nABORTING\n",BoxSize,BoxSize_tmp);
    exit(0);
  }

  // we can safely use L and BoxSize_tmp
  l1dim   = (uint64_t) L;
  BoxSize = BoxSize_tmp;
  
  // allocate space for full(!) Cweb
  // TODO: here we assume that Cweb is given on a *regular* grid!!!
  if(Cweb == NULL) {
    Cweb = (Cweb_t *) calloc(l1dim*l1dim*l1dim, sizeof(Cweb_t));
  }
  CwebNodes += Nnodes;
  
#ifdef VERBOSE
  fprintf(stderr,"o reading %"PRIu64" cells from file (swap=%d,l1dim=%"PRIu64")\n",Nnodes,swap,l1dim);
#endif

  // read in Cweb properties
  for(n=0; n<Nnodes; n++) {
    
    // read node into temporary structure
    ReadFloat(fpin, &(Cweb_tmp.x),           swap);
    ReadFloat(fpin, &(Cweb_tmp.y),           swap);
    ReadFloat(fpin, &(Cweb_tmp.z),           swap);
    ReadFloat(fpin, &(Cweb_tmp.dens),        swap);
    ReadFloat(fpin, &(Cweb_tmp.Vx),          swap);
    ReadFloat(fpin, &(Cweb_tmp.Vy),          swap);
    ReadFloat(fpin, &(Cweb_tmp.Vz),          swap);
    ReadFloat(fpin, &(Cweb_tmp.Wx),          swap);
    ReadFloat(fpin, &(Cweb_tmp.Wy),          swap);
    ReadFloat(fpin, &(Cweb_tmp.Wz),          swap);
    ReadFloat(fpin, &(Cweb_tmp.lambda1),     swap);
    ReadFloat(fpin, &(Cweb_tmp.lambda2),     swap);
    ReadFloat(fpin, &(Cweb_tmp.lambda3),     swap);
    ReadFloat(fpin, &(Cweb_tmp.local_shear[0][0]),     swap);
    ReadFloat(fpin, &(Cweb_tmp.local_shear[1][0]),     swap);
    ReadFloat(fpin, &(Cweb_tmp.local_shear[2][0]),     swap);
    ReadFloat(fpin, &(Cweb_tmp.local_shear[0][1]),     swap);
    ReadFloat(fpin, &(Cweb_tmp.local_shear[1][1]),     swap);
    ReadFloat(fpin, &(Cweb_tmp.local_shear[2][1]),     swap);
    ReadFloat(fpin, &(Cweb_tmp.local_shear[0][2]),     swap);
    ReadFloat(fpin, &(Cweb_tmp.local_shear[1][2]),     swap);
    ReadFloat(fpin, &(Cweb_tmp.local_shear[2][2]),     swap);
#ifdef UonGrid
    ReadFloat(fpin, &(Cweb_tmp.u),         swap);
#endif
    
    // calculate position (i,j,k) in actual Cweb[] array from (x,y,z)
    get_ijk(Cweb_tmp.x,Cweb_tmp.y,Cweb_tmp.z, &i,&j,&k);
    
    // copy over to actual Cweb[] array
    idx = CwebIndex(i,j,k);
    Cweb[idx].x                 = Cweb_tmp.x;
    Cweb[idx].y                 = Cweb_tmp.y;
    Cweb[idx].z                 = Cweb_tmp.z;
    Cweb[idx].dens              = Cweb_tmp.dens;
    Cweb[idx].Vx                = Cweb_tmp.Vx;
    Cweb[idx].Vy                = Cweb_tmp.Vy;
    Cweb[idx].Vz                = Cweb_tmp.Vz;
    Cweb[idx].Wx                = Cweb_tmp.Wx;
    Cweb[idx].Wy                = Cweb_tmp.Wy;
    Cweb[idx].Wz                = Cweb_tmp.Wz;
    Cweb[idx].lambda1           = Cweb_tmp.lambda1;
    Cweb[idx].lambda2           = Cweb_tmp.lambda2;
    Cweb[idx].lambda3           = Cweb_tmp.lambda3;
    Cweb[idx].local_shear[0][0] = Cweb_tmp.local_shear[0][0];
    Cweb[idx].local_shear[1][0] = Cweb_tmp.local_shear[1][0];
    Cweb[idx].local_shear[2][0] = Cweb_tmp.local_shear[2][0];
    Cweb[idx].local_shear[0][1] = Cweb_tmp.local_shear[0][1];
    Cweb[idx].local_shear[1][1] = Cweb_tmp.local_shear[1][1];
    Cweb[idx].local_shear[2][1] = Cweb_tmp.local_shear[2][1];
    Cweb[idx].local_shear[0][2] = Cweb_tmp.local_shear[0][2];
    Cweb[idx].local_shear[1][2] = Cweb_tmp.local_shear[1][2];
    Cweb[idx].local_shear[2][2] = Cweb_tmp.local_shear[2][2];
#ifdef UonGrid
    Cweb[idx].u                 = Cweb_tmp.u;
#endif
    
    /* count web elements */
    if(Cweb_tmp.lambda3 > lambda_threshold) {
      Nknots++;
      Cweb[idx].type = 3;
    }
    if(Cweb_tmp.lambda2 > lambda_threshold && Cweb_tmp.lambda3 < lambda_threshold) {
      Nfilaments++;
      Cweb[idx].type = 2;
    }
    if(Cweb_tmp.lambda1 > lambda_threshold && Cweb_tmp.lambda2 < lambda_threshold) {
      Nsheets++;
      Cweb[idx].type = 1;
    }
    if(Cweb_tmp.lambda1 < lambda_threshold) {
      Nvoids++;
      Cweb[idx].type = 0;
    }
  
  } // for(Nnodes)
}


////////////////////////////////////////////////////////////////////////////////////////


/*
 Read a possibly byte swapped unsigned integer
 */
int ReadUInt32(FILE *fptr,uint32_t *n,int swap)
{
  unsigned char *cptr,tmp;
  
  if(sizeof(int) != 4)
   {
    fprintf(stderr,"ReadInt: sizeof(int)=%ld and not 4\n",sizeof(int));
    exit(0);
   }
  
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
  return(TRUE);
}

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
 * read_halos()
 *=================================================================================================*/
void read_halos(FILE *fpin)
{
  uint64_t numHalos;
  uint32_t numColumns;
  uint64_t i;
  int32_t  one;
  int      swap=0;
  
#ifdef AHF_ASCII
  char     headerstring[2048], buffer[2048];
  fgets(headerstring,2048,fpin);
  i        = Nhalos;
  while(fgets(buffer,2048,fpin) != NULL) {

    // provide space for one more halo
    halo = (halo_t *) realloc(halo, (Nhalos+1)*sizeof(halo_t));

    // read that additional halo
    sscanf(buffer,"%"SCNi64" %"SCNi64" %d %f %d %f %f %f",&(halo[i].ID),&(halo[i].hostHalo),&(halo[i].numSubStruct),&(halo[i].Mvir),&(halo[i].npart),&(halo[i].Xc),&(halo[i].Yc),&(halo[i].Zc));
    
    // re-scale to Cweb units
    halo[i].Xc *= ahf2Cweb;
    halo[i].Yc *= ahf2Cweb;
    halo[i].Zc *= ahf2Cweb;
    
    // move to next halo
    i++;
    Nhalos++;
  }
#else // AHF_ASCII
  // figure out swap status
  fread(&one, sizeof(int32_t), 1, fpin);
  if(one == 1)    swap = 0;
  else            swap = 1;
  
  ReadUInt64(fpin, &numHalos,   swap);
  ReadUInt32(fpin, &numColumns, swap);
  
  // update number of haloes reallocating array to hold them all
  Nhalos += numHalos;
  halo = (halo_t *) realloc(halo, Nhalos*sizeof(halo_t));
  
#ifdef VERBOSE
  fprintf(stderr,"o reading %ld halos from file (swap=%d,numColumns=%d)\n",(long unsigned)numHalos,swap,numColumns);
#endif
  
  // read in halo properties
  for(i=Nhalos-numHalos; i<Nhalos; i++) {
    ReadUInt64(fpin, &(halo[i].ID),           swap);    // ID(1)
    ReadUInt64(fpin, &(halo[i].hostHalo),     swap);    // hostHalo(2)
    ReadUInt32(fpin, &(halo[i].numSubStruct), swap);    // numSubStruct(3)
    ReadFloat(fpin, &(halo[i].Mvir),         swap);    // Mvir(4)
    ReadUInt32(fpin, &(halo[i].npart),        swap);    // npart(5)
    ReadFloat(fpin, &(halo[i].Xc),           swap);    // Xc(6)
    ReadFloat(fpin, &(halo[i].Yc),           swap);    // Yc(7)
    ReadFloat(fpin, &(halo[i].Zc),           swap);    // Zc(8)
    ReadFloat(fpin, &(halo[i].VXc),          swap);    // VXc(9)
    ReadFloat(fpin, &(halo[i].VYc),          swap);    // VYc(10)
    ReadFloat(fpin, &(halo[i].VZc),          swap);    // VZc(11)
    ReadFloat(fpin, &(halo[i].Rvir),         swap);    // Rvir(12)
    ReadFloat(fpin, &(halo[i].Rmax),         swap);    // Rmax(13)
    ReadFloat(fpin, &(halo[i].r2),           swap);    // r2(14)
    ReadFloat(fpin, &(halo[i].mbp_offset),   swap);    // mbp_offset(15)
    ReadFloat(fpin, &(halo[i].com_offset),   swap);    // com_offset(16)
    ReadFloat(fpin, &(halo[i].Vmax),         swap);    // Vmax(17)
    ReadFloat(fpin, &(halo[i].v_esc),        swap);    // v_esc(18)
    ReadFloat(fpin, &(halo[i].sigV),         swap);    // sigV(19)
    ReadFloat(fpin, &(halo[i].lambda),       swap);    // lambda(20)
    ReadFloat(fpin, &(halo[i].lambdaE),      swap);    // lambdaE(21)
    ReadFloat(fpin, &(halo[i].Lx),           swap);    // Lx(22)
    ReadFloat(fpin, &(halo[i].Ly),           swap);    // Ly(23)
    ReadFloat(fpin, &(halo[i].Lz),           swap);    // Lz(24)
    ReadFloat(fpin, &(halo[i].b),            swap);    // b(25)
    ReadFloat(fpin, &(halo[i].c),            swap);    // c(26)
    ReadFloat(fpin, &(halo[i].Eax),          swap);    // Eax(27)
    ReadFloat(fpin, &(halo[i].Eay),          swap);    // Eay(28)
    ReadFloat(fpin, &(halo[i].Eaz),          swap);    // Eaz(29)
    ReadFloat(fpin, &(halo[i].Ebx),          swap);    // Ebx(30)
    ReadFloat(fpin, &(halo[i].Eby),          swap);    // Eby(31)
    ReadFloat(fpin, &(halo[i].Ebz),          swap);    // Ebz(32)
    ReadFloat(fpin, &(halo[i].Ecx),          swap);    // Ecx(33)
    ReadFloat(fpin, &(halo[i].Ecy),          swap);    // Ecy(34)
    ReadFloat(fpin, &(halo[i].Ecz),          swap);    // Ecz(35)
    ReadFloat(fpin, &(halo[i].ovdens),       swap);    // ovdens(36)
    ReadUInt32(fpin, &(halo[i].nbins),        swap);    // nbins(37)
    ReadFloat(fpin, &(halo[i].fMhires),      swap);    // fMhires(38)
    ReadFloat(fpin, &(halo[i].Ekin),         swap);    // Ekin(39)
    ReadFloat(fpin, &(halo[i].Epot),         swap);    // Epot(40)
    ReadFloat(fpin, &(halo[i].SurfP),        swap);    // SurfP(41)
    ReadFloat(fpin, &(halo[i].Phi0),         swap);    // Phi0(42)
    ReadFloat(fpin, &(halo[i].cNFW),         swap);    // cNFW(43)
    if(numColumns > 43) {
      ReadUInt32(fpin, &(halo[i].n_gas),       swap);    // n_gas(44)
      ReadFloat(fpin, &(halo[i].M_gas),       swap);    // M_gas(45)
      ReadFloat(fpin, &(halo[i].lambda_gas),  swap);    // lambda_gas(46)
      ReadFloat(fpin, &(halo[i].lambdaE_gas), swap);    // lambdaE_gas(47)
      ReadFloat(fpin, &(halo[i].Lx_gas),      swap);    // Lx_gas(48)
      ReadFloat(fpin, &(halo[i].Ly_gas),      swap);    // Ly_gas(49)
      ReadFloat(fpin, &(halo[i].Lz_gas),      swap);    // Lz_gas(50)
      ReadFloat(fpin, &(halo[i].b_gas),       swap);    // b_gas(51)
      ReadFloat(fpin, &(halo[i].c_gas),       swap);    // c_gas(52)
      ReadFloat(fpin, &(halo[i].Eax_gas),     swap);    // Eax_gas(53)
      ReadFloat(fpin, &(halo[i].Eay_gas),     swap);    // Eay_gas(54)
      ReadFloat(fpin, &(halo[i].Eaz_gas),     swap);    // Eaz_gas(55)
      ReadFloat(fpin, &(halo[i].Ebx_gas),     swap);    // Ebx_gas(56)
      ReadFloat(fpin, &(halo[i].Eby_gas),     swap);    // Eby_gas(57)
      ReadFloat(fpin, &(halo[i].Ebz_gas),     swap);    // Ebz_gas(58)
      ReadFloat(fpin, &(halo[i].Ecx_gas),     swap);    // Ecx_gas(59)
      ReadFloat(fpin, &(halo[i].Ecy_gas),     swap);    // Ecy_gas(60)
      ReadFloat(fpin, &(halo[i].Ecz_gas),     swap);    // Ecz_gas(61)
      ReadFloat(fpin, &(halo[i].Ekin_gas),    swap);    // Ekin_gas(62)
      ReadFloat(fpin, &(halo[i].Epot_gas),    swap);    // Epot_gas(63)
      ReadUInt32(fpin, &(halo[i].n_star),       swap);    // n_star(64)
      ReadFloat(fpin, &(halo[i].M_star),       swap);    // M_star(65)
      ReadFloat(fpin, &(halo[i].lambda_star),  swap);    // lambda_star(66)
      ReadFloat(fpin, &(halo[i].lambdaE_star), swap);    // lambdaE_star(67)
      ReadFloat(fpin, &(halo[i].Lx_star),      swap);    // Lx_star(68)
      ReadFloat(fpin, &(halo[i].Ly_star),      swap);    // Ly_star(69)
      ReadFloat(fpin, &(halo[i].Lz_star),      swap);    // Lz_star(70)
      ReadFloat(fpin, &(halo[i].b_star),       swap);    // b_star(71)
      ReadFloat(fpin, &(halo[i].c_star),       swap);    // c_star(72)
      ReadFloat(fpin, &(halo[i].Eax_star),     swap);    // Eax_star(73)
      ReadFloat(fpin, &(halo[i].Eay_star),     swap);    // Eay_star(74)
      ReadFloat(fpin, &(halo[i].Eaz_star),     swap);    // Eaz_star(75)
      ReadFloat(fpin, &(halo[i].Ebx_star),     swap);    // Ebx_star(76)
      ReadFloat(fpin, &(halo[i].Eby_star),     swap);    // Eby_star(77)
      ReadFloat(fpin, &(halo[i].Ebz_star),     swap);    // Ebz_star(78)
      ReadFloat(fpin, &(halo[i].Ecx_star),     swap);    // Ecx_star(79)
      ReadFloat(fpin, &(halo[i].Ecy_star),     swap);    // Ecy_star(80)
      ReadFloat(fpin, &(halo[i].Ecz_star),     swap);    // Ecz_star(81)
      ReadFloat(fpin, &(halo[i].Ekin_star),    swap);    // Ekin_star(82)
      ReadFloat(fpin, &(halo[i].Epot_star),    swap);    // Epot_star(83)
    }
    if(numColumns > 83) {
      ReadFloat(fpin, &(halo[i].mean_z_gas),    swap);    // mean_z_gas(84)
      ReadFloat(fpin, &(halo[i].mean_z_star),   swap);    // mean_z_star(85)
    }
    
    halo[i].numColumns = numColumns;
    
  } // for(numHalos)
#endif // AHF_ASCII
}


/*=================================================================================================
 * read_FOFhalos()
 *=================================================================================================*/
void read_FOFhalos(FILE *fpin)
{
  uint64_t i, ID, npart, idummy;
  double   Mvir, Xc, Yc, Zc, VXc, VYc, VZc;
  
  char     headerstring[2048], buffer[2048];
  fgets(headerstring,2048,fpin);
  fgets(headerstring,2048,fpin);
  fgets(headerstring,2048,fpin);
  i = Nhalos;
#ifdef DEBUG
  fprintf(stderr,"%"PRIu64" %"PRIu64"\n",i,Nhalos);
#endif
  while(fgets(buffer,2048,fpin) != NULL) {
    
    // provide space for one more halo
    halo = (halo_t *) realloc(halo, (Nhalos+1)*sizeof(halo_t));
    
    // read that additional halo
    sscanf(buffer,"%"SCNi64" %lf %lf %lf %lf %lf %lf %"SCNi64" %lf",
           &ID,
           &Xc, &Yc, &Zc,
           &VXc, &VYc, &VZc,
           &npart,
           &Mvir             );
    
    // store halo in halo[] array
    halo[i].ID           = ID;
    halo[i].hostHalo     = 0;
    halo[i].numSubStruct = 0;
    halo[i].Mvir         = Mvir;
    halo[i].npart        = npart;
    halo[i].Xc           = Xc*ahf2Cweb;
    halo[i].Yc           = Yc*ahf2Cweb;
    halo[i].Zc           = Zc*ahf2Cweb;
    
#ifdef DEBUG
    fprintf(stderr,"i=%"PRIu64"  Nhalos=%"PRIu64"  %lf %lf %lf\n",i,Nhalos,halo[i].Xc,halo[i].Yc,halo[i].Zc);
#endif
    
    // move to next halo
    i++;
    Nhalos++;
  }
  
#ifdef DEBUG
  for(i=0; i<Nhalos; i++) {
    fprintf(stderr,"read_FOFhalos: %"PRIu64" %f %f %f\n",i,halo[i].Xc,halo[i].Yc,halo[i].Zc);
  }
#endif

}

/*=================================================================================================
 * get_ijk()
 *                  NOTE: makes use of global variable l1dim!
 *=================================================================================================*/
void get_ijk(float x, float y, float z, uint64_t *i, uint64_t *j, uint64_t *k)
{
  *i = (uint64_t) floor(x*((float)l1dim-0.5)/BoxSize);
  *j = (uint64_t) floor(y*((float)l1dim-0.5)/BoxSize);
  *k = (uint64_t) floor(z*((float)l1dim-0.5)/BoxSize);
  
#ifdef DEBUG
  fprintf(stderr,"i=%"PRIu64"  j=%"PRIu64"  k=%"PRIu64"  (x=%f y=%f z=%f, l1dim=%"PRIu64")\n",*i,*j,*k,x,y,z,l1dim);
#endif
}
