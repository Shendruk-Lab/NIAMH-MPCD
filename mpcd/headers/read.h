#ifndef READ_H
#define READ_H

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************ DECLARE FUNCTIONS *********** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/*
   These functions are related reading in the input files
*/
void checkRead( int flag,char failure[],char file[]);
void readin( char fpath[],inputList *in,spec **SP,particleMPC **pSRD,cell ****CL,int *MD_mode );
void readpc( char fpath[],outputFlagsList *out );
void readbc( char fpath[],bc **WALL );
void setPBC( bc *WALL );
void bcin( FILE *fbc,bc *WALL,char fname[] );
void readchckpnt(char fpath[], inputList *in, spec **SP, particleMPC **pSRD, cell ****CL, int *MD_mode, bc **WALL, outputFlagsList *out, int *runtime, int *warmtime, kinTheory *theory, double *AVVEL, double *AVS, double avDIR[_3D], double *S4, double *stdN, double *KBTNOW, double AVV[_3D], double AVNOW[_3D], specSwimmer *specS, swimmer **sw );
void readarg( int argc, char* argv[], char ip[],char op[], int *inMode );
void readJson( char fpath[], inputList *in, spec **SP, particleMPC **pSRD, 
   cell ****CL, int *MDMode, outputFlagsList *out, bc **WALL, 
   specSwimmer *specS, swimmer **sw);

#endif
