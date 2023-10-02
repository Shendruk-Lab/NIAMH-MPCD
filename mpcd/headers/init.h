# include<time.h>

#ifndef INIT_H
#define INIT_H

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************ DECLARE FUNCTIONS *********** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/*
   These functions are used to initialize the program
*/
void zero_bc_var( double *tfrac,double *tdiff,int *g );
void zerocnt( double *KBTNOW,double AVNOW[],double *AVS );
void zeroHISTVEC( int HIST[_3D][BINS] );
void zeroHISTSCALAR( int HIST[BINS] );
void zerocell( cell ***CL );
void zeroparticles( particleMPC *pp );
void zeroPressureColl( cell *CL );

void initvar( unsigned long *seed,time_t *to,clock_t *co,int *runtime,int *warmtime,double *sumM,double AV[_3D],double avDIR[_3D],spec SP[],double *C,double *S,double RA,double *AVVEL,double KBT,bc WALL[],double MAG[_3D],cell ***CL,particleMPC *pp );
void place( double Q[],int PL,FILE *fin );
void replace( particleMPC *p );
void push( double V[],double KBT,int PL,double MASS,FILE *fin );
void orient( double U[],double Q[],int PL );
void setcoord(char dir[], spec SP[], particleMPC *pp, double KBT, double VEL[], bc WALL[], simptr simMD, int MD_mode, int LC );
int checkplaceMPC( int IN,particleMPC *pp,spec SP[],bc WALL[] );
void replacePos_WithCheck( particleMPC *pp,bc WALL[] );
int checkplace( int IN,particleMPC *pp,spec SP[],bc WALL[],simptr simMD,int MDmode );

void theory_trans( double *MFP,double *VISC,double *THERMD,double *SDIFF,double *SPEEDOFSOUND,double RA,double FRICCO,double KBT,double dt,double sumM,int RTECH,int SYNOUT,FILE *fsynopsis );
double ndensity( bc WALL[] );
double mdensity( bc WALL[],double MASS );

void openBasic( FILE **fout,char dir[],char filestring[],char fileextension[] );
void openCheckpoint( FILE **fout,char dir[] );
void opendetails( int i,FILE *fdetail[],char dir[],char fileprefix[],char filesuffix[],char fileextension[] );
void opencoarse( FILE **f,char dir[],char fname[],char ext[] );
void openavvel( FILE **f,char dir[],char fname[],char ext[] );
void openorder( FILE **f,char dir[],char fname[],char ext[] );
void openorderQ( FILE **f,char dir[],char fname[],char ext[] );
void openorderQK( FILE **f,char dir[],char fname[],char ext[] );
void openavs( FILE **f,char dir[],char fname[],char ext[] );
void opendensSTD( FILE **f,char dir[],char fname[],char ext[] );
void openhistVel( FILE **f,char dir[],char fname[],char ext[] );
void openhistSpeed( FILE **f,char dir[],char fname[],char ext[] );
void openhistVort( FILE **f,char dir[],char fname[],char ext[] );
void openhistEnstrophy( FILE **f,char dir[],char fname[],char ext[] );
void openhistDir( FILE **f,char dir[],char fname[],char ext[] );
void openhistDens( FILE **f,char dir[],char fname[],char ext[] );
void openhistS( FILE **f,char dir[],char fname[],char ext[] );
void openavenstrophy( FILE **f,char dir[],char fname[],char ext[] );
void openflow( FILE **f,char dir[],char fname[],char ext[] );
void openvel( FILE **f,char dir[],char fname[],char ext[] );
void opentraj( int bc,FILE *fsolids[],char dir[],char filesolids[],char filesuffix[],char fileextension[] );
void openplace( int i,FILE *fsolids[],char dir[],char filesolids[],char filesuffix[],char fileextension[] );
void openenergy( FILE **f,char dir[],char fname[],char ext[] );
void openenergyfield( FILE **f,char dir[],char fname[],char ext[] );
void openenergyneighbours( FILE **f,char dir[],char fname[],char ext[] );
void opensynopsis( FILE **fsynopsis,char dir[],int firsttime );
void opencorr( FILE **f,char dir[],char fname[],char ext[] );
void openenergyspect( FILE **f,char dir[],char fname[],char ext[] );
void openenstrophyspect( FILE **f,char dir[],char fname[],char ext[] );
void opentopo( FILE **f,char dir[],char fname[],char ext[] );
void opendefect( FILE **f,char dir[],char fname[],char ext[] );
void opendisclin( FILE **f,char dir[],char fname[],char ext[] );
void openmultiphase( FILE **f,char dir[],char fname[],char ext[] );
void openpressure( FILE **f,char dir[],char fname[],char ext[] );
void openbinder( FILE **f,char dir[],char fname[],char ext[],int binSize );
void openswimmer( FILE **f,char dir[],char fname[],char ext[] );
void openswimmerOri( FILE **f,char dir[],char fname[],char ext[] );
void openruntumble( FILE **f,char dir[],char fname[],char ext[] );

void checkSim( FILE *fsynopsis,int SYNOUT,inputList in,spec *SP,bc *WALL,specSwimmer SS );

void initOutput( char op[],outputFlagsList *outFlag,outputFilesList *outFile,inputList in,spec *SP, bc WALL[] );
void initializeSIM(cell ***CL, particleMPC *SRDparticles, spec SP[], bc WALL[], simptr simMD, specSwimmer *specS, swimmer *swimmers, int argc, char* argv[], inputList *in, time_t *to, clock_t *co, int *runtime, int *warmtime, double *AVVEL, kinTheory *theory, double *KBTNOW, double *AVS, double *S4, double *stdN, double AVNOW[_3D], double AVV[_3D], double avDIR[_3D], outputFlagsList outFlags, int MD_mode, FILE *fsynopsis, char ip[] );
void initializeRecovery(cell ***CL, particleMPC *SRDparticles, spec SP[], specSwimmer specS, int RTECH, int LC, int MD_mode, int SYNOUT, FILE *fsynopsis );

#endif
