#ifndef POUT_H
#define POUT_H

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************ DECLARE FUNCTIONS *********** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/*
   These functions are related printing out the
   results of the simulation.
*/

// Printf file headers
void outheader( FILE *fout,int SP );
void coordheader( FILE *fout );
void coarseheader( FILE *fout );
void avvelheader( FILE *fout );
void avvelWithGradVelheader( FILE *fout );
void corrheader( FILE *fout );
void energyspectheader( FILE *fout );
void enstrophyspectheader( FILE *fout );
void topoheader( FILE *fout );
void defectheader( FILE *fout );
void disclinTensorheader( FILE *fout );
void multiphaseheader( FILE *fout );
void pressureheader( FILE *fout );
void avsheader( FILE *fout );
void densheader( FILE *fout );
void avenstrophyheader( FILE *fout );
void orderheader( FILE *fout );
void orderQheader( FILE *fout );
void orderQKheader( FILE *fout );
void flowheader( FILE *fout );
void solidsheader( FILE *fout );
void histVelheader( FILE *fout );
void histVortheader( FILE *fout );
void histDirheader( FILE *fout );
void histSpeedheader( FILE *fout );
void histEnstrheader( FILE *fout );
void histNheader( FILE *fout );
void histSheader( FILE *fout );
void energyheader( FILE *fout );
void energyfieldheader( FILE *fout );
void energyneighboursheader( FILE *fout );
void swimmerheader( FILE *fout );
void swimmeroriheader( FILE *fout );
void runtumbleheader( FILE *fout );

// Print data
void coordout( FILE *fout[10],int pr,double T,particleMPC p[],spec SP[] );
void coarseout( FILE *fout,double t,cell ***CL );
void orderout( FILE *fout,double t,cell ***CL,int LC );
void multiphaseout( FILE *fout,double t,cell ***CL );
void pressureout( FILE *fout,double t,cell ***CL );
void orderQout( FILE *fout,double t,cell ***CL,int LC );
void orderQKout( FILE *fout,double t,particleMPC pMPC[],cell ***CL,int LC );
void histVelout( FILE *fout,int vel[_3D][BINS],double minRange,double maxRange,double t );
void histSpeedout( FILE *fout,int speed[BINS],double minRange,double maxRange,double t );
void histVortout( FILE *fout,int vort[_3D][BINS],double minRange,double maxRange,double t );
void histEnstrout( FILE *fout,int enstrophy[BINS],double minRange,double maxRange,double t );
void histDirout( FILE *fout,int dir[_3D][BINS],double minRange,double maxRange,double t );
void histSout( FILE *fout,int S[BINS],double minRange,double maxRange,double t );
void histNout( FILE *fout,int dens[BINS],double minRange,double maxRange,double t );
void avvelout( FILE *fout,double t,double vel[_3D],double KBT );
void avveloutWithGradVel( FILE *fout,double t,double vel[_3D],double KBT,double gradVel[_3D][_3D] );
void avsout( FILE *fout,double t,double S,double S4,double DIR[] );
void densSTDout( FILE *fout,double t,double stdN );
void avenstrophyout( FILE *fout,double t,double E );
void binderout( FILE *fout,double t,double UL );
void flowout( FILE *fout,cell ***CL,int interval );
void velout( FILE *fout,cell ***CL, double t );
void solidout( FILE *fout,bc WALL,double t );
void topoChargeAndDefectsOut( FILE *ftopo,int TOPOOUT,FILE *fdefect,int DEFECTOUT,double t,cell ***CL,double tolD );
void disclinationTensorOut( FILE *fout,double t,cell ***CL,int LC );
void enout( FILE *fout,particleMPC *pp,spec *pSP,bc WALL[],double t,double KBT,double wmf );
void enfieldout( FILE *fout,cell ***CL,spec *SP,double MFPOT,int LC );
void enneighboursout( FILE *fout,double t,cell ***CL,double MFPOT,int LC );
void corrout( FILE *fout,double corr[],double t );
void spectout( FILE *fout,double spect[],double t );
void binderheader( FILE *fout,int binSize );

//Checkpoint
void checkpoint(FILE *fout, inputList in, spec *SP, particleMPC *pSRD, int MD_mode, bc *WALL, outputFlagsList outFlag, int runtime, int warmtime, double AVVEL, double AVS, double avDIR[_3D], double S4, double stdN, double KBTNOW, double AVV[_3D], double AVNOW[_3D], kinTheory theory, specSwimmer specS, swimmer *sw );
void runCheckpoint(char op[500], time_t *lastCheckpoint, FILE *fout, inputList in, spec *SP, particleMPC *pSRD, int MD_mode, bc *WALL, outputFlagsList outFlag, int runtime, int warmtime, double AVVEL, double AVS, double avDIR[_3D], double S4, double stdN, double KBTNOW, double AVV[_3D], double AVNOW[_3D], kinTheory theory, specSwimmer specS, swimmer *sw );

// Terminal printing
void pcoord( particleMPC p );
void mdcoord( particleMD p );
void bccoord( bc WALL );
void swcoord( swimmer sw );
void pvcoord( double POS[_3D],double VEL[_3D],double ANG[_3D],int dimension );
void pall( particleMPC p[] );
void pvec( double VEC[],int dimension );
void ptens3D( double TENS[][_3D] );
void ptens2D( double TENS[][_2D] );
void ptens( double **TENS,int dimension );

void cellout( cell ***CL );
void listout( cell ***CL );

// Listing Initialization
void stateinput( inputList in,spec SP[],bc WALL[],specSwimmer SS,outputFlagsList out,kinTheory theory,FILE *fsynopsis );
void listinput( inputList in,double AVVEL,spec SP[],kinTheory theory );

// Larger output-control routines
void outputResults(cell ***CL, particleMPC *SRDparticles, spec SP[], bc WALL[], simptr simMD, specSwimmer SS, swimmer swimmers[], double AVNOW[_3D], double AVV[_3D], double avDIR[_3D], int runtime, inputList in, double AVVEL, double KBTNOW, double *AVS, double *S4, double *stdN, int MD_mode, outputFlagsList outFlag, outputFilesList outFiles );
void outputHist( cell ***CL,int runtime, inputList in,outputFlagsList outFlag,outputFilesList outFiles );
void closeOutputFiles( spec *SP,bc WALL[],outputFlagsList outFlag,outputFilesList outFiles );
int writeOutput( int t,outputFlagsList f,int GAL,int zeroNetMom );
int writeHistograms( int t,outputFlagsList f );
void printVersionSummary( );

#endif
