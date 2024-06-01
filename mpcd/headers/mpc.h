#ifndef MPC_H
#define MPC_H

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************ DECLARE FUNCTIONS *********** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/*
   These functions are the simple core of the
   MPC/SRD algorithm.
*/

void timestep(cell ***CL, particleMPC *SRDparticles, spec SP[], bc WALL[], simptr simMD, specSwimmer *SS, swimmer swimmers[], double AVNOW[_3D], double AVV[_3D], double avDIR[_3D], inputList in, double *KBTNOW, double *AVS, int runtime, int MD_mode, outputFlagsList outFlags, outputFilesList outFiles);

double trans( double t,double V,double QOLD );
double acc( double t,double G,double V_OLD );

void stream_all( particleMPC *pp,double t );
void stream_P( particleMPC *p,double t );
void stream_BC( bc *WALL,double t );
void spin_BC( bc *WALL,double t );

double rewind_trans( double t,double V,double P );
double rewind_acc( double t,double G,double V );
void rewind_P( particleMPC *p,double T );
void rewind_BC( bc *WALL,double time );

void acc_BC( bc *WALL,double t,double GRAV[] );
void acc_P( particleMPC *p,double t,double GRAV[] );
void acc_all( particleMPC *pp,double t,double GRAV[] );
void gridShift_all( double SHIFT[],int shiftBack,particleMPC *SRDparticles,bc WALL[],simptr simMD,swimmer swimmers[],int MD_mode );
void rotate_CL( cell CL,spec *SP,double r0[],double n0[],double dw );

void binin( particleMPC p[],cell ***CL );
void bin( cell ***CL,spec *SP,bc WALL[],double KBT,int LC,int shifted );
void bininMD( simptr sim,cell ***CL );
void binMD( cell ***CL );
void addlink( cell *CL,particleMPC *p );
void removelink( particleMPC *current,cell *CL );
void addlinkMD( cell *CL,particleMD *p );
void removelinkMD( particleMD *current,cell *CL );

void rotate( int RT,double Cos,double Sin,double V[3],double L[3],long SIGN,int RAND );
void cellVelForce( cell *CL,double addVel[3] );
void cellVelSet( cell *CL,double vel[3] );

void stochrotMPC( cell *CL,int RTECH,double C,double S,double *CLQ,int outP );
void andersenMPC( cell *CL,spec *SP,specSwimmer SS,double KBT,double *CLQ,int outP );
void andersenROT( cell *CL,spec *SP,specSwimmer SS,double KBT,double *CLQ,int outP );
void andersenMULTIPHASE( cell *CL,spec *SP,specSwimmer SS,double KBT,double *CLQ,int outP );
void langevinMPC( cell *CL,spec *SP,specSwimmer SS,double KBT,double FRICCO,double Step,double *CLQ,int outP );
void vicsek( cell *CL,spec *SP,double *CLQ,int outP );
void activeSRD( cell *CL,spec *SP,int RTECH,double C,double S,double *CLQ,int outP );
void vicsekAndersenMPC( cell *CL,spec *SP,double KBT, double RELAX,double *CLQ,int outP );
void vicsekLangevinMPC( cell *CL,spec *SP,double KBT,double FRICCO,double Step,double RELAX,double *CLQ,int outP );
void chateAndersenMPC( cell *CL,spec *SP,double KBT,double RELAX,double *CLQ,int outP );
void chateLangevinMPC( cell *CL,spec *SP,double KBT,double FRICCO,double Step,double RELAX,double *CLQ,int outP );
void dipoleAndersenMPC( cell *CL,spec *SP,double KBT,double RELAX,double *CLQ,int outP );
void MPCcollision(cell *CL, spec *SP, specSwimmer SS, double KBT, int RTECH, double C, double S, double FRICCO, double TimeStep, int MD_mode, int LC, double RELAX, double *CLQ, int outP );

void incompColl(cell *CL, spec *SP, specSwimmer SS, int INCOMPmode, int MD_mode, double *CLQ, int outP );
void incompAddVirial( cell *CL,double virialCoB, double virialCoC, double virialCoD, spec *SP,specSwimmer SS );
void incompSwap( cell *CL,spec *SP,specSwimmer SS );
void incompSubtractDivergence( cell *CL,spec *SP,specSwimmer SS );

void multiphaseColl(cell *CL, spec *SP, specSwimmer SS, int multiphaseMode, double KBT, int MD_mode, double *CLQ, int outP );
void multiphaseCollPoint(cell *CL, spec *SP, specSwimmer SS, double KBT, int MD_mode, double *CLQ, int outP );

void activeMD(simptr simMD, cell ***CL, spec *SP, inputList in);

void localVCM( double vcm[_3D],cell CL,spec *SP,specSwimmer specS);
void localMPCVCM( double vcm[_3D],cell CL,spec *SP );
double localMASS( cell CL,spec *SP,specSwimmer specS );
double localTEMP( cell CL,spec *SP,specSwimmer specS );
int localPOP( cell CL );
void localPROP( cell ***CL,spec *SP,specSwimmer specS,int RTECH,int LC );
void sumFLOW( cell ***CL );
void sumSWFLOW( cell ***CL, swimmer *sw, specSwimmer *ss);
void localFLOW( cell ***CL,spec *SP );
void localCM( cell *CL,spec *SP,specSwimmer specS );
void localCM_SRD( cell CL,spec *SP,double r_cm[] );
void localMomInertiaTensor( cell *CL,spec *SP,specSwimmer specS );
double localMomInertia_SRD( cell CL,spec *SP,double r0[],double n[] );
void ghostPart( cell ***CL,bc WALL[],double KBT,int LC, spec *SP);

void scramble( particleMPC *p );
void checkEscape_all( particleMPC *pp );

void calcPressureStreaming( cell ***CL,spec *SP );
void normPressureColl( cell *CL,double dt );
void calcPressureColl_preColl( double *relQ,double *dp,particleMPC *,double *CLQ );
void calcPressureColl_postColl( double *relQ,double *dp,double M,double *vel,cell *CL );

void checkParticleNaN(particleMPC p);
void checkAllParticlesNaN(particleMPC *SRDparticles);

#endif
