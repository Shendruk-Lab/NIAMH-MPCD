#ifndef BC_H
#define BC_H

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************ DECLARE FUNCTIONS *********** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/*
   These functions are related to the interaction
   of MPC particles with boundary conditions.
*/

void BC_BCcollision( bc *WALL1,bc *WALL2,double t_step,int *flag );
void MPC_BCcollision( particleMPC *pp,int currentP,bc WALL[],spec *pSP,double KBT,double t_step,int LC,int *bcCNT,int *reCNT,int *rethermCNT, int flagMPConBC );
void BC_MPCcollision(bc WALL[], int BCcurrent, particleMPC *pp, spec *pSP, double KBT, double GRAV[], double t_step, simptr simMD, int MD_mode, int LC, int *bcCNT, int *reCNT, int *rethermCNT);

double calcW( bc WALL,particleMPC P );
double calcWavyW( bc WALL,double POS[], double W );
double calcW_BC( bc movingWall,bc stillWall,int flagCentre );
double calcW_PLANE( bc WALL,particleMPC P );

void shiftBC( double *shift,bc *WALL,particleMPC *pp );
void shiftbackBC( double *shift,bc *WALL );
void MPC_BCrotation( bc *WALL,particleMPC *pp, double sign, int LC );
void rotateBC( bc *WALL,particleMPC *pp, int LC );
void rotatebackBC( bc *WALL,particleMPC *pp, int LC );

void velBC( particleMPC *pp,bc *WALL,double n[3],spec *SP,double KBT );
void posBC( particleMPC *pp,bc WALL,double n[] );
// void BCBCpos( bc *WALL1 ,bc *WALL2,double n[] );
// void BCBCvel( double cp[3],bc *W1,bc *W2,double n[3],double KBT );
double impulse( double n[],double V[],double U[],double Pv[], double Pu[],double Wv[],double Wu[],double invMv,double invMu,double Iv[][3],double Iu[][3],double r[],double E );
void vel_trans( bc *WALL,double VN[],double VT[],double norm[] );

void crosstime( particleMPC p,bc WALL,double *t_pos, double *t_neg,double t_step );
void crosstimeReverse( particleMPC p,bc WALL,double *tc_pos, double *tc_neg,double t_step );
double chooseT( double tstep,double tp,double tn,int p,int *flag );
double secant_time( particleMPC p,bc WALL,double t_step );

// void chooseP( bc WALL,particleMPC *pp,double *time,double *W2,int *chosenP,double t_step,double GRAV[] );
void chooseP( bc WALL,particleMPC *pp,double *chosenW,int *chosenP );
void chooseBC( bc WALL[],int currentP,particleMPC *pp,double *t_min,double *chosenW,int *chosenBC,double time );

double *normal( double *n,bc WALL,double *point,int dimension );
double *normalWavy( double *n,bc WALL,double *point,int dimension );
double *normalNon4foldSymm( double *n,bc WALL,double *point,int dimension );

void rudimentaryPBC( particleMPC *pp,int axis );
void rudimentaryBBBC( particleMPC *pp,int axis );
void rudimentaryPBC_box( particleMPC *pp );
void rudimentaryBBBC_box( particleMPC *pp );
void rudimentaryChannel_x( particleMPC *pp );
void rudimentaryChannel_y( particleMPC *pp );
void rudimentaryChannel_z( particleMPC *pp );

#endif
