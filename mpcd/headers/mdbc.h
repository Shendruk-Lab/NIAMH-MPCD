#ifndef MDBC_H
#define MDBC_H

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************ DECLARE FUNCTIONS *********** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/*
   These functions are for dealing with hard
   boundary conditions for MD particles
*/
//MD BC routines
void MD_BCcollision(particleMD *atom,bc WALL[],double KBT,double t_step);
void shiftBC_MD( double *shift,bc *WALL,particleMD *atom );
double calcW_MD( bc WALL,particleMD *atom );
void chooseBC_MD( bc WALL[],particleMD *atom,double *t_min,double *chosenW,int *chosenBC,double time,double t_step );
void rewind_MD( particleMD *atom,double time );
void stream_MD( particleMD *atom,double t );
void crosstime_MD( particleMD *atom,bc WALL,double *tc_pos, double *tc_neg,double t_step );
double secant_time_MD( particleMD *atom,bc WALL,double t_step );
double *normal_MD( double *n,bc WALL,particleMD *atom,int dimension );
void velBC_MD( particleMD *atom,bc *WALL,double n[3],double KBT );
void posBC_MD( particleMD *atom,bc WALL,double n[3] );
void MD_BCrotation( bc *WALL,particleMD *atom, double sign );
void rotateBC_MD( bc *WALL,particleMD *atom );
void rotatebackBC_MD( bc *WALL,particleMD *atom );

//Swimmer BC routines
void swimmer_BCcollision( smono *atom,bc WALL[],specSwimmer SS,double t_step,int  flagSwimmeronBC);
void BC_Swimmercollision(bc WALL[],int BCcurrent,swimmer swimmers[],specSwimmer SS,double t_step);
void chooseBC_swimmer( bc WALL[],smono *atom,double *t_min,double *chosenW,int *chosenBC,double time,double t_step );
void chooseNS( bc WALL[],int BCcurrent,specSwimmer SS,swimmer swimmers[],double *chosenW,int *chosenNS,int *HM);
void shiftBC_swimmer( double *shift,bc *WALL,smono *atom );
double calcW_swimmer( bc WALL,smono *atom );
void stream_swimmer( smono *atom,double t );
void rewind_swimmer( smono *atom,double time );
void crosstime_swimmer( smono *atom,bc WALL,double *tc_pos, double *tc_neg,double t_step );
double secant_time_swimmer( smono *atom,bc WALL,double t_step );
double *normal_swimmer( double *n,bc WALL,smono *atom,int dimension );
void velBC_swimmer( smono *atom,bc *WALL,specSwimmer SS,double n[_3D] );
void posBC_swimmer( smono *atom,bc WALL,double n[_3D] );
void swimmer_BCrotation( bc *WALL,smono *atom, double sign );
void rotateBC_swimmer( bc *WALL,smono *atom );
void rotatebackBC_swimmer( bc *WALL,smono *atom );

#endif
