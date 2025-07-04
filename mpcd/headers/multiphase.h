#ifndef MULTIPHASE_H
#define MULTIPHASE_H

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************ DECLARE FUNCTIONS *********** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/*
   These functions add a multiphase twist to the NIAMH-MPCD algorithm.
*/

void multiphaseColl(cell *CL, spec *SP, specSwimmer SS, int multiphaseMode, double KBT, int MD_mode, double *CLQ, int outP );
void multiphaseCollPoint(cell *CL, spec *SP, specSwimmer SS, double KBT, int MD_mode, double *CLQ, int outP );

#endif
