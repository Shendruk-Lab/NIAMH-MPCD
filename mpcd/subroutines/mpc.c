///
/// @file
/// 
/// @brief This file contains the core of the MPCD algorithm.
///
/// The file does the majority of MPCD functions. 
/// This includes timestep(), which is iterated over repeatedly in main(), as well as the streaming and collision operators. 
/// It also includes the binning and the calculation of all local (cell-level) fluid properties.
///
/// @see main()
/// @see timestep()
/// @see MPCcollision()
/// @see binin()
///

# include<math.h>
# include<stdio.h>
# include<stdlib.h>

# include "../headers/definitions.h"
# include "../headers/globals.h"
# include "../headers/SRDclss.h"
# include "../headers/mpc.h"
# include "../headers/rand.h"
# include "../headers/mtools.h"
# include "../headers/pout.h"
# include "../headers/init.h"
# include "../headers/bc.h"
# include "../headers/lc.h"
# include "../headers/therm.h"
# include "../headers/swimmers.h"
# include "../headers/mdbc.h"
# include "../headers/ctools.h"

# include "../../md/mdsrd.h"

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ****************** LISTS ***************** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */

/// 
/// @brief This function calculates the properties of all the cells.
///
/// The function calculates properties of all the cell. This includes the cell-level
/// - population (number of particles) 
/// - mass 
/// - centre of mass position
/// - centre of mass velocity
/// - moment of inertia tensor
/// - nematic order parameter tensor
/// - velocity gradient tensor.
/// It calculates all these properties by looping through the linked lists. 
/// @param CL All of the MPCD cells (including the linked list of particles in each cell).
/// @param SP The species-wide information about MPCD particles.
/// @param specS The species-wide information about swimmers.
/// @param RTECH The MPCD collision operator. See `definitions.h` for all options.
/// @param LC Flags whether or not the nematic liquid crystal is turned on.
/// @note localPROP() calculates <b>all</b> local parameters. For single properties, other routines exist. See localMASS() for example
///
void localPROP( cell ***CL,spec *SP,specSwimmer specS,int RTECH,int LC ) {
	int a,b,c,d,id;
	int i;
	double V[_3D],Q[_3D];
	double **S,eigval[_3D];
	double mass;
	particleMPC *pMPC;	//Temporary pointer to MPC particles
	particleMD *pMD;	//Temporary pointer to MD particles
	smono *pSW;			//Temporary pointer to swimmer monomers

	// int flag on whether to compute CM or not
	int computeCM = (RTECH==RAT) || (LC!=ISOF) || (RTECH==DIPOLE_VCM) || (RTECH==DIPOLE_DIR_SUM) || (RTECH==DIPOLE_DIR_AV);

	// Zero
	for( d=0; d<_3D; d++ ) V[d] = 0.0;
	for( d=0; d<_3D; d++ ) eigval[d] = 0.0;

	// Loop over all cells
	for( a=0; a<XYZ_P1[0]; a++ ) for( b=0; b<XYZ_P1[1]; b++ ) for( c=0; c<XYZ_P1[2]; c++ ) {
		// Idea here is to compute POP, MASS, and VCM
		//	CM is only computed if needed

		//Zero everything for recounting
		CL[a][b][c].POP = 0;
		CL[a][b][c].POPSRD = 0;
		CL[a][b][c].POPSW = 0;
		CL[a][b][c].POPMD = 0;
		CL[a][b][c].MASS = 0.0;
		for( d=0; d<DIM; d++ ) CL[a][b][c].VCM[d] = 0.0;
		for( d=0; d<NSPECI; d++ ) CL[a][b][c].SP[d] = 0;
		if (computeCM) for( d=0; d<DIM; d++ ) CL[a][b][c].CM[d] = 0.0;

		//Find local values
		if( CL[a][b][c].pp!=NULL || ( CL[a][b][c].MDpp!=NULL && MDmode==MDinMPC ) || CL[a][b][c].sp!=NULL ) {
			// SRD particles
			if( CL[a][b][c].pp!=NULL ) {
				pMPC = CL[a][b][c].pp;
				while(pMPC!=NULL) {
					CL[a][b][c].POPSRD ++;
					id = pMPC->SPID;
					mass = (SP+id)->MASS;
					CL[a][b][c].MASS += mass;
					CL[a][b][c].SP[id] ++;
					for( d=0; d<DIM; d++ ) CL[a][b][c].VCM[d] +=  pMPC->V[d] * mass;

					if (computeCM){ // add positions if necessary
						for( d=0; d<DIM; d++ ) Q[d] = pMPC->Q[d];
						for( d=0; d<DIM; d++ ) CL[a][b][c].CM[d] += Q[d] * mass;
					}

					//Increment link in list
					pMPC = pMPC->next;
				}
			}
			// MD particles
			if(MDmode==MDinMPC) if( CL[a][b][c].MDpp!=NULL) {
				pMD = CL[a][b][c].MDpp;
				while( pMD!=NULL ) {
					CL[a][b][c].POPMD ++;
					mass = pMD->mass;
					CL[a][b][c].MASS += mass;
					V[0] = pMD->vx;
					if( DIM >= _2D ) V[1] = pMD->vy;
					if( DIM >= _3D ) V[2] = pMD->vz;
					for( d=0; d<DIM; d++ ) CL[a][b][c].VCM[d] += V[d] * mass;

					if (computeCM){
						Q[0] = pMD->rx;
						if( DIM > 1 ) Q[1] = pMD->ry;
						if( DIM > 2 ) Q[2] = pMD->rz;
						for( d=0; d<DIM; d++ ) CL[a][b][c].CM[d] += Q[d] * mass;
					}

					//Increment link in list
					pMD = pMD->nextSRD;
				}
			}
			// Swimmer particles
			if( CL[a][b][c].sp!=NULL) {
				pSW = CL[a][b][c].sp;
				while( pSW!=NULL ) {
					CL[a][b][c].POPSW ++;
					if( pSW->HorM ) mass = (double) specS.middM;
					else mass = (double) specS.headM;
					CL[a][b][c].MASS += mass;
					for( d=0; d<DIM; d++ ) CL[a][b][c].VCM[d] +=  pSW->V[d] * mass;

					if (computeCM){
						for( d=0; d<DIM; d++ ) Q[d] = pSW->Q[d];
						for( d=0; d<DIM; d++ ) CL[a][b][c].CM[d] += Q[d] * mass;
					}

					//Increment link in list
					pSW = pSW->next;
				}
			}
		}
		// Make sums into averages
		if( CL[a][b][c].MASS>0.0 ) for( d=0; d<DIM; d++ ){ 
			CL[a][b][c].VCM[d] /= CL[a][b][c].MASS;
			if (computeCM) CL[a][b][c].CM[d] /= CL[a][b][c].MASS;
		}
		CL[a][b][c].POP = CL[a][b][c].POPSRD + CL[a][b][c].POPSW + CL[a][b][c].POPMD;
	}
	//Calculate moment of inertia
	if( RTECH==RAT || LC!=ISOF ) {
		// moment of inertia to be about the CM of cell
		for( a=0; a<XYZ_P1[0]; a++ ) for( b=0; b<XYZ_P1[1]; b++ ) for( c=0; c<XYZ_P1[2]; c++ ) {
			localMomInertiaTensor( &CL[a][b][c],SP,specS );
		}
	}
	// Find the order parameter tensor, the director and the scalar order parameter
	if( LC!=ISOF || RTECH==CHATE || RTECH==CHATE_MPCAT || RTECH==CHATE_LANG || RTECH==DIPOLE_VCM || RTECH==DIPOLE_DIR_SUM || RTECH==DIPOLE_DIR_AV ) {
		//Allocate memory for S
		S = calloc ( DIM, sizeof( *S ) );
		for( i=0; i<DIM; i++ ) S[i] = calloc ( DIM, sizeof( *S[i] ) );
		for( i=0; i<DIM; i++ ) for( d=0; d<DIM; d++ ) S[i][d] = 0.0;
		// Find the order parameter tensor, the director and the scalar order parameter for each cell
		for( a=0; a<XYZ_P1[0]; a++ ) for( b=0; b<XYZ_P1[1]; b++ ) for( c=0; c<XYZ_P1[2]; c++ ) {
			if( CL[a][b][c].POPSRD > 1 ) {
				// Find the tensor order parameter
				tensOrderParam( &CL[a][b][c],S,LC );				// From the tensor order parameter find eigenvalues and vectors --- S is written over as normalized eigenvectors
				solveEigensystem( S,DIM,eigval );
				//The scalar order parameter is the largest eigenvalue which is given first, ie eigval[0]
				// But can be better approximated (cuz fluctuates about zero) by -2* either of the negative ones (or the average)
				if(DIM==_3D) CL[a][b][c].S = -1.*(eigval[1]+eigval[2]);
				else CL[a][b][c].S=eigval[0];

				if( CL[a][b][c].S<1./(1.-DIM) ){
					#ifdef DBG
					if (DBUG >= DBGRUN){
						printf("Warning: Local scalar order parameter < 1/(1-DIM)\n");
						printf("Cell [%d,%d,%d]\n",a,b,c);
						printf("Eigenvalues=");
						pvec(eigval,DIM);
						printf("Eigenvectors=");
						for( d=0; d<DIM; d++ ) pvec(S[d],DIM);
					}
					#endif
				}

				// The director is the eigenvector corresponding to the largest eigenvalue
				for( i=0; i<DIM; i++ ) CL[a][b][c].DIR[i] = S[0][i];
				if( CL[a][b][c].S>1.0 ) CL[a][b][c].S=1.0;
			}
			//Else just leave as the old value
		}
		for( i=0; i<DIM; i++ ) free( S[i] );
		free( S );
	}
	// Find the velocity gradient tensor
	if( LC!=ISOF ) localVelGrad( CL );
}

/// 
/// @brief This routine calculates the local <b>instantaneous</b> flow in each cell.
///
/// The function simply loops over all cells and calculates the instantaneous center of mass velocity through the routine localMPCVCM(). 
/// @see localMPCVCM()
/// @param CL All of the MPCD cells. 
/// @param SP The species-wide information about MPCD particles.
///
void localFLOW( cell ***CL,spec *SP ) {
	int a,b,c;
	for( a=0; a<XYZ_P1[0]; a++ ) for( b=0; b<XYZ_P1[1]; b++ ) for( c=0; c<XYZ_P1[2]; c++ ) {
		localMPCVCM( CL[a][b][c].VCM,CL[a][b][c],SP ) ;
	}
}

/// 
/// @brief This routine calculates the rolling-average local flow in each cell. 
///
/// The function loops over all cells and adds the current (pre-calculated) local velocity (in `DIM` dimensions) of each cell to a running sum.
/// The sum will become a time-window-averaged flow velocity and will be outputted by flowout().  
/// @param CL All of the MPCD cells.
/// @see flowout()
///
void sumFLOW( cell ***CL ) {
	int a,b,c,d;
	for( a=0; a<XYZ_P1[0]; a++ ) for( b=0; b<XYZ_P1[1]; b++ ) for( c=0; c<XYZ_P1[2]; c++ ) for( d=0; d<DIM; d++ ) {
		CL[a][b][c].FLOW[d] += CL[a][b][c].VCM[d];
	}
}

/// 
/// @brief This routine calculates the rolling-average local flow in each cell, centered around the first swimmer.
///
/// The function loops over all cells and adds the current (pre-calculated) local velocity (in `DIM` dimensions) of each cell to a running sum.
/// The sum will become a time-window-averaged flow velocity and will be outputted by swflowout().  
/// This is done by hijacking the cell structure: the local velocity in question does not correspond to the cell where SWFLOW is stored, but to the cell which holds the same position in a grid centered around the first swimmer.
/// This results in some cells receiving more information than others: this is corrected by SWFLOW[3].
/// @param CL All of the MPCD cells.
/// @param sw The first swimmer's information.
/// @param ss The global information of the swimmers.
/// @see swflowout()
///
void sumSWFLOW( cell ***CL, swimmer *sw , specSwimmer *ss) {
	int a,b,c,d,w;
	double i,j,k,i2,j2,k2;
	int I,J,K;
	double ori[3]={0.0,0.0,0.0},center[3]={0.0,0.0,0.0},finalOrientation[3]={1.0,0.0,0.0};
	double rotMatrix[3][3]={{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};

	// Find the orientation of the first swimmer, and the necessary translation to place its hydrodynamic center at the point (0,0,0)
	for (d=0;d<3; d++) 
	{
		ori[d]=sw[0].H.Q[d]-sw[0].M.Q[d];
		if (ori[d]>0.5*XYZ[d]) ori[d]=ori[d]-XYZ[d];
		else if (ori[d]<-0.5*XYZ[d]) ori[d]=ori[d]+XYZ[d];
		center[d]=0.25*sw[0].M.Q[d]+0.25*sw[0].H.Q[d]-0.5*(sw[0].M.Q[d]-ori[d]*ss->DS);
		center[d]=sw[0].M.Q[d]-ori[d]*ss->DS;
		center[d]=0.5*XYZ[d]-center[d];
	}

	// Find the matrix that places the swimmer's orientation along the x axis
	norm(ori,_3D);
	findRotationMatrix( rotMatrix, ori, finalOrientation);

	// Translates and rewraps, then stores the translated velocity field in the last three slots in SWFLOW (SWFLOW[4],SWFLOW[5],SWFLOW[6]).
	for( a=0; a<XYZ[0]; a++ ) for( b=0; b<XYZ[1]; b++ ) for( c=0; c<XYZ[2]; c++ ) {

		i=a+center[0];j=b+center[1];k=c+center[2];

		if (i>=XYZ[0]) i-=XYZ[0];
		if (i<0) i+=XYZ[0];
		if (j>=XYZ[1]) j-=XYZ[1];
		if (j<0) j+=XYZ[1];
		if (k>=XYZ[2]) k-=XYZ[2];
		if (k<0) k+=XYZ[2];
		I=(int)i;J=(int)j;K=(int)k;

		for (d=0;d<DIM;d++){
			w=d+4;
			CL[I][J][K].SWFLOW[w]=CL[a][b][c].VCM[d];
		}
	}

	// Rotating and rewrapping the cell indices, before rotating the velocity vectors themselves
	for( a=0; a<XYZ[0]; a++ ) for( b=0; b<XYZ[1]; b++ ) for( c=0; c<XYZ[2]; c++ ) {

		for (d=0;d<3; d++) center[d]=XYZ[d]*0.5;

		i2=a-center[0];j2=b-center[1];k2=c-center[2];		
		i=rotMatrix[0][0]*i2+rotMatrix[0][1]*j2+rotMatrix[0][2]*k2+center[0];
		j=rotMatrix[1][0]*i2+rotMatrix[1][1]*j2+rotMatrix[1][2]*k2+center[1];
		k=rotMatrix[2][0]*i2+rotMatrix[2][1]*j2+rotMatrix[2][2]*k2+center[2];

		// Rewrapping
		while (i>=XYZ[0])i=i-XYZ[0];
		while (i<0)i=i+XYZ[0];
		while (j>=XYZ[1])j=j-XYZ[1];
		while (j<0)j=j+XYZ[1];
		while (k>=XYZ[2])k=k-XYZ[2];
		while (k<0)k=k+XYZ[2];

		I=(int)i;J=(int)j;K=(int)k;

		for ( d=0; d<DIM; d++ )
		{
			w=d+4;
			CL[I][J][K].SWFLOW[0] += CL[a][b][c].SWFLOW[w]*rotMatrix[0][d];
			CL[I][J][K].SWFLOW[1] += CL[a][b][c].SWFLOW[w]*rotMatrix[1][d];
			CL[I][J][K].SWFLOW[2] += CL[a][b][c].SWFLOW[w]*rotMatrix[2][d];
		}
		CL[I][J][K].SWFLOW[3] += 1.0; // Counter that keeps track of how many time each grid point corresponds to another grid point.
	}
}

/// 
/// @brief This routine adds the effect of phantom MPCD particles to the centre of mass. 
/// 
/// Then this routine operates (to fill cell up with ghost particles, or apply strong anchoring). 
/// In order to ensure no-slip boundary conditions, ghost particles must be added to all cells that partially overlap with impermeable boundaries. 
/// This is because the MPCD viscosity depends on the particle number density and so cells with partial excluded volume have systematically lower viscosity. 
// To determine if a boundary cuts a cell, a temporary particle assigned to each corner and the routine checks if the corner is inside a boundary. 
/// This is an approximation --- boundaries with sharp corners may be missed. This is done every time step because boundary walls may be mobile. 
/// To ensure no-slip, the centre of mass velocity of the cell is weighted towards zero by the Boltzmann distribution (Gaussian components). 
/// To achieve strong anchoring conditions, the orientation of all the MPCD particles is set equal to the wall normal at that particle's location. 
/// This strengthens anchoring by re-applying orientational boundary conditions (now to the whole cell). 
/// This is so the collision operator can reassign orientations about the director (preferred by the anchoring), with less deviation (as S = 1 for planar, or close to 1 otherwise).
/// @param CL All of the MPCD cells. 
/// @param WALL All of the walls (boundary conditions) that particles might interact with. 
/// @param KBT The thermal energy. 
/// @param LC Flags whether or not the nematic liquid crystal is turned on.
/// @param SP The species-wide information about MPCD particles.
/// @see genrand_gaussMB()
///
void ghostPart( cell ***CL,bc WALL[],double KBT,int LC, spec *SP) {
/*
   This would be oodles better if each BC object had a list of cells
   to worry about and so I wouldn't have to go over every cell all
   the time.
*/
	int a,b,c,d,i,j,k;
	int numCorners = (int) smrtPow(2,DIM);
	double R[DIM];
	double invN;							// Inverse number difference
	particleMPC tp[numCorners];				// Temporary MPC particles for all corners of a square cell
	double W[numCorners];					// W for each of the corners
	double **S, eigval[_3D];
	double n[_3D] = {0.,0.,0.};				// surface normal
	particleMPC *ptMPC;						// temporary pointer to MPC particle
	int setGhostAnch, flagW;
	int numBC;								// Number of walls with anchoring in a given cell
	int wallindex;							// Index of wall with anchoring acting on a cell
	double shift[DIM];
	setGhostAnch = 1; 						// a manual switch to turn on=1 or off=0 the stronger anchoring

	if (setGhostAnch == 1){
		// Allocate memory for S
		S = calloc ( DIM, sizeof( *S ) );
		for( i=0; i<DIM; i++ ) S[i] = calloc ( DIM, sizeof( *S[i] ) );
		for( i=0; i<DIM; i++ ) for( d=0; d<DIM; d++ ) S[i][d] = 0.0;
	}

	// Loop over cells
	for( a=0; a<=XYZ[0]; a++ ) for( b=0; b<=XYZ[1]; b++ ) for( c=0; c<=XYZ[2];c ++ ) {
		// Each cell has a temporary particle assigned to each corner
		// If the boundary cuts a cell (i.e. a corner particle is found inside a boundary)
		// Then this routine operates (to fill cell up with ghost particles, or apply strong anchoring)
		//Position of corners of the cell
		//Origin
		tp[0].Q[0] = (double) a;
		tp[0].Q[1] = (double) b;
		tp[0].Q[2] = (double) c;
		//Shift in x-corner
		tp[1].Q[0] = (double) a+1.0;
		tp[1].Q[1] = (double) b;
		tp[1].Q[2] = (double) c;
		if( DIM >= _2D ) {
			//Shift in y-corners
			tp[2].Q[0] = (double) a;
			tp[2].Q[1] = (double) b+1.0;
			tp[2].Q[2] = (double) c;
			tp[3].Q[0] = (double) a+1.0;
			tp[3].Q[1] = (double) b+1.0;
			tp[3].Q[2] = (double) c;
		}
		if( DIM >= _3D ) {
			//Shift in z-corners
			tp[4].Q[0] = (double) a;
			tp[4].Q[1] = (double) b;
			tp[4].Q[2] = (double) c+1.0;
			tp[5].Q[0] = (double) a+1.0;
			tp[5].Q[1] = (double) b;
			tp[5].Q[2] = (double) c+1.0;
			tp[6].Q[0] = (double) a;
			tp[6].Q[1] = (double) b+1.0;
			tp[6].Q[2] = (double) c+1.0;
			tp[7].Q[0] = (double) a+1.0;
			tp[7].Q[1] = (double) b+1.0;
			tp[7].Q[2] = (double) c+1.0;
		}
		if ( CL[a][b][c].POP > 0 ){

			numBC = 0;
			flagW = 0;
			wallindex = -1; // value not important, just best not between 0 and NBC-1

			// An initial part of the routine to find out how many walls intersect the cell.
			// This is only used for boundaries with anchoring conditions
			// (and strong anchoring turned on) as two conflicting anchoring conditions
			// would take the last one applied (which isn't necessarily what we want).
			// We only want the rest of the ghost particle anchoring routine to continue
			// if 1 anchored wall intersects the cell.
			if ( setGhostAnch == 1 ){
				// Loop over boundary conditions, selecting ones with anchoring
				for( i=0; i<NBC; i++ ) if( WALL[i].PHANTOM && ( feq(WALL[i].MUN,0.0) || feq(WALL[i].MUT,0.0) ) ){

					// Shift moving wall (periodic BC).
					if ( WALL[i].DSPLC ){
						shiftBC( shift, &WALL[i], &tp[0] );
						rotateBC( &WALL[i], &tp[0], LC );
					}

					// Loop over corners. Does the cell get cut by a BC?
					flagW = 0;
					for ( j=0; j<numCorners; j++ ){
						W[j] = calcW( WALL[i], tp[j] );
						if ( W[j] < 0.0 ) { // cell is cut by BC
							flagW = 1;
							wallindex = i; // to identify which wall cut the cell
						}
					}
					if ( flagW == 1 ) numBC += 1;

					// Shift back moving wall
					if ( WALL[i].DSPLC ){
						rotatebackBC( &WALL[i], &tp[0], LC );
						shiftbackBC( shift, &WALL[i] );
					}
				}
			}

			flagW = 0;
			// Boundary loop
			for( i=0; i<NBC; i++ ) if( WALL[i].PHANTOM ){

				// If the colloid is moving through a periodic BC, we want the routine
				// to be applied to its image on the other side of the boundary too.
				// Shift moving wall (periodic BC)
				if ( WALL[i].DSPLC ){
					shiftBC( shift, &WALL[i], &tp[0] );
					rotateBC( &WALL[i], &tp[0], LC );
				}

				// Corner loop. Does the cell get cut by a BC?
				for ( j=0; j<numCorners; j++ ){
					W[j] = calcW( WALL[i], tp[j] );
					if ( W[j] < 0.0 ) flagW = 1; // cell is cut by BC
				}

				// Apply ghost anchoring (see description at top of routine).
				// If turned on, this strengthens the anchoring so the collision operator
				// sees the cell director and scalar order parameter as those preferred
				// by the anchoring conditions.
				if ( LC!=ISOF && setGhostAnch == 1 && flagW && numBC == 1 && wallindex == i){
					if ( feq( WALL[i].MUT, 0.0 ) || feq( WALL[i].MUN, 0.0 ) ){ // if anchored
						// Particle loop
						if (CL[a][b][c].pp!=NULL){
							ptMPC = CL[a][b][c].pp;
							while (ptMPC != NULL) {
								if(WALL[i].INTER[ptMPC->SPID] == BCON) {
									normal(n, WALL[i], ptMPC->Q, DIM);
									norm(n, DIM);
									// apply boundary condition (anchoring) to cell
									// note: this also feeds anchoring torque back to boundary (to add to dV and dL), if mobile
									oriBC(ptMPC, SP, &WALL[i], n);
								}
								ptMPC = ptMPC->next;
							}
						}

						// Re-calculating director and S:
						// Planar wall (manually set S and dir, or eigenvalue calculation fails)
						if ( feq(WALL[i].P[0],1.0) && feq(WALL[i].P[1],1.0) && feq(WALL[i].P[2],1.0) ){
							CL[a][b][c].S = 1.0;
							// Set the director as the orientation of the first particle (as all should be the same)
							for( k=0; k<DIM; k++ ) CL[a][b][c].DIR[k] = CL[a][b][c].pp->U[k];
							norm( CL[a][b][c].DIR, DIM );
						}

						// One particle cell (also manually set)
						else if ( CL[a][b][c].POPSRD == 1 ){
							CL[a][b][c].S = 1.0;
							// Set the director as the orientation of the particle
							for( k=0; k<DIM; k++ ) CL[a][b][c].DIR[k] = CL[a][b][c].pp->U[k];
							norm( CL[a][b][c].DIR, DIM );
						}

						// Curved wall
						else if (CL[a][b][c].POPSRD > 1){
							// Calculate Q tensor
							tensOrderParam( &CL[a][b][c], S, LC );
							// Find S
							solveEigensystem( S,DIM,eigval );
							if(DIM==_3D) CL[a][b][c].S = -1.*(eigval[1]+eigval[2]);
							else CL[a][b][c].S = eigval[0];
							if( CL[a][b][c].S < 1./(1.-DIM) ){
								#ifdef DBG
								if (DBUG >= DBGRUN){
									printf("Warning: Local scalar order parameter < 1/(1-DIM)\n");
									printf("Cell [%d,%d,%d]\n",a,b,c);
									printf("Eigenvalues=");
									pvec(eigval,DIM);
									printf("Eigenvectors=");
									for( d=0; d<DIM; d++ ) pvec(S[d],DIM);
								}
								#endif
							}
							// Cell director (eigenvector corresponding to largest eigenvalue (1st element))
							for( k=0; k<DIM; k++ ) CL[a][b][c].DIR[k] = S[0][k];
							if( CL[a][b][c].S > 1.0 ) CL[a][b][c].S = 1.0;
						}
					}
				} // End ghost anchoring

				// Shift back moving wall
				if ( WALL[i].DSPLC ){
					rotatebackBC( &WALL[i], &tp[0], LC );
					shiftbackBC( shift, &WALL[i] );
				}
			} // End boundary loop

			// Apply ghost particles for centre of mass velocity
			if ( flagW ){
				if( CL[a][b][c].POP < GnDNST && CL[a][b][c].POP > 0 ) {
					invN = 1.0/(GnDNST-(double)CL[a][b][c].POP);
					for( d=0; d<DIM; d++ ) {
						R[d] = genrand_gaussMB( KBT, invN );
						CL[a][b][c].VCM[d] += R[d];
						CL[a][b][c].VCM[d] /= GnDNST;
					}
				}
			}
		}
	} // End cell loop

	// Free memory for S
	if ( setGhostAnch ){
		for( i=0; i<DIM; i++ ) free( S[i] );
		free( S );
	}
}

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* **************** STREAMING *************** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */

/// 
/// @brief This subroutine translates one component of a generic position vector.
/// 
/// This function simply updates one component of a generic position vector based on a constant speed in that direction.
/// @param t The time interval for which the particle translates.
/// @param V The speed with which the particle translates.
/// @param QOLD The initial position of the particle.
/// @return The newly translated component of the position. 
/// @note No acceleration during the time step, which is the philosophy behind the streaming step of MPCD algorithms.
///
double trans( double t,double V, double QOLD ) {
	double QNEW;
	// 	QNEW = QOLD + t * V + 0.5*t*t*G;
	QNEW = QOLD + t * V;
	return QNEW;
}

/// 
/// @brief This subroutine accelerates one component of a generic velocity vector. 
///
/// This function simply updates one component of a generic velocity vector based on a constant acceleration in that direction.
/// @param t The time interval for which the particle accelerates.
/// @param GRAV The acceleration with which the particle speed increases.
/// @param V_OLD The initial speed of the particle.
/// @return The newly accelerated component of the velocity. 
///
double acc( double t,double GRAV,double V_OLD ) {
	double VNEW;
	VNEW = V_OLD + t * GRAV;
	return VNEW;
}

/// 
/// @brief The streaming step of the algorithm translates position of a single wall (boundary condition). 
/// 
/// This function simply updates the position vector (in `DIM` dimensions) of a single wall (boundary condition).
/// @param WALL All of the walls (boundary conditions) that particles might interact with. 
/// @param t The time interval for which the wall translates.
///
void stream_BC( bc *WALL,double t ) {
	int i;
	for( i=0; i<DIM; i++ ) WALL->Q[i] = trans( t,WALL->V[i],WALL->Q[i] );
}

/// 
/// @brief The streaming rotational step of the algorithm rotates the orienation of a single wall (boundary condition). 
/// 
/// This function simply updates the orientation vector of a single wall (boundary condition). 
/// Must be in 3D because the rotation direction in 2D is in the third dimension. 
/// @param WALL All of the walls (boundary conditions) that particles might interact with. 
/// @param t The time interval for which the wall rotates.
///
void spin_BC( bc *WALL,double t ) {
	int i;
	for( i=0; i<_3D	; i++ ) WALL->O[i] += t * WALL->L[i];
}

/// 
/// @brief The streaming step of the algorithm translates the position of a single MPCD particle. 
/// 
/// This function simply updates the position vector (in `DIM` dimensions) of a single MPCD particle.
/// @param t The time interval for which the MPCD particle translates.
/// @param p An MPCD particle. 
///
void stream_P( particleMPC *p,double t ) {
	int i;
	for( i=0; i<DIM; i++ ) p->Q[i] = trans( t,p->V[i],p->Q[i] );
}

/// 
/// @brief Accelerating the velocity of a single wall (boundary condition). 
/// 
/// This function simply updates the velocity vector (in `DIM` dimensions) of a single wall (boundary condition).
/// @param WALL One of the walls (boundary conditions). 
/// @param t The time interval for which the wall (boundary condition) accelerates.
/// @param GRAV The acceleration with which the WALL's speed increases.
///
void acc_BC( bc *WALL,double t,double GRAV[] ) {
	int i;
	for( i=0; i<DIM; i++ ) WALL->V[i] = acc( t,GRAV[i],WALL->V[i] );
}

/// 
/// @brief Accelerating the velocity of a single MPCD particle. 
/// 
/// This function simply updates the velocity vector (in `DIM` dimensions) of a single MPCD particle.
/// @param p An MPCD particle. 
/// @param t The time interval for which the MPCD particle accelerates.
/// @param GRAV The acceleration with which the particle speed increases.
///
void acc_P( particleMPC *p,double t,double GRAV[] ) {

	int i;
	for( i=0; i<DIM; i++ ) p->V[i] = acc( t,GRAV[i],p->V[i] );
}

/// 
/// @brief The streaming step of the algorithm that translates the positions of all MPCD particles.
///
/// This function loops over the global population (`GPOP`) to update all MPCD particle positions.
/// @param pp An MPCD particle. 
/// @param t The time interval for which the MPCD particles translate.
///
void stream_all( particleMPC *pp,double t ) {
	int i;
	for( i=0; i<GPOP; i++ ){
		//Update coordinates --- check if it already streamed
		if( (pp+i)->S_flag ) stream_P( (pp+i),t );
		else (pp+i)->S_flag = STREAM;
	}
}

/// 
/// @brief The accelerating all the MPCD particle velocities.
/// 
/// This function loops over the global population (`GPOP`) to update all MPCD particle velocities.
/// @param pp An MPCD particle. 
/// @param t The time interval for which the MPCD particles accelerate.
/// @param GRAV The acceleration with which the particle speeds increase.
///
void acc_all( particleMPC *pp,double t,double GRAV[] ) {
	int i;
	for( i=0; i<GPOP; i++ ) acc_P( (pp+i),t,GRAV );
}

/// 
/// @brief Shifts the entire system. 
///
/// Since the MPCD algorith happens on a grid, it breaks Galilean invariance. 
/// Therefore, to restore Gallilean invariance, the grid should be randomly shifted. 
/// This was proposed by Ihle and Kroll (https://journals.aps.org/pre/abstract/10.1103/PhysRevE.67.066705). 
/// However, since bin() sorts the particles into cells by truncating the position into cell indices, it is easier to shift <b>everything else</b> and leave the grid in place. 
/// So this routine shifts the entire system by the vector SHIFT. 
/// If `shiftBack` then multiply shift by -1 and shift; otherwise, do the normal shift. 
// It shifts everything by looping over all MPCD particles (global population `GPOP`), all MD particles, all swimmers (`NS`) and all walls (`NBC`).
/// @param SHIFT The random vector everything is shifted by
/// @param shiftBack A flag for whether the routine does the initial shift (==0) or shifts back (==1). 
/// @param SRDparticles All the MPCD particles. 
/// @param WALL All of the walls (boundary conditions) that particles might interact with. 
/// @param simMD A pointer to the entire MD portion of the simulation.
/// @param swimmers All the swimmers, including their head and middle monomers. 
/// @param MD_mode The MD coupling mode. Can be off (noMD), MD particles included in the MPCD collisions (MDinMPC), or MPCD particles included in MD pair interactions (MPCinMD).
///
void gridShift_all( double SHIFT[],int shiftBack,particleMPC *SRDparticles,bc WALL[],simptr simMD,swimmer swimmers[],int MD_mode ) {
	int i,j;							//Counting variables
	double signedSHIFT[_3D];

	if( shiftBack ) for( j=0; j<DIM; j++ ) signedSHIFT[j] = -1.0*SHIFT[j];
	else for( j=0; j<DIM; j++ ) signedSHIFT[j] = SHIFT[j];

	//Shift each particle
	for( i=0; i<GPOP; i++ ) for( j=0; j<DIM; j++ ){
		 SRDparticles[i].Q[j] += signedSHIFT[j]; // perform the shift
		 // The Gallilean shift can (sometimes) force particles outside the boundaries of the system.
		 // These lines will periodically shift those particles to the other side of the domain, while
		 //	keeping those that escaped naturally (ie not due to the shift) in place. It will also return
		 //	them to the correct place when needing to shift back for streaming.
		 if (shiftBack == 0 && SRDparticles[i].Q[j] > XYZ[j] && SRDparticles[i].Q[j] <= XYZ[j] + signedSHIFT[j] )	SRDparticles[i].Q[j] -= XYZ[j];
		 else if (shiftBack == 1 && SRDparticles[i].Q[j] < 0.0 && SRDparticles[i].Q[j] >= signedSHIFT[j] ) SRDparticles[i].Q[j] += XYZ[j];
	 }
	if(MD_mode == MDinMPC) for(i=0; i < (simMD->atom.n); i++ ) {
		(simMD->atom.items+i)->rx += signedSHIFT[0];
		if( DIM>=_2D ) (simMD->atom.items+i)->ry += signedSHIFT[1];
		if( DIM>=_3D ) (simMD->atom.items+i)->rz += signedSHIFT[2];
	}
	#ifdef DBG
		if( DBUG == DBGSWIMMERDEETS && !shiftBack ) {
			printf("\tPre-shift:\n");
			for( i=0; i<NS; i++ ) {
				printf( "\tS%d:\n",i );
				swcoord(swimmers[i]);
			}
		}
	#endif
	for( i=0; i<NS; i++ ) for( j=0; j<DIM; j++ ) {
		swimmers[i].H.Q[j] += signedSHIFT[j];
		swimmers[i].M.Q[j] += signedSHIFT[j];
	}
	#ifdef DBG
		if( DBUG == DBGSWIMMERDEETS && shiftBack ) {
			printf("\tPost-shift:\n");
			for( i=0; i<NS; i++ ) {
				printf( "\tS%d:\n",i );
				swcoord(swimmers[i]);
			}
		}
	#endif
	//Shift each boundary
	for( i=0; i<NBC; i++) for( j=0; j<DIM; j++ ) WALL[i].Q[j] += signedSHIFT[j];
}

/// 
/// @brief This routine applies a solid body rotation to the particles in a cell. 
///
/// Applies a solid body rotation to the particles in a given cell. 
/// This applies a change in angular speed dw to all the MPCD particles in a cell. 
/// The change is about a given point r0 (likely the CM but allowed to be anything) and direction n0. 
/// It updates the particles by looping through the linked list. 
/// @param CL An MPCD cell. 
/// @param SP The species-wide information about MPCD particles.
/// @param r0 The point about which the rotation occurs.
/// @param n0 The axis about which the rotation occurs. 
/// @param dw Change in angular speed. 
///
void rotate_CL( cell CL,spec *SP,double r0[],double n0[],double dw ) {
	int d=0;
	double r[DIM],r_perp[DIM],n_perp[DIM]; 				//Vectors betweem cm and mpcd particles and axis
	double n_v[DIM];									//Direction of velocity change
	double dist,dv;										//Mag of angular velocity on each MPCD particle, distance to line, mag of change in velocity
	particleMPC *pMPC;									//Temporary pointer to MPC particles
	double netMom[DIM],tempAngMom[DIM],netAngMom[DIM];	//Sum momentum and angular momentum
	int cnt;											//Count number of SRD particles

	for( d=0; d<DIM; d++ ) {
		r[d] = 0.0;
		r_perp[d] = 0.0;
		n_perp[d] = 0.0;
		n_v[d] = 0.0;
		netMom[d] = 0.0;
		netAngMom[d] = 0.0;
	}
	cnt = 0;
	if( CL.pp != NULL ) {
		pMPC = CL.pp;
		while(pMPC != NULL) {
			//Vector from head to SRD particle
			for( d=0; d<DIM; d++ ) r[d] = pMPC->Q[d] - r0[d];
			//Projection of r on direction of swimmer
			dist = dotprod( r,n0,DIM );
			//Perpendicular director from the particle to the line of the swimmer's orientation
			for( d=0; d<DIM; d++ ) r_perp[d] = r[d] - dist*n0[d];
			normCopy( r_perp,n_perp,DIM );
			dist = sqrt(dotprod( r_perp,r_perp,DIM ));
			//Direction of angular impulse
			crossprod( n0,n_perp,n_v );
			//Magnitude of the impulse
			dv = dist*dw;
			//Update the velocity
			for( d=0; d<DIM; d++ ) pMPC->V[d] += n_v[d]*dv;
			//Add to net momentum
			for( d=0; d<DIM; d++ ) netMom[d] += n_v[d]*dv;
			//Add to net angular momentum
			crossprod( n_perp,n_v,tempAngMom );
			for( d=0; d<DIM; d++ ) tempAngMom[d] *= dv*dist*SP[pMPC->SPID].MASS;
			for( d=0; d<DIM; d++ ) netAngMom[d] += tempAngMom[d];
			//Increment link in list
			pMPC = pMPC->next;
			cnt++;
		}
	}
	#ifdef DBG
		if( DBUG == DBGSWIMMERTORQUE ) {
			printf( "\t\tNet momentum on cell");
			pvec(netMom,DIM);
			printf( "\t\tMom magnitude %lf\n",length(netMom,DIM) );
			printf( "\t\tMom direction" );
			norm( netMom,DIM );
			pvec(netMom,DIM);
			printf( "\t\tNet ang mom");
			pvec(netAngMom,DIM);
			printf( "\t\tAng mom magnitude %lf\n",length(netAngMom,DIM) );
			printf( "\t\tAng mom direction" );
			norm( netAngMom,DIM );
			pvec(netAngMom,DIM);
			printf( "\t\tNumber of SRD cells %d\n",cnt );
		}
	#endif
}

/// 
/// @brief Rewinds a translation.
/// 
/// This function restores the previous component of a generic position vector based on a constant speed in that direction by subtracting the displacement. 
/// That is to say, it rewinds the particle position to the previous time step. 
/// @param t The time interval for which the object is rewound.
/// @param V The object's velocity component.
/// @param P The object's present position component.
/// @return The position at the previous time step.
/// @note No acceleration during the time step, which is the philosophy behind the streaming step of MPCD algorithms.
///
double rewind_trans( double t,double V,double P ) {
	double QOLD;
	QOLD = P - t*V;
	return QOLD;
}

/// 
/// @brief Rewind an acceleration vector.
/// 
/// This function restores the previous component of a generic velocity vector based on a constant accelearation in that direction by subtracting the displacement. 
/// That is to say, it rewinds the particle velocity to the previous time step. 
/// @param t The time interval for which the object is rewound.
/// @param G The object's acceleration component.
/// @param V The object's present velocity component.
/// @return The velocity at the previous time step.
///
double rewind_acc(double t,double G,double V){
	double V_OLD;
	V_OLD = V - t*G;
	return V_OLD;
}

/// 
/// @brief Bring a given MPCD particle back a time step.
/// 
/// This function restores the previous component of the MPCD position vector based on a constant velocity by subtracting the displacement. 
/// That is to say, it rewinds the particle position to the previous time step. 
/// @param p An MPCD particle. 
/// @param time The time interval for which the MPCD particle is rewound.
///
void rewind_P( particleMPC *p,double time ) {
	int i;
	for( i=0; i<DIM; i++ ) p->Q[i] = rewind_trans(time,p->V[i],p->Q[i]);
}

/// 
/// @brief Bring a given boundary back a time step.
/// 
/// This function restores the previous component of the wall boundary position vector based on a constant velocity by subtracting the displacement. 
/// That is to say, it rewinds the wall position to the previous time step. 
/// @param WALL A moving wall (boundary conditions). 
/// @param time The time interval for which the MPCD particle is rewound.
///
void rewind_BC( bc *WALL,double time ) {
	int i;
	for( i=0; i<DIM; i++ ) WALL->Q[i] = rewind_trans(time,WALL->V[i],WALL->Q[i]);
}

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************ BINNING IN CELLS ************ */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */

/// 
/// @brief This function bins the MPCD particles for the first time for use in the collision steps.
/// 
/// Initial binning loops over all MPCD particles using the global population (`GPOP`) after they have been first initialized. 
/// Particles are placed in cells by truncating their positions into integers, which give the cell array place.
/// This works because the cells <b>must</b> have an cell size of `a=1`. 
/// It is different from bin() in that it uses the actual array of MPCD particles rather than the array of pointers to particles.
/// @param p All of the MPCD particles. 
/// @param CL All of the MPCD cells (including the linked list of particles in each cell). 
/// @see bin()
///
void binin( particleMPC p[],cell ***CL ) {
	int i,a,b,c;
	//Bin Particles
	for( i=0; i<GPOP; i++ ){
		//Truncate the coordinate to see what cell the particleMPC falls in
		a = (int)p[i].Q[0];
		b = (int)p[i].Q[1];
		c = (int)p[i].Q[2];
		addlink( &CL[a][b][c],&p[i] );
	}
}

/// 
/// @brief This function bins the MPCD particles.
/// 
/// This function bins the particleMPCs by placing a pointer to the MPCD particle in the appropriate new list and removing it from it's old list.
/// This does not calculate the local properties, which must be done separately by calling localPROP().
/// @param CL All of the MPCD cells (including the linked list of particles in each cell). 
/// @param SP The species-wide information about MPCD particles.
/// @param WALL All of the walls (boundary conditions) that particles might interact with. 
/// @param KBT The thermal energy. 
/// @param LC Flags whether or not the nematic liquid crystal is turned on.
/// @param shifted Flags if the positions have been randomly shifted (shifted==1), in which case do the PBC wrap.
///
void bin( cell ***CL,spec *SP,bc WALL[],double KBT,int LC,int shifted ) {
	int i,j,k,a,b,c,d;
	double m;
	particleMPC *cp;	//Pointer to current item in list
	particleMPC *tp;	//Temporary pointer
	//Search each cell for particleMPCs that have left the cell
	for( i=0; i<XYZ_P1[0]; i++ ) for( j=0; j<XYZ_P1[1]; j++ ) for( k=0; k<XYZ_P1[2]; k++ ) {
		cp = CL[i][j][k].pp;
		while(cp!=NULL) {
			//Truncate the coordinate to see what cell the particleMPC falls in
			a = (int)cp->Q[0];
			b = (int)cp->Q[1];
			c = (int)cp->Q[2];
			//Save the next address
			tp = cp->next;
			//Apply the periodic BC in binning if the axes have been flagged
			if( shifted ) {
				if( XYZPBC[0] && a==XYZ[0]) a=0;
				if( XYZPBC[1] && b==XYZ[1]) b=0;
				if( XYZPBC[2] && c==XYZ[2]) c=0;
			}
			//Make sure particleMPC didn't escape
			if( a>XYZ[0] || b>XYZ[1] || c>XYZ[2] || a<0 || b<0 || c<0 ){
				#ifdef DBG
					printf( " Warning:\tMPC particle escaped from cell [%i,%i,%i].\n",a,b,c );
					printf( " \t\tSystem size [%i,%i,%i]\n",XYZ[0],XYZ[1],XYZ[2] );
					printf( " \t\tReplacing\n" );
					pcoord( *cp );
				#endif
				replacePos_WithCheck( cp,WALL );
				d = cp->SPID;
				m = (SP+d)->MASS;
				for( d=0; d<DIM; d++ ) cp->V[d] = genrand_gaussMB( KBT,m );
				if( LC>ISOF ) genrand_coneNP( cp->U,pi,DIM );
			}
			//Remove from old list and add to new
			else if( a != i || b != j || c != k ){
				removelink( cp,&CL[i][j][k] );
				addlink( &CL[a][b][c],cp );
			}
			//Return attention to the list under consideration
			cp = tp;
		}
	}
}

/// 
/// @brief This function bins the MD particles for the first time for use in the collision steps.
/// 
/// Initial binning loops over all MD particles after they have been first initialized. Particles are placed in MPCD cells by truncating their positions into integers, which give the cell array place. This works because the MPCD cells <b>must</b> have an cell size of `a=1`. 
/// It is different from binin() only in that it initially bins MD particles.
/// @param sim A pointer to the entire MD portion of the simulation.
/// @param CL All of the MPCD cells (including the linked list of particles in each cell). 
/// @see binin()
///
void bininMD( simptr sim,cell ***CL ) {
	int i,a,b,c;
	int nAtom;				// Number of atoms
	particleMD *atom, *p;

	nAtom = sim->atom.n;
	atom = sim->atom.items;
	//Bin Particles
	for( i=0; i<nAtom; i++ ) {
		p = atom+i;
		//Truncate the coordinate to see what cell the particleMD falls in
		a = (int)p->rx;
		b = (int)p->ry;
		c = (int)p->rz;
		addlinkMD( &CL[a][b][c],p );
	}
}

/// 
/// @brief This function bins the MD particles.
///
/// This function bins the MD particles by placing a pointer to the MD particle in the appropriate new list and removing it from it's old list.
/// This does not calculate the local properties, which must be done separately by calling localPROP().
/// @param CL All of the MPCD cells (including the linked list of particles in each cell). 
/// @see bin()
///
void binMD( cell ***CL ) {
	int i,j,k,a,b,c;
	particleMD *cp;		//Pointer to current item in list
	particleMD *tp;		//Temporary pointer
	//Search each cell for particleMPCs that have left the cell
	for( i=0; i<XYZ_P1[0]; i++ ) for( j=0; j<XYZ_P1[1]; j++ ) for( k=0; k<XYZ_P1[2]; k++ ) {
		cp = CL[i][j][k].MDpp;
		while( cp!=NULL ) {
			//Truncate the coordinate to see what cell the particleMPC falls in
			a = (int)cp->rx;
			b = (int)cp->ry;
			c = (int)cp->rz;

			//Save the next address
			tp = cp->nextSRD;
			//Make sure particleMPC didn't escape
			if( a>XYZ[0] || b>XYZ[1] || c>XYZ[2] || a<0 || b<0 || c<0 ){
				printf( "Error:\tMD particle escaped from cell [%i,%i,%i].\n",a,b,c );
				mdcoord( *cp );
 				exit( 0 );
			}
			//Remove from old list and add to new
			else if( a != i || b != j || c != k ){
				removelinkMD(cp,&CL[i][j][k]);
				addlinkMD(&CL[a][b][c],cp);
			}
			//Return attention to the list under consideration
			cp = tp;
		}
	}
}

/// 
/// @brief This routine adds a link to an MPCD particle to the end of the list. 
/// 
/// This function finds the end of a linked list and adds a new link to a given MPCD particle `p`.
/// @param CL All of the MPCD cells (including the linked list of particles in each cell). 
/// @param p The current MPCD particle in the linked list. 
/// @see bin()
///
void addlink( cell *CL,particleMPC *p ) {
	particleMPC *tp;			//Temporary pointer to particleMPC

	if( CL->pp == NULL ) {
		CL->pp = p;
		p->previous = NULL;
	}
	else {
		tp = CL->pp;
		//Find the end of the list
		while( tp->next!=NULL ) tp = tp->next;
		//Once the end is found, add the particleMPC to the list
		tp->next = p;
		//Point the particleMPC back at the last link
		p->previous = tp;
	}
	//Particle is now at the end of the list so set next to null
	p->next = NULL;
}

/// 
/// @brief This routine removes a link to an MPCD particle from a list and relinks the list.
/// 
/// This function removes the `current` MPCD particle from a linked list and re-stiches the list back together.
/// @param current The current MPCD particle being removed from the list. 
/// @param CL The MPCD cell that the particle is in. 
///
void removelink( particleMPC *current,cell *CL ) {
	//Point the next link back at the previous link (unless last link)
	if( current->next!=NULL ) current->next->previous = current->previous;
	//Point the previous link at the next link (unless first link)
	if( current->previous!=NULL ) current->previous->next = current->next;
	//If the current link is the 1st link link the cell to the new first
	else {
		CL->pp = current->next;
	}
	//Remove the current (this is kinda redundant, no?)
	current->previous = NULL;
	current->next = NULL;
}

/// 
/// @brief This routine adds a link to the end of the list for MD particles.
///
/// This function finds the end of a linked list and adds a new link to a given MD particle `p`.
/// @param CL The MPCD cell that the particle is being added to. 
/// @param p The MPCD particle being added to the list. 
///
void addlinkMD( cell *CL,particleMD *p ) {
	particleMD *tp;					//Temporary pointer to particleMPC

	if( CL->MDpp == NULL ) {
		CL->MDpp = p;
		p->prevSRD = NULL;
	}
	else {
		tp = CL->MDpp;
		//Find the end of the list
		while( tp->nextSRD!=NULL ) tp = tp->nextSRD;
		//Once the end is found, add the particleMD to the list
		tp->nextSRD = p;
		//Point the particleMD back at the last link
		p->prevSRD = tp;
	}
	//Particle is now at the end of the list so set next to null
	p->nextSRD = NULL;
}

/// 
/// @brief This routine removes a link pointing to an MD particle from a list and relinks the list. 
///
/// This function removes the `current` MD particle from a linked list and re-stiches the list back together.
/// @param current The current MD particle being removed from the list. 
/// @param CL The MPCD cell that the particle is in. 
///
void removelinkMD( particleMD *current,cell *CL ) {
	//Point the next link back at the previous link (unless last link)
	if( current->nextSRD!=NULL ) current->nextSRD->prevSRD = current->prevSRD;
	//Point the previous link at the next link (unless first link)
	if( current->prevSRD!=NULL ) current->prevSRD->nextSRD = current->nextSRD;
	//If the current link is the 1st link link the cell to the new first
	else {
		CL->MDpp = current->nextSRD;
	}
	//Remove the current (this is kinda redundant, no?)
	current->prevSRD = NULL;
	current->nextSRD = NULL;
}

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* *************** COLLISIONS *************** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */

/// 
/// @brief Rotates a vector about a given axis. 
///
/// Takes in a vector and rotates it about a given axis. 
/// The rotation axis is either random cartesian axis (if `RT==ORTHAXIS`) or about a random axis L (if `RT==ARBAXIS`).
/// The arbitary axis is chosen randomly for every collision cell and each time step
/// @param RT Rotation technique (collision operator). 
/// For SRD, whether to randomly choose a cartesian axis (`ORTHAXIS`) or generate an arbitrarily random direction (`ARBAXIS`).
/// @param Cos Cosine of the rotation angle.
/// @param Sin Sine of the rotation angle.
/// @param V The vector to be rotated. 
/// @param L The random axis to rotate `V` about.  Used if `RT==ARBAXIS`
/// @param SIGN The sign of the rotation angle `RA` to do the rotation about. 
/// @param RAND Randomly pick one of the `DIM` cartesian axes. Used if `RT==ORTHAXIS`
///
void rotate( int RT,double Cos,double Sin,double V[_3D],double L[_3D],long SIGN,int RAND ) {
	int i,j;
	double R[_3D][_3D][_3D];	//[matrix][row][col]
	double TEMP[_3D]={0.0};
	double C,S;

	C = Cos;
	if( SIGN > 0.5 ) S = Sin;
	else S = -1.*Sin;

	if( DIM==_3D ) {
		if( RT == ARBAXIS ) {
			//We only use one matrix with this technique so we set only one of the three 3X3 matrices in R
			RAND = 0;
			//Set rotation matrix
			R[RAND][0][0] = L[0]*L[0]+(1.-L[0]*L[0])*C;
			R[RAND][0][1] = L[0]*L[1]*(1. - C)-L[2]*S;
			R[RAND][0][2] = L[0]*L[2]*(1. - C)+L[1]*S;
			R[RAND][1][0] = L[0]*L[1]*(1. - C)+L[2]*S;
			R[RAND][1][1] = L[1]*L[1]+(1.-L[1]*L[1])*C;
			R[RAND][1][2] = L[1]*L[2]*(1. - C)-L[0]*S;
			R[RAND][2][0] = L[0]*L[2]*(1. - C)-L[1]*S;
			R[RAND][2][1] = L[1]*L[2]*(1. - C)+L[0]*S;
			R[RAND][2][2] = L[2]*L[2]+(1.-L[2]*L[2])*C;
		}
		else if( RT == ORTHAXIS ) {
			//Rotation about the X-axis
			R[0][0][0] = 1.;
			R[0][0][1] = 0.;
			R[0][0][2] = 0.;
			R[0][1][0] = 0.;
			R[0][1][1] = C;
			R[0][1][2] = S;
			R[0][2][0] = 0.;
			R[0][2][1] = -.1*S;
			R[0][2][2] = C;
			//Rotation about the Y-axis
			R[1][0][0] = C;
			R[1][0][1] = 0.;
			R[1][0][2] = S;
			R[1][1][0] = 0.;
			R[1][1][1] = 1.;
			R[1][1][2] = 0.;
			R[1][2][0] = -1.*S;
			R[1][2][1] = 0.;
			R[1][2][2] = C;
			//Rotation about the Z-axis
			R[2][0][0] = C;
			R[2][0][1] = S;
			R[2][0][2] = 0.;
			R[2][1][0] = -1.*S;
			R[2][1][1] = C;
			R[2][1][2] = 0.;
			R[2][2][0] = 0.;
			R[2][2][1] = 0.;
			R[2][2][2] = 1.;
		}
		else {
			printf( "Error: Rotation technique unacceptable.\nNote: %d D system.\n",DIM );
			exit( 1 );
		}
	}
	else if( DIM == _2D ) {
		R[RAND][0][0] = C;
		R[RAND][0][1] = -1.*S;
		R[RAND][1][0] = S;
		R[RAND][1][1] = C;
	}
	else {
		printf( "Error: Rotation technique unacceptable.\nNote: %d D system.\n",DIM );
		exit( 1 );
	}
	//Multiply the matrix to the vector --- TEMP was inititated to zero
	for( i=0; i<DIM; i++ ) for( j=0; j<DIM; j++ ) TEMP[i] += R[RAND][i][j] * V[j];
	//Save the newly rotated value
	for( i=0; i<DIM; i++ ) V[i] = TEMP[i];
}

/// 
/// @brief Does the stochastic rotation dynamics collision. 
/// 
/// Invented by Malevanets and Kapral (https://doi.org/10.1063/1.478857). 
/// Stochastic rotation dynamics (SRD) version of multi-particle collision dynamics (MPCD). 
/// In SRD, the collision operator rotates the particles through a given rotation angle about a randomly chosen axis.
/// This exchanges momenta between particles but energy and momentum are conserved within the cell. 
/// Does not conserve angular momentum. 
/// It updates the particles by looping through the linked lists. 
/// @param CL An MPCD cell (including the linked list of particles in each cell). 
/// @param RTECH The MPCD collision operator. See `definitions.h` for all options.
/// @param C Cosine of the rotation angle.
/// @param S Sine of the rotation angle.
/// @param CLQ The geometric centre of `CL`, the MPCD cell.
/// @param outP Flag whether or not to output the pressure.
///
void stochrotMPC( cell *CL,int RTECH,double C,double S,double *CLQ,int outP ) {
	int i,j;
	long SIGN;								//SIGN indicates sign of RA
	int CA = 0;								//indicates random cartesian axis chosen
	double RV[_3D] = {0.};					//Random vector
	double V[_3D];
	double dp[CL->POP][DIM],relQ[CL->POP][DIM];	//For pressure
	particleMPC *tmpc;						//Temporary particleMPC
	particleMD *tmd;						//Temporary particleMD
	smono *tsm;								//Temporary swimmer monomer

	//Generate random axis(used if RTECH=ARBAXIS)
	if( RTECH == ARBAXIS ) ranvec( RV,DIM );
	//Randomly pick an axis(used if RTECH=ORTHAXIS)
	if( RTECH == ORTHAXIS ) CA = (int)DIM * genrand_real();
	//Randomly pick a sign for the angle
	SIGN = genrand_real();

	/* ****************************************** */
	/* *************** Collision **************** */
	/* ****************************************** */
	//MPC particles
	i=0;
	tmpc = CL->pp;
	while( tmpc!=NULL ) {
		//Pressure term
		if( outP ) calcPressureColl_preColl( relQ[i],dp[i],tmpc,CLQ );
		for( j=0; j<DIM; j++ ) V[j] = tmpc->V[j] - CL->VCM[j];
		rotate( RTECH,C,S,V,RV,SIGN,CA );
		for( j=0; j<DIM; j++ ) tmpc->V[j] = CL->VCM[j] + V[j];
		//Pressure term
		if( outP ) calcPressureColl_postColl( relQ[i],dp[i],1.0,tmpc->V,CL );
		//Increment link in list
		tmpc = tmpc->next;
		i++;
	}
	//MD particles
	tmd = CL->MDpp;
	while( tmd!=NULL ) {
		V[0] = tmd->vx - CL->VCM[0];
		V[1] = tmd->vy - CL->VCM[1];
		V[2] = tmd->vz - CL->VCM[2];
		rotate( RTECH,C,S,V,RV,SIGN,CA );
		tmd->vx = CL->VCM[0] + V[0];
		tmd->vy = CL->VCM[1] + V[1];
		tmd->vz = CL->VCM[2] + V[2];
		//Increment link in list
		tmd = tmd->nextSRD;
	}
	//Swimmer particles
	tsm = CL->sp;
	while( tsm!=NULL ) {
		for( i=0; i<DIM; i++ ) V[i] = tsm->V[i] - CL->VCM[i];
		rotate( RTECH,C,S,V,RV,SIGN,CA );
		for( i=0; i<DIM; i++ ) tsm->V[i] = CL->VCM[i] + V[i];
		//Increment link in list
		tsm = tsm->next;
	}
}

/// 
/// @brief Does the Andersen-thermostatted collision. 
/// 
/// Invented by Noguchi, Kikuchi and Gompper (https://iopscience.iop.org/article/10.1209/0295-5075/78/10005).
/// Andersen-thermostatted version of multiparticle collision dynamics (MPCD-AT). 
/// MPCD-AT is one particular example of a multi-particle collision dynamics (MPCD). 
/// In MPCD-AT, the collision operator generates random velocity vectors for all the MPCD particles in the cell.
/// The random velocitities are drawn from Gaussian distribtions (genrand_gaussMB()) and so obey Maxwell-Boltzmann statistics.
/// The average of the random velocities is then subtracted from all to conserve momentum.
/// Energy is not conserved; rather the system is thermostatted to `KBT`. 
/// Also does not conserve angular momentum. 
/// @param CL An MPCD cell (including the linked list of particles in each cell). 
/// @param SP The species-wide information about MPCD particles.
/// @param SS The species-wide information about swimmers.
/// @param KBT The thermal energy. 
/// @param CLQ The geometric centre of `CL`, the MPCD cell.
/// @param outP Flag whether or not to output the pressure.
///
void andersenMPC( cell *CL,spec *SP,specSwimmer SS,double KBT,double *CLQ,int outP ) {
	int i,j,id;
	double MASS;
	double RV[CL->POP][DIM];				//Random velocities
	double RS[DIM];							//Sum of random velocities
	double DV[CL->POP][DIM];				//Damping velocity
	double dp[CL->POP][DIM],relQ[CL->POP][DIM];	//For pressure
	particleMPC *tmpc;						//Temporary particleMPC
	particleMD *tmd;						//Temporary particleMD
	smono *tsm;								//Temporary swimmer monomer

	// Zero arrays
	for( i=0; i<DIM; i++ ) RS[i] = 0.;
	for( i=0;i<CL->POP;i++ ) for( j=0;j<DIM;j++ ) {
		RV[i][j] = 0.0;
		DV[i][j] = 0.0;
	}
	
	/* ****************************************** */
	/* ******* Generate random velocities ******* */
	/* ****************************************** */
	i=0;
	// MPC particles
	tmpc = CL->pp;
	while( tmpc!=NULL ) {
		id = tmpc->SPID;
		MASS = (SP+id)->MASS;
		//Pressure term
		if( outP ) calcPressureColl_preColl( relQ[i],dp[i],tmpc,CLQ );
		for( j=0; j<DIM; j++ ) RV[i][j] = genrand_gaussMB( KBT,MASS );
		for( j=0; j<DIM; j++ ) RS[j] += MASS*RV[i][j];
		for( j=0; j<DIM; j++ ) DV[i][j] = ((SP+id)->DAMP)*(CL->VCM[j])/((double)CL->POP);
		tmpc = tmpc->next;
		i++;
	}
	// MD particles
	tmd = CL->MDpp;
	while( tmd!=NULL ) {
		MASS = tmd->mass;
		for( j=0; j<DIM; j++ ) RV[i][j] = genrand_gaussMB( KBT,MASS );
		for( j=0; j<DIM; j++ ) RS[j] += MASS*RV[i][j];
		tmd = tmd->nextSRD;
		i++;
	}
	// Swimmer monomers
	tsm = CL->sp;
	while( tsm!=NULL ) {
		if( tsm->HorM ) MASS = (double) SS.middM;
		else MASS = (double) SS.headM;
		for( j=0; j<DIM; j++ ) RV[i][j] = genrand_gaussMB( KBT,MASS );
		for( j=0; j<DIM; j++ ) RS[j] += MASS*RV[i][j];
		tsm = tsm->next;
		i++;
	}
	// Normalize
	for( j=0; j<DIM; j++ ) RS[j] /= CL->MASS;

	/* ****************************************** */
	/* *************** Collision **************** */
	/* ****************************************** */
	i=0;
	// MPC particles
	tmpc = CL->pp;
	while( tmpc!=NULL ) {
		id = tmpc->SPID;
		MASS = (SP+id)->MASS;
		for( j=0; j<DIM; j++ ) tmpc->V[j] = CL->VCM[j] + RV[i][j] - RS[j] - DV[i][j];
		//Pressure term
		if( outP ) calcPressureColl_postColl( relQ[i],dp[i],MASS,tmpc->V,CL );
		//Increment link in list
		tmpc = tmpc->next;
		i++;
	}
	//MD particles
	tmd = CL->MDpp;
	while( tmd!=NULL ) {
		tmd->vx = CL->VCM[0] + RV[i][0] - RS[0];
		if( DIM > _1D){
			tmd->vy = CL->VCM[1] + RV[i][1] - RS[1];
		}
		if(DIM > _2D){
			tmd->vz = CL->VCM[2] + RV[i][2] - RS[2];
		}
		//Increment link in list
		tmd = tmd->nextSRD;
		i++;
	}
	// Swimmer monomers
	tsm = CL->sp;
	while( tsm!=NULL ) {
		for( j=0; j<DIM; j++ ) tsm->V[j] = CL->VCM[j] + RV[i][j] - RS[j];
		//Increment link in list
		tsm = tsm->next;
		i++;
	}
}

/// 
/// @brief Does the Andersen thermostatted collision that conserves angular momentum. 
/// 
/// Invented by Noguchi, Kikuchi and Gompper (https://iopscience.iop.org/article/10.1209/0295-5075/78/10005).
/// Andersen-thermostatted version of multiparticle collision dynamics (MPCD-AT) that <b>does</b> conserve angular momentum (MPCD-AT+a). 
/// It is just like andersenMPC() but includes a correction to conserve angular momentum.
/// Energy is not conserved; rather the system is thermostatted to `KBT`.
/// @param CL An MPCD cell (including the linked list of particles in each cell). 
/// @param SP The species-wide information about MPCD particles.
/// @param SS The species-wide information about swimmers.
/// @param KBT The thermal energy. 
/// @param CLQ The geometric centre of `CL`, the MPCD cell.
/// @param outP Flag whether or not to output the pressure.
/// @see andersenMPC()
///
void andersenROT( cell *CL,spec *SP,specSwimmer SS,double KBT,double *CLQ,int outP ) {
	int i,j,id;
	double MASS;
	double RV[CL->POP][_3D];					//Random velocities
	double RS[_3D];								//Sum of random velocities
	double DV[CL->POP][DIM];					//Damping velocity
	double relQ[CL->POP][_3D];					//Relative position
	double diffV[_3D];							//Difference in velocity
	double L[_3D];								//Angular momentum
	double angterm[_3D];
	double angmom[_3D];
	double W[_3D];
	double VCM[_3D];
	double II[_3D][_3D];						//Inverse of moment of inertia tensor (3D)
	double dp[CL->POP][DIM],relQP[CL->POP][DIM];	//For pressure
	particleMPC *tmpc;							//Temporary particleMPC
	particleMD *tmd;							//Temporary particleMD
	smono *tsm;									//Temporary swimmer monomer

	for( i=0;i<CL->POP;i++ ) for( j=0;j<_3D;j++ ) {
		RV[i][j] = 0.;
		DV[i][j] = 0.;
		relQ[i][j] = 0.;
	}
	for( j=0;j<_3D;j++ ) {
		RS[j]=0.;
		diffV[j]=0.;
		L[j]=0.;
		angterm[j]=0.;
		angmom[j]=0.;
		W[j]=0.;
		VCM[j] = CL->VCM[j];
		for( i=0;i<_3D;i++ ) II[j][i]=0.;
	}
	/* ****************************************** */
	/* ******* Generate random velocities ******* */
	/* ****************************************** */
	i=0;
	//MPC particles
	tmpc = CL->pp;
	while( tmpc!=NULL ) {
		id = tmpc->SPID;
		MASS = (SP+id)->MASS;
		//Pressure term
		if( outP ) calcPressureColl_preColl( relQP[i],dp[i],tmpc,CLQ );
		for( j=0; j<DIM; j++ ) RV[i][j] = genrand_gaussMB( KBT,MASS );
		for( j=0; j<DIM; j++ ) RS[j] += MASS*RV[i][j];
		for( j=0; j<DIM; j++ ) DV[i][j] = ((SP+id)->DAMP)*(CL->VCM[j]);
		tmpc = tmpc->next;
		i++;
	}
	//MD particles
	tmd = CL->MDpp;
	while( tmd!=NULL ) {
		MASS = tmd->mass;
		for( j=0; j<DIM; j++ ) RV[i][j] = genrand_gaussMB( KBT,MASS );
		for( j=0; j<DIM; j++ ) RS[j] += MASS*RV[i][j];
		tmd = tmd->nextSRD;
		i++;
	}
	//Swimmer particles
	tsm = CL->sp;
	while( tsm!=NULL ) {
		if( tsm->HorM ) MASS = (double) SS.middM;
		else MASS = (double) SS.headM;
		for( j=0; j<DIM; j++ ) RV[i][j] = genrand_gaussMB( KBT,MASS );
		for( j=0; j<DIM; j++ ) RS[j] += MASS*RV[i][j];
		tsm = tsm->next;
		i++;
	}
	for( j=0; j<DIM; j++ ) RS[j] /= CL->MASS;
	/* ****************************************** */
	/* ******** Invert moment of inertia ******** */
	/* ****************************************** */
	for( i=0; i<_3D; i++ ) for( j=0; j<_3D; j++ ) II[i][j] = 0.;
	if( DIM == _3D ) {
		// When the population is 2 (don't even get this far if pop=1) then the moment of inertia tensor is singular
		// In that case, can't invert. Therefore, we must neglect angular momentum conservation in such cell
		// This just makes angterm = [0,0,0]
		if( CL->POP < 3 ) for( i=0; i<_3D; i++ ) for( j=0; j<_3D; j++ ) II[i][j] = 0.;
		else invert3x3(II,CL->I);
	}
	else if( DIM == _2D ) {
		II[2][2] = 1./CL->I[2][2];
		// The rest don't matter and can remain zero
	}
	else {
		printf( "Error: Angular momentum conservation in 1D is nonsequitur. Higher dimensions not done. Change collision technique or dimension.\n" );
		exit( 1 );
	}
	/* ****************************************** */
	/* ******* Find angular velocity term ******* */
	/* ****************************************** */
	for( i=0; i<_3D; i++ ) L[i] = 0.;
	i=0;
	//MPC particles
	tmpc = CL->pp;
	while( tmpc!=NULL ) {
		id = tmpc->SPID;
		MASS = (SP+id)->MASS;
		//Position relative to centre of mass
		for( j=0; j<DIM; j++ ) relQ[i][j] = tmpc->Q[j] - CL->CM[j];
		//Difference in velocity
		for( j=0; j<DIM; j++ ) diffV[j] = MASS * (tmpc->V[j] - RV[i][j]);
		if( DIM < _3D ) {
			relQ[i][2] = 0.;
			diffV[2] = 0.;
		}
		if( DIM < _2D ) {
			relQ[i][1] = 0.;
			diffV[1] = 0.;
		}
		crossprod( relQ[i],diffV,angmom );
		for( j=0; j<_3D; j++ ) L[j] += angmom[j];

		tmpc = tmpc->next;
		i++;
	}
	//MD particles
	tmd = CL->MDpp;
	while( tmd!=NULL ) {
		MASS = tmd->mass;
		//Position relative to centre of mass
		relQ[i][0] = tmd->rx - CL->CM[0];
		relQ[i][1] = tmd->ry - CL->CM[1];
		relQ[i][2] = tmd->rz - CL->CM[2];
		diffV[0] = MASS * (tmd->vx - RV[i][0]);
		diffV[1] = MASS * (tmd->vy - RV[i][1]);
		diffV[2] = MASS * (tmd->vz - RV[i][2]);
		if( DIM < _3D ) {
			relQ[i][2] = 0.;
			diffV[2] = 0.;
		}
		if( DIM < _2D ) {
			relQ[i][1] = 0.;
			diffV[1] = 0.;
		}
		crossprod( relQ[i],diffV,angmom );
		for( j=0; j<_3D; j++ ) L[j] += angmom[j];

		tmd = tmd->nextSRD;
		i++;
	}
	//Swimmer particles
	tsm = CL->sp;
	while( tsm!=NULL ) {
		if( tsm->HorM ) MASS = (double) SS.middM;
		else MASS = (double) SS.headM;
		//Position relative to centre of mass
		for( j=0; j<DIM; j++ ) relQ[i][j] = tsm->Q[j] - CL->CM[j];
		for( j=0; j<DIM; j++ ) diffV[j] = MASS * (tsm->V[j] - RV[i][j]);
		if( DIM < _3D ) {
			relQ[i][2] = 0.;
			diffV[2] = 0.;
		}
		if( DIM < _2D ) {
			relQ[i][1] = 0.;
			diffV[1] = 0.;
		}
		crossprod( relQ[i],diffV,angmom );
		for( j=0; j<_3D; j++ ) L[j] += angmom[j];
		tsm = tsm->next;
		i++;
	}
	// dotprodmat( L,II,W,_3D );
	dotprodMatVec( II,L,W,_3D );

	/* ****************************************** */
	/* *************** Collision **************** */
	/* ****************************************** */
	i=0;
	//MPC particles
	tmpc = CL->pp;
	while( tmpc!=NULL ) {
		id = tmpc->SPID;
		MASS = (SP+id)->MASS;
		crossprod( W,relQ[i],angterm );
		for( j=0; j<DIM; j++ ) tmpc->V[j] = VCM[j] + RV[i][j] - RS[j] -DV[i][j] + angterm[j];
		//Pressure term
		if( outP ) calcPressureColl_postColl( relQP[i],dp[i],MASS,tmpc->V,CL );
		//Increment link in list
		tmpc = tmpc->next;
		i++;
	}
	//MD particles
	tmd = CL->MDpp;
	while( tmd!=NULL ) {
		crossprod( W,relQ[i],angterm );
		tmd->vx = VCM[0] + RV[i][0] - RS[0] + angterm[0];
		tmd->vy = VCM[1] + RV[i][1] - RS[1] + angterm[1];
		tmd->vz = VCM[2] + RV[i][2] - RS[2] + angterm[2];
		//Increment link in list
		tmd = tmd->nextSRD;
		i++;
	}
	//Swimmer particles
	tsm = CL->sp;
	while( tsm!=NULL ) {
		crossprod( W,relQ[i],angterm );
		for( j=0; j<DIM; j++ ) tsm->V[j] = VCM[j] + RV[i][j] - RS[j] -DV[i][j] + angterm[j];
		//Increment link in list
		tsm = tsm->next;
		i++;
	}
}

/// 
/// @brief Does the Langevin thermostatted collision. 
/// 
/// Invented by Noguchi, Kikuchi and Gompper (https://iopscience.iop.org/article/10.1209/0295-5075/78/10005).
/// Langevin-thermostatted version of multiparticle collision dynamics (MPCD). 
/// Energy is not conserved; rather the system is thermostatted to `KBT`.
/// @param CL An MPCD cell (including the linked list of particles in each cell). 
/// @param SP The species-wide information about MPCD particles.
/// @param SS The species-wide information about swimmers.
/// @param KBT The thermal energy. 
/// @param FRICCO Friction coefficient for Langevin thermostat.
/// @param Step The MPCD time step in MPCD units.
/// @param CLQ The geometric centre of `CL`, the MPCD cell.
/// @param outP Flag whether or not to output the pressure.
///
void langevinMPC( cell *CL,spec *SP,specSwimmer SS,double KBT,double FRICCO,double Step,double *CLQ,int outP ) {
	int i,j,id;
	double MASS;
	double WN[CL->POP][DIM];			//White noise
	double WNS[DIM];					//Sum of white noise
	double DV[CL->POP][DIM];			//Damping velocity
	double a,b;
	double VCM[DIM];					//Centre of mass velocity
	double dp[CL->POP][DIM],relQ[CL->POP][DIM];	//For pressure
	particleMPC *tmpc;					//Temporary particleMPC
	particleMD *tmd;					//Temporary particleMD
	smono *tsm;							//Temporary swimmer monomer

	for( i=0; i<DIM; i++ ) VCM[i] = CL->VCM[i];
	
	/* ****************************************** */
	/* ******* Generate random velocities ******* */
	/* ****************************************** */
	for( j=0; j<DIM; j++ ) WNS[j] = 0.;
	for( i=0;i<CL->POP;i++ ) for( j=0;j<DIM;j++ ) DV[i][j] = 0.;
	i=0;
	//MPC particles
	tmpc = CL->pp;
	while( tmpc!=NULL ) {
		id = tmpc->SPID;
		MASS = (SP+id)->MASS;
		//Pressure term
		if( outP ) calcPressureColl_preColl( relQ[i],dp[i],tmpc,CLQ );
		for( j=0; j<DIM; j++ ) {
			WN[i][j] = genrand_gaussMB( KBT,MASS );
			WNS[j] += WN[i][j];
		}
		for( j=0; j<DIM; j++ ) DV[i][j] = ((SP+id)->DAMP)*(CL->VCM[j])/((double)CL->POP);
		tmpc = tmpc->next;
		i++;
	}
	//MD particles
	tmd = CL->MDpp;
	while( tmd!=NULL ) {
		MASS = tmd->mass;
		for( j=0; j<DIM; j++ ) {
			WN[i][j] = genrand_gaussMB( KBT,MASS );
			WNS[j] += WN[i][j];
		}
		tmd = tmd->nextSRD;
		i++;
	}
	//MPC particles
	tsm = CL->sp;
	while( tsm!=NULL ) {
		if( tsm->HorM ) MASS = (double) SS.middM;
		else MASS = (double) SS.headM;
		for( j=0; j<DIM; j++ ) {
			WN[i][j] = genrand_gaussMB( KBT,MASS );
			WNS[j] += WN[i][j];
		}
		tsm = tsm->next;
		i++;
	}
	for( j=0; j<DIM; j++ ) WNS[j] /= (double)(CL->POP);

	/* ****************************************** */
	/* *************** Collision **************** */
	/* ****************************************** */
	i=0;
	//MPC particles
	tmpc = CL->pp;
	while( tmpc!=NULL ) {
		id = tmpc->SPID;
		MASS = (SP+id)->MASS;
		a = (MASS - FRICCO * Step * 0.5) / (MASS + 0.5 * FRICCO * Step);
		b = sqrt( FRICCO*Step ) / ( MASS + 0.5 * FRICCO * Step );
		for( j=0; j<DIM; j++ ) tmpc->V[j] = VCM[j] + a * (tmpc->V[j]-VCM[j]) + b * (WN[i][j] - WNS[j]) - DV[i][j];
		//Pressure term
		if( outP ) calcPressureColl_postColl( relQ[i],dp[i],MASS,tmpc->V,CL );
		//Increment link in list
		tmpc = tmpc->next;
		i++;
	}
	//MD particles
	tmd = CL->MDpp;
	while( tmd!=NULL ) {
		MASS = tmd->mass;
		a = (MASS - FRICCO * Step * 0.5) / (MASS + 0.5 * FRICCO * Step);
		b = sqrt( FRICCO*Step ) / ( MASS + 0.5 * FRICCO * Step );
		tmd->vx = VCM[0] + a * (tmd->vx-VCM[0]) + b * (WN[i][0] - WNS[0]);
		if( DIM > _1D ){
			tmd->vy = VCM[1] + a * (tmd->vy-VCM[1]) + b * (WN[i][1] - WNS[1]);
		}
		if(DIM > _2D){
			tmd->vz = VCM[2] + a * (tmd->vz-VCM[2]) + b * (WN[i][2] - WNS[2]);
		}
		//Increment link in list
		tmd = tmd->nextSRD;
		i++;
	}
	//Swimmer particles
	tsm = CL->sp;
	while( tsm!=NULL ) {
		if( tsm->HorM ) MASS = (double) SS.middM;
		else MASS = (double) SS.headM;
		a = (MASS - FRICCO * Step * 0.5) / (MASS + 0.5 * FRICCO * Step);
		b = sqrt( FRICCO*Step ) / ( MASS + 0.5 * FRICCO * Step );
		for( j=0; j<DIM; j++ ) tsm->V[j] = VCM[j] + a * (tsm->V[j]-VCM[j]) + b * (WN[i][j] - WNS[j]) - DV[i][j];
		//Increment link in list
		tsm = tsm->next;
		i++;
	}
}

/// 
/// @brief Langevin MPC collision that conserves angular momentum (uses Langevin thermostat). 
/// 
/// Invented by Noguchi, Kikuchi and Gompper (https://iopscience.iop.org/article/10.1209/0295-5075/78/10005).
/// Langevin-thermostatted version of multiparticle collision dynamics (MPCD) that <b>does</b> conserver angular momentum. 
/// It is just like langevinMPC() but includes a correction to conserve angular momentum.
/// @param CL An MPCD cell (including the linked list of particles in each cell). 
/// @param SP The species-wide information about MPCD particles.
/// @param SS The species-wide information about swimmers.
/// @param KBT The thermal energy. 
/// @param FRICCO Friction coefficient for Langevin thermostat.
/// @param Step The MPCD time step in MPCD units.
/// @param CLQ The geometric centre of `CL`, the MPCD cell.
/// @param outP Flag whether or not to output the pressure.
/// @see langevinMPC()
///
void langevinROT( cell *CL,spec *SP,specSwimmer SS,double KBT,double FRICCO,double Step,double *CLQ,int outP ) {
	int i,j,id;
	double MASS;
	double WN[CL->POP][DIM];			//White noise
	double WNS[DIM];					//Sum of white noise
	double DV[CL->POP][DIM];			//Damping velocity
	double a,b;
	double VCM[DIM];					//Centre of mass velocity
	double dp[CL->POP][DIM];			//For pressure
	double relQ[CL->POP][_3D];			//Relative position
	double diffV[_3D];					//Difference in velocity
	double L[_3D];						//Angular momentum
	double angterm[_3D];
	double angmom[_3D];
	double W[_3D];
	double II[_3D][_3D];				//Inverse of moment of inertia tensor (3D)
	particleMPC *tmpc;					//Temporary particleMPC
	particleMD *tmd;					//Temporary particleMD
	smono *tsm;							//Temporary swimmer monomer

	for( i=0; i<DIM; i++ ) {
		VCM[i] = CL->VCM[i];
		WNS[i] = 0.0;
	}
	for( i=0;i<CL->POP;i++ ) for( j=0;j<DIM;j++ ) {
		WN[i][j]=0.0;
		DV[i][j] = 0.0;
	}
	for( i=0;i<CL->POP;i++ ) for( j=0;j<_3D;j++ ) {
		relQ[i][j] = 0.0;
		dp[i][j] = 0.0;
	}
	for( j=0;j<_3D;j++ ) {
		diffV[j]=0.;
		L[j]=0.;
		angterm[j]=0.;
		angmom[j]=0.;
		W[j]=0.;
		VCM[j] = CL->VCM[j];
		for( i=0;i<_3D;i++ ) II[j][i]=0.;
	}
	
	/* ****************************************** */
	/* ******* Generate random velocities ******* */
	/* ****************************************** */
	for( j=0; j<DIM; j++ ) WNS[j] = 0.;
	for( i=0;i<CL->POP;i++ ) for( j=0;j<DIM;j++ ) DV[i][j] = 0.;
	i=0;
	//MPC particles
	tmpc = CL->pp;
	while( tmpc!=NULL ) {
		id = tmpc->SPID;
		MASS = (SP+id)->MASS;
		//Pressure term
		if( outP ) calcPressureColl_preColl( relQ[i],dp[i],tmpc,CLQ );
		for( j=0; j<DIM; j++ ) {
			WN[i][j] = genrand_gaussMB( KBT,MASS );
			WNS[j] += WN[i][j];
		}
		for( j=0; j<DIM; j++ ) DV[i][j] = ((SP+id)->DAMP)*(CL->VCM[j])/((double)CL->POP);
		tmpc = tmpc->next;
		i++;
	}
	//MD particles
	tmd = CL->MDpp;
	while( tmd!=NULL ) {
		MASS = tmd->mass;
		for( j=0; j<DIM; j++ ) {
			WN[i][j] = genrand_gaussMB( KBT,MASS );
			WNS[j] += WN[i][j];
		}
		tmd = tmd->nextSRD;
		i++;
	}
	//MPC particles
	tsm = CL->sp;
	while( tsm!=NULL ) {
		if( tsm->HorM ) MASS = (double) SS.middM;
		else MASS = (double) SS.headM;
		for( j=0; j<DIM; j++ ) {
			WN[i][j] = genrand_gaussMB( KBT,MASS );
			WNS[j] += WN[i][j];
		}
		tsm = tsm->next;
		i++;
	}
	for( j=0; j<DIM; j++ ) WNS[j] /= (double)(CL->POP);

	/* ****************************************** */
	/* ******** Invert moment of inertia ******** */
	/* ****************************************** */
	for( i=0; i<_3D; i++ ) for( j=0; j<_3D; j++ ) II[i][j] = 0.;
	if( DIM == _3D ) {
		// When the population is 2 (don't even get this far if pop=1) then the moment of inertia tensor is singular
		// In that case, can't invert. Therefore, we must neglect angular momentum conservation in such cell
		// This just makes angterm = [0,0,0]
		if( CL->POP < 3 ) for( i=0; i<_3D; i++ ) for( j=0; j<_3D; j++ ) II[i][j] = 0.;
		else invert3x3(II,CL->I);
	}
	else if( DIM == _2D ) {
		II[2][2] = 1./CL->I[2][2];
		// The rest don't matter and can remain zero
	}
	else {
		printf( "Error: Angular momentum conservation in 1D is nonsequitur. Higher dimensions not done. Change collision technique or dimension.\n" );
		exit( 1 );
	}

	/* ****************************************** */
	/* ******* Find angular velocity term ******* */
	/* ****************************************** */
	for( i=0; i<_3D; i++ ) L[i] = 0.;
	i=0;
	//MPC particles
	tmpc = CL->pp;
	while( tmpc!=NULL ) {
		id = tmpc->SPID;
		MASS = (SP+id)->MASS;
		//Position relative to centre of mass
		for( j=0; j<DIM; j++ ) relQ[i][j] = tmpc->Q[j] - CL->CM[j];
		//Difference in velocity
		for( j=0; j<DIM; j++ ) diffV[j] = MASS * ( FRICCO*tmpc->V[j] - sqrt(FRICCO) * WN[i][j] );
		if( DIM < _3D ) {
			relQ[i][2] = 0.;
			diffV[2] = 0.;
		}
		if( DIM < _2D ) {
			relQ[i][1] = 0.;
			diffV[1] = 0.;
		}
		crossprod( relQ[i],diffV,angmom );
		for( j=0; j<_3D; j++ ) L[j] += angmom[j];

		tmpc = tmpc->next;
		i++;
	}
	//MD particles
	tmd = CL->MDpp;
	while( tmd!=NULL ) {
		MASS = tmd->mass;
		//Position relative to centre of mass
		relQ[i][0] = tmd->rx - CL->CM[0];
		relQ[i][1] = tmd->ry - CL->CM[1];
		relQ[i][2] = tmd->rz - CL->CM[2];
		diffV[0] = MASS * (FRICCO*tmd->vx - sqrt(FRICCO)*WN[i][0]);
		diffV[1] = MASS * (FRICCO*tmd->vy - sqrt(FRICCO)*WN[i][1]);
		diffV[2] = MASS * (FRICCO*tmd->vz - sqrt(FRICCO)*WN[i][2]);
		if( DIM < _3D ) {
			relQ[i][2] = 0.;
			diffV[2] = 0.;
		}
		if( DIM < _2D ) {
			relQ[i][1] = 0.;
			diffV[1] = 0.;
		}
		crossprod( relQ[i],diffV,angmom );
		for( j=0; j<_3D; j++ ) L[j] += angmom[j];

		tmd = tmd->nextSRD;
		i++;
	}
	//Swimmer particles
	tsm = CL->sp;
	while( tsm!=NULL ) {
		if( tsm->HorM ) MASS = (double) SS.middM;
		else MASS = (double) SS.headM;
		//Position relative to centre of mass
		for( j=0; j<DIM; j++ ) relQ[i][j] = tsm->Q[j] - CL->CM[j];
		for( j=0; j<DIM; j++ ) diffV[j] = MASS * ( FRICCO*tsm->V[j] - sqrt(FRICCO) * WN[i][j] );
		if( DIM < _3D ) {
			relQ[i][2] = 0.;
			diffV[2] = 0.;
		}
		if( DIM < _2D ) {
			relQ[i][1] = 0.;
			diffV[1] = 0.;
		}
		crossprod( relQ[i],diffV,angmom );
		for( j=0; j<_3D; j++ ) L[j] += angmom[j];
		tsm = tsm->next;
		i++;
	}
	dotprodMatVec( II,L,W,_3D );

	/* ****************************************** */
	/* *************** Collision **************** */
	/* ****************************************** */
	i=0;
	//MPC particles
	tmpc = CL->pp;
	while( tmpc!=NULL ) {
		id = tmpc->SPID;
		MASS = (SP+id)->MASS;
		crossprod( W,relQ[i],angterm );
		a = (MASS - FRICCO * Step * 0.5) / (MASS + 0.5 * FRICCO * Step);
		b = sqrt( FRICCO*Step ) / ( MASS + 0.5 * FRICCO * Step );
		for( j=0; j<DIM; j++ ) tmpc->V[j] = VCM[j] + a * (tmpc->V[j]-VCM[j]) + b * (WN[i][j] - WNS[j]) - DV[i][j] + angterm[j];
		//Pressure term
		if( outP ) calcPressureColl_postColl( relQ[i],dp[i],MASS,tmpc->V,CL );
		//Increment link in list
		tmpc = tmpc->next;
		i++;
	}
	//MD particles
	tmd = CL->MDpp;
	while( tmd!=NULL ) {
		MASS = tmd->mass;
		crossprod( W,relQ[i],angterm );
		a = (MASS - FRICCO * Step * 0.5) / (MASS + 0.5 * FRICCO * Step);
		b = sqrt( FRICCO*Step ) / ( MASS + 0.5 * FRICCO * Step );
		tmd->vx = VCM[0] + a * (tmd->vx-VCM[0]) + b * (WN[i][0] - WNS[0]) + angterm[0];
		tmd->vy = VCM[1] + a * (tmd->vy-VCM[1]) + b * (WN[i][1] - WNS[1]) + angterm[1];
		tmd->vz = VCM[2] + a * (tmd->vz-VCM[2]) + b * (WN[i][2] - WNS[2]) + angterm[2];
		//Increment link in list
		tmd = tmd->nextSRD;
		i++;
	}
	//Swimmer particles
	tsm = CL->sp;
	while( tsm!=NULL ) {
		if( tsm->HorM ) MASS = (double) SS.middM;
		else MASS = (double) SS.headM;
		crossprod( W,relQ[i],angterm );
		a = (MASS - FRICCO * Step * 0.5) / (MASS + 0.5 * FRICCO * Step);
		b = sqrt( FRICCO*Step ) / ( MASS + 0.5 * FRICCO * Step );
		for( j=0; j<DIM; j++ ) tsm->V[j] = VCM[j] + a * (tsm->V[j]-VCM[j]) + b * (WN[i][j] - WNS[j]) - DV[i][j] + angterm[j];
		//Increment link in list
		tsm = tsm->next;
		i++;
	}
}

/// 
/// @brief An active stochastic rotation dynamics (SRD) algorithm. 
///
/// @warning Experimental. Proposed by Shendruk but not published or fully characterized yet. 
///
/// An active-version of SRD. 
/// Does the stochastic rotation collision (see stochrotMPC()). 
/// But then rotates the resulting velocities towards the centre of mass velocity injecting momentum but keeping the energy constant.
/// @param CL An MPCD cell (including the linked list of particles in each cell). 
/// @param SP The species-wide information about MPCD particles.
/// @param RTECH The MPCD collision operator. See `definitions.h` for all options.
/// @param C Cosine of the rotation angle.
/// @param S Sine of the rotation angle.
/// @param CLQ The geometric centre of `CL`, the MPCD cell.
/// @param outP Flag whether or not to output the pressure.
/// @see stochrotMPC()
///
void activeSRD( cell *CL,spec *SP,int RTECH,double C,double S,double *CLQ,int outP ) {
	int i;
	double theta;				// Angle between individual velocity and vcm
	double sumACT=0.0;			// Sum of the activity of particles in the cell
	double nv[_3D]={0.0};		// normal vector between vcm and individual vel --- must be _3DD regardless of actual dimension
	double VCM[_3D]={0.0};		// Short-hand CM velocity --- must be 3D regardless of actual dimension
	particleMPC *tmpc;			//Temporary particleMPC

	for( i=0; i<_3D; i++ ) VCM[i] = CL->VCM[i];

	//Do the normal SRD collision
	if( RTECH==ACT_ARBAXIS) stochrotMPC( CL,ARBAXIS,C,S,CLQ,outP );
	else if( RTECH==ACT_ORTHAXIS) stochrotMPC( CL,ORTHAXIS,C,S,CLQ,outP );
	else{
		printf( "Error: Unexpected active collision operator.\n" );
		exit(1);
	}

	//Active part
	//Calculate the activity of the cell
	tmpc = CL->pp;
	while( tmpc!=NULL ) {
		sumACT = (SP+(tmpc->SPID))->ACT;
		//Increment link in list
		tmpc = tmpc->next;
	}
	//Rotate velocities towards the vcm
	//MPC particles
	tmpc = CL->pp;
	while( tmpc!=NULL ) {
		if( fneq( (SP+(tmpc->SPID))->ACT,0.0 ) ) {
			// Find the angle
			theta=absAngle(tmpc->V,VCM,_3D);
			// Find the normal --- don't be sloppy in 2D. You need the sign
			crossprod( tmpc->V,VCM,nv );
			// Rotate the velocity about the normal toward vcm by an amount ACT*theta
			//rotate( ARBAXIS,cos(sumACT*theta),sin(sumACT*theta),tmpc->V,nv,1.0,1.0 );
			// The Rodrigues' rotation formula is faster than using a matrix rotation
			if( fabs(theta)>ROTTOL ) rodriguesRotation( tmpc->V,nv,sumACT*theta );
			//Increment link in list
		}
		tmpc = tmpc->next;
	}
	//MD particles do not participate in the active impulse
}

/// 
/// @brief Active dry polar collision operator based on Vicsek algorithm. 
///
/// @warning Experimental. Proposed by Shendruk but not published or fully characterized yet. 
///
/// This is meant to be an MPCD-version of the Vicsek algorithm.
/// Instead of an interaction radius around every particle (as in traditional Vicsek), the alignment occurs within an MPCD cell. 
///	Here, `ACT` serves as the <b>noise</b> range (typically, \f$\eta \f$ in Vicsek models).
/// - The original Vicsek paper https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.75.1226
/// - A more recent review of the Vicsek model https://link.springer.com/article/10.1140/epjst/e2016-60066-8
/// @param CL An MPCD cell (including the linked list of particles in each cell). 
/// @param SP The species-wide information about MPCD particles.
/// @param CLQ The geometric centre of `CL`, the MPCD cell.
/// @param outP Flag whether or not to output the pressure.
///
void vicsek( cell *CL,spec *SP,double *CLQ,int outP ) {
	int i,j,id;
	double sigma;										//The range on the noise is (1-ACT) and is between [0,1]
	double noiseRange,speed;
	double VCM[_3D],VCMdir[_3D];						//Short-hand CM velocity
	double MASS,dp[CL->POP][DIM],relQ[CL->POP][DIM];	//For pressure
	particleMPC *tmpc;									//Temporary particleMPC

	//Must be 3D because genrand_cone needs 3D vectors even in 2D
	for( i=0; i<DIM; i++ ) VCM[i] = CL->VCM[i];
	for( i=DIM; i<_3D; i++ ) VCM[i] = 0.;
	normCopy( VCM,VCMdir,_3D );

	/* ****************************************** */
	/* ******* Generate random velocities ******* */
	/* ****************************************** */
	//Generate (homogeneously) random velocities on the cone (1-ACT)*pi
	// MPC particles
	tmpc = CL->pp;
	i=0;
	while( tmpc!=NULL ) {
		id = tmpc->SPID;
		MASS = (SP+id)->MASS;
		//Pressure term
		if( outP ) calcPressureColl_preColl( relQ[i],dp[i],tmpc,CLQ );
		sigma = 1.0 - (SP+id)->ACT;
		noiseRange = pi*sigma;
		speed=length(tmpc->V,DIM);
		genrand_cone( VCMdir,tmpc->V,noiseRange,DIM );
		for( j=0; j<DIM; j++ ) tmpc->V[j] *= speed;
		//Pressure term
		if( outP ) calcPressureColl_postColl( relQ[i],dp[i],MASS,tmpc->V,CL );
		tmpc = tmpc->next;
		i++;
	}
}

/// 
/// @brief Active dry nematic collision operator based on Chate algorithm. 
/// 
/// @warning Experimental. Proposed by Shendruk but not published or fully characterized yet. 
///
/// This is meant to be an MPCD-version of Chate's nematic version of the Vicsek model.
/// Instead of an interaction radius around every particle (as in traditional Vicsek), the alignment occurs within an MPCD cell. 
///	Here, `ACT` serves as the <b>noise</b> range (typically, \f$\eta \f$ in Vicsek models).
/// - https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.96.180602
/// It updates the particles by looping through the linked list. 
/// @param CL An MPCD cell (including the linked list of particles in each cell). 
/// @param SP The species-wide information about MPCD particles.
/// @param CLQ The geometric centre of `CL`, the MPCD cell.
/// @param outP Flag whether or not to output the pressure.
/// @see vicsek()
///
void chate( cell *CL,spec *SP,double *CLQ,int outP ) {
	int i,j,id;
	double sigma;											//The range on the noise is (1-ACT) and is between [0,1]
	double noiseRange,speed,pmOne;
	double DIR[_3D];										//Short-hand cell's director
	double MASS,dp[CL->POP][DIM],relQ[CL->POP][DIM];		//For pressure
	particleMPC *tmpc;										//Temporary particleMPC

	//Must be 3D because genrand_cone needs 3D vectors even in 2D
	for( i=0; i<DIM; i++ ) DIR[i] = CL->DIR[i];
	for( i=DIM; i<_3D; i++ ) DIR[i] = 0.;

	/* ****************************************** */
	/* ******* Generate random velocities ******* */
	/* ****************************************** */
	//Generate (homogeneously) random velocities on the cone (1-ACT)*pi
	// MPC particles
	tmpc = CL->pp;
	i=0;
	while( tmpc!=NULL ) {
		id = tmpc->SPID;
		MASS = (SP+id)->MASS;
		//Pressure term
		if( outP ) calcPressureColl_preColl( relQ[i],dp[i],tmpc,CLQ );
		sigma = 1.0 - (SP+id)->ACT;
		noiseRange = 0.5*pi*sigma;
		speed=length(tmpc->V,DIM);
		genrand_cone( DIR,tmpc->V,noiseRange,DIM );
		pmOne = genrand_pmOne();
		for( j=0; j<DIM; j++ ) tmpc->V[j] *= pmOne*speed;
		//Pressure term
		if( outP ) calcPressureColl_postColl( relQ[i],dp[i],MASS,tmpc->V,CL );
		tmpc = tmpc->next;
		i++;
	}
}

/// 
/// @brief Active dry polar collision operator based on Andersen thermostatted collision. 
///  
/// @warning Experimental. Proposed by Shendruk but not published or fully characterized yet. 
///
/// An active-version of Andersen-MPCD. 
/// Like the Andersen-version of MPCD (MPCD-AT) but it relaxes to the speed `ACT` instead of the center of mass speed itself. 
/// The active velocity has speed `ACT` in the direction of the centre of mass velocity. 
/// @param CL An MPCD cell (including the linked list of particles in each cell). 
/// @param SP The species-wide information about MPCD particles.
/// @param KBT The thermal energy. 
/// @param RELAX The relaxation time scale towards the active speed `ACT`.
/// @param CLQ The geometric centre of `CL`, the MPCD cell.
/// @param outP Flag whether or not to output the pressure.
/// @see andersenMPC()
///
void vicsekAndersenMPC( cell *CL,spec *SP,double KBT,double RELAX,double *CLQ,int outP ) {
	int i,j,id;
	double MASS,ACT;
	double RV[CL->POP][DIM];						//Random velocities
	double RS[DIM];									//Sum of random velocities
	double dp[CL->POP][DIM],relQ[CL->POP][DIM];		//For pressure
	particleMPC *tmpc;								//Temporary particleMPC
	particleMD *tmd;								//Temporary particleMD
	double VCM[DIM],VCMdir[DIM];					//Short-hand velocity of CM and its direction

	//The direction of the centre of mass velocity of the cell
	for( i=0; i<DIM; i++ ) VCM[i] = CL->VCM[i];
	normCopy( VCM,VCMdir,DIM );
	// Zero arrays
	for( i=0; i<DIM; i++ ) RS[i] = 0.;
	for( i=0;i<CL->POP;i++ ) for( j=0;j<DIM;j++ ) RV[i][j] = 0.;

	/* ****************************************** */
	/* ******* Generate random velocities ******* */
	/* ****************************************** */
	//Generate random velocities and give the active shift.
	// MPC particles
	tmpc = CL->pp;
	i=0;
	while( tmpc!=NULL ) {
		id = tmpc->SPID;
		MASS = (SP+id)->MASS;
		//Pressure term
		if( outP ) calcPressureColl_preColl( relQ[i],dp[i],tmpc,CLQ );
		for( j=0; j<DIM; j++ ) RV[i][j] = genrand_gaussMB( KBT,MASS );
		for( j=0; j<DIM; j++ ) RS[j] += MASS*RV[i][j];
		tmpc = tmpc->next;
		i++;
	}
	// MD particles
	tmd = CL->MDpp;
	while( tmd!=NULL ) {
		MASS = tmd->mass;
		for( j=0; j<DIM; j++ ) RV[i][j] = genrand_gaussMB( KBT,MASS );
		for( j=0; j<DIM; j++ ) RS[j] += MASS*RV[i][j];
		tmd = tmd->nextSRD;
		i++;
	}
	// Normalize
	for( j=0; j<DIM; j++ ) RS[j] /= CL->MASS;

	/* ****************************************** */
	/* *************** Collision **************** */
	/* ****************************************** */
	// MPC particles
	tmpc = CL->pp;
	i=0;
	while( tmpc!=NULL ) {
		id = tmpc->SPID;
		MASS = (SP+id)->MASS;
		ACT = (SP+id)->ACT;
		//Passive Particles
		if( feq(ACT,0.0) ) for( j=0; j<DIM; j++ ) tmpc->V[j] = VCM[j] + RV[i][j] - RS[j];
		//Active Particles
		else{
			for( j=0; j<DIM; j++ ) tmpc->V[j] = VCM[j] + RELAX*(ACT*VCMdir[j]-tmpc->V[j]) + RV[i][j] - RS[j];
		}
		//Pressure term
		if( outP ) calcPressureColl_postColl( relQ[i],dp[i],MASS,tmpc->V,CL );
		//Increment link in list
		tmpc = tmpc->next;
		i++;
	}
	//MD particles
	tmd = CL->MDpp;
	while( tmd!=NULL ) {
		tmd->vx = VCM[0] + RV[i][0] - RS[0];
		if( DIM > _1D){
			tmd->vy = VCM[1] + RV[i][1] - RS[1];
		}
		if( DIM > _2D){
			tmd->vz = VCM[2] + RV[i][2] - RS[2];
		}
		//Increment link in list
		tmd = tmd->nextSRD;
		i++;
	}
}

/// 
/// @brief Active dry nematic collision operator based on Andersen thermostatted collision. 
/// 
/// @warning Experimental. Proposed by Shendruk but not published or fully characterized yet. 
///
/// Just like the Vicsek-Andersen version of active-MPCD (vicsekAndersenMPC()).
/// The noise is about the director; rather than the centre of mass velocity. 
/// It relaxes to the speed `ACT` in a rondom selection of +/- the director direction. 
/// @param CL An MPCD cell (including the linked list of particles in each cell). 
/// @param SP The species-wide information about MPCD particles.
/// @param KBT The thermal energy. 
/// @param RELAX The temperature relaxation time scale.
/// @param CLQ The geometric centre of `CL`, the MPCD cell.
/// @param outP Flag whether or not to output the pressure.
/// @see vicsekAndersenMPC()
///
void chateAndersenMPC( cell *CL,spec *SP,double KBT,double RELAX,double *CLQ,int outP ) {
	int i,j,id;
	double MASS,ACT,pmOne;
	double RV[CL->POP][DIM];						//Random velocities
	double RS[DIM];									//Sum of random velocities
	double AV[CL->POP][DIM];						//Active velocities
	double AS[DIM];									//Sum of active velocities
	double dp[CL->POP][DIM],relQ[CL->POP][DIM];		//For pressure
	particleMPC *tmpc;								//Temporary particleMPC
	particleMD *tmd;								//Temporary particleMD
	double VCM[DIM],speed;							//Short-hand velocity of CM
	double DIR[DIM];								//Order parameter and director

	//Calculate the director of the cell
	for( i=0; i<DIM; i++ ) {
		VCM[i] = CL->VCM[i];
		DIR[i] = CL->DIR[i];
	}
	// Zero arrays
	for( i=0; i<DIM; i++ ) {
		RS[i] = 0.;
		AS[i] = 0.;
		for( j=0;j<CL->POP;j++ ) {
			RV[j][i] = 0.;
			AV[j][i] = 0.;
		}
	}

	/* ****************************************** */
	/* ******* Generate random velocities ******* */
	/* ****************************************** */
	//Generate random velocities and give the active shift.
	// MPC particles
	tmpc = CL->pp;
	i=0;
	while( tmpc!=NULL ) {
		id = tmpc->SPID;
		MASS = (SP+id)->MASS;
		ACT = (SP+id)->ACT;
		//Pressure term
		if( outP ) calcPressureColl_preColl( relQ[i],dp[i],tmpc,CLQ );
		for( j=0; j<DIM; j++ ) RV[i][j] = genrand_gaussMB( KBT,MASS );
		for( j=0; j<DIM; j++ ) RS[j] += MASS*RV[i][j];
		if( fneq(ACT,0.0) ) {
			speed=length(tmpc->V,DIM);
			pmOne = genrand_pmOne();
			for( j=0; j<DIM; j++ ) AV[i][j] = pmOne*RELAX*(ACT-speed)*DIR[j];
			for( j=0; j<DIM; j++ ) AS[j] += MASS*AV[i][j];
		}
		tmpc = tmpc->next;
		i++;
	}
	// MD particles
	tmd = CL->MDpp;
	while( tmd!=NULL ) {
		MASS = tmd->mass;
		for( j=0; j<DIM; j++ ) RV[i][j] = genrand_gaussMB( KBT,MASS );
		for( j=0; j<DIM; j++ ) RS[j] += MASS*RV[i][j];
		tmd = tmd->nextSRD;
		i++;
	}
	// Normalize
	for( j=0; j<DIM; j++ ) RS[j] /= CL->MASS;
	for( j=0; j<DIM; j++ ) AS[j] /= CL->MASS;

	/* ****************************************** */
	/* *************** Collision **************** */
	/* ****************************************** */
	// MPC particles
	tmpc = CL->pp;
	i=0;
	while( tmpc!=NULL ) {
		id = tmpc->SPID;
		MASS = (SP+id)->MASS;
		for( j=0; j<DIM; j++ ) tmpc->V[j] = VCM[j] + RV[i][j] - RS[j] + AV[i][j] - AS[j];
		//Pressure term
		if( outP ) calcPressureColl_postColl( relQ[i],dp[i],MASS,tmpc->V,CL );
		//Increment link in list
		tmpc = tmpc->next;
		i++;
	}
	//MD particles
	tmd = CL->MDpp;
	while( tmd!=NULL ) {
		tmd->vx = VCM[0] + RV[i][0] - RS[0];
		if( DIM > _1D){
			tmd->vy = VCM[1] + RV[i][1] - RS[1];
		}
		if( DIM > _2D){
			tmd->vz = VCM[2] + RV[i][2] - RS[2];
		}
		//Increment link in list
		tmd = tmd->nextSRD;
		i++;
	}
}

/// 
/// @brief Active dry polar collision operator based on Langevin thermostatted collision. 
/// 
/// @warning Experimental. Proposed by Shendruk but not published or fully characterized yet. 
///
/// Just like the Vicsek-Andersen version of active-MPCD (vicsekAndersenMPC()) but it is Langevin thermostatted. 
/// @param CL An MPCD cell (including the linked list of particles in each cell). 
/// @param SP The species-wide information about MPCD particles.
/// @param KBT The thermal energy. 
/// @param FRICCO Friction coefficient for Langevin thermostat.
/// @param Step The MPCD time step in MPCD units.
/// @param RELAX The temperature relaxation time scale.
/// @param CLQ The geometric centre of `CL`, the MPCD cell.
/// @param outP Flag whether or not to output the pressure.
/// @see vicsekAndersenMPC()
///
void vicsekLangevinMPC( cell *CL,spec *SP,double KBT,double FRICCO,double Step,double RELAX,double *CLQ,int outP ) {
	int i,j,id;
	double MASS,ACT;
	double WN[CL->POP][DIM];						//White noise
	double WNS[DIM];								//Sum of white noise
	double a,b;
	double VCM[DIM],VCMdir[DIM];					//Centre of mass velocity
	double dp[CL->POP][DIM],relQ[CL->POP][DIM];		//For pressure
	particleMPC *tmpc;								//Temporary particleMPC
	particleMD *tmd;								//Temporary particleMD

	for( i=0; i<DIM; i++ ) VCM[i] = CL->VCM[i];
	normCopy( VCM,VCMdir,DIM );

	/* ****************************************** */
	/* ******* Generate random velocities ******* */
	/* ****************************************** */
	//Generate random velocities
	for( j=0; j<DIM; j++ ) WNS[j] = 0.;
	i=0;
	//MPC particles
	tmpc = CL->pp;
	while( tmpc!=NULL ) {
		id = tmpc->SPID;
		MASS = (SP+id)->MASS;
		//Pressure term
		if( outP ) calcPressureColl_preColl( relQ[i],dp[i],tmpc,CLQ );
		for( j=0; j<DIM; j++ ) {
			WN[i][j] = genrand_gaussMB( KBT,MASS );
			WNS[j] += WN[i][j];
		}
		tmpc = tmpc->next;
		i++;
	}
	//MD particles
	tmd = CL->MDpp;
	while( tmd!=NULL ) {
		MASS = tmd->mass;
		for( j=0; j<DIM; j++ ) {
			WN[i][j] = genrand_gaussMB( KBT,MASS );
			WNS[j] += WN[i][j];
		}
		tmd = tmd->nextSRD;
		i++;
	}
	for( j=0; j<DIM; j++ ) WNS[j] /= (double)(CL->POP);

	/* ****************************************** */
	/* *************** Collision **************** */
	/* ****************************************** */
	i=0;
	//MPC particles
	tmpc = CL->pp;
	while( tmpc!=NULL ) {
		id = tmpc->SPID;
		MASS = (SP+id)->MASS;
		ACT = (SP+id)->ACT;
		a = (MASS - FRICCO * Step * 0.5) / (MASS + 0.5 * FRICCO * Step);
		b = sqrt( FRICCO*Step ) / ( MASS + 0.5 * FRICCO * Step );
		//Passive Particles
		if( feq(ACT,0.0) ) for( j=0; j<DIM; j++ ) tmpc->V[j] = VCM[j] + a * (tmpc->V[j]-VCM[j]) + b * (WN[i][j] - WNS[j]);
		//Active Particles
		else{
			//for( j=0; j<DIM; j++ ) tmpc->V[j] = VCM[j] + RELAX*(ACT*VCMdir[j]-VCM[j]) + a * (tmpc->V[j]-VCM[j]) + b * (WN[i][j] - WNS[j]);
			for( j=0; j<DIM; j++ ) tmpc->V[j] = VCM[j] + RELAX*(ACT*VCMdir[j]-tmpc->V[j]) + a * (tmpc->V[j]-VCM[j]) + b * (WN[i][j] - WNS[j]);
		}
		//Pressure term
		if( outP ) calcPressureColl_postColl( relQ[i],dp[i],MASS,tmpc->V,CL );
		//Increment link in list
		tmpc = tmpc->next;
		i++;
	}
	//MD particles
	tmd = CL->MDpp;
	while( tmd!=NULL ) {
		MASS = tmd->mass;
		a = (MASS - FRICCO * Step * 0.5) / (MASS + 0.5 * FRICCO * Step);
		b = sqrt( FRICCO*Step ) / ( MASS + 0.5 * FRICCO * Step );
		tmd->vx = VCM[0] + a * (tmd->vx-VCM[0]) + b * (WN[i][0] - WNS[0]);
		if( DIM > _1D){
			tmd->vy = VCM[1] + a * (tmd->vy-VCM[1]) + b * (WN[i][1] - WNS[1]);
		}
		if( DIM > _2D){
			tmd->vz = VCM[2] + a * (tmd->vz-VCM[2]) + b * (WN[i][2] - WNS[2]);
		}
		//Increment link in list
		tmd = tmd->nextSRD;
		i++;
	}
}

/// 
/// @brief Active dry nematic collision operator based on Langevin thermostatted collision. 
/// 
/// @warning Experimental. Proposed by Shendruk but not published or fully characterized yet. 
///
/// Just like the Chate-Andersen version of active-MPCD (chateAndersenMPC()) but it is Langevin thermostatted. 
/// @param CL An MPCD cell (including the linked list of particles in each cell). 
/// @param SP The species-wide information about MPCD particles.
/// @param KBT The thermal energy. 
/// @param FRICCO Friction coefficient for Langevin thermostat.
/// @param Step The MPCD time step in MPCD units.
/// @param RELAX The temperature relaxation time scale.
/// @param CLQ The geometric centre of `CL`, the MPCD cell.
/// @param outP Flag whether or not to output the pressure.
/// @see chateAndersenMPC()
///
void chateLangevinMPC( cell *CL,spec *SP,double KBT,double FRICCO,double Step,double RELAX,double *CLQ,int outP ) {
	int i,j,id;
	double MASS,ACT,pmOne;
	double WN[CL->POP][DIM];						//White noise
	double WNS[DIM];								//Sum of white noise
	double a,b;
	double dp[CL->POP][DIM],relQ[CL->POP][DIM];		//For pressure
	particleMPC *tmpc;								//Temporary particleMPC
	particleMD *tmd;								//Temporary particleMD
	double VCM[DIM],speed;							//Centre of mass velocity
	double DIR[DIM];								//Order parameter and director

	//Calculate the director of the cell
	for( i=0; i<DIM; i++ ) {
		VCM[i] = CL->VCM[i];
		DIR[i] = CL->DIR[i];
	}
	speed=length(VCM,DIM);

	/* ****************************************** */
	/* ******* Generate random velocities ******* */
	/* ****************************************** */
	//Generate random velocities
	for( j=0; j<DIM; j++ ) WNS[j] = 0.;
	i=0;
	//MPC particles
	tmpc = CL->pp;
	while( tmpc!=NULL ) {
		id = tmpc->SPID;
		MASS = (SP+id)->MASS;
		//Pressure term
		if( outP ) calcPressureColl_preColl( relQ[i],dp[i],tmpc,CLQ );
		for( j=0; j<DIM; j++ ) {
			WN[i][j] = genrand_gaussMB( KBT,MASS );
			WNS[j] += WN[i][j];
		}
		tmpc = tmpc->next;
		i++;
	}
	//MD particles
	tmd = CL->MDpp;
	while( tmd!=NULL ) {
		MASS = tmd->mass;
		for( j=0; j<DIM; j++ ) {
			WN[i][j] = genrand_gaussMB( KBT,MASS );
			WNS[j] += WN[i][j];
		}
		tmd = tmd->nextSRD;
		i++;
	}
	for( j=0; j<DIM; j++ ) WNS[j] /= (double)(CL->POP);

	/* ****************************************** */
	/* *************** Collision **************** */
	/* ****************************************** */
	i=0;
	//MPC particles
	tmpc = CL->pp;
	while( tmpc!=NULL ) {
		id = tmpc->SPID;
		MASS = (SP+id)->MASS;
		ACT = (SP+id)->ACT;
		a = (MASS - FRICCO * Step * 0.5) / (MASS + 0.5 * FRICCO * Step);
		b = sqrt( FRICCO*Step ) / ( MASS + 0.5 * FRICCO * Step );
		if( fneq(ACT,0.0) ) {
			pmOne = genrand_pmOne();
			for( j=0; j<DIM; j++ ) tmpc->V[j] = VCM[j] + a * pmOne*RELAX*(ACT-speed)*DIR[j] + b * (WN[i][j] - WNS[j]);
		}
		else for( j=0; j<DIM; j++ ) tmpc->V[j] = VCM[j] + a * (tmpc->V[j]-VCM[j]) + b * (WN[i][j] - WNS[j]);
		//Pressure term
		if( outP ) calcPressureColl_postColl( relQ[i],dp[i],MASS,tmpc->V,CL );
		//Increment link in list
		tmpc = tmpc->next;
		i++;
	}
	//MD particles
	tmd = CL->MDpp;
	while( tmd!=NULL ) {
		MASS = tmd->mass;
		a = (MASS - FRICCO * Step * 0.5) / (MASS + 0.5 * FRICCO * Step);
		b = sqrt( FRICCO*Step ) / ( MASS + 0.5 * FRICCO * Step );
		tmd->vx = VCM[0] + a * (tmd->vx-VCM[0]) + b * (WN[i][0] - WNS[0]);
		if( DIM > _1D){
			tmd->vy = VCM[1] + a * (tmd->vy-VCM[1]) + b * (WN[i][1] - WNS[1]);
		}
		if( DIM > _2D){
			tmd->vz = VCM[2] + a * (tmd->vz-VCM[2]) + b * (WN[i][2] - WNS[2]);
		}
		//Increment link in list
		tmd = tmd->nextSRD;
		i++;
	}
}

/// 
/// @brief Active wet nematic collision operator, based on Andersen thermostatted collision. 
///
/// Active-nematic version of Andersen-thermostatted version of multiparticle collision dynamics (andersenMPC()). 
/// Andersen thermostat collision, and returns the CM velocity and the local temperature of the cell. 
/// But also applies a force dipole on the cell by applying a kick along the direction of the cell director. 
/// Together, the centre of mass position of the cell and the director of the cell define a plane. 
/// Particles above the plane get a positive impulse (instantaneous change in velocity) along the director direction, while those below get a negative kick.
/// Invented by Kozhukhov and Shendruk (https://www.science.org/doi/full/10.1126/sciadv.abo5788). 
/// @param CL An MPCD cell (including the linked list of particles in each cell). 
/// @param SP The species-wide information about MPCD particles.
/// @param KBT The thermal energy. 
/// @param RELAX The temperature relaxation time scale.
/// @param CLQ The geometric centre of `CL`, the MPCD cell.
/// @param outP Flag whether or not to output the pressure.
/// @see andersenMPC()
///
void dipoleAndersenMPC( cell *CL,spec *SP,double KBT,double RELAX,double *CLQ,int outP ) {
	int i,j,id;
	double MASS,M,ACT,pmOne;
	double RV[CL->POP][DIM];				//Random velocities
	double RS[DIM];							//Sum of random velocities
	double DV[CL->POP][DIM];				//Damping velocities
	double AV[CL->POP][DIM];				//Active velocities
	double AS[DIM];							//Sum of active velocities
	double dp[CL->POP][DIM],relQ[CL->POP][DIM];	//For pressure
	particleMPC *tmpc;						//Temporary particleMPC
	particleMD *tmd;						//Temporary particleMD
	bc PLANE;								//The plane that cuts the cell in half
	double W;								//The particle's W for passing the plane
	double speed;							//Speed of particle

	//Define the plane normal to the centre of mass velocity at the centre of mass position
	for( i=0;i<4;i++ ) PLANE.P[i]=1;
	PLANE.INV=0;
	PLANE.ABS=0;
	PLANE.R=0.0;
	PLANE.ROTSYMM[0]=4.0;
	PLANE.ROTSYMM[1]=4.0;
	//Normal
	for( i=0; i<DIM; i++ ) PLANE.A[i] = CL->DIR[i];
	//Position
	for( i=0; i<DIM; i++ ) PLANE.Q[i] = CL->CM[i];
	MASS = CL->MASS;
	// Zero arrays
	for( i=0; i<DIM; i++ ) {
		RS[i] = 0.;
		AS[i] = 0.;
		for( j=0;j<CL->POP;j++ ) {
			RV[j][i] = 0.;
			DV[j][i] = 0.;
			AV[j][i] = 0.;
		}
	}
	//Calculate total activity of cell
	tmpc = CL->pp;
	ACT=0.;
	while( tmpc!=NULL ) {
		id = tmpc->SPID;
		ACT += (SP+id)->ACT;
		tmpc = tmpc->next;
	}

	/* ****************************************** */
	/* ******* Generate random velocities ******* */
	/* ****************************************** */
	// MPC particles
	tmpc = CL->pp;
	i=0;
	while( tmpc!=NULL ) {
		id = tmpc->SPID;
		M = (SP+id)->MASS;
		//Pressure term
		if( outP ) calcPressureColl_preColl( relQ[i],dp[i],tmpc,CLQ );
		//Random perturbation
		for( j=0; j<DIM; j++ ) RV[i][j] = genrand_gaussMB( KBT,M );
		for( j=0; j<DIM; j++ ) RS[j] += M*RV[i][j];
		for( j=0; j<DIM; j++ ) DV[i][j] = ((SP+id)->DAMP)*(CL->VCM[j])/((double)CL->POP);
		if( fneq(ACT,0.0) ) {
			speed=length(tmpc->V,DIM);
			//Check which side of the plane
			W = calcW( PLANE,*tmpc );
			if( W<=0 ) pmOne=-1.;
			else pmOne=1.;
			for( j=0; j<DIM; j++ ) AV[i][j] = pmOne*RELAX*(ACT-speed)*PLANE.A[j];
			for( j=0; j<DIM; j++ ) AS[j] += M*AV[i][j];
		}
		tmpc = tmpc->next;
		i++;
	}
	// MD particles
	tmd = CL->MDpp;
	while( tmd!=NULL ) {
		M = tmd->mass;
		for( j=0; j<DIM; j++ ) RV[i][j] = genrand_gaussMB( KBT,M );
		for( j=0; j<DIM; j++ ) RS[j] += M*RV[i][j];
		tmd = tmd->nextSRD;
		i++;
	}
	// Normalize
	for( j=0; j<DIM; j++ ) RS[j] /= MASS;
	for( j=0; j<DIM; j++ ) AS[j] /= MASS;

	/* ****************************************** */
	/* *************** Collision **************** */
	/* ****************************************** */
	// MPC particles
	tmpc = CL->pp;
	i=0;
	while( tmpc!=NULL ) {
		id = tmpc->SPID;
		M = (SP+id)->MASS;
		for( j=0; j<DIM; j++ ) tmpc->V[j] = CL->VCM[j] + RV[i][j] - RS[j] - DV[i][j] + AV[i][j] - AS[j];
		//Pressure term
		if( outP ) calcPressureColl_postColl( relQ[i],dp[i],M,tmpc->V,CL );
		//Increment link in list
		tmpc = tmpc->next;
		i++;
	}
	//MD particles
	tmd = CL->MDpp;
	while( tmd!=NULL ) {
		tmd->vx = CL->VCM[0] + RV[i][0] - RS[0];
		if( DIM > _1D){
			tmd->vy = CL->VCM[1] + RV[i][1] - RS[1];
		}
		if( DIM > _2D){
			tmd->vz = CL->VCM[2] + RV[i][2] - RS[2];
		}
		//Increment link in list
		tmd = tmd->nextSRD;
		i++;
	}
}

/// 
/// @brief The MPC collision operation. 
///
/// The collision operator --- this is the heart of the MPCD algorithm. 
/// This function routes the code to the user-selected version of the collision operator. 
/// Each of these collision operators are cell-based coarse-grained multi-particle collision events. 
/// In each of these, the momentum of the MPCD particles in the cell are stochastically exchanged. 
/// The goals of each of these particle-based mesoscale simulation techniques is to incorporate thermal fluctuations and hydrodynamic interactions. 
/// @param CL An MPCD cell. 
/// @param SP The species-wide information about MPCD particles.
/// @param SS The species-wide information about swimmers.
/// @param KBT The thermal energy. 
/// @param RTECH The MPCD collision operator. See `definitions.h` for all options.
/// @param C Cosine of the rotation angle.
/// @param S Sine of the rotation angle.
/// @param FRICCO Friction coefficient for Langevin thermostat.
/// @param TimeStep The interval of each MPCD time step in MPCD time units.
/// @param MD_mode The MD coupling mode. Can be off (noMD), MD particles included in the MPCD collisions (MDinMPC), or MPCD particles included in MD pair interactions (MPCinMD).
/// @param LC Flags whether or not the nematic liquid crystal is turned on.
/// @param RELAX The temperature relaxation time scale.
/// @param CLQ The geometric centre of `CL`, the MPCD cell.
/// @param outP Flag whether or not to output the pressure.
///
/// For momentum collision operators:
/// @see andersenMPC()
/// @see andersenROT()
/// @see langevinMPC()
/// @see langevinROT()
/// @see stochrotMPC()
/// 
/// For liquid crystal collision operator:
/// @see andersenROT_LC()
///
/// For active matter collision operators:
/// @see dipoleAndersenROT_LC()
/// @see dipoleAndersenMPC()
/// @see vicsek()
/// @see chate()
/// @see activeSRD()
/// @see vicsekAndersenMPC()
/// @see vicsekLangevinMPC()
/// @see chateAndersenMPC()
/// @see chateLangevinMPC()
///
void MPCcollision(cell *CL, spec *SP, specSwimmer SS, double KBT, int RTECH, double C, double S, double FRICCO, double TimeStep, int MD_mode, int LC, double RELAX, double *CLQ, int outP ) {
	particleMD *tmd;	//Temporary particleMD

	tmd = CL->MDpp;
	if (MD_mode == MPCinMD) CL->MDpp=NULL;
	if( outP ) zeroPressureColl( CL );
	//Liquid Crystal
	if( LC!=ISOF ) {
		if( RTECH==DIPOLE_DIR_SUM || RTECH==DIPOLE_DIR_AV || RTECH==DIPOLE_DIR_SIG || RTECH==DIPOLE_DIR_SIG_SUM ) dipoleAndersenROT_LC( CL,SP,SS,KBT,RELAX,TimeStep,RTECH,CLQ,outP );
		else andersenROT_LC( CL,SP,SS,KBT,TimeStep,CLQ,outP );
	}
	//Newtonian Fluid
	else{
		//Passive MPC operators
		if( RTECH == MPCAT ) andersenMPC( CL,SP,SS,KBT,CLQ,outP );
		else if( RTECH == RAT ) andersenROT( CL,SP,SS,KBT,CLQ,outP );
		else if( RTECH == LANG ) langevinMPC( CL,SP,SS,KBT,FRICCO,TimeStep,CLQ,outP );
		else if( RTECH == RLANG ) langevinROT( CL,SP,SS,KBT,FRICCO,TimeStep,CLQ,outP );
		else if( RTECH==ORTHAXIS || RTECH==ARBAXIS ) stochrotMPC( CL,RTECH,C,S,CLQ,outP );
		//Active MPC operators
		else if( RTECH==VICSEK ) vicsek( CL,SP,CLQ,outP );
		else if( RTECH==CHATE ) chate( CL,SP,CLQ,outP );
		else if( RTECH==ACT_ARBAXIS || RTECH==ACT_ORTHAXIS ) activeSRD( CL,SP,RTECH,C,S,CLQ,outP );
		else if( RTECH==VICSEK_MPCAT ) vicsekAndersenMPC( CL,SP,KBT,RELAX,CLQ,outP );
		else if( RTECH==VICSEK_LANG) vicsekLangevinMPC( CL,SP,KBT,FRICCO,TimeStep,RELAX,CLQ,outP );
		else if( RTECH==CHATE_MPCAT ) chateAndersenMPC( CL,SP,KBT,RELAX,CLQ,outP );
		else if( RTECH==CHATE_LANG ) chateLangevinMPC( CL,SP,KBT,FRICCO,TimeStep,RELAX,CLQ,outP );
		else if( RTECH==DIPOLE_VCM ) dipoleAndersenMPC( CL,SP,KBT,RELAX,CLQ,outP );
		else{
			printf( "Error: Collision technique unacceptable.\n" );
			exit( 1 );
		}
	}
	if( outP ) normPressureColl( CL,TimeStep );
	if (MD_mode == MPCinMD) CL->MDpp=tmd;
}

/// 
/// @brief Collision operation for phase separating fluids. 
/// 
/// This is a supplement to the collision operator that allows multiphase fluids to phase separate. 
/// In theory, it works equally well with any of the collision operators in MPCcollision(). 
/// Modifies the collision operation to allow fluid particles of different species to interact. 
/// Phase separation requires estimating the gradient between species --- there are different ways to approximate the gradient. 
/// This function routes the code to the user-selected version of the multiphase collision operator for different ways of estimating the gradients. 
/// @param CL An MPCD cell. 
/// @param SP The species-wide information about MPCD particles.
/// @param SS The species-wide information about swimmers.
/// @param multiphaseMode The interactions between different species that allows phase segregation. 
/// @param KBT The thermal energy. 
/// @param MD_mode The MD coupling mode. Can be off (noMD), MD particles included in the MPCD collisions (MDinMPC), or MPCD particles included in MD pair interactions (MPCinMD).
/// @param CLQ The geometric centre of `CL`, the MPCD cell.
/// @param outP Flag whether or not to output the pressure.
///
void multiphaseColl(cell *CL, spec *SP, specSwimmer SS, int multiphaseMode, double KBT, int MD_mode, double *CLQ, int outP ) {
/*
    THIS IS WHERE KIRA SHOULD ADD HER NEW BIT!!!
*/
	if( multiphaseMode==MPHSURF ) {
		printf("Error: Multiphase interaction not yet implemented.\nTo be implemented by Kira");
		exit( 1 );
	}
	else if( multiphaseMode==MPHPOINT ) multiphaseCollPoint(CL, SP, SS, KBT, MD_mode, CLQ, outP );
	else {
		printf( "Error: Multiphase interaction  technique unacceptable.\n" );
		exit( 1 );
	}
}

/// 
/// @brief Collision operation for phase separating fluids that estimates gradients by a point-particle method. 
/// 
/// This routine supplements to the collision operator to allow different species of particles to interact. 
/// This can produce multiphase fluids to phase separate. 
/// This version uses the point particle positions to approximate the phase gradient. 
/// @param CL An MPCD cell (including the linked list of particles in each cell). 
/// @param SP The species-wide information about MPCD particles.
/// @param SS The species-wide information about swimmers.
/// @param KBT The thermal energy. 
/// @param MD_mode The MD coupling mode. Can be off (noMD), MD particles included in the MPCD collisions (MDinMPC), or MPCD particles included in MD pair interactions (MPCinMD).
/// @param CLQ The geometric centre of `CL`, the MPCD cell.
/// @param outP Flag whether or not to output the pressure.
///
void multiphaseCollPoint(cell *CL, spec *SP, specSwimmer SS, double KBT, int MD_mode, double *CLQ, int outP ) {
	int i,j,k,id;
	int mixedCell=0;
	double N,NSP[NSPECI];			//Number of each type
	particleMPC *tmpc;              //Temporary particleMPC
	particleMD *tmd;                //Temporary particleMD
	smono *tsm;                     //Temporary swimmer monomer
	double relQ[DIM];               //Relative position
	double VMUtot[DIM];             //Velocity due to chemical potential
	double gradSP[NSPECI][DIM];     //Directional gradient of each species
	double thisGrad;				//A temporary gradient contribution component
	double VMU[CL->POP][DIM];       //Grad. chemical potential  velocity of type A (B is negative this)

	// Zero arrays
	for( i=0; i<DIM; i++ ) {
		relQ[i] = 0.0;
		for( j=0; j<NSPECI; j++ ) gradSP[j][i] = 0.0;
	}
	for( i=0;i<CL->POP;i++ ) for( j=0;j<DIM;j++ ) {
		VMU[i][j] = 0.0;
		VMUtot[j] = 0.;
	}
	for( j=0; j<NSPECI; j++ ) NSP[j]=0.0;

	//Calculate the number of each type
	//MPCD particles
	tmpc = CL->pp;
	while( tmpc!=NULL ) {
		id = tmpc->SPID;
		NSP[id] += 1.0;
		//Increment link in list
		tmpc = tmpc->next;
	}
	//Swimmer monomers
	tsm = CL->sp;
	while( tsm!=NULL ) {
		if( tsm->HorM ) id = SS.MSPid;
		else id = SS.HSPid;
		NSP[id] += 1.0;
		//Increment link in list
		tsm = tsm->next;
	}
	//MD particles --- ALWAYS type 0
	id=0;
	tmd = CL->MDpp;
	while( tmd!=NULL ) {
		NSP[id] += 1.0;
		//Increment link in list
		tmd = tmd->nextSRD;
	}
	N=0.0;
	for( j=0; j<NSPECI; j++ ) N += NSP[j];

	//Generate separation velocities
	mixedCell=0;
	for( j=0; j<NSPECI; j++ ) if( NSP[j]>0.0 ) mixedCell+=1;
	if( mixedCell>1 ) {
		// Calculate the gradient of the different species
		// MPC particles
		tmpc = CL->pp;
		while( tmpc!=NULL ) {
			id = tmpc->SPID;
			//Particle-based gradient of this species
			for( j=0; j<DIM; j++ ) relQ[j] = tmpc->Q[j] - CLQ[j];
			for( j=0; j<DIM; j++ ) {
				thisGrad = 8.0*relQ[j]*relQ[j]-3.0;
				for( k=0; k<DIM; k++ ) thisGrad += 6.0*relQ[k]*relQ[k];
				thisGrad *= 30.0*relQ[j];
				gradSP[id][j] += thisGrad;
			}
			//Increment link in list
			tmpc = tmpc->next;
		}
		//Swimmer monomers
		tsm = CL->sp;
		while( tsm!=NULL ) {
			if( tsm->HorM ) id = SS.MSPid;
			else id = SS.HSPid;
			//Particle-based gradient of this species
			for( j=0; j<DIM; j++ ) relQ[j] = tsm->Q[j] - CLQ[j];
			for( j=0; j<DIM; j++ ) {
				thisGrad = 8.0*relQ[j]*relQ[j]-3.0;
				for( k=0; k<DIM; k++ ) thisGrad += 6.0*relQ[k]*relQ[k];
				thisGrad *= 30.0*relQ[j];
				gradSP[id][j] += thisGrad;
			}
			//Increment link in list
			tsm = tsm->next;
		}
		//MD particles --- ALWAYS type 0
		id=0;
		tmd = CL->MDpp;
		while( tmd!=NULL ) {
			if( DIM>=_1D ) relQ[0] = tmd->rx - CLQ[0];
			if( DIM>=_2D ) relQ[1] = tmd->ry - CLQ[1];
			if( DIM>=_3D ) relQ[2] = tmd->rz - CLQ[2];
			for( j=0; j<DIM; j++ ) {
				thisGrad = 8.0*relQ[j]*relQ[j]-3.0;
				for( k=0; k<DIM; k++ ) thisGrad += 6.0*relQ[k]*relQ[k];
				thisGrad *= 30.0*relQ[j];
				gradSP[id][j] += thisGrad;
			}
			//Increment link in list
			tmd = tmd->nextSRD;
		}

		//Calculate the velocities due to the cell's chemical potential
		i=0;
		//MPCD particles
		tmpc = CL->pp;
		while( tmpc!=NULL ) {
			id = tmpc->SPID;
			for( j=0; j<DIM; j++ ) {
				for( k=0; k<NSPECI; k++ ) VMU[i][j] += gradSP[k][j]*((SP+id)->M[k]) / NSP[id];
				VMUtot[j] += VMU[i][j];
			}
			//Increment link in list
			tmpc = tmpc->next;
			i++;
		}
		//Swimmer monomers
		tsm = CL->sp;
		while( tsm!=NULL ) {
			if( tsm->HorM ) id = SS.MSPid;
			else id = SS.HSPid;
			for( j=0; j<DIM; j++ ) {
				for( k=0; k<NSPECI; k++ ) VMU[i][j] += gradSP[k][j]*((SP+id)->M[k]) / NSP[id];
				VMUtot[j] += VMU[i][j];
			}
			//Increment link in list
			tsm = tsm->next;
			i++;
		}
	}
	//MD particles --- ALWAYS type 0
	id=0;
	tmd = CL->MDpp;
	while( tmd!=NULL ) {
		for( j=0; j<DIM; j++ ) {
			for( k=0; k<NSPECI; k++ ) VMU[i][j] += gradSP[k][j]*((SP+id)->M[k]) / NSP[id];
			VMUtot[j] += VMU[i][j];
		}
		//Increment link in list
		tmd = tmd->nextSRD;
		i++;
	}
	//Turn sums into averages
	for( j=0; j<DIM; j++ ) VMUtot[j] /= N;

	/* ****************************************** */
	/* *************** Collision **************** */
	/* ****************************************** */
	i=0;
	// MPC particles
	tmpc = CL->pp;
	while( tmpc!=NULL ) {
		id = tmpc->SPID;
		for( j=0; j<DIM; j++ ) tmpc->V[j] += VMU[i][j] - VMUtot[j];
		//Increment link in list
		tmpc = tmpc->next;
		i++;
	}
	// Swimmer monomers
	tsm = CL->sp;
	while( tsm!=NULL ) {
		if( tsm->HorM ) id = SS.MSPid;
		else id = SS.HSPid;
		for( j=0; j<DIM; j++ ) tsm->V[j] += VMU[i][j] - VMUtot[j];
		//Increment link in list
		tsm = tsm->next;
		i++;
	}
	//MD particles --- ALWAYS type 0
	id=0;
	tmd = CL->MDpp;
	while( tmd!=NULL ) {
		if( DIM>=_1D ) {
			j=0;
			tmd->vx += VMU[i][j] - VMUtot[j];
		}
		if( DIM>=_2D ) {
			j=1;
			tmd->vy += VMU[i][j] - VMUtot[j];
		}
		if( DIM>=_3D ) {
			j=2;
			tmd->vz += VMU[i][j] - VMUtot[j];
		}
		//Increment link in list
		tmd = tmd->nextSRD;
		i++;
	}
}

/// 
/// @brief Applies a correction to the collision operation to give the fluid a non-ideal equation of state. 
///
/// @warning Experimental. Proposed by Shendruk but not published or fully characterized yet. 
///
/// This is a supplement to the collision operator that attempts to give the fluid a non-ideal equation of state. 
/// The goal is to make the MPCD fluid less compressible; however, these are all currently experimental. 
/// In theory, it works equally well with any of the collision operators in MPCcollision(). 
/// Inspired by J. Chem. Phys. 154, 024105 (2021); https://doi.org/10.1063/5.0037934
/// This subroutine simply calls the possible ways to make the fluid non-ideal. 
/// @see incompSwap()
/// @see incompAddVirial()
/// @see incompSubtractDivergence()
/// @param CL An MPCD cell. 
/// @param SP The species-wide information about MPCD particles.
/// @param SS The species-wide information about swimmers.
/// @param INCOMPmode The technique for making the fluid a non-ideal gas.
/// @param MD_mode The MD coupling mode. Can be off (noMD), MD particles included in the MPCD collisions (MDinMPC), or MPCD particles included in MD pair interactions (MPCinMD).
/// @param CLQ The geometric centre of `CL`, the MPCD cell.
/// @param outP Flag whether or not to output the pressure.
///
void incompColl(cell *CL, spec *SP, specSwimmer SS, int INCOMPmode, int MD_mode, double *CLQ, int outP ) {
	if(INCOMPmode==INCOMPSWAP) incompSwap( CL,SP,SS );
	else if(INCOMPmode==INCOMPVIRIAL) incompAddVirial( CL,7.0,49.0,343.0,SP,SS );
	else if(INCOMPmode==INCOMPSUB) incompSubtractDivergence( CL,SP,SS );
	else{
		printf("Incompressibility correction not yet implemented.\n");
		exit( 1 );
	}
}

/// 
/// @brief MPCD operation to attempt to produce a non-ideal equation of state. 
///
/// @warning Experimental. Proposed by Shendruk but not published or fully characterized yet. 
///
/// An attempt to give the fluid a non-ideal equation of state to make the MPCD fluid less compressible. 
/// MPC collision that applies a correction radial kicks to try to keep the density constant. 
/// If this works, virialCoB,C,D will need to be an input parameters. 
/// @param CL An MPCD cell (including the linked list of particles in each cell). 
/// @param virialCoB First virial coefficient
/// @param virialCoC Second virial coefficient
/// @param virialCoD Third virial coefficient
/// @param SP The species-wide information about MPCD particles.
/// @param SS The species-wide information about swimmers.
/// @note Experimental and might not work.
///
void incompAddVirial( cell *CL,double virialCoB, double virialCoC, double virialCoD, spec *SP,specSwimmer SS ) {
	int i,j,id;
	double MASS;
	double RV[CL->POP][DIM];	//Radial velocities
	double RS[DIM];				//Sum of radial velocities
	double relQ[DIM];			//Relative position
	double dN;					//Relative density
	particleMPC *tmpc;			//Temporary particleMPC
	particleMD *tmd;			//Temporary particleMD
	smono *tsm;					//Temporary swimmer monomer

	if( CL->POP > GnDNST ) {
		dN = ((double)CL->POP)/GnDNST - 1.0;
		// Zero arrays
		for( i=0;i<CL->POP;i++ ) for( j=0;j<DIM;j++ ) RV[i][j] = 0.0;
		for( i=0; i<DIM; i++ ) {
			RS[i] = 0.;
			relQ[i] = 0.;
		}

		/* ****************************************** */
		/* ******* Generate radial velocities ******* */
		/* ****************************************** */
		i=0;
		// MPC particles
		tmpc = CL->pp;
		while( tmpc!=NULL ) {
			id = tmpc->SPID;
			MASS = (SP+id)->MASS;
			for( j=0; j<DIM; j++ ) relQ[j] = tmpc->Q[j] - CL->CM[j];
			norm( relQ,DIM );
			for( j=0; j<DIM; j++ ) RV[i][j] = (virialCoB+virialCoC*dN+virialCoD*dN*dN)*dN*relQ[j];
			for( j=0; j<DIM; j++ ) RS[j] += MASS*RV[i][j];
			tmpc = tmpc->next;
			i++;
		}
		// MD particles
		tmd = CL->MDpp;
		while( tmd!=NULL ) {
			MASS = tmd->mass;
			relQ[0] = tmd->rx - CL->CM[0];
			if(DIM>=_2D) relQ[1] = tmd->ry - CL->CM[1];
			if(DIM>=_3D) relQ[2] = tmd->rz - CL->CM[2];
			norm( relQ,DIM );
			for( j=0; j<DIM; j++ ) RV[i][j] = (virialCoB+virialCoC*dN+virialCoD*dN*dN)*dN*relQ[j];
			for( j=0; j<DIM; j++ ) RS[j] += MASS*RV[i][j];
			tmd = tmd->nextSRD;
			i++;
		}
		// Swimmer monomers
		tsm = CL->sp;
		while( tsm!=NULL ) {
			if( tsm->HorM ) MASS = (double) SS.middM;
			else MASS = (double) SS.headM;
			for( j=0; j<DIM; j++ ) relQ[j] = tsm->Q[j] - CL->CM[j];
			norm( relQ,DIM );
			for( j=0; j<DIM; j++ ) RV[i][j] = (virialCoB+virialCoC*dN+virialCoD*dN*dN)*dN*relQ[j];
			for( j=0; j<DIM; j++ ) RS[j] += MASS*RV[i][j];
			tsm = tsm->next;
			i++;
		}
		// Normalize
		for( j=0; j<DIM; j++ ) RS[j] /= CL->MASS;

		/* ****************************************** */
		/* *************** Collision **************** */
		/* ****************************************** */
		i=0;
		// MPC particles
		tmpc = CL->pp;
		while( tmpc!=NULL ) {
			id = tmpc->SPID;
			for( j=0; j<DIM; j++ ) tmpc->V[j] += RV[i][j] - RS[j];
			//Increment link in list
			tmpc = tmpc->next;
			i++;
		}
		//MD particles
		tmd = CL->MDpp;
		while( tmd!=NULL ) {
			tmd->vx += RV[i][0] - RS[0];
			if(DIM>=_2D) tmd->vy += RV[i][1] - RS[1];
			if(DIM>=_3D) tmd->vz += RV[i][2] - RS[2];
			//Increment link in list
			tmd = tmd->nextSRD;
			i++;
		}
		// Swimmer monomers
		tsm = CL->sp;
		while( tsm!=NULL ) {
			for( j=0; j<DIM; j++ ) tsm->V[j] += RV[i][j] - RS[j];
			//Increment link in list
			tsm = tsm->next;
			i++;
		}
	}
}

/// 
/// @brief MPCD operation to attempt to constrain the divergence of momentum density to be zero. 
///
/// @warning Experimental. Proposed by Shendruk but not published or fully characterized yet. 
///
/// An attempt to give the fluid a non-ideal equation of state to make the MPCD fluid less compressible. 
/// MPC collision that applies a correction to constrain the divergence of momentum density to be zero by swapping velocities. 
/// Only swaps velocities of MPCD fluid particles (not swimmers or MD particles).
/// @param CL An MPCD cell (including the linked list of particles in each cell). 
/// @param SP The species-wide information about MPCD particles.
/// @param SS The species-wide information about swimmers.
/// @note Experimental and might not work.
///
void incompSwap( cell *CL,spec *SP,specSwimmer SS ) {
	int i,j,k,id,id2,POP;
	double M,M2,MCM,MC;						//Masses
	double relQ[CL->POP][_3D];				//Relative position
	double DIV[CL->POP];					//Divergence contribution from each particle
	double V[_3D],V2[_3D];					//Velocities
	double VCM[_3D],RCM[_3D],PCM[_3D];
	double DIV1,DIV2,DIVC,DIVC_swapped;
	//Checks for debugging
	double DIVCI,DIVCF;						//Divergence before and after (I=initial; F=final). 
	// double EI,EF;							//Energy before and after (I=initial; F=final). 
	double diffP,PMI,PMF;					//Momenumt magnitude before and after (I=initial; F=final). 
	double PI[_3D],PF[_3D];					//Momentum before and after
	//Linked list variables
	particleMPC *tmpc,*tmpc2;				//Temporary particleMPC
	particleMD *tmd;						//Temporary particleMD
	smono *tsm;								//Temporary swimmer monomer

	for( i=0;i<CL->POP;i++ ) {
		DIV[i] = 0.;
		for( j=0;j<_3D;j++ ) relQ[i][j] = 0.;
	}
	for( j=0;j<_3D;j++ ) {
		V[j]=0.;
		V2[j]=0.;
		VCM[j] = 0.;
		RCM[j] = 0.;
		PCM[j] = 0.;
	}
	DIVC = 0.;
	POP=CL->POP;
	MCM=CL->MASS;
	MC=MCM/(double)POP;
	for( j=0;j<DIM;j++ ) {
		VCM[j] = CL->VCM[j];
		RCM[j] = CL->CM[j];
		PCM[j] = MC*VCM[j];
	}
	// Debug checks
	#ifdef DBG
		if ( DBUG == DBGINCOMP ) {
			for( j=0;j<_3D;j++ ) {
				PI[j] = 0.;
				PF[j] = 0.;
			}
			DIVCI=0.;
			DIVCF=0.;
			// EI=0.;
			// EF=0.;
		}
	#endif

	/* ****************************************** */
	/* ********** Calculate divergence ********** */
	/* ****************************************** */
	i=0;
	//MPC particles
	tmpc = CL->pp;
	while( tmpc!=NULL ) {
		id = tmpc->SPID;
		M = (SP+id)->MASS;
		for( j=0; j<DIM; j++ ) V[j] = tmpc->V[j];
		//Position relative to centre of mass
		for( j=0; j<DIM; j++ ) relQ[i][j] = tmpc->Q[j] - RCM[j];
		//Initial divergence
		for( j=0; j<DIM; j++ ) DIV[i] += (PCM[j]-M*V[j])/relQ[i][j];
		tmpc = tmpc->next;
		i++;
	}
	//MD particles
	tmd = CL->MDpp;
	while( tmd!=NULL ) {
		M = tmd->mass;
		V[0] = tmd->vx;
		V[1] = tmd->vy;
		V[2] = tmd->vz;
		//Position relative to centre of mass
		relQ[i][0] = tmd->rx - CL->CM[0];
		relQ[i][1] = tmd->ry - CL->CM[1];
		relQ[i][2] = tmd->rz - CL->CM[2];
		//Initial divergence
		for( j=0; j<DIM; j++ ) DIV[i] += (PCM[j]-M*V[j])/relQ[i][j];
		tmd = tmd->nextSRD;
		i++;
	}
	//Swimmer particles
	tsm = CL->sp;
	while( tsm!=NULL ) {
		if( tsm->HorM ) M = (double) SS.middM;
		else M = (double) SS.headM;
		for( j=0; j<DIM; j++ ) V[j] = tsm->V[j];
		//Position relative to centre of mass
		for( j=0; j<DIM; j++ ) relQ[i][j] = tsm->Q[j] - CL->CM[j];
		//Initial divergence
		for( j=0; j<DIM; j++ ) DIV[i] += (PCM[j]-M*V[j])/relQ[i][j];
		tsm = tsm->next;
		i++;
	}
	for( i=0; i<POP; i++ ) DIVC += DIV[i];
	DIVC/=(double)POP;
	// Debug checks
	#ifdef DBG
		if ( DBUG == DBGINCOMP ) {
			i=0;
			//MPC particles
			tmpc = CL->pp;
			while( tmpc!=NULL ) {
				id = tmpc->SPID;
				M = (SP+id)->MASS;
				for( j=0; j<DIM; j++ ) V[j] = tmpc->V[j];
				// for( j=0; j<DIM; j++ ) EI += M*V[j]*V[j];
				for( j=0; j<DIM; j++ ) PI[j] += M*V[j];
				tmpc = tmpc->next;
				i++;
			}
			//MD particles
			tmd = CL->MDpp;
			while( tmd!=NULL ) {
				M = tmd->mass;
				V[0] = tmd->vx;
				V[1] = tmd->vy;
				V[2] = tmd->vz;
				// for( j=0; j<DIM; j++ ) EI += M*V[j]*V[j];
				for( j=0; j<DIM; j++ ) PI[j] += M*V[j];
				tmd = tmd->nextSRD;
				i++;
			}
			//Swimmer particles
			tsm = CL->sp;
			while( tsm!=NULL ) {
				if( tsm->HorM ) M = (double) SS.middM;
				else M = (double) SS.headM;
				for( j=0; j<DIM; j++ ) V[j] = tsm->V[j];
				// for( j=0; j<DIM; j++ ) EI += M*tsm->V[j]*V[j];
				for( j=0; j<DIM; j++ ) PI[j] += M*V[j];
				tsm = tsm->next;
				i++;
			}
			// EI*=0.5;
			DIVCI=DIVC;
		}
	#endif

	/* ****************************************** */
	/* ********** Swapping Collision ************ */
	/* ****************************************** */
	i=0;
	//MPC particles --- only allow swapping of MPC particles (not MD or swimmers)
	tmpc = CL->pp;
	while( tmpc!=NULL ) {
		id = tmpc->SPID;
		M = (SP+id)->MASS;

		//Loop through all other particles
		tmpc2 = tmpc->next;
		k=i+1;
		while( tmpc2!=NULL ) {
			id2 = tmpc2->SPID;
			M2 = (SP+id2)->MASS;
			//Swap velocities
			for( j=0; j<DIM; j++ ) {
				V[j] = (M2/M)*tmpc2->V[j];
				V2[j] = (M/M2)*tmpc->V[j];
			}
			// Calculate divergence if swap is kept
			DIV1=0.0;
			DIV2=0.0;
			for( j=0; j<DIM; j++ ) {
				DIV1 += (PCM[j]-M*V[j])/relQ[i][j];
				DIV2 += (PCM[j]-M2*V2[j])/relQ[k][j];
			}
			DIVC_swapped = DIVC - DIV[i] + DIV1 - DIV[k] + DIV2;
			//If divergence is smaller then save the swap. Otherwise, don't do anything
			if( fabs(DIVC_swapped)<fabs(DIVC) ) {
				DIVC=DIVC_swapped;
				DIV[i] = DIV1;
				DIV[k] = DIV2;
				for( j=0; j<DIM; j++ ) {
					tmpc->V[j] = V[j];
					tmpc2->V[j] = V2[j];
				}
			}
			//Increment link in list
			tmpc2 = tmpc2->next;
			k++;
		}
		//Increment link in list
		tmpc = tmpc->next;
		i++;
	}
	// Debug checks
	#ifdef DBG
		if ( DBUG == DBGINCOMP ) {
			i=0;
			//MPC particles
			tmpc = CL->pp;
			while( tmpc!=NULL ) {
				id = tmpc->SPID;
				M = (SP+id)->MASS;
				for( j=0; j<DIM; j++ ) V[j] = tmpc->V[j];
				// for( j=0; j<DIM; j++ ) EF += M*V[j]*V[j];
				for( j=0; j<DIM; j++ ) PF[j] += M*V[j];
				tmpc = tmpc->next;
				i++;
			}
			//MD particles
			tmd = CL->MDpp;
			while( tmd!=NULL ) {
				M = tmd->mass;
				V[0] = tmd->vx;
				V[1] = tmd->vy;
				V[2] = tmd->vz;
				// for( j=0; j<DIM; j++ ) EF += M*V[j]*V[j];
				for( j=0; j<DIM; j++ ) PF[j] += M*V[j];
				tmd = tmd->nextSRD;
				i++;
			}
			//Swimmer particles
			tsm = CL->sp;
			while( tsm!=NULL ) {
				if( tsm->HorM ) M = (double) SS.middM;
				else M = (double) SS.headM;
				for( j=0; j<DIM; j++ ) V[j] = tsm->V[j];
				// for( j=0; j<DIM; j++ ) EF += M*tsm->V[j]*V[j];
				for( j=0; j<DIM; j++ ) PF[j] += M*V[j];
				tsm = tsm->next;
				i++;
			}
			// EF*=0.5;
			DIVCF=DIVC;
			diffP=0.;
			PMI=0.;
			PMF=0.;
			for( j=0; j<DIM; j++ ) {
				PMI += PI[j]*PI[j];
				PMF += PF[j]*PF[j];
				diffP += (PI[j]-PF[j])*(PI[j]-PF[j]);
			}
			PMI=sqrt(PMI);
			PMF=sqrt(PMF);
			diffP=sqrt(diffP);
			// printf( "Divergence:\tInit=%e\tFin=%e\tperFracDiff=%f%%\n",DIVCI,DIVCF,100.0*DIVCF/DIVCI );
			// printf( "Energy:\t\tInit=%e\tFin=%e\tperDiff=%e%%\n",EI,EF,100.0*(EI-EF)/EI );
			// printf( "Momentum:\tInit=%e\tFin=%e\tperDiff=%e%%\n\n",PMI,PMF,100.0*diffP/PMI );
			printf( "%e\n",100.0*DIVCF/DIVCI );
		}
	#endif
}

/// 
/// @brief MPCD operation to attempt to constrain the divergence of momentum density to be zero. 
///
/// @warning Experimental. Proposed by Shendruk but not published or fully characterized yet. 
///
/// An attempt to give the fluid a non-ideal equation of state to make the MPCD fluid less compressible. 
/// Collision that applies a correction to constrain the divergence to be zero. 
/// @param CL An MPCD cell (including the linked list of particles in each cell). 
/// @param SP The species-wide information about MPCD particles.
/// @param SS The species-wide information about swimmers.
/// @note Experimental and might not work.
///
void incompSubtractDivergence( cell *CL,spec *SP,specSwimmer SS ) {
	int i,j,id;
	double MASS,MCM,MC;
	double CV[CL->POP][_3D];				//Correction velocities
	double CS[_3D];							//Sum of correction velocities
	double relQ[CL->POP][_3D];				//Relative position
	double V[_3D],L[_3D];					//Velocity and angular momentum
	double VCM[_3D],RCM[_3D],PCM[_3D];
	//Checks for debugging
	double DIV0,DIV1,E0,E1,diffP,PM0,PM1;	//Divergence before and after. Energy before and after
	double DIVAV[_3D];						//The average of the divergence term applied to each particle in the cells
	double P0[_3D],P1[_3D];					//Momentum before and after
	//Linked list variables
	particleMPC *tmpc;						//Temporary particleMPC
	particleMD *tmd;						//Temporary particleMD
	smono *tsm;								//Temporary swimmer monomer

	for( i=0;i<CL->POP;i++ ) for( j=0;j<_3D;j++ ) {
		CV[i][j] = 0.;
		relQ[i][j] = 0.;
	}
	for( j=0;j<_3D;j++ ) {
		CS[j]=0.;
		L[j]=0.;
		V[j]=0.;
		VCM[j] = 0.;
		PCM[j] = 0.;
		RCM[j] = 0.;
		DIVAV[j] = 0.;
		P0[j] = 0.;
		P1[j] = 0.;
	}
	DIV0=0.;
	DIV1=0.;
	E0=0.;
	E1=0.;
	MCM=CL->MASS;
	MC=MCM/(double)CL->POP;
	for( j=0;j<DIM;j++ ) {
		VCM[j] = CL->VCM[j];
		RCM[j] = CL->CM[j];
		PCM[j] = MC*VCM[j];
	}

	/* ****************************************** */
	/* ********** Calculate divergence ********** */
	/* ****************************************** */
	i=0;
	//MPC particles
	tmpc = CL->pp;
	while( tmpc!=NULL ) {
		id = tmpc->SPID;
		MASS = (SP+id)->MASS;
		for( j=0; j<DIM; j++ ) V[j] = tmpc->V[j];
		//Position relative to centre of mass
		for( j=0; j<DIM; j++ ) relQ[i][j] = tmpc->Q[j] - RCM[j];
		//Initial divergence, kinetic energy and momentum
		for( j=0; j<DIM; j++ ) DIV0 += (PCM[j]-MASS*V[j])/relQ[i][j];
		for( j=0; j<DIM; j++ ) DIVAV[j] += relQ[i][j]/MASS;
		for( j=0; j<DIM; j++ ) E0 += MASS*V[j]*V[j];
		for( j=0; j<DIM; j++ ) P0[j] += MASS*V[j];
		tmpc = tmpc->next;
		i++;
	}
	//MD particles
	tmd = CL->MDpp;
	while( tmd!=NULL ) {
		MASS = tmd->mass;
		V[0] = tmd->vx;
		V[1] = tmd->vy;
		V[2] = tmd->vz;
		//Position relative to centre of mass
		relQ[i][0] = tmd->rx - CL->CM[0];
		relQ[i][1] = tmd->ry - CL->CM[1];
		relQ[i][2] = tmd->rz - CL->CM[2];
		//Initial divergence, kinetic energy and momentum
		for( j=0; j<DIM; j++ ) DIV0 += (PCM[j]-MASS*V[j])/relQ[i][j];
		for( j=0; j<DIM; j++ ) DIVAV[j] += relQ[i][j]/MASS;
		for( j=0; j<DIM; j++ ) E0 += MASS*V[j]*V[j];
		for( j=0; j<DIM; j++ ) P0[j] += MASS*V[j];
		tmd = tmd->nextSRD;
		i++;
	}
	//Swimmer particles
	tsm = CL->sp;
	while( tsm!=NULL ) {
		if( tsm->HorM ) MASS = (double) SS.middM;
		else MASS = (double) SS.headM;
		for( j=0; j<DIM; j++ ) V[j] = tsm->V[j];
		//Position relative to centre of mass
		for( j=0; j<DIM; j++ ) relQ[i][j] = tsm->Q[j] - CL->CM[j];
		//Initial divergence, kinetic energy and momentum
		for( j=0; j<DIM; j++ ) DIV0 += (PCM[j]-MASS*V[j])/relQ[i][j];
		for( j=0; j<DIM; j++ ) DIVAV[j] += relQ[i][j]/MASS;
		for( j=0; j<DIM; j++ ) E0 += MASS*tsm->V[j]*V[j];
		for( j=0; j<DIM; j++ ) P0[j] += MASS*V[j];
		tsm = tsm->next;
		i++;
	}
	DIV0/=(double)(DIM*(CL->POP));			//DIM always appears in denominator so just include it
	for( j=0; j<DIM; j++ ) DIVAV[j]*=(DIV0/((double)CL->POP));
	E0*=0.5;
	/* ****************************************** */
	/* *************** Collision **************** */
	/* ****************************************** */
	i=0;
	//MPC particles
	tmpc = CL->pp;
	while( tmpc!=NULL ) {
		id = tmpc->SPID;
		MASS = (SP+id)->MASS;
		//Perfectly remove divergence but doesn't conserve momentum
		// for( j=0; j<DIM; j++ ) tmpc->V[j] = tmpc->V[j] + relQ[i][j]*DIV0/MASS;
		// //Perfectly conserver momentum at cost of removing divergence
		for( j=0; j<DIM; j++ ) tmpc->V[j] = tmpc->V[j] + relQ[i][j]*DIV0/MASS - DIVAV[j];
		//Switched sign
		// for( j=0; j<DIM; j++ ) tmpc->V[j] = tmpc->V[j] - relQ[i][j]*DIV0/MASS + DIVAV[j];
		//Just average
		// for( j=0; j<DIM; j++ ) tmpc->V[j] = tmpc->V[j] + DIVAV[j];
		//Final divergence, kinetic energy and momentum
		for( j=0; j<DIM; j++ ) V[j] = tmpc->V[j];
		for( j=0; j<DIM; j++ ) DIV1 += (PCM[j]-MASS*V[j])/relQ[i][j];
		for( j=0; j<DIM; j++ ) E1 += MASS*V[j]*V[j];
		for( j=0; j<DIM; j++ ) P1[j] += MASS*V[j];
		//Increment link in list
		tmpc = tmpc->next;
		i++;
	}
	//MD particles
	tmd = CL->MDpp;
	while( tmd!=NULL ) {
		MASS = tmd->mass;
		// tmd->vx = tmd->vx + relQ[i][0]*DIV0/MASS;
		// tmd->vy = tmd->vy + relQ[i][1]*DIV0/MASS;
		// tmd->vz = tmd->vz + relQ[i][2]*DIV0/MASS;
		tmd->vx = tmd->vx + relQ[i][0]*DIV0/MASS - DIVAV[0];
		tmd->vy = tmd->vy + relQ[i][1]*DIV0/MASS - DIVAV[1];
		tmd->vz = tmd->vz + relQ[i][2]*DIV0/MASS - DIVAV[2];
		// tmd->vx = tmd->vx - relQ[i][0]*DIV0/MASS + DIVAV[0];
		// tmd->vy = tmd->vy - relQ[i][1]*DIV0/MASS + DIVAV[1];
		// tmd->vz = tmd->vz - relQ[i][2]*DIV0/MASS + DIVAV[2];
		//Final divergence, kinetic energy and momentum
		V[0] = tmd->vx;
		V[1] = tmd->vy;
		V[2] = tmd->vz;
		for( j=0; j<DIM; j++ ) DIV1 += (PCM[j]-MASS*V[j])/relQ[i][j];
		for( j=0; j<DIM; j++ ) E1 += MASS*V[j]*V[j];
		for( j=0; j<DIM; j++ ) P1[j] += MASS*V[j];
		//Increment link in list
		tmd = tmd->nextSRD;
		i++;
	}
	//Swimmer particles
	tsm = CL->sp;
	while( tsm!=NULL ) {
		if( tsm->HorM ) MASS = (double) SS.middM;
		else MASS = (double) SS.headM;
		// for( j=0; j<DIM; j++ ) tsm->V[j] = tsm->V[j] + relQ[i][j]*DIV0/MASS;
		for( j=0; j<DIM; j++ ) tsm->V[j] = tsm->V[j] + relQ[i][j]*DIV0/MASS - DIVAV[j];
		// for( j=0; j<DIM; j++ ) tsm->V[j] = tsm->V[j] - relQ[i][j]*DIV0/MASS + DIVAV[j];
		//Final divergence, kinetic energy and momentum
		for( j=0; j<DIM; j++ ) V[j] = tsm->V[j];
		for( j=0; j<DIM; j++ ) DIV1 += (PCM[j]-MASS*V[j])/relQ[i][j];
		for( j=0; j<DIM; j++ ) E1 += MASS*V[j]*V[j];
		for( j=0; j<DIM; j++ ) P1[j] += MASS*V[j];
		//Increment link in list
		tsm = tsm->next;
		i++;
	}
	DIV1/=(double)CL->POP;
	E1*=0.5;
	#ifdef DBG
		if ( DBUG == DBGINCOMP ) {
			printf( "Divergence:\tInit=%e\tFin=%e\tperFracDiff=%e%%\n",DIV0,DIV1,100.0*DIV1/DIV0 );
			printf( "Energy:\t\tInit=%e\tFin=%e\tperDiff=%e%%\n",E0,E1,100.0*(E0-E1)/E0 );
			diffP=0.;
			PM0=0.;
			PM1=0.;
			for( j=0; j<DIM; j++ ) {
				PM0 += P0[j]*P0[j];
				PM1 += P1[j]*P1[j];
				diffP += (P0[j]-P1[j])*(P0[j]-P1[j]);
			}
			PM0=sqrt(PM0);
			PM1=sqrt(PM1);
			diffP=sqrt(diffP);
			printf( "Momentum:\tInit=%e\tFin=%e\tperDiff=%e%%\n\n",PM0,PM1,100.0*diffP/PM0 );
		}
	#endif
}

/// 
/// @brief This routine calculates the centre of mass velocity of a single cell.
/// 
/// It loops through the linked lists, to calculate the centre of mass velocity of a given cell. 
/// It includes MPCD, MD and swimmer particles. 
/// @param vcm The centre of mass velocity vector of cell `CL`. Velocity is returned through this variable.
/// @param CL An MPCD cell (including the linked list of particles in each cell). 
/// @param SP The species-wide information about MPCD particles.
/// @param specS The species-wide information about swimmers.
///
void localVCM( double vcm[_3D],cell CL,spec *SP,specSwimmer specS ) {
	int id,i;
	double summ = 0.0;
	double mass;
	particleMPC *tmpc;
	particleMD *tmd;
	smono *tsm;

	// Zero the centre of mass velocity
	for( i=0; i<DIM; i++ ) vcm[i] = 0.;
	// MPC particles
	if( CL.pp!=NULL ) {
		tmpc = CL.pp;
		while( tmpc!=NULL ) {
			id = tmpc->SPID;
			mass=(SP+id)->MASS;
			summ += mass;
			for( i=0; i<DIM; i++ ) vcm[i] += tmpc->V[i] * mass;
			//Increment link in list
			tmpc = tmpc->next;
		}
	}
	// MD particles
	if( CL.MDpp!=NULL ) {
		tmd = CL.MDpp;
		while( tmd!=NULL ) {
			mass = tmd->mass;
			summ += mass;
			vcm[0] += tmd->vx * mass;
			vcm[1] += tmd->vy * mass;
			vcm[2] += tmd->vz * mass;
			//Increment link in list
			tmd = tmd->nextSRD;
		}
	}
	// Swimmer particles
	if( CL.sp!=NULL ) {
		tsm = CL.sp;
		while( tsm!=NULL ) {
			if( tsm->HorM ) mass = (double) specS.middM;
			else mass = (double) specS.headM;
			summ += mass;
			for( i=0; i<DIM; i++ ) vcm[i] += tsm->V[i] * mass;
			//Increment link in list
			tsm = tsm->next;
		}
	}
	for( i=0; i<DIM; i++ ) vcm[i] /= (summ ? summ : 1.);
}

///
/// @brief This routine applies force dipoles from active MD particles onto the MPCD particles.
///
/// The function loops over all MD particles and finds their location and backbone tangent. It
/// then finds the plane and applies a kick to all the MPCD particles in the cell either away from
/// or towards the plane, depending on if extensile or contractile.
/// @param simMD A pointer to the entire MD portion of the simulation.
/// @param CL ALL cells. 
/// @param SP The species-wide information about MPCD particles.
///
void activeMD(simptr simMD, cell ***CL, spec *SP, inputList in) {
	int i, j, k, nAtom, cx, cy, cz, N_cl, id;
	double pmOne;
	double MDloc[_3D], prevloc[DIM], nextloc[DIM], vecprev[DIM], vecnext[DIM], tangent[DIM];  // stuff to calc plane and velocity
	double mpcdCM[DIM];
	double pW;  //The particle's pW for passing the plane
	double force; // the dipole 'force'(=dt*d/m) on one mpcd particle, will be a constant for each cell. delta_v = 'force'/mass
	double dt, deltaV;
	float dipole;
	bc PLANE;  //The plane that cuts the cell in half
	particleMPC *tmpc;
	particleMD	*atoms;
	
	atoms  = simMD->atom.items;
	nAtom = simMD->atom.n;
	dt = in.dt;

	for (i=0; i<nAtom; i++) {
		// zero things
		for (j=0; j<_3D; j++) {
			MDloc[j] = 0.;	// this has to be 3D always because CL is 3d always
		}

		for (j=0; j<DIM; j++) {
			prevloc[j] = 0.;
			nextloc[j] = 0.;
			vecprev[j] = 0.;
			vecnext[j] = 0.;
			tangent[j] = 0.;
		}
		N_cl = 0;
		cx = 0.;
		cy = 0.;
		cz = 0.;
		force = 0.;
		// find location of MD particle(s) in question
		// there is probably a more concise way of doing this
		MDloc[0] = atoms[i].rx;
		MDloc[1] = atoms[i].ry;
		MDloc[2] = atoms[i].rz;
		
		// get tangent of backbone at MD particle, and normalise it - different if first or last
		
		if (i==0) {
			nextloc[0] = atoms[i+1].rx;
			nextloc[1] = atoms[i+1].ry;
			if (DIM>2){
				nextloc[2] = atoms[i+1].rz;  // index 2 only exists if 3d
			}
			for (j=0; j<DIM; j++) {
			tangent[j] = nextloc[j]-MDloc[j];
			}
		}
		else if (i==nAtom) {
			prevloc[0] = atoms[i-1].rx;
			prevloc[1] = atoms[i-1].ry;
			if (DIM>2){
				prevloc[2] = atoms[i-1].rz;
			}
			for (j=0; j<DIM; j++) {
			tangent[j] = MDloc[j]-prevloc[j];
			}
		}
		else {
			prevloc[0] = atoms[i-1].rx;
			prevloc[1] = atoms[i-1].ry;
			nextloc[0] = atoms[i+1].rx;
			nextloc[1] = atoms[i+1].ry;
			if (DIM>2){
				prevloc[2] = atoms[i-1].rz;
				nextloc[2] = atoms[i+1].rz;
			}
			for (j=0; j<DIM; j++) {
				vecprev[j] = MDloc[j]-prevloc[j];
				vecnext[j] = nextloc[j]-MDloc[j];
			}
			norm(vecprev, DIM); //normprev = sqrt(pow(vecprev[0],2) + pow(vecprev[1],2)+pow(vecprev[2],2));
			norm(vecnext, DIM); //normnext = sqrt(pow(vecnext[0],2) + pow(vecnext[1],2)+pow(vecnext[2],2));
			for (j=0; j<DIM; j++) {
				//vecprev[j] = vecprev[j]/normprev;
				//vecnext[j] = vecnext[j]/normnext;
				tangent[j] = vecprev[j]+vecnext[j];
			}
		}
		norm(tangent, DIM);//normtan = sqrt(pow(tangent[0],2) + pow(tangent[1],2)+pow(tangent[2],2));
		//for (j=0; j<_3D; j++) {
		//	tangent[j] = tangent[j]/normtan;
		//}

		// check
		//printf("tangent for atom %d ", i+1);
		//printf("is %f, ", tangent[0]);
		//printf("%f, ", tangent[1]);
		//if (DIM>2){
		//	printf("%f\n", tangent[2]);
		//}

		// identify which MPCD cell it's in
		cx = (int)MDloc[0];
		cy = (int)MDloc[1];
		cz = (int)MDloc[2];	// gotta be 3d always
		// then use like CL[cx][cy][cz].pp

		// find const 'force' on each MPCD particle from the MD in this cell
		dipole = atoms[i].dipole;
		N_cl = CL[cx][cy][cz].POPSRD;
		force = dt*dipole/N_cl;

		// Think these need introduced here rather than at the very start because looking at one cell rather than all
		//double AV[CL[cx][cy][cz].POP][DIM];	//Active velocities
		double AS[DIM];			//Sum of active velocities
		// Zero arrays
		//for( k=0;k<CL[cx][cy][cz].POP;k++ ) for( j=0;j<_3D;j++ ) {
		//	AV[k][j] = 0.;
		//}
		for( j=0;j<DIM;j++ ) {
			AS[j]=0.;
			mpcdCM[j]=0.;
		}

		// find plane - line 1775 in lc.c
		// Define the plane normal to the backbone tangent at the MD particle position
		for( k=0;k<4;k++ ) PLANE.P[k]=1;
		PLANE.INV=0;
		PLANE.ABS=0;
		PLANE.R=0.0;
		PLANE.ROTSYMM[0]=4.0;
		PLANE.ROTSYMM[1]=4.0;
		// Normal is TANGENT
		for( k=0; k<DIM; k++ ) PLANE.A[k] = tangent[k];
		// Position is cell COM
		// ONLY TAKE INTO ACCOUNT MPCD PARTICLES NOT MD AS WELL

		localCM_SRD(CL[cx][cy][cz],SP,mpcdCM); //i think that's how this works??
		for( k=0; k<DIM; k++ ) PLANE.Q[k] = mpcdCM[k];
		
		if(CL[cx][cy][cz].pp!=NULL) {
			tmpc = CL[cx][cy][cz].pp;
			// loop through MPCD particles
			//m = 0;// counter for momentum conservation bits
			while( tmpc!=NULL ) {
				// zero things
				pW = 0.;
				id = 0;
				deltaV = 0.;
				// give a kick according to which side of plane - dipole AndersenROT_LC
				// Check which side of the plane
				pW = calcW( PLANE,*tmpc );
				if( pW<=0 ) pmOne=-1.;
				else pmOne=1.;

				// change in velocity to be applied according to species mass
				id = tmpc->SPID;
				deltaV = force/(double)(SP+id)->MASS; //think this is correct way of getting mass?
				//printf("vel kick: %f\n",deltaV);
				// apply the change (hopefully). In direction of tangent, away from (or towards) plane.
				for( k=0; k<DIM; k++ ){
					tmpc->V[k] += tangent[k]*deltaV*pmOne;
				}
				
				// stuff for momentum conservation
				//for( j=0; j<DIM; j++ ) AV[m][j] = tmpc->V[j];
				//for( j=0; j<DIM; j++ ) AS[j] += AV[m][j]*(double)(SP+id)->MASS;
				for( j=0; j<DIM; j++ ) AS[j] += tangent[k]*force*pmOne;
				//m++;
				// Increment link in list
				tmpc = tmpc->next;
			}
		}
		// find net momentum - divided by N because then shared between all particles
		for( j=0; j<DIM; j++ ) AS[j] = AS[j]/N_cl; // is this correct?

		// subtract net momentum from all mpcd if non-zero - new loop
		if(CL[cx][cy][cz].pp!=NULL) {
			tmpc = CL[cx][cy][cz].pp;
			while( tmpc!=NULL ) {
				for (j=0; j<DIM; j++){
					tmpc->V[j] -= AS[j]/(double)(SP+id)->MASS; // ;
				}
				// Increment link in list
				tmpc = tmpc->next;
			}
		}	
			
	}

}

/// 
/// @brief This routine calculates the centre of mass velocity of the MPCD particles in a single cell or bin.
///
/// The function loops over all MPCD particles within a given cell and calculates the center of mass velocity.
/// It differs from localVCM(), in that it only considers MPCD particles and not all particles. 
/// It updates the particles by looping through the linked list. 
/// @param vcm The centre of mass velocity vector of MPCD particles in cell `CL`. Velocity is returned through this variable.
/// @param CL One specific MPCD cell.
/// @param SP The species-wide information about MPCD particles.
/// @see localFLOW()
/// @see localVCM()
///
void localMPCVCM( double vcm[_3D],cell CL,spec *SP ) {
	int id,i;
	double summ = 0.0;
	double mass;
	particleMPC *tmpc;

	// Zero the centre of mass velocity
	for( i=0; i<DIM; i++ ) vcm[i] = 0.;
	// MPC particles
	if( CL.pp!=NULL ) {
		tmpc = CL.pp;
		while( tmpc!=NULL ) {
			id = tmpc->SPID;
			mass=(SP+id)->MASS;
			summ += mass;
			for( i=0; i<DIM; i++ ) vcm[i] += tmpc->V[i] * mass;
			//Increment link in list
			tmpc = tmpc->next;
		}
	}
	for( i=0; i<DIM; i++ ) vcm[i] /= (summ ? summ : 1.);
}

/// 
/// @brief This routine calculates the local mass of a cell.
/// 
/// It loops through the linked lists, to calculate the centre of mass position of a given cell. 
/// It includes MPCD, MD and swimmer particles. 
/// @param CL An MPCD cell (including the linked list of particles in each cell). 
/// @param SP The species-wide information about MPCD particles.
/// @param specS The species-wide information about swimmers.
/// @return The local mass of a cell (including MPCD, MD and swimmer particles).
///
double localMASS( cell CL,spec *SP,specSwimmer specS ) {
	int id;
	double M = 0.0;
	particleMPC *tmpc;
	particleMD  *tmd;
	smono *tsm;

	// MPC particles
	if( CL.pp!=NULL ) {
		tmpc = CL.pp;
		while( tmpc!=NULL ) {
			id = tmpc->SPID;
			M += (SP+id)->MASS;
			//Increment link in list
			tmpc = tmpc->next;
		}
	}
	// MD particles
	if( CL.MDpp!=NULL ) {
		tmd = CL.MDpp;
		while( tmd!=NULL ) {
			M += tmd->mass;
			//Increment link in list
			tmd = tmd->nextSRD;
		}
	}
	// swimmer particles
	if( CL.sp!=NULL ) {
		tsm = CL.sp;
		while( tsm!=NULL ) {
			if( tsm->HorM ) M += (double) specS.middM;
			else M += (double) specS.headM;
			//Increment link in list
			tsm = tsm->next;
		}
	}
	return M;
}

/// 
/// @brief This routine calculates the local temperature in the cell, via the equipartition function.
/// 
/// It loops through the linked lists, to calculate the local temperature of a given cell. 
/// It includes MPCD, MD and swimmer particles. 
/// @param CL An MPCD cell (including the linked list of particles in each cell). 
/// @param SP The species-wide information about MPCD particles.
/// @param specS The species-wide information about swimmers.
/// @return The local thermal energy of a cell (via equipartition theorem). 
///
double localTEMP( cell CL,spec *SP,specSwimmer specS ) {
	int d,id,p = 0;
	double V[_3D],mass,KBT = 0.0;
	particleMPC *tmpc;
	particleMD *tmd;
	smono *tsm;

	// Zero
	for( d=0; d<_3D; d++ ) V[d] = 0.0;
	// MPC particles
	if( CL.pp!=NULL ) {
		tmpc = CL.pp;
		while( tmpc!=NULL ) {
			p++;
			id = tmpc->SPID;
			mass = (SP+id)->MASS;
			for( d=0; d<DIM; d++ ) V[d] = tmpc->V[d];
			KBT += mass*(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]);
			//Increment link in list
			tmpc = tmpc->next;
		}
	}
	// Md particles
	if( CL.MDpp!=NULL ) {
		tmd = CL.MDpp;
		while( tmd!=NULL ) {
			p++;
			mass=(double) (tmd->mass);
			V[0] = tmd->vx;
			V[1] = tmd->vy;
			V[2] = tmd->vz;
			KBT += mass*(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]);
			//Increment link in list
			tmd = tmd->nextSRD;
		}
	}
	// Swimmer particles
	if( CL.sp!=NULL ) {
		tsm = CL.sp;
		while( tsm!=NULL ) {
			p++;
			if( tsm->HorM ) mass = (double) specS.middM;
			else mass = (double) specS.headM;
			for( d=0; d<DIM; d++ ) V[d] = tsm->V[d];
			KBT += mass*(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]);
			//Increment link in list
			tsm = tsm->next;
		}
	}
	//By equipartition function
	KBT /= (double)(DIM*p);
	return KBT;
}

/// 
/// @brief This routine calculates the local population (number of particles in the cell). 
/// 
/// It loops through the linked lists, to sum the total number of particles. 
/// It includes MPCD, MD and swimmer particles. 
/// @param CL An MPCD cell (including the linked list of particles in each cell). 
/// @return Total number of particles in this cell.
///
int localPOP( cell CL ) {
	int i = 0;
	particleMPC *tmpc;
	particleMD *tmd;
	smono *tsm;

	// MPC particles
	if( CL.pp!=NULL ) {
		tmpc = CL.pp;
		while( tmpc!=NULL ) {
			i++;
			//Increment link in list
			tmpc = tmpc->next;
		}
	}
	// MD particles
	if( CL.MDpp!=NULL ) {
		tmd = CL.MDpp;
		while( tmd!=NULL ) {
			i++;
			//Increment link in list
			tmd = tmd->nextSRD;
		}
	}
	// MPC particles
	if( CL.sp!=NULL ) {
		tsm = CL.sp;
		while( tsm!=NULL ) {
			i++;
			//Increment link in list
			tsm = tsm->next;
		}
	}
	return i;
}

/// 
/// @brief This routine randomly `scrambles` the MPCD particle velocities. 
///
/// The function is used by when noHI is set to `HIOFF`. 
/// It creates a computationally expensive Brownian thermostat, with no hydrodynamic interactions, but with the same material properties. 
/// It does this by "scrambling the velocities". 
/// It does this by randomly exchanges particle velocities.
/// It loops over the global population (`GPOP`) and randomly choses another MPCD particle. 
/// The two velocities are swapped. 
/// @param p All the MPCD particles. 
/// @note It is slow because it switches some velocities more than once. 
/// But compared to copying all the velocities and managing/organizing the array every time, this is faster.
/// Plus, it is expected that noHI will only be seldom used, only to check the effects of hydrodynamic interactions. 
///
void scramble( particleMPC *p ) {
	int i,j,k;
	double temp;

	for( i=0; i<GPOP; i++ ) {
		//Randomly chose a particle from those left
		k = rand_particle( GPOP );
		//Give the velocity to the current particle
		for( j=0; j<DIM; j++ ) {
			temp = p[i].V[j];
			p[i].V[j] = p[k].V[j];
			p[k].V[j] = temp;
		}
	}
}

/// 
/// @brief Calculates the centre of mass position for a given cell. 
///
/// This routine finds the centre of mass position of a cell, including MPCD, MD and swimmer particles. 
/// This is now a legacy method after being refactored into localPROP(). 
/// @param CL An MPCD cell (including the linked list of particles in each cell). 
/// @param SP The species-wide information about MPCD particles.
/// @param specS The species-wide information about swimmers.
/// @note It assumes the cell mass has already been computed.
///
void localCM( cell *CL,spec *SP,specSwimmer specS ) {
	int id,d;
	double mass,cellMass,Q[_3D];
	particleMPC *pMPC;	//Temporary pointer to MPC particles
	particleMD *pMD;	//Temporary pointer to MD particles
	smono *pSW;			//Temporary pointer to swimmer monomers

	//Zero everything for recounting
	for( d=0; d<DIM; d++ ) CL->CM[d] = 0.;
	cellMass = 0.0;
	//Find local values
	if( CL->pp!=NULL || ( CL->MDpp!=NULL && MDmode==MDinMPC ) || CL->sp!=NULL ) {
		// SRD particles
		if( CL->pp!=NULL ) {
			pMPC = CL->pp;
			while(pMPC!=NULL) {
				id = pMPC->SPID;
				mass = (SP+id)->MASS;
				cellMass+=mass;
				for( d=0; d<DIM; d++ ) Q[d] = pMPC->Q[d];
				for( d=0; d<DIM; d++ ) CL->CM[d] += Q[d] * mass;
				//Increment link in list
				pMPC = pMPC->next;
			}
		}
		// MD particles
		if(MDmode==MDinMPC) if( CL->MDpp!=NULL) {
			pMD = CL->MDpp;
			while( pMD!=NULL ) {
				mass = pMD->mass;
				cellMass+=mass;
				Q[0] = pMD->rx;
				if( DIM > 1 ) Q[1] = pMD->ry;
				if( DIM > 2 ) Q[2] = pMD->rz;
				for( d=0; d<DIM; d++ ) CL->CM[d] += Q[d] * mass;
				//Increment link in list
				pMD = pMD->nextSRD;
			}
		}
		// Swimmer particles
		if( CL->sp!=NULL ) {
			pSW = CL->sp;
			while(pSW!=NULL) {
				if( pSW->HorM ) mass = (double) specS.middM;
				else mass = (double) specS.headM;
				cellMass+=mass;
				for( d=0; d<DIM; d++ ) Q[d] = pSW->Q[d];
				for( d=0; d<DIM; d++ ) CL->CM[d] += Q[d] * mass;
				//Increment link in list
				pSW = pSW->next;
			}
		}
		// Make sums into averages
		CL->MASS = cellMass;
		for( d=0; d<DIM; d++ ) CL->CM[d] /= cellMass;
	}
}

/// 
/// @brief Calculates the centre of mass position for a given cell's MPCD (SRD) particles. 
///
/// This routine finds the centre of mass position of only the MPCD particles in a cell. 
/// @param CL An MPCD cell (including the linked list of particles in each cell). 
/// @param SP The species-wide information about MPCD particles.
/// @param r_cm The centre of mass position vector of MPCD particles in cell `CL`. Position is returned through this variable.
///
void localCM_SRD( cell CL,spec *SP,double r_cm[] ) {
	int id,d;
	double mass,cellMass;
	particleMPC *pMPC;	//Temporary pointer to MPC particles

	//Zero everything for recounting
	for( d=0; d<DIM; d++ ) r_cm[d] = 0.;
	cellMass = 0.0;
	// SRD particles
	if( CL.pp!=NULL ) {
		pMPC = CL.pp;
		while(pMPC!=NULL) {
			id = pMPC->SPID;
			mass = (SP+id)->MASS;
			cellMass+=mass;
			for( d=0; d<DIM; d++ ) r_cm[d] += pMPC->Q[d] * mass;
			//Increment link in list
			pMPC = pMPC->next;
		}
	}
	// Make sums into averages
	for( d=0; d<DIM; d++ ) r_cm[d] /= cellMass;
}

/// 
/// @brief Calculates the moment of inertia tensor for a given cell. 
///
/// This routine calculates the moment of inertia tensor for a given cell, including MPCD, MD and swimmer particles. 
/// It updates the particles by looping through the linked lists. 
/// @param CL An MPCD cell (including the linked list of particles in each cell). 
/// @param SP The species-wide information about MPCD particles.
/// @param specS The species-wide information about swimmers.
/// @note It assumes the cell mass has already been computed.
///
void localMomInertiaTensor( cell *CL,spec *SP,specSwimmer specS ) {
	int i,j,id,d;
	double mass,Q[_3D];
	particleMPC *pMPC;	//Temporary pointer to MPC particles
	particleMD *pMD;	//Temporary pointer to MD particles
	smono *pSW;			//Temporary pointer to swimmer monomers

	// Zero
	for( i=0; i<_3D; i++ ) for( j=0; j<_3D; j++ ) CL->I[i][j] = 0.0;
	for( i=0; i<_3D; i++ ) Q[i] = 0.0;
	// SRD particles
	if( CL->pp!=NULL ) {
		pMPC = CL->pp;
		while( pMPC!=NULL ) {
			id = pMPC->SPID;
			mass = (SP+id)->MASS;
			for( d=0; d<DIM; d++ ) Q[d] = pMPC->Q[d] - CL->CM[d];
			// Calculate moment of inertia
			// I11
			CL->I[0][0] += mass * (Q[1]*Q[1] + Q[2]*Q[2]);
			// I22
			CL->I[1][1] += mass * (Q[0]*Q[0] + Q[2]*Q[2]);
			// I33
			CL->I[2][2] += mass * (Q[0]*Q[0] + Q[1]*Q[1]);
			// I12
			CL->I[0][1] += mass * Q[0] * Q[1];
			// I13
			CL->I[0][2] += mass * Q[0] * Q[2];
			// I23
			CL->I[1][2] += mass * Q[1] * Q[2];
			//Increment link in list
			pMPC = pMPC->next;
		}
	}
	// MD particles
	if( CL->MDpp!=NULL ) {
		pMD = CL->MDpp;
		while( pMD!=NULL ) {
			mass = pMD->mass;
			Q[0] = pMD->rx - CL->CM[0];
			if( DIM > 1 ) Q[1] = pMD->ry - CL->CM[1];
			if( DIM > 2 ) Q[2] = pMD->rz - CL->CM[2];
			// Calculate moment of inertia
			// I11
			CL->I[0][0] += mass * (Q[1]*Q[1] + Q[2]*Q[2]);
			// I22
			CL->I[1][1] += mass * (Q[0]*Q[0] + Q[2]*Q[2]);
			// I33
			CL->I[2][2] += mass * (Q[0]*Q[0] + Q[1]*Q[1]);
			// I12
			CL->I[0][1] += mass * Q[0] * Q[1];
			// I13
			CL->I[0][2] += mass * Q[0] * Q[2];
			// I23
			CL->I[1][2] += mass * Q[1] * Q[2];
			//Increment link in list
			pMD = pMD->nextSRD;
		}
	}
	// Swimmer particles
	if( CL->sp!=NULL ) {
		pSW = CL->sp;
		while( pSW!=NULL ) {
			if( pSW->HorM ) mass = (double) specS.middM;
			else mass = (double) specS.headM;
			for( d=0; d<DIM; d++ ) Q[d] = pSW->Q[d] - CL->CM[d];
			// Calculate moment of inertia
			// I11
			CL->I[0][0] += mass * (Q[1]*Q[1] + Q[2]*Q[2]);
			// I22
			CL->I[1][1] += mass * (Q[0]*Q[0] + Q[2]*Q[2]);
			// I33
			CL->I[2][2] += mass * (Q[0]*Q[0] + Q[1]*Q[1]);
			// I12
			CL->I[0][1] += mass * Q[0] * Q[1];
			// I13
			CL->I[0][2] += mass * Q[0] * Q[2];
			// I23
			CL->I[1][2] += mass * Q[1] * Q[2];
			//Increment link in list
			pSW = pSW->next;
		}
	}
	// Symmetry of the moment of inertia
	CL->I[0][1] *= -1.;
	CL->I[0][2] *= -1.;
	CL->I[1][2] *= -1.;
	CL->I[1][0] = CL->I[0][1];
	CL->I[2][0] = CL->I[0][2];
	CL->I[2][1] = CL->I[1][2];
}

/// 
/// @brief Calculates the moment of inertia value of the MPCD particles (here called SRD) for a given cell. 
/// 
/// This routine calculates the moment of inertia of only the MPCD particles in a cell. 
/// The value returned is the moment of inertia about a given position r0 and axis n (assumed normalized). 
/// @param CL An MPCD cell (including the linked list of particles in each cell). 
/// @param SP The species-wide information about MPCD particles.
/// @param r0 The point about which the rotation occurs.
/// @param n The axis about which the rotation occurs. 
/// @return Magnitude of the moment of inertia about the given position and axis. 
///
double localMomInertia_SRD( cell CL,spec *SP,double r0[],double n[] ) {
	int id,d;
	double mass,momI,d2,r[DIM],r_perp[DIM];
	particleMPC *pMPC;	//Temporary pointer to MPC particles

	//Zero
	momI = 0.0;
	// SRD particles
	if( CL.pp!=NULL ) {
		pMPC = CL.pp;
		while(pMPC!=NULL) {
			id = pMPC->SPID;
			mass = (SP+id)->MASS;
			for( d=0; d<DIM; d++ ) r[d] = pMPC->Q[d] - r0[d];
			//Projection of r on direction (d2 used as temporary value)
			d2 = dotprod( r,n,DIM );
			//Perpendicular director from the particle to the line of the swimmer's orientation
			for( d=0; d<DIM; d++ ) r_perp[d] = r[d] - d2*n[d];
			//Distance from axis squared
			d2 = dotprod( r_perp,r_perp,DIM );
			momI += d2*mass;
			//Increment link in list
			pMPC = pMPC->next;
		}
	}
	return momI;
}

/// 
/// @brief This is a very coarse check to make sure that none of the MPCD particles have escaped the simulation. 
/// 
/// This is a strong-armed check to ensure that none of the MPCD particles have escaped the control volume. 
/// This function loops over the global population (`GPOP`) to update all MPCD particle positions. 
/// It is a legacy function that was used for debugging --- it might be valuable again in future developments. 
/// @param pp All the MPCD particles. 
///
void checkEscape_all( particleMPC *pp ) {
	int i,d;
	double Q[_3D];
	for( i=0; i<GPOP; i++ ){
		for( d=0; d<_3D; d++ ) Q[d] = (pp+i)->Q[d];
		for( d=0; d<_3D; d++ ) if( Q[d]<0. || Q[d]>XYZ_P1[d] )  {
			#ifdef DBG
				printf( "Warning: Particle %d escaped control volume",i );
				pvec( Q,_3D );
			#endif
		}
	}
}

/// 
/// @brief Applies a change in velocity to every MPCD particle in a cell. 
///
/// This routine loops through the linked list of a given cell, adding a constant velocity to each particle velocity. 
/// It adds to every particle in the cell, including MPCD, MD and swimmer particles. 
/// @param CL An MPCD cell (including the linked list of particles in each cell). 
/// @param addVel The vector that is being added to the velocity of every MPCD particle in cell `CL`.
///
void cellVelForce( cell *CL,double addVel[3] ) {
	// int i;
	int j;
	particleMPC *tmpc;	//Temporary particleMPC
	particleMD *tmd;	//Temporary particleMD
	smono *pSW;			//Temporary pointer to swimmer monomers

	//Give particles a kick
	// MPC particles
	tmpc = CL->pp;
	// i=0;
	while( tmpc!=NULL ) {
		for( j=0; j<DIM; j++ ) tmpc->V[j] += addVel[j];
		//Increment link in list
		tmpc = tmpc->next;
		// i++;
	}
	//MD particles
	tmd = CL->MDpp;
	while( tmd!=NULL ) {
		tmd->vx += addVel[0];
		tmd->vy += addVel[1];
		tmd->vz += addVel[2];
		//Increment link in list
		tmd = tmd->nextSRD;
		// i++;
	}
	// Swimmer particles
	pSW = CL->sp;
	while( pSW!=NULL ) {
		for( j=0; j<DIM; j++ ) pSW->V[j] += addVel[j];
		//Increment link in list
		pSW = pSW->next;
		// i++;
	}
}

/// 
/// @brief Overrides the velocity of every MPCD particle in a cell. 
/// 
/// This routine loops through the linked list of a given cell, overwriting every particle velocity to a set value. 
/// It sets the velocity of every particle in the cell, including MPCD, MD and swimmer particles. 
/// @param CL An MPCD cell (including the linked list of particles in each cell). 
/// @param vel The vector that every particle's velocity is set to within the cell `CL`. 
///
void cellVelSet( cell *CL,double vel[3] ) {
	// int i;
	int j;
	particleMPC *tmpc;	//Temporary particleMPC
	particleMD *tmd;	//Temporary particleMD
	smono *pSW;			//Temporary pointer to swimmer monomers

	//Give particles a kick
	// MPC particles
	tmpc = CL->pp;
	// i=0;
	while( tmpc!=NULL ) {
		for( j=0; j<DIM; j++ ) tmpc->V[j] = vel[j];
		//Increment link in list
		tmpc = tmpc->next;
		// i++;
	}
	//MD particles
	tmd = CL->MDpp;
	while( tmd!=NULL ) {
		tmd->vx += vel[0];
		tmd->vy += vel[1];
		tmd->vz += vel[2];
		//Increment link in list
		tmd = tmd->nextSRD;
		// i++;
	}
	// Swimmer particles
	pSW = CL->sp;
	while( pSW!=NULL ) {
		for( j=0; j<DIM; j++ ) pSW->V[j] += vel[j];
		//Increment link in list
		pSW = pSW->next;
		// i++;
	}
}

/// 
/// @brief The timestep routine contains all the routines that happen in each time iteration. 
///
/// This routine includes all the major aspects of the MPCD code. 
/// Everything in a time step, except writing output and checkpointing, is included within this function. 
/// The outline of everything done in a single time step is:
/// - Integrate molecular dynamics type particles (integrateMD() and integrateSwimmers()).
/// - Grid shift (ranshift() and gridShift_all()).
/// - Bin (bin(), binSwimmers() and binMD()).
/// - Accelerate the particles (acc_all()).
/// - Calculate local properties (localPROP()).
/// - Apply ghost particles at surfaces (ghostPart()). 
/// - Collision operation
///   * Liquid crystal collision operation (LCcollision(), magTorque_all() and jefferysTorque()).
///   * Velocity collision operation (MPCcollision(), multiphaseColl(), incompColl() and scramble())
/// - Temperature scaling (scaleT()).
/// - Grid shift back (gridShift_all()).
/// - Swimmer forces on fluid (swimmerDipole()).
/// - Streaming (stream_all()).
/// - BCs
///   + MPCD particle collision with BCs (MPC_BCcollision()).
///   + Stream/rotate the BCs (stream_BC() and spin_BC()).
///   + BC collision with other BCs (BC_BCcollision()).
///   + Moving BC collision with MPC particles (BC_MPCcollision()).
/// - Re-Bin (bin(), binSwimmers() and binMD()).
/// - Re-Local properties (localPROP()).
/// @param CL All of the MPCD cells. 
/// @param SRDparticles All the MPCD particles. 
/// @param SP The species-wide information about MPCD particles.
/// @param WALL All of the walls (boundary conditions) that particles might interact with. 
/// @param simMD A pointer to the entire MD portion of the simulation.
/// @param SS The species-wide information about swimmers.
/// @param swimmers All the swimmers, including their head and middle monomers. 
/// @param AVNOW The average velocity vector of the entire system after the collision operations are applied to the fluid. 
/// @param AVV The average flow velocity that is subtracted off when thermostatting the temperature. 
/// @param avDIR The average global director from the nematic Q-tensor.
/// @param in The complete list of all input parameters. 
/// @param KBTNOW The current thermal energy (not the thermostats target temperature). 
/// @param AVS  The average global scalar order parameter from the nematic Q-tensor.
/// @param runtime The number of iterations that timestep will ultimately be looped over. 
/// @param MD_mode The MD coupling mode. Can be off (noMD), MD particles included in the MPCD collisions (MDinMPC), or MPCD particles included in MD pair interactions (MPCinMD).
/// @param outFlags The complete list of what is being outputted as results. 
/// @param outFiles The complete list of pointers to output files. 
///
void timestep(cell ***CL, particleMPC *SRDparticles, spec SP[], bc WALL[], simptr simMD, specSwimmer *SS, swimmer swimmers[], double AVNOW[_3D], double AVV[_3D], double avDIR[_3D], inputList in, double *KBTNOW, double *AVS, int runtime, int MD_mode, outputFlagsList outFlags, outputFilesList outFiles) {

	int i,j,k,l;					//Counting variables
	double RSHIFT[_3D];				//Random vector to positively shift all components of the simulation
	double CLQ[_3D];				//Position of the cell since calculating pressure sucks
	int BC_FLAG;					//Flags if the BC moved in this time step
	int outPressure=0;				//Whether to make the pressure calculations (never used just outputted)
	int bcCNT,reCNT,rethermCNT;		//Count if any particles had problems with the BCs
    // Active layer variables
    int zeroMFPot;                	// Flag to switch to zero effective potential
    double savedAct[NSPECI];        // Storage to save the activity of species, to make it height dependent

	#ifdef DBG
		if ( DBUG >= DBGSTEPS ) {
			if( in.warmupSteps ) printf( "\nBegin warmup time step %i. Simulation time = %lf\n",runtime,runtime*in.dt );
			else printf( "\nBegin time step %i. Simulation time = %lf\n",runtime,runtime*in.dt );
		}
	#endif
	//Zero counters
	if( outFlags.PRESOUT>=OUT && runtime%outFlags.PRESOUT==0 ) outPressure=1;
	zerocnt( KBTNOW,AVNOW,AVS );
	// Zero impulse on BCs
	// NOTE: Louise thinks this should be fine being here (and did check),
	// This was moved when editing ghostPart to increase anchoring strength
	// If something looks bad related to mobile walls, maybe start looking here.
	for( i=0; i<NBC; i++ ) {
		zerovec(WALL[i].dV,DIM);
		zerovec(WALL[i].dL,_3D);
	}

	/* ****************************************** */
	/* ************* INTEGRATE MD *************** */
	/* ****************************************** */
	#ifdef DBG
		if( DBUG >= DBGTITLE && MD_mode != noMD ) printf("Integrate MD.\n" );
	#endif
	if(MD_mode){
		if(simMD->warmupMD == FREE_WARMUP || simMD->warmupMD == POS_WARMUP){
			integrateMD(simMD, MD_mode, in.stepsMD, SRDparticles, WALL, SP, GPOP, NBC, CL);
		}else if(simMD->warmupMD == PINNED_WARMUP){
			integrateMD_Pinned(simMD, MD_mode, in.stepsMD, SRDparticles, WALL, SP, GPOP, NBC, CL);
		}
	}
	/* ****************************************** */
	/* *********** INTEGRATE SWIMMER ************ */
	/* ****************************************** */
	if( NS>0 ) {
		#ifdef DBG
			if( DBUG >= DBGTITLE ) printf( "Integrate swimmer.\n" );
		#endif
		//Run/Tumble --- if runTime or tumbleTime==0 then don't do it (if either is zero both are zeroed in initialization)
		if( SS->runTime>TOL ) {
			#ifdef DBG
				if( DBUG == DBGSWIMMER || DBUG == DBGSWIMMERDEETS ) printf( "\tApply run/tumble dynamics.\n" );
			#endif
			runTumbleDynamics( SS,swimmers,WALL,in.stepsMD,in.MAG,in.dt,outFlags.RTOUT,outFiles.fruntumble );
		}
		//MD integration
		#ifdef DBG
			if( DBUG == DBGSWIMMER || DBUG == DBGSWIMMERDEETS ) printf( "\tMD integration.\n" );
		#endif
		integrateSwimmers( *SS,swimmers,WALL,in.stepsMD,in.dt,in.MAG,HOOKESPRING );
	}
	// //Apply swimmer dipole (both force and torque dipoles)
	// #ifdef DBG
	// 	if( DBUG == DBGSWIMMER || DBUG == DBGSWIMMERDEETS ) printf( "\tApply swimmer dipole.\n" );
	// #endif
	// swimmerDipole( *SS,swimmers,CL,SP,in.dt,SRDparticles,WALL,simMD );
	/* ****************************************** */
	/* *** GRID SHIFT FOR GALILEAN INVARIANCE *** */
	/* ****************************************** */
	if( in.GALINV ) {
		#ifdef DBG
			if( DBUG >= DBGTITLE ) printf( "Shift Grid.\n" );
		#endif
		//Generate random shift
		ranshift( RSHIFT,in.GALINV,DIM );
		//Shift entire system by RSHIFT
		gridShift_all(RSHIFT, 0, SRDparticles, WALL, simMD, swimmers, MD_mode );
	}

	/* ******************************************************/
	/* * SUBTRACTING THE ORIENTATION PART FROM THE VELOCITY */
	/* ******************************************************/
	//For now only one specie is supported
	if (in.LC == BCT) {
		for (int i = 0; i < GPOP; i++) {
			for( j=0; j<DIM; j++ ) {
				(SRDparticles+i)->V[j] -= (SRDparticles+i)->U[j]*((SP+0)->BS);
			}
		}
	}

	/* ****************************************** */
	/* ******************* BIN ****************** */
	/* ****************************************** */
	#ifdef DBG
		if( DBUG >= DBGTITLE ) printf( "Bin Particles.\n" );
	#endif
	// Bin SRD particles
	bin( CL,SP,WALL,in.KBT,in.LC,1 );
	// Bin swimmer monomers
	binSwimmers( CL,1 );
	// Bin MD particles
	if( MD_mode ) binMD(CL );
	/* ****************************************** */
	/* ************** PRE-COLLISION ************* */
	/* ****************************************** */
	#ifdef DBG
		if (DBUG == DBGTHERM) {
			*KBTNOW = TEMP( SRDparticles,SP,WALL,AVV );
			printf( "\tSystem Temperature pre-collision: %09e\n",*KBTNOW );
		}
	#endif
	/* ****************************************** */
	/* **************** PRESSURE **************** */
	/* ****************************************** */
	if( outPressure ) calcPressureStreaming( CL,SP );
	/* ****************************************** */
	/* ************** ACCELERATION ************** */
	/* ****************************************** */
	// The acceleration is really part of the collision
	#ifdef DBG
		if( DBUG >= DBGTITLE && in.GRAV_FLAG) printf( "Apply Acceleration Operator.\n" );
	#endif
	if( in.GRAV_FLAG ) acc_all( SRDparticles,in.dt,in.GRAV );
	/* ****************************************** */
	/* **************** COLLISION *************** */
	/* ****************************************** */
	#ifdef DBG
		if( DBUG >= DBGTITLE ) printf( "Collisions.\n" );
	#endif
	//Calculate the local properties of each cell (VCM,KBT,POPulation,Mass)
	//Do this AFTER acceleration so that use accelerated VCM in collision
	localPROP( CL,SP,*SS,in.RTECH,in.LC );
	/* ****************************************** */
	/* *********** ADD GHOST PARTICLES ********** */
	/* ****************************************** */
	ghostPart( CL,WALL,in.KBT,in.LC,SP );

	/* ****************************************** */
	/* *********LIQUID CRYSTAL COLLISION ******** */
	/* ****************************************** */
	//Liquid Crystal collision operations
	if( in.LC!=ISOF ) {
		//Calculate the average scalar order parameter
		if(in.LC==LCG) *AVS = avOrderParam( SRDparticles,in.LC,avDIR );
		#ifdef DBG
			if( DBUG >= DBGTITLE ) printf( "Orientation Collision Step.\n" );
		#endif
		for( i=0; i<XYZ_P1[0]; i++ ) for( j=0; j<XYZ_P1[1]; j++ ) for( k=0; k<XYZ_P1[2]; k++ ) {
            if (in.MFPLAYERH > 0) { // if MFPLAYERH set then perform the "active layer" hack
                // Active layer
                // Check if cell is within the layer, if not then run with an effective MFP of 0
                if (j <= in.MFPLAYERH) zeroMFPot=0;
                else zeroMFPot = 1;
                if( CL[i][j][k].POPSRD > 1 ) LCcollision( &CL[i][j][k],SP,in.KBT,zeroMFPot,in.dt,*AVS,in.LC );
            } else {
                // otherwise perform typical LC collision
				zeroMFPot=0;
                // LC collision algorithm (no collision if only 1 particle in cell)
                if( CL[i][j][k].POPSRD > 1 ) LCcollision( &CL[i][j][k],SP,in.KBT,zeroMFPot,in.dt,*AVS,in.LC );
            }
		}
		// Magnetic alignment is really part of the collision
		#ifdef DBG
			if( DBUG >= DBGTITLE && in.MAG_FLAG) printf( "Apply Magnetic Torque Operator.\n" );
		#endif
		if( in.MAG_FLAG ) magTorque_all( SRDparticles,SP,in.dt,in.MAG );
		#ifdef DBG
			if( DBUG >= DBGTITLE ) printf( "Orientation Shear Alignment.\n" );
		#endif
		for( i=0; i<XYZ_P1[0]; i++ ) for( j=0; j<XYZ_P1[1]; j++ ) for( k=0; k<XYZ_P1[2]; k++ ) {
			//Coupling shear to orientation
			if( CL[i][j][k].POPSRD > 1 ) jefferysTorque( &CL[i][j][k],SP,in.dt );
		}
		#ifdef DBG
			if (DBUG == DBGTHERM) {
				*KBTNOW = TEMP( SRDparticles,SP,WALL,AVV );
				printf( "\tSystem Temperature post-LCcollision: %09e\n",*KBTNOW );
			}
		#endif
	}
	/* ****************************************** */
	/* *********** VELOCITY COLLISION *********** */
	/* ****************************************** */
	//Apply linear momentum collision operator
	#ifdef DBG
		if( DBUG >= DBGTITLE ) printf( "Velocity Collision Step.\n" );
	#endif
	#ifdef DBG
		if( DBUG == DBGSWIMMERDEETS ) {
			printf("\tPre-collision:\n");
			for( i=0; i<NS; i++ ) {
				printf( "\tS%d:\n",i );
				swcoord(swimmers[i]);
			}
		}
    #endif

    for (i = 0; i < NSPECI; i++) { // Temporarily store the intended activity for each species
        savedAct[i] = (double)(SP+i)->ACT;
    }
	for( i=0; i<XYZ_P1[0]; i++ ) for( j=0; j<XYZ_P1[1]; j++ ) for( k=0; k<XYZ_P1[2]; k++ ) {
		//MPC/SRD collision algorithm (no collision if only 1 particle in cell)
		CLQ[0]=i+0.5;
		CLQ[1]=j+0.5;
		CLQ[2]=k+0.5;

        // handle active layer logic
        if (in.MFPLAYERH > 0) {
            for (l = 0; l < NSPECI; l++) { // if not inside a layer, set activity to 0
                if (j <= in.MFPLAYERH) (SP+l)->ACT = savedAct[l];
                else (SP+l)->ACT = 0.0;
            }
        }

        if( CL[i][j][k].POP > 1 ) MPCcollision(&CL[i][j][k], SP, *SS, in.KBT, in.RTECH, in.C, in.S, in.FRICCO, in.dt, MD_mode, in.LC, in.TAU, CLQ, outPressure );
	}

	// Apply the multiphase interactions
	if( in.MULTIPHASE != MPHOFF && NSPECI>1 ) {
		for( i=0; i<XYZ_P1[0]; i++ ) for( j=0; j<XYZ_P1[1]; j++ ) for( k=0; k<XYZ_P1[2]; k++ ) if( CL[i][j][k].POP > 1 ) multiphaseColl(&CL[i][j][k], SP, *SS, in.MULTIPHASE, in.KBT, MD_mode, CLQ, outPressure );
	}
	// Apply the incompressibility correction
	if( in.inCOMP != INCOMPOFF ) {
		for( i=0; i<XYZ_P1[0]; i++ ) for( j=0; j<XYZ_P1[1]; j++ ) for( k=0; k<XYZ_P1[2]; k++ ) if( CL[i][j][k].POP > 1 ) incompColl(&CL[i][j][k], SP, *SS, in.inCOMP, MD_mode, CLQ, outPressure );
	}

    for (i = 0; i < NSPECI; i++) { // revert activity
        (SP+i)->ACT = savedAct[i];
    }


	// Apply the MD active dipoles
	if( MDmode && fabs(simMD->dStrength)>=TOL) {
		activeMD(simMD, CL, SP, in );
	}

	// Brownian thermostat (no hydrodynamic interactions -scramble velocities)
	if( in.noHI == HIOFF ) scramble( SRDparticles );
	//Calculate average
	avVel( CL,AVNOW );
	#ifdef DBG
		if( DBUG == DBGSWIMMERDEETS ) {
			printf("\tPost-collision:\n");
			for( i=0; i<NS; i++ ) {
				printf( "\tS%d:\n",i );
				swcoord(swimmers[i]);
			}
		}
	#endif
	/* ****************************************** */
	/* *********** TEMPERATURE SCALING ********** */
	/* ****************************************** */
	//Update post collision values
	*KBTNOW = TEMP( SRDparticles,SP,WALL,AVV );
	#ifdef DBG
		if (DBUG == DBGTHERM) printf( "\tSystem Temperature post-collision: %09e\n",*KBTNOW );
		if (DBUG >= DBGTITLE) printf( "Apply thermostat.\n" );
	#endif
	//Scale the termperature
	if( in.TSTECH != NOTHERM ) {
		scaleT( in.KBT,*KBTNOW,in.dt,in.TAU,AVV,AVNOW,in.TSTECH,SP,in.LC,WALL,SRDparticles,CL );
		*KBTNOW = TEMP( SRDparticles,SP,WALL,AVV );
		#ifdef DBG
			if (DBUG == DBGTHERM) printf( "\tSystem Temperature post-scaling: %09e\n",*KBTNOW );
		#endif
	}
	//Check for simulations with high activity at EARLY times activity since there is an early time tempreature spike
	if( in.RTECH >= VICSEK && in.RTECH <= DIPOLE_DIR_AV ) {
		//Check if temperature is VERY large
		if( *KBTNOW>50.0*in.KBT ) {
			#ifdef DBG
				if (DBUG >= DBGWARN) printf( "\tWarning active system too energetic: %09e\n",*KBTNOW );
			#endif
			scaleT( in.KBT,*KBTNOW,in.dt,in.TAU,AVV,AVNOW,VSC,SP,in.LC,WALL,SRDparticles,CL );
			*KBTNOW = TEMP( SRDparticles,SP,WALL,AVV );
			#ifdef DBG
				if (DBUG >= DBGWARN) printf( "\tVelocities rescaled: %09e\n",*KBTNOW );
			#endif
		}
	}

	/* *******************************************/
	/* * ADDING ORIENTATION PART TO THE VELOCITY */
	/* *******************************************/
	if (in.LC == BCT) {
		for (int i = 0; i < GPOP; i++) {
			for( j=0; j<DIM; j++ ) {
				(SRDparticles+i)->V[j] += (SRDparticles+i)->U[j]*((SP+0)->BS);
			}
		}
	}
	/* ****************************************** */
	/* ************ GRID SHIFT BACK ************* */
	/* ****************************************** */
	if( in.GALINV ) {
		#ifdef DBG
			if( DBUG >= DBGTITLE ) printf( "Shift Grid Back.\n" );
		#endif
		//Shift entire system by back
		gridShift_all( RSHIFT,1,SRDparticles,WALL,simMD,swimmers,MDmode );
	}
	/* ****************************************** */
	/* ********** APPLY SWIMMER DIPOLE ********** */
	/* ****************************************** */
	#ifdef DBG
		if( DBUG == DBGSWIMMER || DBUG == DBGSWIMMERDEETS ) printf( "\tApply swimmer dipole.\n" );
	#endif
	//Apply swimmer dipole (both force and torque dipoles)
	swimmerDipole( *SS,swimmers,CL,SP,in.dt,SRDparticles,WALL,simMD );
	//Allow streaming and boundary conditions before re-binning
	/* ****************************************** */
	/* ************ SRD TRANSLATION ************* */
	/* ****************************************** */
	#ifdef DBG
		if( DBUG >= DBGTITLE ) printf( "Translate MPC particles.\n" );
	#endif
	if( MDmode != MPCinMD ) {
		stream_all( SRDparticles,in.dt );
	}

	/* ****************************************** */
	/* ******************* BCs ****************** */
	/* ****************************************** */
	#ifdef DBG
		if( DBUG >= DBGTITLE ) printf( "Check MPCs Against BCs.\n" );
	#endif
	//Zero impulse on BCs --> moved on 27/01/21
	//for( i=0; i<NBC; i++ ) {
	//	zerovec(WALL[i].dV,DIM);
	//	zerovec(WALL[i].dL,_3D);
	//}
	//Check each particle
	bcCNT=0;
	reCNT=0;
	rethermCNT=0;
	for( i=0; i<GPOP; i++ ) MPC_BCcollision( SRDparticles,i,WALL,SP,in.KBT,in.dt,in.LC,&bcCNT,&reCNT,&rethermCNT,1 );
	#ifdef DBG
		if( DBUG == DBGBCCNT ) if(bcCNT>0) printf( "\t%d particles had difficulty with the BCs (%d rewind events; %d rethermalization events).\n",bcCNT,reCNT,rethermCNT );
	#endif
	/* ****************************************** */
	/* ************* APPLY IMPULSE ************** */
	/* ****************************************** */
	#ifdef DBG
		if( DBUG >= DBGTITLE ) printf( "Impulse on BCs from MPCD collisions.\n" );
	#endif
	//Apply impulse from MPC_BCcollision()

	if(!in.warmupSteps){
		for( i=0; i<NBC; i++ ) if( (WALL+i)->DSPLC ) {
			for( j=0; j<DIM; j++ ) (WALL+i)->V[j] += (WALL+i)->dV[j];
			for( j=0; j<_3D; j++ ) (WALL+i)->L[j] += (WALL+i)->dL[j];
			zerovec(WALL[i].dV,DIM);
			zerovec(WALL[i].dL,_3D);
		}
	/* ****************************************** */
	/* ************* TRANSLATE BCs ************** */
	/* ****************************************** */
	#ifdef DBG
		if( DBUG >= DBGTITLE ) printf( "Translate BCs.\n" );
	#endif
	//Save the old position in case a BC-BC collision occurs
	for( i=0; i<NBC; i++ ) if( (WALL+i)->DSPLC ) for( j=0; j<DIM; j++ ) {
		(WALL+i)->Q_old[j] = (WALL+i)->Q[j];
		(WALL+i)->O_old[j] = (WALL+i)->O[j];
	}
	//Translate each of the BCs --- using velocity from BEFORE MPC_BCcollision()
	for( i=0; i<NBC; i++ ) if( (WALL+i)->DSPLC ) stream_BC( (WALL+i),in.dt );
	/* ****************************************** */
	/* *************** ROTATE BCs *************** */
	/* ****************************************** */
	#ifdef DBG
		if( DBUG >= DBGTITLE ) printf( "Spin BCs.\n" );
	#endif
	for( i=0; i<NBC; i++ ) if( (WALL+i)->DSPLC ) spin_BC( (WALL+i),in.dt );
	/* ****************************************** */
	/* ***************** BC-BC ****************** */
	/* ****************************************** */
	#ifdef DBG
		if( DBUG >= DBGTITLE ) printf( "Check BCs Against BCs.\n" );
	#endif
	//Check each BC
	for( i=0; i<NBC; i++ ) if( (WALL+i)->DSPLC ) {
		BC_FLAG = 0;
		for( j=0; j<NBC; j++ ) if( j != i ) {
			//Check BC number i for collisions other BCs
			#ifdef DBG
				if( DBUG == DBGBCBC ) printf( "BC%d BC%d\n",i,j );
			#endif
			BC_BCcollision( WALL+i,WALL+j,in.dt,&BC_FLAG );
		}
	}
	/* ****************************************** */
	/* ***************** BC-MPCD **************** */
	/* ****************************************** */
	// if( BC_FLAG ) {
		#ifdef DBG
			if( DBUG >= DBGTITLE ) printf( "Check BCs Against MPCs after BC-BC collisions.\n" );
		#endif
		bcCNT=0;
		reCNT=0;
		rethermCNT=0;
		// Check each BC for collisions MPC particles
		for( i=0; i<NBC; i++ ) if( (WALL+i)->DSPLC ) {
			BC_MPCcollision( WALL,i,SRDparticles,SP,in.KBT,in.GRAV,in.dt,simMD,MDmode,in.LC,&bcCNT,&reCNT,&rethermCNT );
		}
		#ifdef DBG
			if( DBUG == DBGBCCNT ) if( bcCNT>0 ) printf( "\t%d particles had difficulty with the BCs when the BCs moved (%d rewind events; %d rethermalization events).\n",bcCNT,reCNT,rethermCNT );
		#endif
	// }
	/* ****************************************** */
	/* ************* APPLY IMPULSE ************** */
	/* ****************************************** */
	#ifdef DBG
		if( DBUG >= DBGTITLE ) printf( "Impulse on BCs from BC-translations.\n" );
	#endif
	//Apply impulse from BC_MPCcollision()
	for( i=0; i<NBC; i++ ) if( (WALL+i)->DSPLC ) {
		for( j=0; j<DIM; j++ ) (WALL+i)->V[j] += (WALL+i)->dV[j];
		//THERE SHOULD BE NO dL since BC_MPCcollision() ignores ang mom
		for( j=0; j<_3D; j++ ) (WALL+i)->L[j] += (WALL+i)->dL[j];
	}
	/* ****************************************** */
	/* ************* ACCELERATE BCs ************* */
	/* ****************************************** */
	#ifdef DBG
		if( DBUG >= DBGTITLE ) printf( "Accelerate BCs.\n" );
	#endif
	//Accelerate each of the BCs
	if( in.GRAV_FLAG ) for( i=0; i<NBC; i++ ) if( (WALL+i)->DSPLC ) acc_BC( (WALL+i),in.dt,(WALL+i)->G );
}
	/* ****************************************** */
	/* ***************** RE-BIN ***************** */
	/* ****************************************** */
	#ifdef DBG
		if( DBUG >= DBGTITLE ) printf( "Re-bin Particles.\n" );
	#endif
	// Bin SRD particles
	bin( CL,SP,WALL,in.KBT,in.LC,0 );
	// Bin swimmer monomers
	binSwimmers( CL,0 );
	// Bin MD particles
	if( MDmode ) binMD( CL );
	//Recalculate localPROP to ensure updated cell properties
	localPROP( CL,SP,*SS,in.RTECH,in.LC );
	/* ****************************************** */
	/* ********** SAVE SOME PROPERTIES ********** */
	/* ****************************************** */
	//Build the flow profile by summing over everytime step and averaging after FLOWOUT or SWFLOWOUT iterations
	if( (outFlags.FLOWOUT>=OUT || outFlags.SWFLOWOUT>=OUT) && !in.warmupSteps ) {
		localFLOW( CL,SP );
		//The 'lab frame' flowout
		if( outFlags.FLOWOUT>=OUT ) sumFLOW( CL );
		//The first swimmer 'frame' averaged swflowout
		if( outFlags.SWFLOWOUT>=OUT && NS>0 ) sumSWFLOW( CL, swimmers,SS);
	}
}

/// 
/// @brief This function calculates the pre-collisional part (ballistic/streaming part) of the pressure tensor.
///
/// The routine only calculates the contribution due to the MPCD fluid --- not MD or swimmers
/// The volume of the MPCD cell is 1 (in MPCD units).
/// It updates the particles by looping through the linked lists. 
/// @param CL All of the MPCD cells (including the linked list of particles in each cell). 
/// @param SP The species-wide information about MPCD particles.
///
void calcPressureStreaming( cell ***CL,spec *SP ) {
	int a,b,c,i,j,id;
	double V[DIM];
	double mass;
	particleMPC *pMPC;	//Temporary pointer to MPC particles

	// Calculate POP, MASS and VCM (don't calculate CM. Only if needed - see below)
	for( a=0; a<XYZ_P1[0]; a++ ) for( b=0; b<XYZ_P1[1]; b++ ) for( c=0; c<XYZ_P1[2]; c++ ) {
		//Zero everything for recounting
		for( i=0; i<DIM; i++ ) for( j=0; j<DIM; j++ ) CL[a][b][c].Ps[i][j] = 0.0;
		// SRD particles
		if( CL[a][b][c].pp!=NULL ) {
			pMPC = CL[a][b][c].pp;
			while(pMPC!=NULL) {
				id = pMPC->SPID;
				mass = (SP+id)->MASS;
				for( i=0; i<DIM; i++ ) V[i] =  pMPC->V[i];
				//Stress tensory
				for( i=0; i<DIM; i++ ) for( j=0; j<DIM; j++ ) CL[a][b][c].Ps[i][j] += mass*V[i]*V[j];
				//Increment link in list
				pMPC = pMPC->next;
			}
			//Stress is actually the negative (and divided by cell volume [unity]) of sum above
			for( i=0; i<DIM; i++ ) for( j=0; j<DIM; j++ ) CL[a][b][c].Ps[i][j] *= (-1.0);
		}
	}
}

/// 
/// @brief Zero-collisional pressure.
///
/// The zero-collisional term in the pressure, i.e. the kinetic contribution. 
/// It is divided by volume a^DIM but cell volume=1, so this is not included. 
/// Dividing by volume and time make the changes in momentum into pressure
/// - \f$P_{ij} = \frac{1}{V} \left\langle \delta r_i F_j \right\rangle = \frac{1}{V \ \delta t} \left\langle \delta r_i \Delta p_j \right\rangle \f$. 
/// - https://link.springer.com/chapter/10.1007/978-3-540-87706-6_1
/// @param CL An MPCD cell. 
/// @param dt The MPCD time step.
/// @see calcPressureColl_preColl()
/// @see calcPressureColl_postColl()
///
void normPressureColl( cell *CL,double dt ) {
	int i,j;
	for( i=0; i<DIM; i++ ) for( j=0; j<DIM; j++ ) CL->Pc[i][j] /= (-dt);
}

/// 
/// @brief The <b>pre</b>-collision calculations needed to calculate the collisional pressure term.
///
/// To calculate the pressure, need to know change in momentum. 
/// So record initial velocity before the MPCD collision and relative position from the centre of the cell. 
/// - \f$P_{ij} = \frac{1}{V} \left\langle \delta r_i F_j \right\rangle = \frac{1}{V \ \delta t} \left\langle \delta r_i \Delta p_j \right\rangle \f$. 
/// - https://link.springer.com/chapter/10.1007/978-3-540-87706-6_1
/// @param relQ The MPCD particle position relative to the geometric centre of the cell `CLQ`. The relative position is returned through this variable.
/// @param dp The change in momentum, which is being set to the initial velocity here in order to later find the difference. The change in momentum is returned through this variable.
/// @param p An MPCD particle. 
/// @param CLQ The geometric centre of `CL`, the MPCD cell.
/// @see calcPressureColl_postColl()
/// @see normPressureColl()
///
void calcPressureColl_preColl( double *relQ,double *dp,particleMPC *p,double *CLQ ) {
	int d;
	for( d=0; d<DIM; d++ ) {
		//Calculate the relative position from the cell centre
		relQ[d] = p->Q[d] - CLQ[d];
		//Prepare the change in momentum
		dp[d] = p->V[d];
	}
}

/// 
/// @brief The <b>post</b>-collision calculations needed to calculate the collisional pressure term.
///
/// To calculate the pressure, need to find in momentum after the collision operation. 
/// So record difference in velocity. 
/// Impulse is equivalent to force over the MPCD time step, where the time step is divided in normPressureColl(). 
/// - \f$P_{ij} = \frac{1}{V} \left\langle \delta r_i F_j \right\rangle = \frac{1}{V \ \delta t} \left\langle \delta r_i \Delta p_j \right\rangle \f$. 
/// - https://link.springer.com/chapter/10.1007/978-3-540-87706-6_1
/// @param relQ The MPCD particle position relative to the geometric centre of the cell `CLQ`. 
/// @param dp The change in momentum, which was previously initialized in calcPressureColl_preColl(). The change in momentum is returned through this variable.
/// @param M The MPCD particle mass.
/// @param vel The MPCD particle velocity.
/// @param CL An MPCD cell.  
///
void calcPressureColl_postColl( double *relQ,double *dp,double M,double *vel,cell *CL ) {
	int i,j;
	//Finish calculating the change in momentum
	for( i=0; i<DIM; i++ ) {
		dp[i] -= vel[i];
		dp[i] *= -M;
	}
	//Calculate the collisional contribution to the stress
	for( i=0; i<DIM; i++ ) for( j=0; j<DIM; j++ ) CL->Pc[i][j] += dp[i]*relQ[j];
}

///
/// @brief Check if a particular particle has any NaN values.
///
/// Checks is particle position, velocity or orientation are NaN. Primarily for debugging purposes. Further checks can
/// be implemented here if needed.
///
/// @param p A particle to check for NaNs.
///
void checkParticleNaN(particleMPC p) {
    // Check position
    if (isNaNs(p.Q, 3)) {
        printf("Particle position found to be NaN\n");
    }
    // Check velocity
    if (isNaNs(p.V, 3)) {
        printf("Particle velocity found to be NaN\n");
    }
    // Check orientation
    if (isNaNs(p.U, 3)) {
        printf("Particle orientation found to be NaN\n");
    }
}

///
/// @brief Check all particles to see if any particular particle has NaN values.
///
/// Checks all particles for NaN values. Calls on checkParticleNaN() to do the actual checking. Primarily for debugging
/// purposes.
///
/// @param p Array of particles to check for NaNs. Assumed to be the SRDparticles array.
/// @see checkParticleNaN()
///
void checkAllParticlesNaN(particleMPC *SRDparticles) {
    int i;
    for (i=0; i < GPOP; i++) {
        checkParticleNaN(SRDparticles[i]);
    }
}
