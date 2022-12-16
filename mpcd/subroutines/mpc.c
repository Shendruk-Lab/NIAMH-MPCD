# include<math.h>
# include<stdio.h>
# include<stdlib.h>
#ifdef _OPENMP
# include<omp.h>
#endif

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

# include "../../md/mdsrd.h"

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ****************** LISTS ***************** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
void localPROP( cell ***CL,spec *SP,specSwimmer specS,int RTECH,int LC ) {
/*
   This routine finds the local
   properties of all the cells.
*/
	int a,b,c,d,id;
	int i;
	double V[_3D],Q[_3D];
	double **S,eigval[_3D];
	double mass;
	particleMPC *pMPC;	//Temporary pointer to MPC particles
	particleMD *pMD;	//Temporary pointer to MD particles
	smono *pSW;			//Temporary pointer to swimmer monomers

	// int flag on whether to compute CM or not
	int computeCM = (RTECH==RAT) || (LC!=ISOF) || (RTECH==DIPOLE_VCM) || (RTECH==DIPOLE_DIR_SUM) || (RTECH==DIPOLE_DIR_AV) || (RTECH==MULTIPHASE);

	// Zero
	for( d=0; d<_3D; d++ ) V[d] = 0.0;
	for( d=0; d<_3D; d++ ) eigval[d] = 0.0;

	// Loop over all cells
	for( a=0; a<XYZ_P1[0]; a++ ) for( b=0; b<XYZ_P1[1]; b++ ) for( c=0; c<XYZ_P1[2]; c++ ) {
		// Idea here is to compute POP, MASS, and VCM
		//	CM is only computed if needed

		//Zero everything for recounting
		CL[a][b][c].POP = 0;
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
					CL[a][b][c].POP ++;
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
					CL[a][b][c].POP ++;
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
					CL[a][b][c].POP ++;
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
		S = malloc ( DIM * sizeof( *S ) );
		for( i=0; i<DIM; i++ ) S[i] = malloc ( DIM * sizeof( *S[i] ) );
		for( i=0; i<DIM; i++ ) for( d=0; d<DIM; d++ ) S[i][d] = 0.0;
		// Find the order parameter tensor, the director and the scalar order parameter for each cell
		for( a=0; a<XYZ_P1[0]; a++ ) for( b=0; b<XYZ_P1[1]; b++ ) for( c=0; c<XYZ_P1[2]; c++ ) {
			if( CL[a][b][c].POP > 1 ) {
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
void localFLOW( cell ***CL,spec *SP ) {
/*
   This routine calculates the local flow in each cell
*/
	int a,b,c;
	for( a=0; a<XYZ_P1[0]; a++ ) for( b=0; b<XYZ_P1[1]; b++ ) for( c=0; c<XYZ_P1[2]; c++ ) {
		localMPCVCM( CL[a][b][c].VCM,CL[a][b][c],SP ) ;
	}
}
void sumFLOW( cell ***CL ) {
/*
   This routine sums the local flow in each cell --- to be turned into an average in flowout()
*/
	int a,b,c,d;
	for( a=0; a<XYZ_P1[0]; a++ ) for( b=0; b<XYZ_P1[1]; b++ ) for( c=0; c<XYZ_P1[2]; c++ ) for( d=0; d<DIM; d++ ) {
		CL[a][b][c].FLOW[d] += CL[a][b][c].VCM[d];
	}
}

void ghostPart( cell ***CL,bc WALL[],double KBT,int LC, spec *SP) {
/*
   This routine adds the effect of phantom particleMPCs to the CM
   This would be oodles better if each BC object had a list of cells
   to worry about and so I wouldn't have to go over every cell all
   the time
	 This routine also strengthens anchoring by re-applying orientational boundary
	 conditions (now to the whole cell).  This is so the collision operator can
	 reassign orientations about the director (preferred by the anchoring), with
	 less deviation (as S = 1 for planar, or close to 1 otherwise).
*/

	int a,b,c,d,i,j,k;
	int numCorners = (int) smrtPow(2,DIM);
	double R[DIM];
	double invN;									// Inverse number difference
	particleMPC tp[numCorners];		// Temporary MPC particles for all corners of a square cell
	double W[numCorners];					// W for each of the corners
	double **S, eigval[_3D];
	double n[_3D] = {0.,0.,0.};		// surface normal
	particleMPC *ptMPC;						// temporary pointer to MPC particle
	int setGhostAnch, flagW;
	int numBC;										// Number of walls with anchoring in a given cell
	int wallindex;								// Index of wall with anchoring acting on a cell
	double shift[DIM];
	setGhostAnch = 1; 						// a manual switch to turn on=1 or off=0 the stronger anchoring


	if (setGhostAnch == 1){
		// Allocate memory for S
		S = malloc ( DIM * sizeof( *S ) );
		for( i=0; i<DIM; i++ ) S[i] = malloc ( DIM * sizeof( *S[i] ) );
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
							while (ptMPC != NULL){
								normal(n, WALL[i], ptMPC->Q, DIM);
								norm(n, DIM);
								// apply boundary condition (anchoring) to cell
								// note: this also feeds anchoring torque back to boundary (to add to dV and dL), if mobile
								oriBC(ptMPC, SP, &WALL[i], n);
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
						else if ( CL[a][b][c].POP == 1 ){
							CL[a][b][c].S = 1.0;
							// Set the director as the orientation of the particle
							for( k=0; k<DIM; k++ ) CL[a][b][c].DIR[k] = CL[a][b][c].pp->U[k];
							norm( CL[a][b][c].DIR, DIM );
						}

						// Curved wall
						else if (CL[a][b][c].POP > 1){
							// Calculate Q tensor
							tensOrderParam( &CL[a][b][c], S, LC );
							// Find S
							solveEigensystem( S, DIM, eigval );
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
							if( CL[a][b][c].S > 1.0 ){
								CL[a][b][c].S = 1.0;
							}
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
				if( CL[a][b][c].POP < nDNST && CL[a][b][c].POP > 0 ) {
					invN = 1.0/(nDNST-(double)CL[a][b][c].POP);
					for( d=0; d<DIM; d++ ) {
						R[d] = genrand_gaussMB( KBT, invN );
						CL[a][b][c].VCM[d] += R[d];
						CL[a][b][c].VCM[d] /= nDNST;
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

					// // Same as above but doing EACH phantom particle separately
					// for( d=0; d<DIM; d++ ) CL[a][b][c].VCM[d] *= (double)CL[a][b][c].POP;
					// for( j=0; j<(int)(nDNST-(double)CL[a][b][c].POP); j++ ) {
					// 	for( d=0; d<DIM; d++ ) {
					// 		//Generate random ghost velocity
					// 		R[d] = genrand_gaussMB( KBT,1.0 );
					// 		CL[a][b][c].VCM[d] += R[d];
					// 	}
					// }
					// for( d=0; d<DIM; d++ ) CL[a][b][c].VCM[d] /= nDNST;

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* **************** STREAMING *************** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
double trans( double t,double V, double QOLD ) {
/*
   This subroutine translates the particleMPCs'
   positions
   NO ACCELERATION DURING TIME STEP! See papers on SRD
*/
	double QNEW;
// 	QNEW = QOLD + t * V + 0.5*t*t*G;
	QNEW = QOLD + t * V;
	return QNEW;
}
double acc( double t,double GRAV,double VOLD ) {
/*
   This subroutine accelerates the particleMPCs'
   velocity
*/
	double VNEW;
	VNEW = VOLD + t * GRAV;
	return VNEW;
}
void stream_BC( bc *WALL,double t ) {
/*
    The streaming step of the algorithm translates
    position
*/
	int i;
	for( i=0; i<DIM; i++ ) WALL->Q[i] = trans( t,WALL->V[i],WALL->Q[i] );
}
void spin_BC( bc *WALL,double t ) {
/*
    The streaming rotational step of the algorithm rotates
    orienation and accelerates the velocity
*/
	int i;
	for( i=0; i<_3D	; i++ ) WALL->O[i] += t * WALL->L[i];
}
void stream_P( particleMPC *p,double t ) {
/*
    The streaming step of the algorithm translates position
*/
	int i;
	for( i=0; i<DIM; i++ ) p->Q[i] = trans( t,p->V[i],p->Q[i] );
}
void acc_BC( bc *WALL,double t,double GRAV[] ) {
/*
    The streaming step of the algorithm translates
    position and accelerates the velocity
*/
	int i;
	for( i=0; i<DIM; i++ ) WALL->V[i] = acc( t,GRAV[i],WALL->V[i] );
}
void acc_P( particleMPC *p,double t,double GRAV[] ) {
/*
    The streaming step of the algorithm translates
    position and accelerates the velocity
*/
	int i;
	for( i=0; i<DIM; i++ ) p->V[i] = acc( t,GRAV[i],p->V[i] );
}
void stream_all( particleMPC *pp,double t ) {
/*
    The streaming step of the algorithm translates the
    positions and accelerates the velocities
*/
	int i;
    ///TODO: OMP here
	for( i=0; i<GPOP; i++ ){
		//Update coordinates --- check if it already streamed
		if( (pp+i)->S_flag ) stream_P( (pp+i),t );
		else (pp+i)->S_flag = STREAM;
	}
}
void acc_all( particleMPC *pp,double t,double GRAV[] ) {
/*
    The streaming step of the algorithm translates the
    positions and accelerates the velocities
*/
	int i;
	for( i=0; i<GPOP; i++ ) acc_P( (pp+i),t,GRAV );
}
void gridShift_all( double SHIFT[],int shiftBack,particleMPC *SRDparticles,bc WALL[],simptr simMD,swimmer swimmers[],int MDmode ) {
	/*
	    Shifts the entire system by the vector SHIFT
			If shiftBack then multiply shift by -1 and shift. Else do the normal shift
	*/
	int i,j;							//Counting variables
	double signedSHIFT[_3D];		//

	if( shiftBack ) for( j=0; j<DIM; j++ ) signedSHIFT[j] = -1.0*SHIFT[j];
	else for( j=0; j<DIM; j++ ) signedSHIFT[j] = SHIFT[j];

	//Shift each particle
	for( i=0; i<GPOP; i++ ) for( j=0; j<DIM; j++ ){
		 SRDparticles[i].Q[j] += signedSHIFT[j]; // perform the shift

		 // The Gallilean shift can (sometimes) force particles outside the boundaries of the system.
		 // 	These lines will periodically shift those particles to the other side of the domain, while
		 //	keeping those that escaped naturally (ie not due to the shift) in place. It will also return
		 //	them to the correct place when needing to shift back for streaming.
		 if (shiftBack == 0 && SRDparticles[i].Q[j] > XYZ[j] && SRDparticles[i].Q[j] <= XYZ[j] + signedSHIFT[j] )	SRDparticles[i].Q[j] -= XYZ[j];
		 else if (shiftBack == 1 && SRDparticles[i].Q[j] < 0.0 && SRDparticles[i].Q[j] >= signedSHIFT[j] ) SRDparticles[i].Q[j] += XYZ[j];
	 }
	if( MDmode == MDinMPC) for( i=0; i<(simMD->atom.n); i++ ) {
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
void rotate_CL( cell CL,spec *SP,double r0[],double n0[],double dw ) {
/*
   This routine applies a solid body rotation to the cell i.e.
	 change in angular speed dw to all the SRD particles in a cell
	 about a given point r0 (likely the CM but allowed to be anything) and direction n0
*/
	int d=0;
	double r[DIM],r_perp[DIM],n_perp[DIM]; 		//Vectors betweem cm and mpcd particles and axis
	double n_v[DIM];													//Direction of velocity change
	double dist,dv;														//Mag of angular velocity on each MPCD particle, distance to line, mag of change in velocity
	particleMPC *pMPC;												//Temporary pointer to MPC particles
	double netMom[DIM],tempAngMom[DIM],netAngMom[DIM];				//Sum momentum and angular momentum
	int cnt;																	//Count number of SRD particles

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
double rewind_trans( double t,double V,double P ) {
/*
     Rewinds a translation
     NO ACCELERATION DURING TIME STEP! See papers on SRD
*/
	double QOLD;
	QOLD = P - t*V;
	return QOLD;
}
double rewind_acc(double t,double G,double V){
/*
     Rewinds the acceleration
*/
	double VOLD;
	VOLD = V - t*G;
	return VOLD;
}
void rewind_P( particleMPC *p,double time ) {
/*
     Bring the particleMPC back in time step.
*/
	int i;
	for( i=0; i<DIM; i++ ) p->Q[i] = rewind_trans(time,p->V[i],p->Q[i]);
}
void rewind_BC( bc *WALL,double time ) {
/*
     Bring the BC back in time step.
*/
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
void binin( particleMPC p[],cell ***CL ) {
/*
   This function does the initial binning of the
   particleMPCs after they have been first initialized.
   It is different from bin in that it uses the
   actual array of particleMPCs rather than the array
   of pointers to particleMPCs.
*/
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
void bin( cell ***CL,spec *SP,bc WALL[],double KBT,int LC,int shifted ) {
/*
   This function bins the particleMPCs i.e. it places
   a pointer to the particleMPC in the appropriate new
   list and removes it from it's old list.
   This does not calculate the local
   properties. That must be done separately
   by calling localPROP().
	 If the positions have been randomly shifted (shifted=1) then
	 do PBC wrap, otherwise (shifted=0) don't.
*/
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
void bininMD( simptr sim,cell ***CL ) {
/*
   This function bins the MD particles for use in the collision steps
*/
	int i,a,b,c;
	int nAtom;		// Number of atoms
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
void binMD( cell ***CL ) {
/*
   This function bins the particleMDs i.e. it places
   a pointer to the particleMD in the appropriate new
   list and removes it from it's old list
*/
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
void addlink( cell *CL,particleMPC *p ) {
/*
	This routine adds a link to the end of the list
*/
	particleMPC *tp;	//Temporary pointer to particleMPC

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
void removelink( particleMPC *current,cell *CL ) {
/*
	This routine removes a link from
	a list and relinks the list
*/
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
void addlinkMD( cell *CL,particleMD *p ) {
/*
	This routine adds a link to the end of the list
*/
	particleMD *tp;		//Temporary pointer to particleMPC

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
void removelinkMD( particleMD *current,cell *CL ) {
/*
	This routine removes a link from
	a list and relinks the list
*/
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
void rotate( int RT,double Cos,double Sin,double V[_3D],double L[_3D],long SIGN,int RAND ) {
/*
   Takes in a vector V and rotates it about
   a random cartesian axis (if RT=ORTHAXIS)
   or about a random axis L (if RT=ARBAXIS).
   The arbitary axis is chosen randomly for
   every collision cell and each time step
*/
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
void stochrotMPC( cell *CL,int RTECH,double C,double S,double *CLQ,int outP ) {
/*
    Does the stochastic rotation collision, and returns the
    CM velocity and the local temperature of the cell.
*/
	int i,j;
	long SIGN;		//SIGN indicates sign of RA
	int CA = 0;		//indicates random cartesian axis chosen
	double RV[_3D] = {0.};	//Random vector
	double V[_3D];
	double dp[CL->POP][DIM],relQ[CL->POP][DIM];		//For pressure
	particleMPC *tmpc;	//Temporary particleMPC
	particleMD *tmd;	//Temporary particleMD
	smono *tsm;			//Temporary swimmer monomer

	//Generate random axis(used if RTECH=ARBAXIS)
	if( RTECH == ARBAXIS ) ranvec( RV,DIM );
	//Randomly pick an axis(used if RTECH=ORTHAXIS)
	if( RTECH == ORTHAXIS || RTECH == NOHI_ARBAXIS ) CA = (int)DIM * genrand_real();

	//Randomly pick a sign for the angle
	SIGN = genrand_real();

	//Collision
	//MPC particles
	i=0;
	tmpc = CL->pp;
	while( tmpc!=NULL ) {
		//Pressure term
		if( outP ) calcPressureColl_preColl( relQ[i],dp[i],tmpc,CLQ );
		for( j=0; j<DIM; j++ ) V[j] = tmpc->V[j] - CL->VCM[j];
		if( RTECH == NOHI_ARBAXIS ) rotate( ORTHAXIS,C,S,V,RV,SIGN,CA );
		else rotate( RTECH,C,S,V,RV,SIGN,CA );
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
		if( RTECH == NOHI_ARBAXIS ) rotate( ORTHAXIS,C,S,V,RV,SIGN,CA );
		else rotate( RTECH,C,S,V,RV,SIGN,CA );
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
		if( RTECH == NOHI_ARBAXIS ) rotate( ORTHAXIS,C,S,V,RV,SIGN,CA );
		else rotate( RTECH,C,S,V,RV,SIGN,CA );
		for( i=0; i<DIM; i++ ) tsm->V[i] = CL->VCM[i] + V[i];
		//Increment link in list
		tsm = tsm->next;
	}
}
void andersenMPC( cell *CL,spec *SP,specSwimmer SS,double KBT,double *CLQ,int outP ) {
/*
    Does the Andersen thermostat collision, and returns the
    CM velocity and the local temperature of the cell.
*/
	int i,j,id;
	double MASS;
	double RV[CL->POP][DIM];	//Random velocities
	double RS[DIM];		//Sum of random velocities
	double DV[CL->POP][DIM];	//Damping velocity
	double dp[CL->POP][DIM],relQ[CL->POP][DIM];		//For pressure
	particleMPC *tmpc;	//Temporary particleMPC
	particleMD *tmd;	//Temporary particleMD
	smono *tsm;			//Temporary swimmer monomer

	// Zero arrays
	for( i=0; i<DIM; i++ ) RS[i] = 0.;
	for( i=0;i<CL->POP;i++ ) for( j=0;j<DIM;j++ ) {
		RV[i][j] = 0.0;
		DV[i][j] = 0.0;
	}

	//Generate random velocities
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

	//Collision
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
		tmd->vy = CL->VCM[1] + RV[i][1] - RS[1];
		tmd->vz = CL->VCM[2] + RV[i][2] - RS[2];
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
void andersenMULTIPHASE( cell *CL,spec *SP,specSwimmer SS,double KBT,double *CLQ,int outP ) {
/*
    Does the Andersen thermostat collision, and returns the
    CM velocity and the local temperature of the cell.
*/
  int i,j,k,id;
	int mixedCell=0;
  double MASS,MASSsum;
  double RV[CL->POP][DIM];        //Random velocities
	double N,NSP[NSPECI];          //Number of each type
  double RVSPsum[NSPECI][DIM]; //Sum of random velocities of each type
  double RVsum[DIM];              //Sum of random velocities that aren't A or B-type (monomers/swimmers)

  double DV[CL->POP][DIM];        //Damping velocity
  particleMPC *tmpc;              //Temporary particleMPC
  particleMD *tmd;                //Temporary particleMD
  smono *tsm;                     //Temporary swimmer monomer

  double relQ[DIM];               //Relative position
  double VMUtot[DIM];             //Relative position
  double gradSP[NSPECI][DIM];     //Directional gradient of each species
	double thisGrad;								//A temporary gradient contribution component
  double VMU[CL->POP][DIM];       //Grad. chemical potential  velocity of type A (B is negative this)

  // Zero arrays
  for( i=0; i<DIM; i++ ) {
    RVsum[i] = 0.0;
    relQ[i] = 0.0;
		for( j=0; j<NSPECI; j++ ) {
			RVSPsum[j][i] = 0.0;
	    gradSP[j][i] = 0.0;
		}
  }
  for( i=0;i<CL->POP;i++ ) for( j=0;j<DIM;j++ ) {
    RV[i][j] = 0.0;
    VMU[i][j] = 0.0;
    VMUtot[j] = 0.;             //Relative position
    DV[i][j] = 0.0;
  }
	for( j=0; j<NSPECI; j++ ) NSP[j]=0.0;
	MASSsum=0.0;

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
      // MASS = (double) (SP+id)->MASS;
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
			if( tsm->HorM ) {
				id = SS.MSPid;
				// MASS = (double) SS.middM;
			}
			else {
				id = SS.HSPid;
				// MASS = (double) SS.headM;
			}
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
		// MASS = (double) tmd->mass;
		for( j=0; j<DIM; j++ ) {
			for( k=0; k<NSPECI; k++ ) VMU[i][j] += gradSP[k][j]*((SP+id)->M[k]) / NSP[id];
			VMUtot[j] += VMU[i][j];
		}
		//Increment link in list
		tmd = tmd->nextSRD;
		i++;
	}

  //Generate random velocities
  i=0;
  //MPC particles
  tmpc = CL->pp;
  while( tmpc!=NULL ) {
    id = tmpc->SPID;
    MASS = (double) (SP+id)->MASS;
		MASSsum+=MASS;
    for( j=0; j<DIM; j++ ) {
      RV[i][j] = genrand_gaussMB( KBT,MASS );
      // RVsum[j] += MASS*RV[i][j];
			RVsum[j] += RV[i][j];
			for( k=0; k<NSPECI; k++ ) RVSPsum[k][j] += RV[i][j];
      DV[i][j] = ((SP+id)->DAMP)*(CL->VCM[j])/((double)CL->POP);
    }
		//Increment link in list
    tmpc = tmpc->next;
    i++;
  }
  //Swimmer monomers
	tsm = CL->sp;
	while( tsm!=NULL ) {
		if( tsm->HorM ) {
			id = SS.MSPid;
			MASS = (double) SS.middM;
		}
		else {
			id = SS.HSPid;
			MASS = (double) SS.headM;
		}
		MASSsum+=MASS;
		for( j=0; j<DIM; j++ ) {
      RV[i][j] = genrand_gaussMB( KBT,MASS );
      // RVsum[j] += MASS*RV[i][j];
			RVsum[j] += RV[i][j];
			for( k=0; k<NSPECI; k++ ) RVSPsum[k][j] += RV[i][j];
			//No dampening on the swimmer itself
      DV[i][j] = 0.0;
    }
		//Increment link in list
		tsm = tsm->next;
		i++;
	}
  //MD particles --- ALWAYS type 0
	id=0;
  tmd = CL->MDpp;
  while( tmd!=NULL ) {
    MASS = tmd->mass;
    MASSsum+=MASS;
		for( j=0; j<DIM; j++ ) {
			RV[i][j] = genrand_gaussMB( KBT,MASS );
			// RVsum[j] += MASS*RV[i][j];
			RVsum[j] += RV[i][j];
			for( k=0; k<NSPECI; k++ ) RVSPsum[k][j] += RV[i][j];
			//No dampening on the swimmer itself
			DV[i][j] = 0.0;
		}
		//Increment link in list
    tmd = tmd->nextSRD;
    i++;
  }
  //Turn sums into averages
  // if( MASSsum>0.0 ) for( j=0; j<DIM; j++ ) RVsum[j] /= MASSsum;
	if( MASSsum>0.0 ) for( j=0; j<DIM; j++ ) RVsum[j] /= N;
	for( k=0; k<NSPECI; k++ ) if( NSP[k]>0.0 ) for( j=0; j<DIM; j++ ) RVSPsum[k][j] /= NSP[k];
	for( j=0; j<DIM; j++ ) VMUtot[j] /= N;

  //Collision
  i=0;
  // MPC particles
  tmpc = CL->pp;
  while( tmpc!=NULL ) {
    id = tmpc->SPID;
    // MASS = (double)(SP+id)->MASS;
    // for( j=0; j<DIM; j++ ) tmpc->V[j] = CL->VCM[j] + RV[i][j] - RVsum[j] - DV[i][j] - RVSPsum[id][j] + VMU[i][j] - VMUtot[j];
		for( j=0; j<DIM; j++ ) tmpc->V[j] = CL->VCM[j] + RV[i][j] - DV[i][j] - RVSPsum[id][j] + VMU[i][j] - VMUtot[j];
    //Increment link in list
    tmpc = tmpc->next;
    i++;
  }
  // Swimmer monomers
  tsm = CL->sp;
  while( tsm!=NULL ) {
		if( tsm->HorM ) {
			id = SS.MSPid;
			// MASS = (double) SS.middM;
		}
		else {
			id = SS.HSPid;
			// MASS = (double) SS.headM;
		}
		// for( j=0; j<DIM; j++ ) tsm->V[j] = CL->VCM[j] + RV[i][j] - RVsum[j] - DV[i][j] - RVSPsum[id][j] + VMU[i][j] - VMUtot[j];
		for( j=0; j<DIM; j++ ) tsm->V[j] = CL->VCM[j] + RV[i][j] - DV[i][j] - RVSPsum[id][j] + VMU[i][j] - VMUtot[j];
    //Increment link in list
    tsm = tsm->next;
    i++;
  }
  //MD particles --- ALWAYS type 0
	id=0;
  tmd = CL->MDpp;
  while( tmd!=NULL ) {
		// MASS = tmd->mass;
		if( DIM>=_1D ) {
			j=0;
			// tmd->vx = CL->VCM[j] + RV[i][j] - RVsum[j] - DV[i][j] - RVSPsum[id][j] + VMU[i][j] - VMUtot[j];
			tmd->vx = CL->VCM[j] + RV[i][j] - DV[i][j] - RVSPsum[id][j] + VMU[i][j] - VMUtot[j];
		}
		if( DIM>=_2D ) {
			j=1;
			// tmd->vy = CL->VCM[j] + RV[i][j] - RVsum[j] - DV[i][j] - RVSPsum[id][j] + VMU[i][j] - VMUtot[j];
			tmd->vy = CL->VCM[j] + RV[i][j] - DV[i][j] - RVSPsum[id][j] + VMU[i][j] - VMUtot[j];
		}
		if( DIM>=_3D ) {
			j=2;
			// tmd->vz = CL->VCM[j] + RV[i][j] - RVsum[j] - DV[i][j] - RVSPsum[id][j] + VMU[i][j] - VMUtot[j];
			tmd->vz = CL->VCM[j] + RV[i][j] - DV[i][j] - RVSPsum[id][j] + VMU[i][j] - VMUtot[j];
		}
    //Increment link in list
    tmd = tmd->nextSRD;
    i++;
  }
}
void andersenROT( cell *CL,spec *SP,specSwimmer SS,double KBT,double *CLQ,int outP ) {
/*
    MPC collision that conserves angular momentum (uses andersen thermostat), and returns the
    CM velocity and the local temperature of the cell.
*/
	int i,j,id;
	double MASS;
	double RV[CL->POP][_3D];	//Random velocities
	double RS[_3D];		//Sum of random velocities
	double DV[CL->POP][DIM];	//Damping velocity
	double relQ[CL->POP][_3D];	//Relative position
	double diffV[_3D];	//Difference in velocity
	double L[_3D];		//Angular momentum
	double angterm[_3D];
	double angmom[_3D];
	double W[_3D];
	double VCM[_3D];
	double II[_3D][_3D];	//Inverse of moment of inertia tensor (3D)
	double dp[CL->POP][DIM],relQP[CL->POP][DIM];		//For pressure
	particleMPC *tmpc;	//Temporary particleMPC
	particleMD *tmd;	//Temporary particleMD
	smono *tsm;			//Temporary swimmer monomer

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
		for( j=0; j<DIM; j++ ) DV[i][j] = ((SP+id)->DAMP)*(CL->VCM[j])/((double)CL->POP);
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
// 	dotprodmat( L,II,W,_3D );
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
void langevinMPC( cell *CL,spec *SP,specSwimmer SS,double KBT,double FRICCO,double Step,double *CLQ,int outP ) {
/*
    Does the Langevin thermostat collision, and returns the
    CM velocity and the local temperature of the cell.
*/
	int i,j,id;
	double MASS;
	double WN[CL->POP][DIM];	//White noise
	double WNS[DIM];	//Sum of white noise
	double DV[CL->POP][DIM];	//Damping velocity
	double a,b;
	double VCM[DIM];	//Centre of mass velocity
	double dp[CL->POP][DIM],relQ[CL->POP][DIM];		//For pressure
	particleMPC *tmpc;	//Temporary particleMPC
	particleMD *tmd;	//Temporary particleMD
	smono *tsm;			//Temporary swimmer monomer

	for( i=0; i<DIM; i++ ) VCM[i] = CL->VCM[i];

	//Generate random velocities
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

	//Collision
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
		tmd->vy = VCM[1] + a * (tmd->vy-VCM[1]) + b * (WN[i][1] - WNS[1]);
		tmd->vz = VCM[2] + a * (tmd->vz-VCM[2]) + b * (WN[i][2] - WNS[2]);
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
void langevinROT( cell *CL,spec *SP,specSwimmer SS,double KBT,double FRICCO,double Step,double *CLQ,int outP ) {
/*
	Langevin MPC collision that conserves angular momentum (uses andersen thermostat), and returns the
	CM velocity and the local temperature of the cell.
*/
	int i,j,id;
	double MASS;
	double WN[CL->POP][DIM];	//White noise
	double WNS[DIM];	//Sum of white noise
	double DV[CL->POP][DIM];	//Damping velocity
	double a,b;
	double VCM[DIM];	//Centre of mass velocity
	double dp[CL->POP][DIM];		//For pressure
	double relQ[CL->POP][_3D];	//Relative position
	double diffV[_3D];	//Difference in velocity
	double L[_3D];		//Angular momentum
	double angterm[_3D];
	double angmom[_3D];
	double W[_3D];
	double II[_3D][_3D];	//Inverse of moment of inertia tensor (3D)
	particleMPC *tmpc;	//Temporary particleMPC
	particleMD *tmd;	//Temporary particleMD
	smono *tsm;			//Temporary swimmer monomer

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
		for( j=0; j<DIM; j++ ) tmpc->V[j] = VCM[j] + a * (tmpc->V[j]-VCM[j]) + b * (WN[i][j] - WNS[j]) - DV[i][j] + angterm[j];;
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
void activeSRD( cell *CL,spec *SP,int RTECH,double C,double S,double *CLQ,int outP ) {
	/*
	 Does the stochastic rotation collision but then rotates the
	 resulting velocities towards the centre of mass velocity
	 injecting momentum but keeping the energy constant.
	 */
	int i;
	double theta;			// Angle between individual velocity and vcm
	double sumACT=0.0;		// Sum of the activity of particles in the cell
	double nv[_3D]={0.0};		// normal vector between vcm and individual vel --- must be _3DD regardless of actual dimension
	double VCM[_3D]={0.0};	// Short-hand CM velocity --- must be 3D regardless of actual dimension
	particleMPC *tmpc;		//Temporary particleMPC

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
void vicsek( cell *CL,spec *SP,double *CLQ,int outP ) {
	/*
	This is meant to be the Vicsek algorithm BUT instead of an
	interaction radius the alignment occurs within an MPCD cell
	Here ACT serves as the NOISE range (eta in Vicsek)
	*/
	int i,j,id;
	double sigma;			//The range on the noise is (1-ACT) and is between [0,1]
	double noiseRange,speed;
	double VCM[_3D],VCMdir[_3D];	// Short-hand CM velocity
	double MASS,dp[CL->POP][DIM],relQ[CL->POP][DIM];		//For pressure
	particleMPC *tmpc;		//Temporary particleMPC

	//Must be 3D because genrand_cone needs 3D vectors even in 2D
	for( i=0; i<DIM; i++ ) VCM[i] = CL->VCM[i];
	for( i=DIM; i<_3D; i++ ) VCM[i] = 0.;
	normCopy( VCM,VCMdir,_3D );

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
void chate( cell *CL,spec *SP,double *CLQ,int outP ) {
	/*
	This is meant to be the Vicsek algorithm BUT instead of an
	interaction radius the alignment occurs within an MPCD cell
	Here ACT serves as the NOISE range (eta in Vicsek)
	*/
	int i,j,id;
	double sigma;			//The range on the noise is (1-ACT) and is between [0,1]
	double noiseRange,speed,pmOne;
	double DIR[_3D];		// Short-hand cell's director
	double MASS,dp[CL->POP][DIM],relQ[CL->POP][DIM];		//For pressure
	particleMPC *tmpc;		//Temporary particleMPC

	//Must be 3D because genrand_cone needs 3D vectors even in 2D
	for( i=0; i<DIM; i++ ) DIR[i] = CL->DIR[i];
	for( i=DIM; i<_3D; i++ ) DIR[i] = 0.;

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
void vicsekAndersenMPC( cell *CL,spec *SP,double KBT,double RELAX,double *CLQ,int outP ) {
	/*
	 Just like the andersen version of MPCD but it relaxes to the speed ACT
	 in the direction of the centre of mass velocity instead of the CMV.
	 */
	int i,j,id;
	double MASS,ACT;
	double RV[CL->POP][DIM];	//Random velocities
	double RS[DIM];			//Sum of random velocities
	double dp[CL->POP][DIM],relQ[CL->POP][DIM];		//For pressure
	particleMPC *tmpc;		//Temporary particleMPC
	particleMD *tmd;		//Temporary particleMD
	double VCM[DIM],VCMdir[DIM];	// Short-hand velocity of CM and its direction

	//The direction of the centre of mass velocity of the cell
	for( i=0; i<DIM; i++ ) VCM[i] = CL->VCM[i];
	normCopy( VCM,VCMdir,DIM );

	// Zero arrays
	for( i=0; i<DIM; i++ ) RS[i] = 0.;
	for( i=0;i<CL->POP;i++ ) for( j=0;j<DIM;j++ ) RV[i][j] = 0.;

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

	//Collision
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
		tmd->vy = VCM[1] + RV[i][1] - RS[1];
		tmd->vz = VCM[2] + RV[i][2] - RS[2];
		//Increment link in list
		tmd = tmd->nextSRD;
		i++;
	}
}
void chateAndersenMPC( cell *CL,spec *SP,double KBT,double RELAX,double *CLQ,int outP ) {
	/*
	 Just like the andersen version of MPCD but the noise is about the director (time activity)
	 instead of the centre of mass velocity.
	 */
	int i,j,id;
	double MASS,ACT,pmOne;
	double RV[CL->POP][DIM];	//Random velocities
	double RS[DIM];			//Sum of random velocities
	double AV[CL->POP][DIM];	//Active velocities
	double AS[DIM];			//Sum of active velocities
	double dp[CL->POP][DIM],relQ[CL->POP][DIM];		//For pressure
	particleMPC *tmpc;		//Temporary particleMPC
	particleMD *tmd;		//Temporary particleMD
	double VCM[DIM],speed;		// Short-hand velocity of CM
	double DIR[DIM];		//Order parameter and director

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

	//Collision
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
		tmd->vy = VCM[1] + RV[i][1] - RS[1];
		tmd->vz = VCM[2] + RV[i][2] - RS[2];
		//Increment link in list
		tmd = tmd->nextSRD;
		i++;
	}
}
void vicsekLangevinMPC( cell *CL,spec *SP,double KBT,double FRICCO,double Step,double RELAX,double *CLQ,int outP ) {
	/*
	 Does the Langevin thermostat collision, and returns the
	 CM velocity and the local temperature of the cell.
	 */
	int i,j,id;
	double MASS,ACT;
	double WN[CL->POP][DIM];	//White noise
	double WNS[DIM];	//Sum of white noise
	double a,b;
	double VCM[DIM],VCMdir[DIM];	//Centre of mass velocity
	double dp[CL->POP][DIM],relQ[CL->POP][DIM];		//For pressure
	particleMPC *tmpc;	//Temporary particleMPC
	particleMD *tmd;	//Temporary particleMD

	for( i=0; i<DIM; i++ ) VCM[i] = CL->VCM[i];
	normCopy( VCM,VCMdir,DIM );

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

	//Collision
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
		tmd->vy = VCM[1] + a * (tmd->vy-VCM[1]) + b * (WN[i][1] - WNS[1]);
		tmd->vz = VCM[2] + a * (tmd->vz-VCM[2]) + b * (WN[i][2] - WNS[2]);
		//Increment link in list
		tmd = tmd->nextSRD;
		i++;
	}
}
void chateLangevinMPC( cell *CL,spec *SP,double KBT,double FRICCO,double Step,double RELAX,double *CLQ,int outP ) {
	/*
	 Does the Langevin thermostat collision, and returns the
	 CM velocity and the local temperature of the cell.
	 */
	int i,j,id;
	double MASS,ACT,pmOne;
	double WN[CL->POP][DIM];	//White noise
	double WNS[DIM];		//Sum of white noise
	double a,b;
	double dp[CL->POP][DIM],relQ[CL->POP][DIM];		//For pressure
	particleMPC *tmpc;		//Temporary particleMPC
	particleMD *tmd;		//Temporary particleMD
	double VCM[DIM],speed;		//Centre of mass velocity
	double DIR[DIM];		//Order parameter and director

	//Calculate the director of the cell
	for( i=0; i<DIM; i++ ) {
		VCM[i] = CL->VCM[i];
		DIR[i] = CL->DIR[i];
	}
	speed=length(VCM,DIM);

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

	//Collision
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
		tmd->vy = VCM[1] + a * (tmd->vy-VCM[1]) + b * (WN[i][1] - WNS[1]);
		tmd->vz = VCM[2] + a * (tmd->vz-VCM[2]) + b * (WN[i][2] - WNS[2]);
		//Increment link in list
		tmd = tmd->nextSRD;
		i++;
	}
}
void dipoleAndersenMPC( cell *CL,spec *SP,double KBT,double RELAX,double *CLQ,int outP ) {
/*
    Does the Andersen thermostat collision, and returns the
    CM velocity and the local temperature of the cell.
    But also applies a force dipole on the cell by applying
    a kick along the direction of centre of mass velocity VCM
    (or director DIR) to all particles. Those above the plane
    defined by this direction and the cnetre of mass point CM
    get a positive kick, while those below get a negative kick.
*/
	int i,j,id;
	double MASS,M,ACT,pmOne;
	double RV[CL->POP][DIM];	//Random velocities
	double RS[DIM];			//Sum of random velocities
	double DV[CL->POP][DIM];	//Damping velocities
	double AV[CL->POP][DIM];	//Active velocities
	double AS[DIM];			//Sum of active velocities
	double dp[CL->POP][DIM],relQ[CL->POP][DIM];		//For pressure
	particleMPC *tmpc;		//Temporary particleMPC
	particleMD *tmd;		//Temporary particleMD
	bc PLANE;			//The plane that cuts the cell in half
	double W;			//The particle's W for passing the plane
	double speed;			//Speed of particle

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

	//Collision
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
		tmd->vy = CL->VCM[1] + RV[i][1] - RS[1];
		tmd->vz = CL->VCM[2] + RV[i][2] - RS[2];
		//Increment link in list
		tmd = tmd->nextSRD;
		i++;
	}
}
void MPCcollision( cell *CL,spec *SP,specSwimmer SS,double KBT,int RTECH,double C,double S,double FRICCO,double TimeStep,int MDmode,int LC,double RELAX,double *CLQ,int outP ) {
/*
    Does the MPC collision, and returns the
    CM velocity and the local temperature of the cell.
*/
	particleMD *tmd;	//Temporary particleMD

	tmd = CL->MDpp;
	if (MDmode==MPCinMD) CL->MDpp=NULL;
	if( outP ) zeroPressureColl( CL );

	//Liquid Crystal
	if( LC!=ISOF ) {
		if( RTECH==DIPOLE_DIR_SUM || RTECH==DIPOLE_DIR_AV || RTECH==DIPOLE_DIR_SIG || RTECH==DIPOLE_DIR_SIG_SUM ) dipoleAndersenROT_LC( CL,SP,SS,KBT,RELAX,TimeStep,RTECH,CLQ,outP );
		else andersenROT_LC( CL,SP,SS,KBT,TimeStep,CLQ,outP );
	}
	//Newtonian Fluid
	else{
		//Passive MPC operators
		if( RTECH == MPCAT || RTECH == NOHI_MPCAT) andersenMPC( CL,SP,SS,KBT,CLQ,outP );
		else if( RTECH == RAT ) andersenROT( CL,SP,SS,KBT,CLQ,outP );
		else if( RTECH == LANG ) langevinMPC( CL,SP,SS,KBT,FRICCO,TimeStep,CLQ,outP );
		else if( RTECH == RLANG ) langevinROT( CL,SP,SS,KBT,FRICCO,TimeStep,CLQ,outP );
		else if( RTECH==ORTHAXIS || RTECH==ARBAXIS || RTECH==NOHI_ARBAXIS) stochrotMPC( CL,RTECH,C,S,CLQ,outP );
		else if( RTECH==MULTIPHASE ) andersenMULTIPHASE( CL,SP,SS,KBT,CLQ,outP );
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
	if (MDmode==MPCinMD) CL->MDpp=tmd;
}
void localVCM( double vcm[_3D],cell CL,spec *SP,specSwimmer specS ) {
/*
   This routine calculates the centre of mass
   velocity of a single cell or bin.
*/
	int id,i;
	double summ = 0.0;
	double mass;
	particleMPC *tmpc;
	particleMD *tmd;
	smono *tsm;

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
void localMPCVCM( double vcm[_3D],cell CL,spec *SP ) {
/*
   This routine calculates the centre of mass
   velocity of the MPCD particles in a single cell or bin.
*/
	int id,i;
	double summ = 0.0;
	double mass;
	particleMPC *tmpc;

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
double localMASS( cell CL,spec *SP,specSwimmer specS ) {
/*
   This routine calculates the local
   mass from the particleMPCs listed in
   the linked list (i.e. assuming
   localPROP hasn't been called).
*/
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
double localTEMP( cell CL,spec *SP,specSwimmer specS ) {
/*
   This routine calculates the local
   mass from the particleMPCs listed in
   the linked list (i.e. assuming
   localPROP hasn't been called).
*/
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
int localPOP( cell CL ) {
/*
   This routine calculates the local
   mass from the particleMPCs listed in
   the linked list (i.e. assuming
   localPROP hasn't been called).
*/
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
void scramble( particleMPC *p ) {
/*
   This routine randomly scrambles the velocities.
   Notice that it is slow because it switches some
   velocities more than once. But compared to copying
   all the velocities and managing/organizing the
   array every time, this is much much faster.
*/
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
void localCM( cell *CL,spec *SP,specSwimmer specS ) {
/*
	This is now a legacy method after being refactored into localPROP.
   This routine just calculates the centre of mass position for a given cell
*/
	int id,d;
	double mass,cellMass,Q[_3D];
	particleMPC *pMPC;	//Temporary pointer to MPC particles
	particleMD *pMD;	//Temporary pointer to MD particles
	smono *pSW;			//Temporary pointer to swimmer monomers

	//mass has already been computed

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
void localCM_SRD( cell CL,spec *SP,double r_cm[] ) {
/*
   This routine just calculates the centre of mass position for a given cell's SRD particles
*/
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
void localMomInertiaTensor( cell *CL,spec *SP,specSwimmer specS ) {
/*
   This routine calculates the moment of inertia tensor for a given cell
   --- Need to have already calculated the centre of mass position
*/
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
double localMomInertia_SRD( cell CL,spec *SP,double r0[],double n[] ) {
/*
   This routine calculates the moment of inertia value about a given position r0 and axis n (assumed normalized)
*/
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
void checkEscape_all( particleMPC *pp ) {
/*
    This is a very coarse check to make sure that none of the
    MPCD particles have escaped the control volume
*/
	int i,d;
	double Q[_3D];
	for( i=0; i<GPOP; i++ ){
		for( d=0; d<_3D; d++ ) Q[d] = (pp+i)->Q[d];
		for( d=0; d<_3D; d++ ) if( Q[d]<0. )  {
			#ifdef DBG
				printf( "Warning: Particle %d escaped control volume",i );
				pvec( Q,_3D );
			#endif
		}
	}
}
void cellVelForce( cell *CL,double addVel[3] ) {
/*
    Does the Andersen thermostat collision, and returns the
    CM velocity and the local temperature of the cell.
*/
	int i,j;
	particleMPC *tmpc;	//Temporary particleMPC
	particleMD *tmd;		//Temporary particleMD
	smono *pSW;					//Temporary pointer to swimmer monomers

	//Give particles a kick
	// MPC particles
	tmpc = CL->pp;
	i=0;
	while( tmpc!=NULL ) {
		for( j=0; j<DIM; j++ ) tmpc->V[j] += addVel[j];
		//Increment link in list
		tmpc = tmpc->next;
		i++;
	}
	//MD particles
	tmd = CL->MDpp;
	while( tmd!=NULL ) {
		tmd->vx += addVel[0];
		tmd->vy += addVel[1];
		tmd->vz += addVel[2];
		//Increment link in list
		tmd = tmd->nextSRD;
		i++;
	}
	// Swimmer particles
	pSW = CL->sp;
	while( pSW!=NULL ) {
		for( j=0; j<DIM; j++ ) pSW->V[j] += addVel[j];
		//Increment link in list
		pSW = pSW->next;
		i++;
	}
}
void cellVelSet( cell *CL,double vel[3] ) {
/*
    Does the Andersen thermostat collision, and returns the
    CM velocity and the local temperature of the cell.
*/
	int i,j;
	particleMPC *tmpc;	//Temporary particleMPC
	particleMD *tmd;		//Temporary particleMD
	smono *pSW;					//Temporary pointer to swimmer monomers

	//Give particles a kick
	// MPC particles
	tmpc = CL->pp;
	i=0;
	while( tmpc!=NULL ) {
		for( j=0; j<DIM; j++ ) tmpc->V[j] = vel[j];
		//Increment link in list
		tmpc = tmpc->next;
		i++;
	}
	//MD particles
	tmd = CL->MDpp;
	while( tmd!=NULL ) {
		tmd->vx += vel[0];
		tmd->vy += vel[1];
		tmd->vz += vel[2];
		//Increment link in list
		tmd = tmd->nextSRD;
		i++;
	}
	// Swimmer particles
	pSW = CL->sp;
	while( pSW!=NULL ) {
		for( j=0; j<DIM; j++ ) pSW->V[j] += vel[j];
		//Increment link in list
		pSW = pSW->next;
		i++;
	}
}

void timestep( cell ***CL,particleMPC *SRDparticles,spec SP[],bc WALL[],simptr simMD,specSwimmer *SS,swimmer swimmers[],double AVNOW[_3D],double AVV[_3D],double avDIR[_3D],inputList in,double *KBTNOW, double *AVS,int runtime,int MDmode,outputFlagsList outFlags,outputFilesList outFiles ) {

	int i,j,k;						//Counting variables
	double RSHIFT[_3D];		//Random vector to positively shift all components of the simulation
	double CLQ[_3D];					//Position of the cell since calculating pressure sucks
	int BC_FLAG;					//Flags if the BC moved in this time step
	int outPressure=0;		//Whether to make the pressure calculations (never used just outputted)
	int bcCNT,reCNT,rethermCNT;					//Count if any particles had problems with the BCs

	#ifdef DBG
		if ( DBUG >= DBGSTEPS ) {
			if( in.warmupSteps ) printf( "\nBegin warmup time step %i. Simulation time = %lf\n",runtime,runtime*in.dt );
			else printf( "\nBegin time step %i. Simulation time = %lf\n",runtime,runtime*in.dt );
		}
	#endif
	if( outFlags.PRESOUT>=OUT && runtime%outFlags.PRESOUT==0 ) outPressure=1;
	//Zero counters
	zerocnt( KBTNOW,AVNOW,AVS );

	// Zero impulse on BCs
	// NOTE: Louise thinks this should be fine being here (and did check),
	// This was moved when editing ghostPart to increase anchoring strength
	// (as oriBC is now called in ghostPart too).
	// If something looks bad related to mobile walls, maybe start looking here.
	for( i=0; i<NBC; i++ ) {
		zerovec(WALL[i].dV,DIM);
		zerovec(WALL[i].dL,_3D);
	}

	/* ****************************************** */
	/* ************* INTEGRATE MD *************** */
	/* ****************************************** */
	#ifdef DBG
		if( DBUG >= DBGTITLE && MDmode != noMD ) printf( "Integrate MD.\n" );
	#endif
	if( MDmode ) integrateMD(simMD,MDmode,in.stepsMD,SRDparticles,WALL,SP,GPOP,NBC,CL);
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
		gridShift_all( RSHIFT,0,SRDparticles,WALL,simMD,swimmers,MDmode );
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
	if( MDmode ) binMD( CL );
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
        ///TODO: OMP here
		for( i=0; i<XYZ_P1[0]; i++ ) for( j=0; j<XYZ_P1[1]; j++ ) for( k=0; k<XYZ_P1[2]; k++ ) {
			//LC collision algorithm (no collision if only 1 particle in cell)
			if( CL[i][j][k].POP > 1 ) LCcollision( &CL[i][j][k],SP,in.KBT,in.MFPOT,in.dt,*AVS,in.LC );
		}
		// Magnetic alignment is really part of the collision
		#ifdef DBG
			if( DBUG >= DBGTITLE && in.MAG_FLAG) printf( "Apply Magnetic Torque Operator.\n" );
		#endif
		if( in.MAG_FLAG ) magTorque_all( SRDparticles,SP,in.dt,in.MAG );
		#ifdef DBG
			if( DBUG >= DBGTITLE ) printf( "Orientation Shear Alignment.\n" );
        #endif
        ///TODO: OMP here
		for( i=0; i<XYZ_P1[0]; i++ ) for( j=0; j<XYZ_P1[1]; j++ ) for( k=0; k<XYZ_P1[2]; k++ ) {
			//Coupling shear to orientation
			if( CL[i][j][k].POP > 1 ) jefferysTorque( &CL[i][j][k],SP,in.dt );
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
#ifdef _OPENMP
#pragma omp parallel for private(i,j,k) default(shared)
#endif
	for( i=0; i<XYZ_P1[0]; i++ ) for( j=0; j<XYZ_P1[1]; j++ ) for( k=0; k<XYZ_P1[2]; k++ ) {
		//MPC/SRD collision algorithm (no collision if only 1 particle in cell)
		CLQ[0]=i+0.5;
		CLQ[1]=j+0.5;
		CLQ[2]=k+0.5;
		if( CL[i][j][k].POP > 1 ) MPCcollision( &CL[i][j][k],SP,*SS,in.KBT,in.RTECH,in.C,in.S,in.FRICCO,in.dt,MDmode,in.LC,in.TAU,CLQ,outPressure );
	}
	// Brownian thermostat (no hydrodynamic interactions -scramble velocities)
	if( in.RTECH == NOHI_ARBAXIS || in.RTECH == NOHI_MPCAT ) scramble( SRDparticles );
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
	if( MDmode != MPCinMD ) stream_all( SRDparticles,in.dt );
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
    ///TODO: OMP here
	for( i=0; i<GPOP; i++ ) MPC_BCcollision( SRDparticles,i,WALL,SP,in.KBT,in.dt,in.LC,&bcCNT,&reCNT,&rethermCNT,1 );
	// XYZPBC[0]=1;
	// XYZPBC[1]=1;
	// if(DIM>=_3D) XYZPBC[2]=1;
	// for( i=0; i<GPOP; i++ ) rudimentaryPBC_box( (SRDparticles+i) );
	// XYZPBC[0]=0
	// XYZPBC[1]=0;
	// if(DIM>=_3D) XYZPBC[2]=0;
	// for( i=0; i<GPOP; i++ ) rudimentaryBBBC_box( (SRDparticles+i) );
	// XYZPBC[0]=1;
	// XYZPBC[1]=0;
	// if(DIM>=_3D) XYZPBC[2]=1;
	// for( i=0; i<GPOP; i++ ) rudimentaryChannel_y( (SRDparticles+i) );
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
        ///TODO: OMP here
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
		//Build the flow profile by summing over everytime step and averaging after FLOWOUT iterations
	if( outFlags.FLOWOUT>=OUT && !in.warmupSteps) {
		localFLOW( CL,SP );
		sumFLOW( CL );
	}
}
void calcPressureStreaming( cell ***CL,spec *SP ) {
/*
   This function calculates the pre-collisional part (ballistic/streaming part) of the pressure tensor
	 It only calculates the contribution due to the MPCD fluid --- not MD or swimmers
	 The volume of the MPCD cell is 1
*/
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
void normPressureColl( cell *CL,double dt ) {
/*
    Zero collisional pressure term --- also divided by volume but cell volume=1
*/
	int i,j;
	for( i=0; i<DIM; i++ ) for( j=0; j<DIM; j++ ) CL->Pc[i][j] /= (-dt);
}
void calcPressureColl_preColl( double *relQ,double *dp,particleMPC *p,double *CLQ ) {
/*
    The pre-collision calculations needed to calculate the collisional pressure term
*/
	int d;
	for( d=0; d<DIM; d++ ) {
		//Calculate the relative position from the cell centre
		relQ[d] = p->Q[d] - CLQ[d];
		//Prepare the change in momentum
		dp[d] = p->V[d];
	}
}
void calcPressureColl_postColl( double *relQ,double *dp,double M,double *vel,cell *CL ) {
/*
    The post-collision calculations needed to calculate the collisional pressure term
*/
	int i,j;
	//Finish calculating the change in momentum
	for( i=0; i<DIM; i++ ) {
		dp[i] -= vel[i];
		dp[i] *= -M;
	}
	//Calculate the collisional contribution to the stress
	for( i=0; i<DIM; i++ ) for( j=0; j<DIM; j++ ) CL->Pc[i][j] += dp[i]*relQ[j];
}
