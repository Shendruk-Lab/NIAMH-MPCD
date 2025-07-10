///
/// @file
///
/// @brief This file contains methods to implement thermostats, compute averages and energies.
///
/// This file contains methods to implement thermostats, compute averages and energies.
///

# include<stdio.h>
# include<math.h>
# include<stdlib.h>

# include "../headers/SRDclss.h"
# include "../headers/definitions.h"
# include "../headers/globals.h"
# include "../headers/mpc.h"
# include "../headers/rand.h"
# include "../headers/mtools.h"
# include "../headers/therm.h"

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ***************** THERMO ***************** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */

///
/// @brief Function that calculates the velocity scaling factor g.
///
/// This routine calculates the velocity scaling factor g. It can do so
/// in 4 different ways, depending on which thermostat is being requested
/// as specified by the input parameter `TSTECH`.
/// i) No thermostat, at which case g = 1.
/// ii) Using the velocity scaling as the thermometer.
/// iii) Using the Berendsen thermostat.
/// iv) Finally, although not a proper thermostat, this routine can yield
///     a g that maximizes the velocity vector. Used in combination with
///     the rotation technique `RTECH = MPCAT`.
///
/// @param KBT The target temperature.
/// @param KBTNOW The current temperature.
/// @param t MPCD time-step.
/// @param RELAX Temperature relaxation time scale.
/// @param TSTECH Integer specifying the thermostat being employed.
///        0 for no thermostat, 1 for velocity scaling as thermostat,
///        2 for Berendsen thermostat and 4 for maximum velocity vector.
/// @see scaleT()
/// @return The velocity scaling factor g.
///
double thermostat( double KBT,double KBTNOW,double t,double RELAX,int TSTECH ) {
	double g;
	//No thermostat used
	if( TSTECH == NOTHERM ) g = 1.;
	//Velocity rescaling
	else if( TSTECH==VSC ) g = sqrt( KBT/KBTNOW );
	//Berendsen rescaling
	else if( TSTECH==BEREND ) g = sqrt( 1.+(KBT/KBTNOW-1.)*t/RELAX );
	else if( TSTECH==MAXV ) {
		g = sqrt( 1.+(KBT/KBTNOW-1.)*t/RELAX );	//Same as BEREND
	}
	else{
			printf( "Error:\tThermostat unacceptable.\n" );
			exit( 1 );
		}
	return g;
}

///
/// @brief Function that calls the thermostat and scales the velocities.
///
/// This routine calls the thermostat being requested and scales the velocities of
/// all MPCD particles and movable boundary conditions (obstacles). There are 5 possible
/// thermostats. These are:
/// i) No thermostat, at which case g = 1.
/// ii) Using the velocity scaling as the thermometer.
/// iii) Using the Berendsen thermostat, see H.J.C. Berendsen et. al.,
/// Molecular dynamics with coupling to an external bath. The Journal of Chemical
/// Physics, 81(8):3684–3690, 1984.
/// iv) Using the Heyes thermostat.
/// v) Finally, although not a proper thermostat, this routine can yield
///     a g that maximizes the velocity vector. Used in combination with
///     the rotation technique `RTECH = MPCAT`.
///
/// @param KBT The target temperature.
/// @param KBTNOW The current temperature.
/// @param t MPCD Time-step.
/// @param RELAX Temperature relaxation time scale.
/// @param VEL Average velocity.
/// @param VELNOW Current velocity.
/// @param TSTECH Integer specifying the thermostat being employed.
///        0 for no thermostat, 1 for velocity scaling as thermostat,
///        2 for Berendsen thermostat, 3 for Heyes thermostat and
///        4 for maximum velocity vector.
/// @param SP Array of all species.
/// @param LC Integer specifying LC mode being employed.
/// @param WALL Return pointer for array of all boundary conditions (obstacles).
/// @param p Return pointer to first element in array of all MPCD particles.
/// @param CL Array of all cells.
/// @see thermostat()
/// @see heyes_cell()
/// @note Heyes thermostat is a local one and implemented thorugh heyes_cell(), the rest are global
/// ones and are called through thermostat().
///
void scaleT( double KBT,double KBTNOW,double t,double RELAX,double VEL[],double VELNOW[],int TSTECH,spec SP[],int LC,bc *WALL,particleMPC *p,cell ***CL ) {
	int i,j,k;
	double TSC;	//The velocity scaling factor (commonly lambda)

	//Local thermostat
	if( TSTECH==HEYES ) {
		bin( CL,SP,WALL,KBT,LC,0 );
		//bin( CL );
		for( i=0;i<XYZ_P1[0];i++ ) for( j=0;j<XYZ_P1[1];j++ ) for( k=0;k<XYZ_P1[2];k++ ) heyes_cell( CL[i][j][k],KBT,KBTNOW,RELAX );
	}
// 	Global thermostats
	else if( TSTECH!=NOTHERM ) {
		if( TSTECH==MAXV ) {
			//Relaxes to a given velocity (actually it relaxes to a given kinetic energy rather than thermal energy (or shifted thermal energy))
			TSC = thermostat( KBT,KBTNOW,t,RELAX,TSTECH );
			for( j=0; j<DIM; j++ ) {
				if( fneq(VEL[j],0.0) ) {
// 					TSC = thermostat( VEL[j],VELNOW[j],t,RELAX,TSTECH );
					for( i=0; i<GPOP; i++ ) (p+i)->V[j] = VEL[j] + TSC * ((p+i)->V[j] -VEL[j]);
					for( i=0; i<NBC; i++ ) if( WALL[i].DSPLC ) {
						(WALL+i)->V[j] *= TSC;
						(WALL+i)->L[j] *= TSC;
					}
				}
			}
		}
		else {
			//Relaxes to zero
			TSC = thermostat( KBT,KBTNOW,t,RELAX,TSTECH );
			//Scale the velocities
			for( i=0; i<GPOP; i++ ) for( j=0; j<DIM; j++ ) (p+i)->V[j] *= TSC;
			for( i=0; i<NBC; i++ ) if( WALL[i].DSPLC ) for( j=0; j<DIM; j++ ) {
				(WALL+i)->V[j] *= TSC;
				(WALL+i)->L[j] *= TSC;
			}
		}
	}
}

///
/// @brief Function that rescales velocities when employing a Heyes thermostat.
///
/// This routines rescales the velocities of the particles within an MPCD cell
/// for a Heyes thermostat. See G. Gompper, T. Ihle, D.M. Kroll, and R.G. Winkler.
/// Multi-particle collision dynamics: A particle-based mesoscale simulation
/// approach to the hydrodynamics of complex fluids. Advances in Polymer Science,
/// pages 1–87. Springer Berlin Heidelberg, 2008.
///
/// @param CL An MPCD cell.
/// @param KBT Target temperature.
/// @param KBTNOW Current temperature.
/// @param RELAX Relaxation time scale.
/// @see scaleT()
///
void heyes_cell( cell CL,double KBT,double KBTNOW,double RELAX ) {
	int d;
	double sc_fctr;	//Scaling factor
	double prob;		//Probability of rescaling
	particleMPC *pMPC;	//Temporary pointer to MPC particles

	sc_fctr= 1.+RELAX*genrand_real( );
	if( genrand_real( ) < 0.5 ) sc_fctr = 1./sc_fctr;
	prob = 0;

	//Create probability of rescaling
	if( CL.pp != NULL ) {
		pMPC = CL.pp;
		while(pMPC != NULL) {
			for( d=0;d<DIM;d++ ) prob += (pMPC->V[d]-CL.VCM[d])*(pMPC->V[d]-CL.VCM[d]);
			//Increment link in list
			pMPC = pMPC->next;
		}
		prob *= -0.5*(sc_fctr*sc_fctr-1.)*CL.MASS;
		prob /= KBT;
		prob = exp( prob );
		prob *= smrtPow( sc_fctr, DIM*(CL.POP - 1) );
		if( prob>1. ) prob = 1.;
	}

	//Randomly check against this probability of acceptance. Scale if passes
	if( genrand_real() < prob ) if( CL.pp != NULL ) {
		pMPC = CL.pp;
		while(pMPC != NULL) {
			for( d=0;d<DIM;d++ ) pMPC->V[d] *= sqrt( KBT/KBTNOW );
			//Increment link in list
			pMPC = pMPC->next;
		}
	}
}

///
/// @brief Function that calculates the temperature out of the total energy.
///
/// This routine calculates the temperature by adding up the total energy,
/// from particles and mobile boundary conditions and employing the
/// equipartition theorem.
///
/// @param pp Pointer to the first element in the array of all MPCD particles.
/// @param SP Array of all species.
/// @param WALL Array of all boundary conditions (obstacles).
/// @param VEL Average velocity.
/// @return Temperature.
///
double TEMP( particleMPC *pp,spec SP[],bc WALL[],double VEL[] ) {
	double KBT = 0.;
	double temp = 0.;
	int i,j,k,c = 0;
	//Fluid particleMPC contribution
	for( i=0; i<GPOP; i++ ) {
		temp = ((pp+i)->V[0]-VEL[0])*((pp+i)->V[0]-VEL[0]);
		temp += ((pp+i)->V[1]-VEL[1])*((pp+i)->V[1]-VEL[1]);
		temp += ((pp+i)->V[2]-VEL[2])*((pp+i)->V[2]-VEL[2]);
		KBT +=  SP[(pp+i)->SPID].MASS * temp;
// 		KBT += SP[(pp+i)->SPID].M * ( (pp+i)->V[0]*(pp+i)->V[0]+(pp+i)->V[1]*(pp+i)->V[1]+(pp+i)->V[2]*(pp+i)->V[2] );
	}
	for( i=0; i<NBC; i++ ) if( WALL[i].DSPLC ) {
		//BC object's kinetic contribution
		temp = (WALL[i].V[0]-VEL[0])*(WALL[i].V[0]-VEL[0]);
		temp += (WALL[i].V[1]-VEL[1])*(WALL[i].V[1]-VEL[1]);
		temp += (WALL[i].V[2]-VEL[2])*(WALL[i].V[2]-VEL[2]);
		KBT += WALL[i].MASS * temp;
// 		KBT += WALL[i].M * (WALL[i].V[0]*WALL[i].V[0] + WALL[i].V[1]*WALL[i].V[1] + WALL[i].V[2]*WALL[i].V[2]);
		//BC object's rotational contribution
		for( j=0; j<_3D; j++ ) for( k=0; k<_3D; k++ ) KBT += WALL[i].L[j] * WALL[i].I[j][k] * WALL[i].L[k];
		c++;
	}
	//Equipartition function ( sum(mv^2/2) = DIM*N*kBT/2 )
	KBT /= ( (double)((GPOP+c)*DIM) );
	return KBT;
}

///
/// @brief Function that calculates the total energy of the system.
///
/// This function calculates the total energy of the system. It considers
/// the kinetic energy of all MPCD particles as well as the translational and rotational
/// energy of all mobile boundary conditions (obstacles).
///
/// @param p Pointer to first element in the array of all MPCD particles.
/// @param pSP Pointer to first element in the array of all species.
/// @param WALL Array of all boundary conditions (obstacles).
/// @return The system's total energy.
///
double calcE( particleMPC *p,spec *pSP,bc WALL[] ) {
	double E,TE = 0.;
	int i,j,k;
	for( i=0; i<GPOP; i++ ) {
		E = 0.;
		for( j=0; j<DIM; j++ ) E += (p+i)->V[j] * (p+i)->V[j];
		E *= 0.5 * (pSP+(p+i)->SPID)->MASS;
		TE += E;
	}
	for( i=0; i<NBC; i++ ) if( WALL[i].DSPLC ) {
		//Kinetic
		E = 0.;
		for( j=0; j<DIM; j++ ) E += WALL[i].V[j]*WALL[i].V[j];
		E *= 0.5 * WALL[i].MASS;
		TE += E;
		//Rotational
		E = 0.;
		for( j=0; j<_3D; j++ ) for( k=0; k<_3D; k++ ) E += WALL[i].L[j] * WALL[i].I[j][k] * WALL[i].L[j];
		E *= 0.5;
		TE += E;
	}
	return TE;
}

///
/// @brief Function that calculates a particle's kinetic energy.
///
/// This function returns the kinetic energy of a single MPCD particle.
///
/// @param p Pointer to the MPCD particle whose energy is being calculated.
/// @param pSP Pointer to the species of the MPCD particle whose energy is being calculated.
/// @return The particle's kinetic energy.
///
double calcE_MPC( particleMPC *p,spec *pSP ) {
	double E = 0.;
	int j;

	for( j=0; j<DIM; j++ ) E += p->V[j] * p->V[j];
	E *= 0.5 * (pSP+p->SPID)->MASS;
	return E;
}

///
/// @brief Function that computes the kinetic energy of a boundary condition.
///
/// This function returns the kinetic energy (translational + rotational).
/// of a single mobile boundary condition (obstacle).
///
/// @param WALL Pointer to the boundary condition whose energy is being calulated.
/// @return Energy of the boundary condition (obstacle).
///
double calcE_BC( bc *WALL ) {
	double E=0., TE=0.;
	int j,k;

	if( WALL->DSPLC ) {
		//Translational
		for( j=0; j<DIM; j++ ) E += WALL->V[j] * WALL->V[j];
		E *= 0.5 * WALL->MASS;
		TE += E;
		//Rotational
		E = 0.;
		for( j=0; j<_3D; j++ ) for( k=0; k<_3D; k++ ) E += WALL->L[j] * WALL->I[j][k] * WALL->L[k];
		E *= 0.5;
		TE += E;
	}

	return TE;
}

///
/// @brief Function that computes the total potential orientational energy of all nematic particles.
///
/// This function returns the sum of the potential orientational energy of all nematic particles.
/// Since this is a potential energy, this is a negative number.
///
/// @param CL Array of all cells.
/// @param LC Integer specifying type of liquid crystal/
/// @param pSP Pointer to the species of the MPCD particle whose energy is being calculated.
/// @return Total potential orientational energy.
///
double calcE_LC( cell ***CL,int LC,spec *pSP ) {
	double wmf=0.;
	double S,un,DIR[_3D],u[_3D];
	particleMPC *tmpc;	//Temporary particleMPC
	int a,b,c,i,id;
	//double invdim=1./((double)DIM);

	if( LC ) for( a=0; a<XYZ[0]; a++ ) for( b=0; b<XYZ[1]; b++ ) for( c=0; c<XYZ[2]; c++ ) if( CL[a][b][c].POP > 1 ) {
		S = CL[a][b][c].S;
		for( i=0; i<DIM; i++ ) DIR[i] = CL[a][b][c].DIR[i];
		tmpc = CL[a][b][c].pp;
		while( tmpc != NULL ) {
			id = tmpc->SPID;
			for( i=0; i<DIM; i++ ) u[i] = tmpc->U[i];
			un = dotprod( u,DIR,DIM );
			wmf += ( S*un*un )*( (pSP+id)->sMFPOT );
			//wmf += (1.-S)*invdim;		//Don't include constant (wrt u.n) term
			tmpc = tmpc->next;
		}
	}
	return wmf;
}

///
/// @brief Function that calculates the average global velocity.
///
/// This function computes the average velocity of the system by averaging
/// over the center of mass velocities of all cells.
///
/// @param CL Array of all cells.
/// @param AVVEL Return pointer for the average global velocity.
///
void avVel( cell ***CL,double AVVEL[] ) {
	int a,b,c,d;
	for( d=0; d<DIM; d++ ) AVVEL[d]=0.;
	for( a=0; a<XYZ_P1[0]; a++ ) for( b=0; b<XYZ_P1[1]; b++ ) for( c=0; c<XYZ_P1[2]; c++ ) for( d=0; d<DIM; d++ ) AVVEL[d] += CL[a][b][c].VCM[d];
	for( d=0; d<DIM; d++ ) AVVEL[d] /= VOL;
}

///
/// @brief Function that calculates the average global orientation.
///
/// This function computes the average orientation of the system by averaging all the orientation of particles
///
/// @param p Return pointer to first element in array of all MPCD particles.
/// @param AVVEL Return pointer for the average global orientation.
///
void avOri( particleMPC *p,double AVORI[] ) {
	int j,i,d;
	for( d=0; d<DIM; d++ ) AVORI[d]=0.;
	for( j=0; j<DIM; j++ ) {
		for( i=0; i<GPOP; i++ ) AVORI[j]+=(p+i)->U[j];
		AVORI[j] /= GPOP;
	}
}


///
/// @brief Function that calculates the global average enstrophy.
///
/// This function returns the global average enstrophy `E` (mean squared vorticity)
/// by averaging it over all cells.
///
/// @param CL Array of all cells.
/// @return Global average enstrophy.
///
double avEnstrophy( cell ***CL ) {
	int a,b,c;
	double w[_3D];
	double E;

	E=0.;
	for( a=0; a<XYZ_P1[0]; a++ ) for( b=0; b<XYZ_P1[1]; b++ ) for( c=0; c<XYZ_P1[2]; c++ ) {
		w[0]=(CL[a][b][c].E[2][1] - CL[a][b][c].E[1][2]);
		w[1]=(CL[a][b][c].E[0][2] - CL[a][b][c].E[2][0]);
		w[2]=(CL[a][b][c].E[1][0] - CL[a][b][c].E[0][1]);
		E += dotprod( w,w, _3D );
	}
	E /= VOL;
	E *= 0.5;
	return E;
}
