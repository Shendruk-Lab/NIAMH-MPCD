
//================================================================================
//
// name:   md.c
// author: ftessier
// date:   2005-05-03 @ 11:04:36
//
//================================================================================
///
/// @file md.c
/// Core Molecular Dynamics algorithm. The functions provided in this file
/// constitute the essence of the md simulation: md stepping, force calculation
/// and neighbor lists updates. Inspired in part from "The Art of Molecular Dynamics"
/// by D. C. Rapaport.
///
/// This code has been optimized somewhat, but not at the expense of clarity.
/// As Knuth says, it is much easier to optimize correct code than to correct
/// optimized code!
///
/// All the information about the simulation is encapsulated in a local sim
/// structure declared in the main function. It is then easy to pass this
/// information via a pointer to sim. Function often use local copies of the
/// variables in sim for clarity and portability.
///
/// Some basic OpenMP directives are included for some loops, but are only
/// effective when the code is compiled with OpenMP (e.g., -xopenmp on hpcvl).
///
//================================================================================


#include <stdlib.h>
#include <math.h>
#include "mdtypes.h"
#include "mderror.h"
#include "mdfiles.h"
#include "mdsetup.h"
#include "mdmeasure.h"
#include "mdmemory.h"
#include "mdutil.h"
#include "mdthermostat.h"
#include "mdsrd.h"

#include "../mpcd/headers/bc.h"
#include "../mpcd/headers/lc.h"
#include "../mpcd/headers/mdbc.h"
#include "../mpcd/headers/definitions.h"
#include "../mpcd/headers/globals.h"

#ifdef _OPENMP
#include <omp.h>
#endif


/// The standard program entry point. The sim variables holds all the
/// information about the simulation. It is allocated and initialized within
/// SetupSimulation function///
/// @param		sim md simulation pointer

//================================================================================
simptr launchMD (int argc, char *argv[])
//================================================================================
{
	// setup simulation
	return SetupSimulation (argc, argv);
}


/// The standard program entry point. The sim variables holds all the
/// information about the simulation. The function simply iterates the MD step until done.
///
/// @param		sim md simulation pointer
/// @param		steps number of md steps to integrate

//================================================================================
void integrateMD (simptr sim, int MDmode, int steps,struct particleMPC *pSRD,struct bc *pBC,struct spec *SP,int GPOP,int NBC,struct cell ***CL)
//================================================================================
{
	int	i;

	sim->drTotMax = 0.0;	//reset maximum displacement to zero since particles were rebinned

	// Add the following line if you need electostatics
	if (sim->monoCharge[0]!=0 || sim->rCutCoul!=0.) ComputeElectrostaticForcesSRD(sim,pSRD,SP,GPOP,steps,CL);
	// main simulation loop
	for(i=0;i<steps;i++) {
		MDstep (sim,MDmode,pSRD,pBC,SP,GPOP,NBC,CL);
	}
	PeriodicBoundaries (sim);
}


/// The standard program entry point dealing with MD warmup during MPCD warmup. The sim variables holds all 
/// the information about the simulation. The function simply iterates the MD steps and pinns the central 
/// monomer in its initial place until done.
///
/// @param		sim md simulation pointer
/// @param		steps number of md steps to integrate

//================================================================================
void integrateMD_Pinned (simptr sim, int MDmode, int steps,struct particleMPC *pSRD,struct bc *pBC,struct spec *SP,int GPOP,int NBC,struct cell ***CL)
//================================================================================
{
	int	i,nAtom;
	particleMD	*p ;
	double RX=0.0,RY=0.0,RZ=0.0,WX=0.0,WY=0.0,WZ=0.0;

	sim->drTotMax = 0.0;	//reset maximum displacement to zero since particles were rebinned

	// Add the following line if you need electostatics
	if (sim->monoCharge[0]!=0 || sim->rCutCoul!=0.) ComputeElectrostaticForcesSRD(sim,pSRD,SP,GPOP,steps,CL);

	nAtom	= sim->atom.n;
	int nRod = (int)(nAtom/2) + TRANS_SHIFT;  // offset +TRANS_SHIFT

	// the central monomer 
	p	= (sim->atom.items) + nRod ;
	RX = p->rx;
	RY = p->ry;
	RZ = p->rz;
	WX = p->wx;
	WY = p->wy;
	WZ = p->wz;

	
	// main simulation loop
	for(i=0;i<steps;i++) {
		MDstep (sim,MDmode,pSRD,pBC,SP,GPOP,NBC,CL);
		// the central monomer 
		p	= (sim->atom.items) + nRod ;
		// Put it back to its initial position
		p->rx = RX;
		p->ry = RY;
		p->rz = RZ;
		// Put it back to its initial position
		p->wx = WX;
		p->wy = WY;
		p->wz = WZ;
		// Kept pinned
		p->vx = 0.0;
		p->vy = 0.0;
		p->vz = 0.0;
		// Kept pinned
		p->ax = 0.0;
		p->ay = 0.0;
		p->az = 0.0;
	}
	PeriodicBoundaries (sim);
}


//================================================================================
void ComputeElectrostaticForcesSRD (simptr sim,struct particleMPC *pSRD,struct spec *SP,int GPOP, int steps, struct cell ***CL)
//================================================================================
{
	int			i, j, n;
	particleMD	*p1;
	particleMPC	*p2;
	item1STD	*charge;
	real		*box;
	real		dx, dy, dz;
	real		lambda_D=0, rCutCoul=0;
	real		potE=0, coulE=0, Efield=0, dt;
	real		r, *weight;
	int		xbox,ybox,zbox;
	int 		xc,yc,zc;
	int		a,b,c,pbcx1,pbcy1,pbcz1,pbcx2,pbcy2,pbcz2;
	real		drTotMax;
	real		condenseCriteria,q0;

	condenseCriteria=sim->condenseCriteria;
	q0 = sim->monoCharge[0];

	// calculate coulomb interaction variables (here for clarity)
	if (sim->charge.n > 0) {
		rCutCoul = sim->rCutCoul;
	} else return;

	// local sim variables
	charge  	= sim->charge.items;
	box		= sim->box;
	lambda_D	= sim->lambda_D;
	Efield  	= sim->Efield[sim->phase];
	dt		= sim->dt*steps;
	drTotMax	= sim->drTotMax;

	// compute coulomb forces for ALL charge pairs
	n = sim->charge.n;
	weight = mycalloc(n, sizeof(real));

	#ifdef _OPENMP
	#pragma omp parallel  for \
		schedule  (static) \
		default   (none) \
		shared	(n, charge, rCutCoul2, bjerrumkT, sim) \
		private   (i, j, p1, p2, dx, dy, dz, E, Efield, Eforce) \
		reduction (+: coulE, potE)
	#endif

	for (j=0; j<GPOP; j++) {
		(pSRD+j)->q=0;
	}

	for (i=0; i<n; i++) {
		// extract first particle pointer
		p1 = charge[i].p1;
		p1->q = q0/fabs(q0);
		a = p1->rx;
		b = p1->ry;
		c = p1->rz;
		weight[i] = 0;

		//note since xc, yc, and zc are ints the floor is implicit for positive values
		//this does the MPC-MD bead interactions
		//we have to go one extra box on the plus side since a,b,c is can be one lower than the MD bead's SRD box
		pbcx1 = 0.0;//these pbc variables account for the fact that box 10 is box 0 meaning we need to add 1 to the halo
		if (a-rCutCoul-drTotMax<1) pbcx1 = 1.0;
		pbcx2 = 0.0;
		if (a+rCutCoul+drTotMax>box[x_]-1) pbcx2 = 1.0;
		for (xc=floor(a-rCutCoul-drTotMax-pbcx1);xc<a+rCutCoul+drTotMax+pbcx2+2;xc++) {
			xbox=xc;
			if (xbox>box[x_]) {
				if (sim->pbcond & PBC_COND_x) while (xbox>box[x_]) xbox = xbox-box[x_]-1;
				else continue;
			}
			if (xbox<0) {
				if (sim->pbcond & PBC_COND_x) while (xbox<0) xbox = xbox+box[x_]+1;
				else continue;
			}
			pbcy1 = 0.0;
			if (b-rCutCoul-drTotMax<1) pbcy1 = 1.0;
			pbcy2 = 0.0;
			if (b+rCutCoul+drTotMax>box[y_]-1) pbcy2 = 1.0;
			for (yc=floor(b-rCutCoul-drTotMax-pbcy1);yc<b+rCutCoul+drTotMax+pbcy2+2;yc++) {
				ybox=yc;
				if (ybox>box[y_]) {
					if (sim->pbcond & PBC_COND_y) while (ybox>box[y_]) ybox = ybox-box[y_]-1;
					else continue;
				}
				if (ybox<0) {
					if (sim->pbcond & PBC_COND_y) while (ybox<0) ybox = ybox+box[y_]+1;
					else continue;
				}
				pbcz1 = 0.0;
				if (c-rCutCoul-drTotMax<1) pbcz1 = 1.0;
				pbcz2 = 0.0;
				if (c+rCutCoul+drTotMax>box[z_]-1) pbcz2 = 1.0;
				for (zc=floor(c-rCutCoul-drTotMax-pbcz1);zc<c+rCutCoul+drTotMax+pbcz2+2;zc++) {
					zbox=zc;
					if (zbox>box[z_]) {
						if (sim->pbcond & PBC_COND_z) while (zbox>box[z_]) zbox = zbox-box[z_]-1;
						else continue;
					}
					if (zbox<0) {
						if (sim->pbcond & PBC_COND_z) while (zbox<0) zbox = zbox+box[z_]+1;
						else continue;
					}
					p2 = CL[xbox][ybox][zbox].pp;
					while (p2!=NULL) {
						// calculate dr
						dx = ( p2->Q[0] - p1->rx );
						dy = ( p2->Q[1] - p1->ry );
						dz = ( p2->Q[2] - p1->rz );
						// apply periodic boundary conditions
						ApplyPBC_dr (sim,&dx,&dy,&dz);
						// calculate DH weighting
						r = sqrt(dx*dx+dy*dy+dz*dz);
						if (r<rCutCoul) weight[i] += exp(-r/lambda_D)/r;
						p2=p2->next;
					}
				}
			}
		}

		//note since xc, yc, and zc are ints the floor is implicit for positive values
		//this does the MPC-MD bead interactions
		//we have to go one extra box on the plus side since a,b,c is can be one lower than the MD bead's SRD box
		pbcx1 = 0.0;//these pbc variables account for the fact that box 10 is box 0 meaning we need to add 1 to the halo
		if (a-rCutCoul-drTotMax<1) pbcx1 = 1.0;
		pbcx2 = 0.0;
		if (a+rCutCoul+drTotMax>box[x_]-1) pbcx2 = 1.0;
		for (xc=floor(a-rCutCoul-drTotMax-pbcx1);xc<a+rCutCoul+drTotMax+pbcx2+2;xc++) {
			xbox=xc;
			if (xbox>box[x_]) {
				if (sim->pbcond & PBC_COND_x) while (xbox>box[x_]) xbox = xbox-box[x_]-1;
				else continue;
			}
			if (xbox<0) {
				if (sim->pbcond & PBC_COND_x) while (xbox<0) xbox = xbox+box[x_]+1;
				else continue;
			}
			pbcy1 = 0.0;
			if (b-rCutCoul-drTotMax<1) pbcy1 = 1.0;
			pbcy2 = 0.0;
			if (b+rCutCoul+drTotMax>box[y_]-1) pbcy2 = 1.0;
			for (yc=floor(b-rCutCoul-drTotMax-pbcy1);yc<b+rCutCoul+drTotMax+pbcy2+2;yc++) {
				ybox=yc;
				if (ybox>box[y_]) {
					if (sim->pbcond & PBC_COND_y) while (ybox>box[y_]) ybox = ybox-box[y_]-1;
					else continue;
				}
				if (ybox<0) {
					if (sim->pbcond & PBC_COND_y) while (ybox<0) ybox = ybox+box[y_]+1;
					else continue;
				}
				pbcz1 = 0.0;
				if (c-rCutCoul-drTotMax<1) pbcz1 = 1.0;
				pbcz2 = 0.0;
				if (c+rCutCoul+drTotMax>box[z_]-1) pbcz2 = 1.0;
				for (zc=floor(c-rCutCoul-drTotMax-pbcz1);zc<c+rCutCoul+drTotMax+pbcz2+2;zc++) {
					zbox=zc;
					if (zbox>box[z_]) {
						if (sim->pbcond & PBC_COND_z) while (zbox>box[z_]) zbox = zbox-box[z_]-1;
						else continue;
					}
					if (zbox<0) {
						if (sim->pbcond & PBC_COND_z) while (zbox<0) zbox = zbox+box[z_]+1;
						else continue;
					}
					p2 = CL[xbox][ybox][zbox].pp;
					while (p2!=NULL) {
						// calculate dr
						dx = ( p2->Q[0] - p1->rx );
						dy = ( p2->Q[1] - p1->ry );
						dz = ( p2->Q[2] - p1->rz );
						// apply periodic boundary conditions
						ApplyPBC_dr (sim,&dx,&dy,&dz);
						// calculate DH weighting
						r = sqrt(dx*dx+dy*dy+dz*dz);
						if (r<rCutCoul) p2->q -= p1->q*exp(-r/lambda_D)/weight[i]/r;
						p2=p2->next;
					}
				}
			}
		}
	}

	if (q0<0) {
		for (j=0; j<GPOP; j++) {
			// flip charge to show condensation
			if ((pSRD+j)->q>condenseCriteria) (pSRD+j)->q = -(pSRD+j)->q;
		}
		for (i=0; i<n; i++) { //here we "condense charges...
			// extract first particle pointer
			p1 = charge[i].p1;
			a = p1->rx;
			b = p1->ry;
			c = p1->rz;

			//note since xc, yc, and zc are ints the floor is implicit for positive values
			//this does the MPC-MD bead interactions
			//we have to go one extra box on the plus side since a,b,c is can be one lower than the MD bead's SRD box
			pbcx1 = 0.0;//these pbc variables account for the fact that box 10 is box 0 meaning we need to add 1 to the halo
			if (a-rCutCoul-drTotMax<1) pbcx1 = 1.0;
			pbcx2 = 0.0;
			if (a+rCutCoul+drTotMax>box[x_]-1) pbcx2 = 1.0;
			for (xc=floor(a-rCutCoul-drTotMax-pbcx1);xc<a+rCutCoul+drTotMax+pbcx2+2;xc++) {
				xbox=xc;
				if (xbox>box[x_]) {
					if (sim->pbcond & PBC_COND_x) while (xbox>box[x_]) xbox = xbox-box[x_]-1;
					else continue;
				}
				if (xbox<0) {
					if (sim->pbcond & PBC_COND_x) while (xbox<0) xbox = xbox+box[x_]+1;
					else continue;
				}
				pbcy1 = 0.0;
				if (b-rCutCoul-drTotMax<1) pbcy1 = 1.0;
				pbcy2 = 0.0;
				if (b+rCutCoul+drTotMax>box[y_]-1) pbcy2 = 1.0;
				for (yc=floor(b-rCutCoul-drTotMax-pbcy1);yc<b+rCutCoul+drTotMax+pbcy2+2;yc++) {
					ybox=yc;
					if (ybox>box[y_]) {
						if (sim->pbcond & PBC_COND_y) while (ybox>box[y_]) ybox = ybox-box[y_]-1;
						else continue;
					}
					if (ybox<0) {
						if (sim->pbcond & PBC_COND_y) while (ybox<0) ybox = ybox+box[y_]+1;
						else continue;
					}
					pbcz1 = 0.0;
					if (c-rCutCoul-drTotMax<1) pbcz1 = 1.0;
					pbcz2 = 0.0;
					if (c+rCutCoul+drTotMax>box[z_]-1) pbcz2 = 1.0;
					for (zc=floor(c-rCutCoul-drTotMax-pbcz1);zc<c+rCutCoul+drTotMax+pbcz2+2;zc++) {
						zbox=zc;
						if (zbox>box[z_]) {
							if (sim->pbcond & PBC_COND_z) while (zbox>box[z_]) zbox = zbox-box[z_]-1;
							else continue;
						}
						if (zbox<0) {
							if (sim->pbcond & PBC_COND_z) while (zbox<0) zbox = zbox+box[z_]+1;
							else continue;
						}
						p2 = CL[xbox][ybox][zbox].pp;
						while (p2!=NULL) {
							// calculate dr
							dx = ( p2->Q[0] - p1->rx );
							dy = ( p2->Q[1] - p1->ry );
							dz = ( p2->Q[2] - p1->rz );
							// apply periodic boundary conditions
							ApplyPBC_dr (sim,&dx,&dy,&dz);
							// calculate DH weighting
							r = sqrt(dx*dx+dy*dy+dz*dz);
							if (r<rCutCoul && p2->q/q0>0) p1->q -= q0*exp(-r/lambda_D)/weight[i]/r;
							p2=p2->next;
						}
					}
				}
			}
		}
	}

	for (j=0; j<GPOP; j++) {
		// adjust counter ion speed
		if ( (pSRD+j)->q/q0<0) {
			(pSRD+j)->V[0] += Efield*(pSRD+j)->q*dt; //N.B. dt is not sim->dt but sim->dt*mdSteps
		}
	}
	// add potential energies in sim structure
	sim->potE  += potE;
	sim->coulE  += coulE;
	free(weight);
}


/// Performs a single MD step. This is a wrapper function which simply calls other
/// functions to:
/// <ol>
///		<li> perform one molecular dynamics step in phase space
///		<li> increments step counters
///		<li> update the neighbor list if needed
///		<li> perform physical measurements
///		<li> update the simulation phase.
///	</ol>
///
/// @param		sim	a pointer to a simulation structure
/// @return		the status of the simualtion, DONE or NOT_DONE

//================================================================================
// inline int MDstep (simptr sim,int MDmode,struct particleMPC *pSRD,struct bc *pBC,struct spec *SP,int GPOP,int NBC,struct cell ***CL)
int MDstep (simptr sim,int MDmode,struct particleMPC *pSRD,struct bc *pBC,struct spec *SP,int GPOP,int NBC,struct cell ***CL)
//================================================================================
{
	int done=0;

	// integrate motion for one time step
// 	PeriodicBoundaries (sim);
	VelocityVerletStep (sim,MDmode,pSRD,pBC,SP,GPOP,NBC,CL);

	// increment step counters and simulation time
	IncrementStepCounters (sim);

	// update neighbor list if needed
//	CheckNebrList (sim);
// 	PeriodicBoundaries (sim);

	// perform measurements
	Measure (sim);
// 	Measure ( sim,*(pBC+0) );	//Hard coded that the BC of interest is the first in the list (bc.inp)

	// check simulation phase progress
	done = SimulationPhaseCheck (sim);

	// return to main simulation loop after each MD step
	return done;
}


/// Integrate one time step of Newton's equations of motion. We are using the
/// Velocity Verlet algorithm (Allen and Tilsdeley, page 81). According to Allen,
/// this is the only algorithm that makes sense for very large systems on
/// modern hardware. In his words, "predictor-corrector" is dead!
///
/// @param		sim	a pointer to a simulation structure
/// @return 	void

//================================================================================
void VelocityVerletStep (simptr sim,int MDmode,struct particleMPC *pSRD,struct bc *pBC,struct spec *SP,int GPOP,int NBC,struct cell ***CL)
//================================================================================
{
	int	 	i,j, nAtom, groupThermRescale, bcCNT, reCNT, rethermCNT;
	particleMD	*atom, *p;
	real		dt, dtHalf, kinE, kinETherm, v2, v2max;
	real		kT;

	// local sim variables
	atom	= sim->atom.items;
	nAtom	= sim->atom.n;
	dt	= sim->dt;
	dtHalf	= 0.5 * dt;
	groupThermRescale = sim->groupThermRescale[sim->phase];
	kT = sim->kT[sim->phase];

	// update positions and mid-step velocities
	#ifdef _OPENMP
	#pragma omp parallel for		 	 \
			schedule (static)			 \
			private  (n, p, dx, dy, dz)  \
			shared   (atom, nAtom)
	#endif
	kinE = kinETherm = 0;
	v2 = v2max = 0;
	for (i=0; i<nAtom; i++) {
		p = atom+i;
		// increment velocities
		p->vx += p->ax*dtHalf;
		p->vy += p->ay*dtHalf;
		p->vz += p->az*dtHalf;
		// increment positions
		p->rx += p->vx*dt;
		p->ry += p->vy*dt;
		p->rz += p->vz*dt;
		// increment world (true) positions
		p->wx += p->vx*dt;
		p->wy += p->vy*dt;
		p->wz += p->vz*dt;
		v2 = (p->vx*p->vx) + (p->vy*p->vy) + (p->vz*p->vz);
		kinE += p->mass*v2;
		if (v2 > v2max)
			v2max = v2;
	}
	for (i=0; i<nAtom; i++) {
		p = atom+i;
		MD_BCcollision(p,pBC,kT,dt);
	}

	if(MDmode==MPCinMD) for (j=0; j<GPOP; j++) {
		//First half of Verlet
		(pSRD+j)->Q[0] += dtHalf*(pSRD+j)->V[0];
		(pSRD+j)->Q[1] += dtHalf*(pSRD+j)->V[1];
		(pSRD+j)->Q[2] += dtHalf*(pSRD+j)->V[2];
		v2 = ((pSRD+j)->V[0]*(pSRD+j)->V[0]) + ((pSRD+j)->V[1]*(pSRD+j)->V[1]) + ((pSRD+j)->V[2]*(pSRD+j)->V[2]);
		if (v2 > v2max)
			v2max = v2;
	}

	bcCNT=0;
	reCNT=0;
	rethermCNT=0;
	if(MDmode==MPCinMD) for( i=0; i<GPOP; i++ ) {
		MPC_BCcollision( pSRD,i,pBC,SP,kT,dtHalf,0,&bcCNT,&reCNT,&rethermCNT,1 );
	}
	sim->drTotMax  += sqrt(v2max)*dt;
	// calculate the net force on all particles
	ComputeForcesSRD (sim,MDmode,pSRD,SP,GPOP,CL);

	if(MDmode==MPCinMD) for (j=0; j<GPOP; j++) {
		//Second half of Verlet
		(pSRD+j)->Q[0] += dtHalf*(pSRD+j)->V[0];
		(pSRD+j)->Q[1] += dtHalf*(pSRD+j)->V[1];
		(pSRD+j)->Q[2] += dtHalf*(pSRD+j)->V[2];
		(pSRD+j)->S_flag = NO_STREAM;
	}
	if(MDmode==MPCinMD) for( i=0; i<GPOP; i++ ) {
		MPC_BCcollision( pSRD,i,pBC,SP,sim->kT[sim->phase],dtHalf,0,&bcCNT,&reCNT,&rethermCNT,1 );
	}

	// divide by mass to get acceleration
	for (i=0; i<nAtom; i++) {
		p = atom+i;
		if (p->mass != 1.0) {
			p->ax /= p->mass;
			p->ay /= p->mass;
			p->az /= p->mass;
		}
	}

	// update velocities, maximum velocity and kinetic energies
	for (i=0; i<nAtom; i++) {
		p = atom+i;
		p->vx += p->ax*dtHalf;
		p->vy += p->ay*dtHalf;
		p->vz += p->az*dtHalf;
		if (p->group & groupThermRescale)
			kinETherm += p->mass*v2;
	}

	// update sim variables
	sim->kinE	    = 0.5*kinE;
	sim->kinETherm  = 0.5*kinETherm;
	if (sim->drTotMax>0.5) {
		printf ("%E drTotMax\n", sim->drTotMax);
		fflush(stdout);
		// exit(1);
	}
}


/// Calculation of all forces acting on the atoms. First, the energy, accelearations
/// and force variables are all reset to 0, and then separate functions are called to
/// handle each type of interaction.
///
/// @param		sim	a pointer to a simulation structure
/// @return 	void

//================================================================================
void ComputeForcesSRD (simptr sim,int MDmode,struct particleMPC *pSRD,struct spec *SP,int GPOP,struct cell ***CL)
//================================================================================
{
	int		i, nAtom;
	particleMD	*atom;

	// local sim variables
	atom  = sim->atom.items;
	nAtom = sim->atom.n;

	// reset all potential energies
	sim->potE  = 0;
	sim->ljE   = 0;
	sim->harmE = 0;
	sim->coulE = 0;
	sim->feneE = 0;
	sim->bendE = 0;
	sim->nemE = 0;

	// reset accelerations, forces, etc.
	for (i=0; i<nAtom; i++) {
		atom[i].ax	= 0;
		atom[i].ay	= 0;
		atom[i].az	= 0;
		#ifdef TEMPERATURE_CONF
		// configurational temperature
		atom[i].Tfx   = 0;
		atom[i].Tfy   = 0;
		atom[i].Tfz   = 0;
		atom[i].Tdivf = 0;
		#endif
	}

	// compute forces
	// ComputeDispersionForcesSRDCell (sim,MDmode,pSRD,SP,GPOP,CL);
	if (MDmode == MPCinMD) ComputeDispersionForcesSRDCell (sim,MDmode,pSRD,SP,GPOP,CL);
	else ComputeDispersionForces(sim);
	ComputeElectrostaticForces(sim);
	ComputeAnchorForces(sim);
	ComputeFeneForces(sim);
	ComputeBendForces(sim);
	ComputeDihedralForces(sim);
	ComputeNemForces(sim,SP,CL);
	ComputeSqueezeForces(sim);
}

/// Calculation of all forces acting on the atoms. First, the energy, accelearations
/// and force variables are all reset to 0, and then separate functions are called to
/// handle each type of interaction. Note that this is identical to above but no electrostatics and a cap on the LJ
///
/// @param		sim	a pointer to a simulation structure
/// @return 	void

//================================================================================
void ComputeForces (simptr sim)
//================================================================================
{
	int		i, nAtom;
	particleMD	*atom;

	// local sim variables
	atom  = sim->atom.items;
	nAtom = sim->atom.n;

	// reset all potential energies
	sim->potE  = 0;
	sim->ljE   = 0;
	sim->harmE = 0;
	sim->coulE = 0;
	sim->feneE = 0;
	sim->bendE = 0;
	sim->nemE = 0;

	// reset accelerations, forces, etc.
	for (i=0; i<nAtom; i++) {
		atom[i].ax	= 0;
		atom[i].ay	= 0;
		atom[i].az	= 0;
		#ifdef TEMPERATURE_CONF
		// configurational temperature
		atom[i].Tfx   = 0;
		atom[i].Tfy   = 0;
		atom[i].Tfz   = 0;
		atom[i].Tdivf = 0;
		#endif
	}

	// compute forces
	ComputeDispersionForces(sim);
	ComputeElectrostaticForces(sim);
	ComputeAnchorForces(sim);
	ComputeFeneForces(sim);
	ComputeBendForces(sim);
	ComputeDihedralForces(sim);
	ComputeSqueezeForces(sim);
}

/// Calculation of all forces acting on the atoms. First, the energy, accelearations
/// and force variables are all reset to 0, and then separate functions are called to
/// handle each type of interaction. Note that this is identical to above but no electrostatics and a cap on the LJ
///
/// @param		sim	a pointer to a simulation structure
/// @return 	void

//================================================================================
void ComputeCapForces (simptr sim)
//================================================================================
{
	int			i, nAtom;
	particleMD	*atom;

	// local sim variables
	atom  = sim->atom.items;
	nAtom = sim->atom.n;

	// reset all potential energies
	sim->potE  = 0;
	sim->ljE   = 0;
	sim->harmE = 0;
	sim->coulE = 0;
	sim->feneE = 0;
	sim->bendE = 0;
	sim->nemE = 0;

	// reset accelerations, forces, etc.
	for (i=0; i<nAtom; i++) {
		atom[i].ax	= 0;
		atom[i].ay	= 0;
		atom[i].az	= 0;
		#ifdef TEMPERATURE_CONF
		// configurational temperature
		atom[i].Tfx   = 0;
		atom[i].Tfy   = 0;
		atom[i].Tfz   = 0;
		atom[i].Tdivf = 0;
		#endif
	}

	// compute forces
	ComputeCapDispersionForces(sim);
	ComputeAnchorForces(sim);
	ComputeFeneForces(sim);
	ComputeBendForces(sim);
	ComputeDihedralForces(sim);
	ComputeSqueezeForces(sim);
}


/// Calculates the net short-range dispersion force (LJ) acting on each atom. Note
/// that the neighbors are stored as pointers to particle structures to avoid
/// double indirections (Steve Guillouzic, 2001). Separate atom lists are
/// used for STD neighbors and PBC neighbors, in order to avoid costly conditional
/// branches in the main loop (Frederic Tessier, 2002). We also calculate here
/// the DPD thermostat contribution (Soddemann, Dunweg, Kremer, PRE 68, 046702, 2003),
/// to avoid traversing the neigbour lists twice.
///
/// @param		sim	a pointer to a simulation structure
/// @return 	void
/// @note		The sigma variable includes the dt-dependent correction factor
///				to account for timestep artifacts in the DPD thermostat
///				(Peters, Europhys. Lett. 66, 311, 2004).

//================================================================================
void ComputeDispersionForces (simptr sim)
//================================================================================
{
	int	  		i, j, n;
	particleMD	*atom, *p1, *p2;
	real		dx, dy, dz;
	real		rCut2, ljShift;
	real		E=0, potE=0, ljE=0;
	real		kT, dtSqrti, eta, sigma;
	int		groupThermDPD;
	real		sigma_lj;

	// local sim variables
	sigma_lj		= sim->sigma_lj;
	ljShift	 		= sim->ljShift;
	rCut2	 		= sim->rCut*sim->rCut;
	kT		 		= sim->kT[sim->phase];
	dtSqrti	 		= 1.0/sqrt(sim->dt);
	eta		 		= sim->eta[sim->phase];
	sigma    		= sqrt (2*kT*eta*(1-eta*sim->dt/2.0));
	groupThermDPD 	= sim->groupThermDPD[sim->phase];
	atom	 = sim->atom.items;
	n 	 = sim->atom.n;

	// compute LJ forces for STD pairs
	#ifdef _OPENMP
	#pragma omp parallel  for \
		schedule  (static) \
		default   (none) \
		shared	  (n, nebrSTD, rCut2, ljShift) \
		private   (i, p1, p2, dx, dy, dz, E) \
		reduction (+: ljE, potE)
	#endif
	for (i=0; i<n; i++) {
		// Calculate force between all MD-MD particles
		for (j=0; j<i; j++) {
			p1 = atom+i;
			p2 = atom+j;
			// calculate dr
			dx = (p2->rx - p1->rx);
			dy = (p2->ry - p1->ry);
			dz = (p2->rz - p1->rz);
			// apply periodic boundary conditions
			ApplyPBC_dr (sim,&dx,&dy,&dz);
			// calculate LJ interaction
			dx /= sigma_lj;
			dy /= sigma_lj;
			dz /= sigma_lj;
			E = LennardJones (p1, p2, dx, dy, dz, rCut2) + ljShift;
			ljE  += E;
			potE += E;
			#ifdef THERMOSTAT_DPD
				// add fluctuation-dissipation forces (DPD thermostat)
	 			ThermostatDPD (p1, p2, dx, dy, dz, rCut2, eta, sigma, dtSqrti, groupThermDPD);
			#endif
		}
	}
	// update energies in sim structure
	sim->potE += potE;
	sim->ljE   = ljE;
}

/// Calculates the net short-range dispersion force (LJ) acting on each atom. Note
/// that the neighbors are stored as pointers to particle structures to avoid
/// double indirections (Steve Guillouzic, 2001). Separate atom lists are
/// used for STD neighbors and PBC neighbors, in order to avoid costly conditional
/// branches in the main loop (Frederic Tessier, 2002). We also calculate here
/// the DPD thermostat contribution (Soddemann, Dunweg, Kremer, PRE 68, 046702, 2003),
/// to avoid traversing the neigbour lists twice.
/// Note that in this implementation the LJ force with SRD particles is capped.
///
/// @param		sim	a pointer to a simulation structure
/// @return 	void
/// @note		The sigma variable includes the dt-dependent correction factor
///				to account for timestep artifacts in the DPD thermostat
///				(Peters, Europhys. Lett. 66, 311, 2004).

//================================================================================
void ComputeDispersionForcesSRD (simptr sim,particleMPC *pSRD,spec *SP,int GPOP)
//================================================================================
{
	int	  		i, j, n;
	particleMD	*atom, *p1, *p2;
	real		dx, dy, dz;
	real		rCut2, ljShift;
	real		E=0., potE=0., ljE=0.;
	real		kT, dtSqrti, eta, sigma;
	int		groupThermDPD;
	real		sigma_lj;

	// local sim variables
	sigma_lj		= sim->sigma_lj;
	ljShift	 		= sim->ljShift;
	rCut2	 		= sim->rCut*sim->rCut;
	kT		 		= sim->kT[sim->phase];
	dtSqrti	 		= 1.0/sqrt(sim->dt);
	eta		 		= sim->eta[sim->phase];
	sigma    		= sqrt (2*kT*eta*(1-eta*sim->dt/2.0));
	groupThermDPD 	= sim->groupThermDPD[sim->phase];
	atom	 = sim->atom.items;
	n 	 = sim->atom.n;

	// compute LJ forces for STD pairs
	#ifdef _OPENMP
	#pragma omp parallel  for \
		schedule  (static) \
		default   (none) \
		shared	  (n, nebrSTD, rCut2, ljShift) \
		private   (i, p1, p2, dx, dy, dz, E) \
		reduction (+: ljE, potE)
	#endif

	for (i=0; i<n; i++) {
		// Calculate force between all MD-MD particles
		for (j=0; j<i; j++) {
			p1 = atom+i;
			p2 = atom+j;
			// calculate dr
			dx = (p2->rx - p1->rx);
			dy = (p2->ry - p1->ry);
			dz = (p2->rz - p1->rz);
			// apply periodic boundary conditions
			ApplyPBC_dr (sim,&dx,&dy,&dz);
			// calculate LJ interaction
			dx /= sigma_lj;
			dy /= sigma_lj;
			dz /= sigma_lj;
			E = LennardJones (p1, p2, dx, dy, dz, rCut2) + ljShift;
			ljE  += E;
			potE += E;
			#ifdef THERMOSTAT_DPD
				// add fluctuation-dissipation forces (DPD thermostat)
				ThermostatDPD (p1, p2, dx, dy, dz, rCut2, eta, sigma, dtSqrti, groupThermDPD);
			#endif
		}
		// Calculate force between all MD-SRD particles
		for (j=0; j<GPOP; j++) {
			p1 = atom+i;
			// calculate dr
			dx = ( (pSRD+j)->Q[0] - p1->rx );
			dy = ( (pSRD+j)->Q[1] - p1->ry );
			dz = ( (pSRD+j)->Q[2] - p1->rz );
			// apply periodic boundary conditions
			ApplyPBC_dr (sim,&dx,&dy,&dz);
			// calculate LJ interaction
			dx /= sigma_lj;
			dy /= sigma_lj;
			dz /= sigma_lj;
			E = LennardJonesSRD (p1, &pSRD[j], SP, dx, dy, dz, sim->dt, rCut2) + ljShift;
			ljE  += E;
			potE += E;
		}
	}
	// update energies in sim structure
	sim->potE += potE;
	sim->ljE   = ljE;
}

/// Calculates the net short-range dispersion force (LJ) acting on each atom.
/// Note the neighbors are accesed via the SRD cell list.
/// This does NOT calculate forces between neighbours more than one SRD cell away
/// We also calculate here the DPD thermostat contribution (Soddemann, Dunweg, Kremer, PRE 68, 046702, 2003),
/// to avoid traversing the neigbour lists twice.
/// Note that in this implementation the LJ force with SRD particles is capped.
///
/// @param		sim	a pointer to a simulation structure
/// @return 	void
/// @note		The sigma variable includes the dt-dependent correction factor
///				to account for timestep artifacts in the DPD thermostat
///				(Peters, Europhys. Lett. 66, 311, 2004).

//================================================================================
void ComputeDispersionForcesSRDCell (simptr sim,int MDmode,struct particleMPC *pSRD,struct spec *SP,int GPOP,struct cell ***CL)
//================================================================================
{
	int	  	a,b,c,xbox,ybox,zbox,xc,yc,zc;
	particleMD	*p1, *p3;
	particleMPC	*p2;
	real		*box;
	real		dx, dy, dz;
	real		rCut2, ljShift;
	real		E=0., potE=0., ljE=0.;
	real		kT, dtSqrti, eta, sigma;
	int		groupThermDPD;
	real		sigma_lj, drTotMax;
	int		pbcx1, pbcy1,pbcz1,pbcx2,pbcy2,pbcz2;

	// local sim variables
	sigma_lj		= sim->sigma_lj;
	ljShift	 	= sim->ljShift;
	rCut2	 		= sim->rCut*sim->rCut;
	kT			= sim->kT[sim->phase];
	dtSqrti		= 1.0/sqrt(sim->dt);
	eta			= sim->eta[sim->phase];
	sigma    		= sqrt (2*kT*eta*(1-eta*sim->dt/2.0));
	groupThermDPD 	= sim->groupThermDPD[sim->phase];
	box			= sim->box;
	drTotMax = sim->drTotMax;

	#ifdef _OPENMP
	#pragma omp parallel  for \
		schedule  (static) \
		default   (none) \
		shared	  (n, nebrSTD, rCut2, ljShift) \
		private   (i, p1, p2, dx, dy, dz, E) \
		reduction (+: ljE, potE)
	#endif

	for( a=0; a<=box[x_]; a++ ) for( b=0; b<=box[y_]; b++ ) for( c=0; c<=box[z_]; c++ ) {
		// MD-MD bead interactions
		if (CL[a][b][c].MDpp!=NULL) {
			p1 = CL[a][b][c].MDpp;
			while (p1!=NULL) {
				p3 = p1->nextSRD;
				while (p3!=NULL) {
					// calculate dr
					dx = ( p3->rx - p1->rx );
					dy = ( p3->ry - p1->ry );
					dz = ( p3->rz - p1->rz );
					// apply periodic boundary conditions
					ApplyPBC_dr (sim,&dx,&dy,&dz);
					// calculate LJ interaction
					dx /= sigma_lj;
					dy /= sigma_lj;
					dz /= sigma_lj;
					E = LennardJones (p1, p3, dx, dy, dz, rCut2) + ljShift;
					E = LennardJones (p1, p3, dx, dy, dz, rCut2) + ljShift;
					ljE  += E;
					potE += E;
					#ifdef THERMOSTAT_DPD
						// add fluctuation-dissipation forces (DPD thermostat)
						ThermostatDPD (p1, p3, dx, dy, dz, rCut2, eta, sigma, dtSqrti, groupThermDPD);
					#endif
					p3=p3->nextSRD;
				}

				// Account for periodic boundaries
				//note since xc, yc, and zc are ints the floor is implicit for positive values
				//X-PBCs
				pbcx1 = 0.0;//these pbc variables account for the fact that box 10 is box 0 meaning we need to add 1 to the halo
				if (a-sigma_lj-drTotMax<1) pbcx1 = 1.0;
				pbcx2 = 0.0;
				if (a+sigma_lj+drTotMax>box[x_]-1) pbcx2 = 1.0;
				for (xc=floor(a-sigma_lj-drTotMax-pbcx1); xc<a+sigma_lj+drTotMax+pbcx2+2; xc++) {
					xbox=xc;
					if (xbox>box[x_]) {
						if (sim->pbcond & PBC_COND_x) while (xbox>box[x_]) xbox = xbox-box[x_]-1;
						//Otherwise skip this xc
						else continue;
					}
					if (xbox<0) {
						if (sim->pbcond & PBC_COND_x) while (xbox<0) xbox = xbox+box[x_]+1;
						//Otherwise skip this xc
						else continue;
					}
					//Y-PBCs
					pbcy1 = 0.0;
					if (b-sigma_lj-drTotMax<1) pbcy1 = 1.0;
					pbcy2 = 0.0;
					if (b+sigma_lj+drTotMax>box[y_]-1) pbcy2 = 1.0;
					for (yc=floor(b-sigma_lj-drTotMax-pbcy1);yc<b+sigma_lj+drTotMax+pbcy2+2;yc++) {
						ybox=yc;
						if (ybox>box[y_]) {
							if (sim->pbcond & PBC_COND_y) while (ybox>box[y_]) ybox = ybox-box[y_]-1;
							//Otherwise skip this yc
							else continue;
						}
						if (ybox<0) {
							if (sim->pbcond & PBC_COND_y) while (ybox<0) ybox = ybox+box[y_]+1;
							//Otherwise skip this yc
							else continue;
						}
						//Z-PBCs
						pbcz1 = 0.0;
						if (c-sigma_lj-drTotMax<1) pbcz1 = 1.0;
						pbcz2 = 0.0;
						if (c+sigma_lj+drTotMax>box[z_]-1) pbcz2 = 1.0;
						for (zc=floor(c-sigma_lj-drTotMax-pbcz1);zc<c+sigma_lj+drTotMax+pbcz2+2;zc++) {
							zbox=zc;
							if (zbox>box[z_]) {
								if (sim->pbcond & PBC_COND_z) while (zbox>box[z_]) zbox = zbox-box[z_]-1;
								//Otherwise skip this zc
								else continue;
							}
							if (zbox<0) {
								if (sim->pbcond & PBC_COND_z) while (zbox<0) zbox = zbox+box[z_]+1;
								//Otherwise skip this zc
								else continue;
							}
							//Only go forward if a change was made
							if (xbox==a && ybox==b && zbox==c) continue;
							//Some periodic BCs were triggered so keep going
							p3 = CL[xbox][ybox][zbox].MDpp;
							while (p3!=NULL) {
								// calculate dr
								dx = ( p3->rx - p1->rx );
								dy = ( p3->ry - p1->ry );
								dz = ( p3->rz - p1->rz );
								// apply periodic boundary conditions
								ApplyPBC_dr (sim,&dx,&dy,&dz);
								// calculate LJ interaction
								dx /= sigma_lj;
								dy /= sigma_lj;
								dz /= sigma_lj;
								E = LennardJones (p1, p3, dx, dy, dz, rCut2) + ljShift;
								ljE  += E;
								potE += E;
								#ifdef THERMOSTAT_DPD
								// add fluctuation-dissipation forces (DPD thermostat)
								ThermostatDPD (p1, p3, dx, dy, dz, rCut2, eta, sigma, dtSqrti, groupThermDPD);
								#endif
								p3=p3->nextSRD;
							}
						}
					}
				}

				//MPC-MD bead interactions
				if (MDmode == MPCinMD) {
					pbcx1 = 0.0;//these pbc variables account for the fact that box 10 is box 0 meaning we need to add 1 to the halo
					if (a-sigma_lj-drTotMax<1) pbcx1 = 1.0;
					pbcx2 = 0.0;
					if (a+sigma_lj+drTotMax>box[x_]-1) pbcx2 = 1.0;
					for (xc=floor(a-sigma_lj-drTotMax-pbcx1);xc<a+sigma_lj+drTotMax+pbcx2+2;xc++) {
						xbox=xc;
						if (xbox>box[x_]) {
							if (sim->pbcond & PBC_COND_x) while (xbox>box[x_]) xbox = xbox-box[x_]-1;
							else continue;
						}
						if (xbox<0) {
							if (sim->pbcond & PBC_COND_x) while (xbox<0) xbox = xbox+box[x_]+1;
							else continue;
						}
						pbcy1 = 0.0;
						if (b-sigma_lj-drTotMax<1) pbcy1 = 1.0;
						pbcy2 = 0.0;
						if (b+sigma_lj+drTotMax>box[y_]-1) pbcy2 = 1.0;
						for (yc=floor(b-sigma_lj-drTotMax-pbcy1);yc<b+sigma_lj+drTotMax+pbcy2+2;yc++) {
							ybox=yc;
							if (ybox>box[y_]) {
								if (sim->pbcond & PBC_COND_y) while (ybox>box[y_]) ybox = ybox-box[y_]-1;
								else continue;
							}
							if (ybox<0) {
								if (sim->pbcond & PBC_COND_y) while (ybox<0) ybox = ybox+box[y_]+1;
								else continue;
							}
							pbcz1 = 0.0;
							if (c-sigma_lj-drTotMax<1) pbcz1 = 1.0;
							pbcz2 = 0.0;
							if (c+sigma_lj+drTotMax>box[z_]-1) pbcz2 = 1.0;
							for (zc=floor(c-sigma_lj-drTotMax-pbcz1);zc<c+sigma_lj+drTotMax+pbcz2+2;zc++) {
								zbox=zc;
								if (zbox>box[z_]) {
									if (sim->pbcond & PBC_COND_z) while (zbox>box[z_]) zbox = zbox-box[z_]-1;
									else continue;
								}
								if (zbox<0) {
									if (sim->pbcond & PBC_COND_z) while (zbox<0) zbox = zbox+box[z_]+1;
									else continue;
								}
								p2 = CL[xbox][ybox][zbox].pp;
								while (p2!=NULL) {
									// calculate dr
									dx = ( p2->Q[0] - p1->rx );
									dy = ( p2->Q[1] - p1->ry );
									dz = ( p2->Q[2] - p1->rz );
									// apply periodic boundary conditions
									ApplyPBC_dr (sim,&dx,&dy,&dz);
									// calculate LJ interaction
									dx /= sigma_lj;
									dy /= sigma_lj;
									dz /= sigma_lj;
									E = LennardJonesSRD (p1, p2, SP, dx, dy, dz, sim->dt, rCut2) + ljShift;
									ljE  += E;
									potE += E;
									p2=p2->next;
								}
							}
						}
					}
				}
				p1 = p1->nextSRD;
			}
		}
	}
	// update energies in sim structure
	sim->potE += potE;
	sim->ljE   = ljE;
}


/// Calculates the net short-range dispersion force (LJ) acting on each atom. Note
/// that the neighbors are stored as pointers to particle structures to avoid
/// double indirections (Steve Guillouzic, 2001). Separate atom lists are
/// used for STD neighbors and PBC neighbors, in order to avoid costly conditional
/// branches in the main loop (Frederic Tessier, 2002). We also calculate here
/// the DPD thermostat contribution (Soddemann, Dunweg, Kremer, PRE 68, 046702, 2003),
/// to avoid traversing the neigbour lists twice.
/// Note that in this implementation the LJ force is capped.
///
/// @param		sim	a pointer to a simulation structure
/// @return 	void
/// @note		The sigma variable includes the dt-dependent correction factor
///				to account for timestep artifacts in the DPD thermostat
///				(Peters, Europhys. Lett. 66, 311, 2004).

//================================================================================
void ComputeCapDispersionForces (simptr sim)
//================================================================================
{
	int	  		i,j,n;
	particleMD	*atom, *p1, *p2;
	real		dx, dy, dz;
	real		rCut2, ljShift;
	real		E=0, potE=0;
	real		r2min=10;
	real		sigma_lj;

	// local sim variables
	sigma_lj		= sim->sigma_lj;
	ljShift	 	= sim->ljShift;
	rCut2	 		= sim->rCut*sim->rCut;
	atom	 = sim->atom.items;
	n 	 = sim->atom.n;

	// compute LJ forces for STD pairs
	#ifdef _OPENMP
	#pragma omp parallel  for \
		schedule  (static) \
		default   (none) \
		shared	  (n, nebrSTD, rCut2, ljShift) \
		private   (i, p1, p2, dx, dy, dz, E) \
		reduction (+: r2min, potE)
	#endif

	for (i=0; i<n; i++) {
		for (j=0; j<i; j++) {
			p1 = atom+i;
			p2 = atom+j;
			// calculate dr
			dx = (p2->rx - p1->rx);
			dy = (p2->ry - p1->ry);
			dz = (p2->rz - p1->rz);
			// apply periodic boundary conditions
			ApplyPBC_dr (sim,&dx,&dy,&dz);
			// calculate LJ interaction
			dx /= sigma_lj;
			dy /= sigma_lj;
			dz /= sigma_lj;
			E = LennardJonesCap (p1, p2, dx, dy, dz, rCut2, sim) + ljShift;
			if (E<r2min) r2min = E;
			potE += E;
		}
	}
	// update energies in sim structure
	sim->potE += potE;
	sim->ljE   = r2min;
}

/// Computes all the forces due to electrostatic interactions between
/// charges. See notes in ComputeLJForces above. The coulomb interaction
/// variables are calculated here for clarity. The function involves a sum
/// over all charge pairs, and thus scales as nCharge^2. We use a simple
/// interaction cutoff at a distance rCutCoul, which is ok in 1-dimensional
/// systems only (see Allen and Tildesley).
///
/// @param		sim	a pointer to a simulation structure
/// @return 	void
/// @warning    Don't include electric force in sum(f2) for configurational
///				temperature (Denis J. Evans, by email: it's a dissipative force. But
///				there is still confusion regarding Delhomelle's comments, 2003). But
///				if the field is small enough this is negligible.
/// @note		We don't call the ApplyPBC and rather calculate it inline
///				to save some time.
/// @see		ComputeDispersionForces for an explanation of the PBC and STD
///				list handling.

//================================================================================
void ComputeElectrostaticForces (simptr sim)
//================================================================================
{
	int			i, j, n;
	particleMD	*p1, *p2;
	item1STD	*charge;
// 	real		*box, *boxHalf;
	real		dx, dy, dz;
	real		rCutCoul2=0, bjerrumkT=0, lambda_D=0;
	real		E=0, potE=0, coulE=0, Efield=0, Eforce=0;

	// calculate coulomb interaction variables (here for clarity)
	if (sim->charge.n > 0) {
		rCutCoul2  = sim->rCutCoul * sim->rCutCoul;
	} else return;

	// local sim variables
	charge  	= sim->charge.items;
	bjerrumkT   = sim->bjerrum * sim->kT[sim->phase];
	lambda_D	= sim->lambda_D;
	Efield  	= sim->Efield[sim->phase];
	// compute coulomb forces for ALL charge pairs
	n = sim->charge.n;

	#ifdef _OPENMP
	#pragma omp parallel  for \
		schedule  (static) \
		default   (none) \
		shared	(n, charge, rCutCoul2, bjerrumkT, sim) \
		private   (i, j, p1, p2, dx, dy, dz, E, Efield, Eforce) \
		reduction (+: coulE, potE)
	#endif

	for (i=0; i<n; i++) {
		// extract first particle pointer
		p1 = charge[i].p1;
		for (j=0; j<i; j++) {
			// extract second particle pointer
			p2 = charge[j].p1;
			// calculate dr
			dx = (p2->rx - p1->rx);
			dy = (p2->ry - p1->ry);
			dz = (p2->rz - p1->rz);
			// apply periodic boundary conditions
			ApplyPBC_dr (sim,&dx,&dy,&dz);
			// calculate Coulomb interaction
			E = Coulomb (p1, p2, dx, dy, dz, rCutCoul2, bjerrumkT, lambda_D);
			coulE += E;
			potE  += E;
		}
		// driving electric field (along x)
		Eforce  = p1->q * Efield;
		p1->ax += Eforce;
	}
	// add potential energies in sim structure
	sim->potE  += potE;
	sim->coulE  = coulE;
}


/// Computes the restoring forces acting on anchored atoms. We do not need to
/// separate STD and PBC anchored atoms because we use the anchor point x0, y0,
/// z0 and keep track of the real world position of the atom wx, wy, wz.
///
/// @param		sim	a pointer to a simulation structure
/// @return		void

//================================================================================
void ComputeAnchorForces (simptr sim)
//================================================================================
{
	int			i, n;
	particleMD	*p1;
	item1STD	*anchor;
	real		E=0, potE=0, harmE=0;
	real		dx, dy, dz;

	// local sim variables
	anchor = sim->anchor.items;

	// compute forces for anchor atoms
	n = sim->anchor.n;

	#ifdef _OPENMP
	#pragma omp parallel  for \
		schedule  (static) \
		default   (none) \
		shared	(n, anchor) \
		private   (i, p1, dx, dy, dz, E) \
		reduction (+: harmE, potE)
	#endif

	for (i=0; i<n; i++) {
		// extract particle pointer
		p1 = anchor[i].p1;
		// compute dr
		dx = p1->wx - p1->x0;
		dy = p1->wy - p1->y0;
		dz = p1->wz - p1->z0;
		// calculate harmonic interaction
		E = HarmonicSpring (p1, dx, dy, dz, p1->kspring);
		harmE += E;
		potE  += E;
	}
	// add potential energies in sim structure
	sim->potE  += potE;
	sim->harmE  = harmE;
}


/// Computes all the forces between FENE-bonded atom pairs. We compute the dr from
/// the WORLD positions, therefore we do not ever have to worry about the PBC
/// here! But the world coordinates of the FENE pairs MUST be initialized
/// accordingly.
///
/// @param		sim	a pointer to a simulation structure
/// @return 	void
/// @warning	Real-world coordinates of the fene pairs MUST be initialized correctly
///				because we don't consider the PBC in the FENE calculation.

//================================================================================
void ComputeFeneForces (simptr sim)
//================================================================================
{
	int	  		i, nFene;
	particleMD	*p1, *p2;
	item2STD	*fene;
	real		E=0, potE=0, feneE=0, kFene, r0Fene;
	real		dx, dy, dz;

	// local sim variables
	fene	= sim->fene.items;
	nFene   = sim->fene.n;
	kFene   = sim->kFene;
	r0Fene  = sim->r0Fene;

	// loop over fene pairs
	for (i=0; i<nFene; i++) {
		// extract pair pointers
		p1 = fene[i].p1;
		p2 = fene[i].p2;
		// compute dr (using WORLD positions)
		dx = p2->wx - p1->wx;
		dy = p2->wy - p1->wy;
		dz = p2->wz - p1->wz;
		// calculate FENE interaction
		E = FENE (p1, p2, dx, dy, dz, kFene, r0Fene);
		feneE += E;
		potE  += E;
	}
	// add potential energies in sim structure
	sim->potE  += potE;
	sim->feneE  = feneE;
}

/// Computes all the forces between BEND-bonded atom pairs. We compute the dr from
/// the WORLD positions, therefore we do not ever have to worry about the PBC
/// here! But the world coordinates of the FENE pairs MUST be initialized
/// accordingly.
///
/// @param		sim	a pointer to a simulation structure
/// @param		sim	a pointer to a simulation structure
/// @return 	void
/// @warning	Real-world coordinates of the fene pairs MUST be initialized correctly
///				because we don't consider the PBC in the BEND calculation.

//================================================================================
void ComputeBendForces (simptr sim)
//================================================================================
{
	int	  		i, nBend;
	particleMD	*p1, *p2, *p3;
	item3STD	*bend;
	real		E=0, potE=0, bendE=0, kBend=0, theta0=0;
	real		dx12, dy12, dz12, dx23, dy23, dz23;

	// local sim variables
	bend	= sim->bend.items;
	nBend   = sim->bend.n;
	kBend   = sim->kBend;
	theta0	= sim->theta0Bend;

	if(kBend>TOL) {
		// loop over bend pairs
		for (i=0; i<nBend; i++) {
			// extract pair pointers
			p1 = bend[i].p1;
			p2 = bend[i].p2;
			p3 = bend[i].p3;
			// compute dr (using WORLD positions)
			dx12 = p1->wx - p2->wx;
			dy12 = p1->wy - p2->wy;
			dz12 = p1->wz - p2->wz;
			dx23 = p2->wx - p3->wx;
			dy23 = p2->wy - p3->wy;
			dz23 = p2->wz - p3->wz;
			// calculate harmonic three-particle bend interaction
			E = bendHarmonic (p1, p2, p3, dx12, dy12, dz12, dx23, dy23, dz23, kBend, theta0);
			bendE += E;
			potE  += E;
		}
	}

	// add potential energies in sim structure
	sim->potE  += potE;
	sim->bendE  = bendE;
}

/// Computes all the forces between DIHEDRAL-bonded atom pairs. We compute the dr from
/// the WORLD positions, therefore we do not ever have to worry about the PBC
/// here! But the world coordinates of the FENE pairs MUST be initialized
/// accordingly.
///
/// @param		sim	a pointer to a simulation structure
/// @param		sim	a pointer to a simulation structure
/// @return 	void
/// @warning	Real-world coordinates of the fene pairs MUST be initialized correctly
///				because we don't consider the PBC in the BEND calculation.

//================================================================================
void ComputeDihedralForces (simptr sim)
//================================================================================
{
	int	  		i, nDihedral;
	particleMD	*p1, *p2, *p3, *p4;
	item4STD	*dihedral;
	real		E=0, potE=0, dihedralE=0, kDihedral=0, phi0=0;
	real		dx12, dy12, dz12, dx23, dy23, dz23, dx34, dy34, dz34;

	// local sim variables
	dihedral	= sim->dihedral.items;
	nDihedral   = sim->dihedral.n;
	kDihedral   = sim->kDihedral;
	phi0		= sim->phi0Dihedral;
	if(kDihedral>=TOL) {
		// loop over dihedral pairs
		for (i=0; i<nDihedral; i++) {
			// extract pair pointers
			p1 = dihedral[i].p1;
			p2 = dihedral[i].p2;
			p3 = dihedral[i].p3;
			p4 = dihedral[i].p4;
			// compute dr (using WORLD positions)
			dx12 = p2->wx - p1->wx;
			dy12 = p2->wy - p1->wy;
			dz12 = p2->wz - p1->wz;
			dx23 = p3->wx - p2->wx;
			dy23 = p3->wy - p2->wy;
			dz23 = p3->wz - p2->wz;
			dx34 = p4->wx - p3->wx;
			dy34 = p4->wy - p3->wy;
			dz34 = p4->wz - p3->wz;
			// calculate harmonic four-particle dihedral interaction
			E = dihedralHarmonic (p1, p2, p3, p4, dx12, dy12, dz12, dx23, dy23, dz23, dx34, dy34, dz34, kDihedral, phi0);
			dihedralE += E;
			potE  += E;
		}
	}

	// add potential energies in sim structure
	sim->potE  += potE;
	sim->dihedralE  = dihedralE;
}

/// Computes all the forces between the bond defined by FENE-bonded atom pairs
/// and the background nematic solvent. We compute the dr from
/// the WORLD positions, therefore we do not ever have to worry about the PBC
/// here! But the world coordinates of the FENE pairs MUST be initialized
/// accordingly.
///
/// @param		sim	a pointer to a simulation structure
/// @return 	void
/// @warning	Real-world coordinates of the fene pairs MUST be initialized correctly
///				because we don't consider the PBC in the FENE calculation.

//================================================================================
void ComputeNemForces (simptr sim,struct spec *SP,struct cell ***CL)
//================================================================================
{
	int	  		i,j,k,l, nBend;
	particleMD	*p1, *p2;
	item2STD	*bend;
	real		E=0, potE=0, nemE=0, kNem, dt;
	real		dxMD, dyMD, dzMD, dxMPCD, dyMPCD, dzMPCD, S;

	// local sim variables. Use fene as an already established convenient set of nearest-neighbors
	bend	= sim->fene.items;
	nBend   = sim->fene.n;
	kNem   = sim->kNemMPC;
	dt	= sim->dt;

	if(kNem>0.0) {
		for (i=0; i<nBend; i++) {
		// loop over bend pairs
			// extract pair pointers
			p1 = bend[i].p1;
			p2 = bend[i].p2;
			// compute dr (using WORLD positions)
			dxMD = p1->wx - p2->wx;
			dyMD = p1->wy - p2->wy;
			dzMD = p1->wz - p2->wz;
			//Find what MPCD cell MD bead i occupies
			j = (int)p1->rx;
			k = (int)p1->ry;
			l = (int)p1->rz;
			//Set the bond vector in that directions
			dxMPCD = CL[j][k][l].DIR[0];
			dyMPCD = CL[j][k][l].DIR[1];
			dzMPCD = CL[j][k][l].DIR[2];
			S = CL[j][k][l].S;
			// calculate BEND interaction
			if (fabs(kNem*S)>=TOL && (dxMPCD*dxMPCD+dyMPCD*dyMPCD+dzMPCD*dzMPCD)>=TOL) {
				E = bendNematic (p1, p2, dxMD, dyMD, dzMD, dxMPCD, dyMPCD, dzMPCD, kNem, S, dt, SP, &CL[j][k][l]);
				nemE += E;
				potE  += E;
			}
		}
	}

	// add potential energies in sim structure
	sim->potE  += potE;
	sim->nemE  += nemE;
}

/// Compute radial squeezing force to force polymer into longitudinal config.
//================================================================================
void ComputeSqueezeForces (simptr sim)
//================================================================================

{
    itemPoly    *polymer;
    particleMD    *p1;
    real        E=0, potE=0, kSqu;
    real        dy, dz;
    real        *boxHalf;

    // local sim variables
    kSqu = sim->kSqu;
    polymer  = sim->polymer.items;
    p1 = polymer[0].p1;
    boxHalf = sim->boxHalf;
    // loop over monomers
    while (p1) {
        dy = p1->wy - boxHalf[y_];		//The well is centred on the middle of the box
        dz = p1->wz - boxHalf[z_];
				// Each of the different tube potential types
	 			E = HarmonicSpringTube(p1, dy, dz, kSqu);
				// E = QuadraticSpringTube(p1, dy, dz, kSqu);
				// E = SexticSpringTube(p1, dy, dz, kSqu);
				// E = LennardJonesTube(p1, dy, dz, kSqu, sim);
				potE += E;
				p1 = p1->next;
    }
    // add potential energies in sim structure
    sim->potE  += potE;
}

/// Compute the LJ interaction between an atom pair. The ljShift is not added
/// here, but in the calling function. Divergence of the LJ force is also calculated
/// for the configurational temperature calculation (- div f is stored in each
/// particle's Tdivf variable) if it's compiled in.
///
/// @param		p1 pointer to the first atom
/// @param		p2 pointer to the second atom
/// @param		dx,dy,dz position difference between the two atoms
/// @param		rCut2 square of the LJ cutoff distance
/// @return		the value of the LJ interaction energy
/// @warning	The ljShift must be added in the parent function.

//================================================================================
// inline real LennardJones (particleMD *p1, particleMD *p2, real dx, real dy, real dz,
// 						  real rCut2)
real LennardJones (particleMD *p1, particleMD *p2, real dx, real dy, real dz,
						  real rCut2)
//================================================================================
{
	real r2, r2i, r6i, fMag, fx, fy, fz;
	real potE=0;
	#ifdef TEMPERATURE_CONF
		real Tdivf;
	#endif

	// calculate the distance squared
	r2 = dx*dx + dy*dy + dz*dz;
	// compute the force and the energy if we are inside the cutoff range
	if (r2 < rCut2) {
		r2i   = 1/r2;
		r6i   = r2i*r2i*r2i;
		potE  = 4*r6i*(r6i-1);  
		fMag  = 48 * r2i * r6i;
		#ifdef TEMPERATURE_CONF
			Tdivf = fMag * (11*r6i - 2.5);
		#endif
		fMag  = fMag * (r6i - 0.5);
		fx = fMag * dx;
		fy = fMag * dy;
		fz = fMag * dz;
		p1->ax -= fx;
		p1->ay -= fy;
		p1->az -= fz;
		p2->ax += fx;
		p2->ay += fy;
		p2->az += fz;
		#ifdef TEMPERATURE_CONF
			// configurational temperature
			p1->Tfx -= fx;
			p1->Tfy -= fy;
			p1->Tfz -= fz;
			p2->Tfx += fx;
			p2->Tfy += fy;
			p2->Tfz += fz;
			p1->Tdivf += Tdivf;
			p2->Tdivf += Tdivf;
		#endif
	}

	return potE;
}

//================================================================================
real LennardJonesSRD (particleMD *p1, struct particleMPC *p2,struct spec *SP, real dx, real dy, real dz, real dt,
						  real rCut2)
//================================================================================
{
	real r2, r2i, r6i, fMag, fx, fy, fz;
	real potE=0;
	real mass = ( SP + (p2->SPID) )->MASS;
	#ifdef TEMPERATURE_CONF
		real Tdivf;
	#endif

	// calculate the distance squared
	r2 = dx*dx + dy*dy + dz*dz;
	// compute the force and the energy if we are inside the cutoff range
	if (r2 < rCut2) {
		r2i   = 1/r2;
		r6i   = r2i*r2i*r2i;
		potE  = 4*r6i*(r6i-1);
		fMag  = 48 * r2i * r6i;
		#ifdef TEMPERATURE_CONF
			Tdivf = fMag * (11*r6i - 2.5);
		#endif
		fMag  = fMag * (r6i - 0.5);
		fx = fMag * dx;
		fy = fMag * dy;
		fz = fMag * dz;
		p1->ax -= fx;
		p1->ay -= fy;
		p1->az -= fz;
		p2->V[0] += dt*fx/mass;
		p2->V[1] += dt*fy/mass;
		p2->V[2] += dt*fz/mass;
		#ifdef TEMPERATURE_CONF
		// configurational temperature
		p1->Tfx -= fx;
		p1->Tfy -= fy;
		p1->Tfz -= fz;
		p1->Tdivf += Tdivf;
		#endif
	}

	return potE;
}

/// Compute the LJ interaction between an atom pair. The ljShift is not added
/// here, but in the calling function. Divergence of the LJ force is also calculated
/// for the configurational temperature calculation (- div f is stored in each
/// particle's Tdivf variable) if it's compiled in.
/// Note that in this implementation the force is capped according to relaxStep
///
/// @param		p1 pointer to the first atom
/// @param		p2 pointer to the second atom
/// @param		dx,dy,dz position difference between the two atoms
/// @param		rCut2 square of the LJ cutoff distance
/// @return		the value of the LJ interaction energy
/// @warning	The ljShift must be added in the parent function.

//================================================================================
// inline real LennardJonesCap (particleMD *p1, particleMD *p2, real dx, real dy, real dz,
// 						  real rCut2, simptr sim)
real LennardJonesCap (particleMD *p1, particleMD *p2, real dx, real dy, real dz,
						  real rCut2, simptr sim)
//================================================================================
{
	int  stepRelax;
	real r2, r2i, r6i, fMag, fx, fy, fz;
	real forceCap;
// 	real potE=0;
	#ifdef TEMPERATURE_CONF
		real Tdivf;
	#endif

	stepRelax = sim->stepRelax;
	// forceCap = 10000./stepRelax;
	forceCap = (1010-stepRelax);
	if (stepRelax == 1) forceCap = 10000000.0;
	// calculate the distance squared
	r2 = dx*dx + dy*dy + dz*dz;

	// compute the force and the energy if we are inside the cutoff range
	if (r2 < rCut2) {
		r2i   = 1/r2;
		r6i   = r2i*r2i*r2i;
		// potE  = 4*r6i*(r6i-1);
		fMag  = 48 * r2i * r6i;
		#ifdef TEMPERATURE_CONF
			Tdivf = fMag * (11*r6i - 2.5);
		#endif
		fMag  = fMag * (r6i - 0.5);
		if (fMag>forceCap && p1->type != TYPE_WALL && p2->type != TYPE_WALL) {
			fMag = forceCap;
		}
		if (fMag>forceCap && (p1->type == TYPE_WALL || p2->type == TYPE_WALL)) {
			if (fMag>1000. && fMag>forceCap) fMag = 1000.;
		}
		fx = fMag * dx;
		fy = fMag * dy;
		fz = fMag * dz;
		p1->ax -= fx;
		p1->ay -= fy;
		p1->az -= fz;
		p2->ax += fx;
		p2->ay += fy;
		p2->az += fz;
		#ifdef TEMPERATURE_CONF
			// configurational temperature
			p1->Tfx -= fx;
			p1->Tfy -= fy;
			p1->Tfz -= fz;
			p2->Tfx += fx;
			p2->Tfy += fy;
			p2->Tfz += fz;
			p1->Tdivf += Tdivf;
			p2->Tdivf += Tdivf;
		#endif
	}
	if (r2 < 0.0005) {
		p1->ax +=0.5;
		p1->ay +=0.5;
		p1->az +=0.5;
	}

	return r2;
}


/// Compute the coulomb interaction between two charged atoms. Normally it is
/// incorrect in 3D to cutoff the Coulomb interaction, but for 1D problems
/// (e.g., a capillary), the interaction is short-ranged and we can use a cutoff
/// (See Allen and Tildesley). Divergence of the coulomb force is 0, so we don't
/// need to consider Tdivf for the configurational temperature calculation.
///
/// @param		p1 pointer to the first atom
/// @param		p2 pointer to the second atom
/// @param		dx,dy,dz position difference between the two atoms
/// @param		rCutCoul2 square of the Coulomb cutoff distance
/// @param		bjerrumkT the bjerrum length mulitplied by the thermal energy kT
/// @return		the value of the coulomb interaction energy
/// @warning	Don't use this function with 2D or 3D periodic boundaries.
/// @todo		Implement Ewald type methods to handle 2D and 3D systems correctly.

//================================================================================
// inline real Coulomb (particleMD *p1, particleMD *p2,  real dx, real dy, real dz,
// 					 real rCutCoul2, real bjerrumkT, real lambda_D)
real Coulomb (particleMD *p1, particleMD *p2,  real dx, real dy, real dz,
					 real rCutCoul2, real bjerrumkT, real lambda_D)
//================================================================================
{
	real r, r2;
	real q1q2, fx, fy, fz, fMag;
	real potE=0;

	// calculate the distance
	r2 = dx*dx + dy*dy + dz*dz;
	// compute the force and the energy below cutoff
	if (r2 <= rCutCoul2) {
		r = sqrt(r2);
		q1q2 = bjerrumkT * p1->q * p2->q;
		potE = q1q2 * (exp(-r/lambda_D)/r);
		fMag = q1q2 * (exp(-r/lambda_D))/ (r2*r) + q1q2 * (exp(-r/lambda_D))/ (r2*lambda_D);
		fx = fMag * dx;
		fy = fMag * dy;
		fz = fMag * dz;
		p1->ax -= fx;
		p1->ay -= fy;
		p1->az -= fz;
		p2->ax += fx;
		p2->ay += fy;
		p2->az += fz;
		#ifdef TEMPERATURE_CONF
			// configurational temperature
			p1->Tfx -= fx;
			p1->Tfy -= fy;
			p1->Tfz -= fz;
			p2->Tfx += fx;
			p2->Tfy += fy;
			p2->Tfz += fz;
		#endif
	}

	return potE;
}


/// Compute the harmonic force on an atom, given its diplacement from its rest
/// position. Divergence of the harmonic force is -3*kspring, for the configurational
/// temperature calculation.
///
/// @param		p1 pointer to the atom
/// @param		dx,dy,dz position difference between atom and the centre of the tube
/// @param		k the harmonic spring constant
///
/// @return		the value of the harmonic interaction energy

//================================================================================
// inline real HarmonicSpring (particleMD *p1, real dx, real dy, real dz, real k)
real HarmonicSpring (particleMD *p1, real dx, real dy, real dz, real k)
//================================================================================
{
	real r2, fx, fy, fz;
	real potE=0;

	// calculate the distance squared
	r2 = dx*dx + dy*dy + dz*dz;
	// compute the force and the energy
	potE = (0.5*k) * r2;
	fx = k * dx;
	fy = k * dy;
	fz = k * dz;
	p1->ax -= fx;
	p1->ay -= fy;
	p1->az -= fz;
	#ifdef TEMPERATURE_CONF
		// configurational temperature
		p1->Tfx -= fx;
		p1->Tfy -= fy;
		p1->Tfz -= fz;
		p1->Tdivf += 3 * k;
	#endif

	return potE;
}

/// Compute a 2D (cylindrical) harmonic force on an atom, given its diplacement from its rest
/// axis. Divergence of the harmonic force is -3*kspring, for the configurational
/// temperature calculation.
///
/// @param		p1 pointer to the atom
/// @param		dx,dy,dz position difference between atom and the centre of the tube
/// @param		k the harmonic spring constant
///
/// @return		the value of the harmonic interaction energy

//================================================================================
// inline real HarmonicSpringTube (particleMD *p1, real dy, real dz, real k)
real HarmonicSpringTube (particleMD *p1, real dy, real dz, real k)
//================================================================================
{
    real r2, fy, fz;
    real potE=0;

    // calculate the distance squared
    r2 = dy*dy + dz*dz;
    // compute the force and the energy
    potE = (0.5*k) * r2;
    fy = k * dy;
    fz = k * dz;
    p1->ay -= fy;
    p1->az -= fz;
    #ifdef TEMPERATURE_CONF
	    // configurational temperature
	    p1->Tfy -= fy;
	    p1->Tfz -= fz;
	    p1->Tdivf += 2 * k;
    #endif

    return potE;
}

/// Compute a 2D (cylindrical) quartic force on an atom, given its diplacement from its rest
/// axis. Divergence of the harmonic force is -3*kspring, for the configurational
/// temperature calculation.
///
/// @param		p1 pointer to the atom
/// @param		dx,dy,dz position difference between atom and the centre of the tube
/// @param		k the harmonic spring constant
///
/// @return		the value of the harmonic interaction energy

//================================================================================
// inline real QuadraticSpringTube (particleMD *p1, real dy, real dz, real k)
real QuadraticSpringTube (particleMD *p1, real dy, real dz, real k)
//================================================================================
{
    real r2,r4, fy, fz;
    real potE=0;

    // calculate the distance squared
    r2 = dy*dy + dz*dz;
    r4 = r2*r2;
    // compute the force and the energy
    potE = (0.25*k) * r4;
    fy = k * dy*dy*dy;
    fz = k * dz*dz*dz;
    p1->ay -= fy;
    p1->az -= fz;
    #ifdef TEMPERATURE_CONF
	    // configurational temperature
	    p1->Tfy -= fy;
	    p1->Tfz -= fz;
	    p1->Tdivf += 2 * k;
    #endif

    return potE;
}

/// Compute a 2D (cylindrical) quintic force on an atom, given its diplacement from its rest
/// axis. Divergence of the harmonic force is -3*kspring, for the configurational
/// temperature calculation.
///
/// @param		p1 pointer to the atom
/// @param		dx,dy,dz position difference between atom and the centre of the tube
/// @param		k the harmonic spring constant
///
/// @return		the value of the harmonic interaction energy

//================================================================================
// inline real SexticSpringTube (particleMD *p1, real dy, real dz, real k)
real SexticSpringTube (particleMD *p1, real dy, real dz, real k)
//================================================================================
{
    real r2,r6, fy, fz;
    real potE=0;

    // calculate the distance squared
    r2 = dy*dy + dz*dz;
    r6 = r2*r2*r2;
    // compute the force and the energy
    potE = (1.0/6.0) * k * r6;
    fy = k * dy*dy*dy*dy*dy;
    fz = k * dz*dz*dz*dz*dz;
    p1->ay -= fy;
    p1->az -= fz;
    #ifdef TEMPERATURE_CONF
	    // configurational temperature
	    p1->Tfy -= fy;
	    p1->Tfz -= fz;
	    p1->Tdivf += 2 * k;
    #endif

    return potE;
}

/// Compute a 2D (cylindrical) force on an atom, given its diplacement from its rest
/// axis. The tube is a Lennard Jones Tube that only applies to the MD particles.
///
/// @param		p1 pointer to the atom
/// @param		dx,dy,dz difference between atom and the centre of the tube
/// @param		k the harmonic spring constant
///
/// @return		the value of the harmonic interaction energy

//================================================================================
// inline real LennardJonesTube (particleMD *p1, real dy, real dz, real k, simptr sim)
real LennardJonesTube (particleMD *p1, real dy, real dz, real k, simptr sim)
//================================================================================
{
    if( dy==0.0 && dz==0.0 ) return 0.0;
    else {
        // Initialize
        particleMD p_wall;			// Imaginary wall particle
        real kT = sim->kT[sim->phase];
				real ljShift = sim->ljShift;
        real tubeRad = sqrt( k/kT );
        real p1Rad = sqrt(  dy*dy+dz*dz );
        real rCut2 = sim->rCut*sim->rCut;
        real dY=0.0, dZ=0.0;
        real potE = 0.0;

        // y and z components of the "wall particle"
        dY = tubeRad*dy/p1Rad;
        dZ = tubeRad*dz/p1Rad;
        p_wall = *p1;
        p_wall.ry = dY;
        p_wall.rz = dZ;
        p_wall.wy = dY;
        p_wall.wz = dZ;
				p_wall.ax = 0.0;
        p_wall.ay = 0.0;
				p_wall.az = 0.0;
				p_wall.vx = 0.0;
				p_wall.vy = 0.0;
				p_wall.vz = 0.0;

				potE = LennardJones (p1,&p_wall,0.0,(dY-dy),(dZ-dz),rCut2) + ljShift;
        return potE;
    }
}

/// Compute the FENE (finite-extensible non-linear elastic) interaction between
/// two atoms (see Warner). The position difference, the maximum extension and
/// the strenght of the FENE potential are passed as parameters. Divergence of
/// the force is also calculated for the configurational temperature.
///
/// @param		p1 pointer to the first atom
/// @param		p2 pointer to the second atom
/// @param		dx,dy,dz position difference between the two atoms
/// @param		k FENE spring constant
/// @param		r0 maximum FENE spring extension
/// @return		the value of the FENE interaction energy
/// @warning	If the finite integration increase dr over r0, the bond is broken
///				(and the program will likely crash).

//================================================================================
// inline real FENE (particleMD *p1, particleMD *p2, real dx, real dy, real dz,
// 				  real k, real r0)
real FENE (particleMD *p1, particleMD *p2, real dx, real dy, real dz,
				  real k, real r0)
//================================================================================
{
	real r2, r02, fene, fMag, fx, fy, fz;
	real potE=0;
	#ifdef TEMPERATURE_CONF
		real Tdivf;
	#endif

	// calculate the distance squared
	r2 = dx*dx + dy*dy + dz*dz;

	// compute the force and the energy
	r02   = r0*r0;
	fene  = 1-r2/r02;
	potE  = -0.5*k*r02*log(fene);
	fMag  = -k/fene;
	#ifdef TEMPERATURE_CONF
		Tdivf = k*(2+fene)/(fene*fene);
	#endif
	fx	= fMag * dx;
	fy	= fMag * dy;
	fz	= fMag * dz;
	p1->ax -= fx;
	p1->ay -= fy;
	p1->az -= fz;
	p2->ax += fx;
	p2->ay += fy;
	p2->az += fz;

	#ifdef TEMPERATURE_CONF
		// configurational temperature
		p1->Tfx -= fx;
		p1->Tfy -= fy;
		p1->Tfz -= fz;
		p2->Tfx += fx;
		p2->Tfy += fy;
		p2->Tfz += fz;
		p1->Tdivf += Tdivf;
		p2->Tdivf += Tdivf;
	#endif

	return potE;
}

/// Compute the bendHarmonic (simple harmonic angle) interaction between
/// two bond vectors (see Allen). The two bond vectors and k
/// the strength of the BEND potential are passed as parameters. Divergence of
/// the force is also calculated for the configurational temperature.
///
/// @param		p1 pointer to the first atom
/// @param		p2 pointer to the second atom
/// @param		p3 pointer to the third atom
/// @param		dx12,dy12,dz12 position difference between atoms 1 and 2
/// @param		dx23,dy23,dz23 position difference between atoms 2 and 3
/// @param		k bend spring constant
/// @param		equi equilibrium angle
/// @return		the value of the harmonic angle interaction energy

//================================================================================
real bendHarmonic (particleMD *p1, particleMD *p2, particleMD *p3,
				  real dx12, real dy12, real dz12, real dx23, real dy23, real dz23, real k, real equi)
//================================================================================
{
	real r12r23, ir12, ir23;
	real theta, c, c0, isinc, dcx1, dcy1, dcz1, dcx3, dcy3, dcz3;
	real fx1=0, fy1=0, fz1=0, fx3=0, fy3=0, fz3=0;
	real potE=0;
	int bendStyle=0; //0==angularHarmonic 1==cosineHarmonic 2==cosineExpansion

	// calculate the distances squared and make sure these aren't too close to zero to explode
	r12r23 = dx12*dx23 + dy12*dy23 + dz12*dz23;
	ir12 = sqrt(dx12*dx12 + dy12*dy12 + dz12*dz12);
	ir23 = sqrt(dx23*dx23 + dy23*dy23 + dz23*dz23);
	if( ir12>TOL && ir23>TOL ) {
		// Invert
		ir12 = 1.0/ir12;
		ir23 = 1.0/ir23;
		// Trig
		c=r12r23*ir12*ir23;
		if( c>0.99999 && c<1.00001 ) theta=0.0;
		else if( c<-0.999999 && c>-1.00001 ) theta=0.0;
		else theta=acos(c);
		//Shift by equibrium angle
		c0=cos(equi);

		// compute the derivative of the cosine for the first and third particle *k
		// first
		dcx1=k*ir12*( dx23*ir23 - c*dx12*ir12 );
		dcy1=k*ir12*( dy23*ir23 - c*dy12*ir12 );
		dcz1=k*ir12*( dz23*ir23 - c*dz12*ir12 );

		// last
		dcx3=k*ir23*( c*dx23*ir23 - dx12*ir12 );
		dcy3=k*ir23*( c*dy23*ir23 - dy12*ir12 );
		dcz3=k*ir23*( c*dz23*ir23 - dz12*ir12 );

		if (bendStyle==0) { //angular harmonic
			// computes 1/sinc(theta)
			if(fabs(theta)<0.00001) isinc=1.0;
			else if(fabs(fabs(theta)-M_PI)<0.00001) isinc=0.0;
			else isinc=(theta-equi)/sin(theta);

			//computes the forces
			fx1 = isinc*dcx1;
			fy1 = isinc*dcy1;
			fz1 = isinc*dcz1;
			fx3 = isinc*dcx3;
			fy3 = isinc*dcy3;
			fz3 = isinc*dcz3;

			// compute the energy
			potE  = 0.5*k*theta*theta;
		}
		else if (bendStyle==1) { //cosine harmonic
			//computes the forces
			fx1 = -2.0*(c-c0)*dcx1;
			fy1 = -2.0*(c-c0)*dcy1;
			fz1 = -2.0*(c-c0)*dcz1;
			fx3 = -2.0*(c-c0)*dcx3;
			fy3 = -2.0*(c-c0)*dcy3;
			fz3 = -2.0*(c-c0)*dcz3;

			// compute the energy
			potE  = k*(c-c0)*(c-c0);
		}
		else if (bendStyle==2) { // cosine expansion
			printf("Error: BendStyle=2 (cosine expansion) needs to be fixed.\n");
			exit(1);
			//computes the forces
			fx1 = dcx1;
			fy1 = dcy1;
			fz1 = dcz1;
			fx3 = dcx3;
			fy3 = dcy3;
			fz3 = dcz3;

			// compute the energy
			potE  = k*(1.0-c);
		}
		// Apply the force p1, p2 and p3
		p1->ax += fx1;
		p1->ay += fy1;
		p1->az += fz1;
		p3->ax += fx3;
		p3->ay += fy3;
		p3->az += fz3;
		p2->ax -= (fx1+fx3);
		p2->ay -= (fy1+fy3);
		p2->az -= (fz1+fz3);

	}

	return potE;
}


/// Compute the dihedralHarmonic (simple harmonic angle) interaction between
/// four bond vectors . 
/// This is written the same as origami dihedralforce in origami code written by Tyler Shendruk and the Origami Senior Honours Group.
///The four bond vectors and k
/// the strength of the DIHEDRAL potential are passed as parameters. Divergence of
/// the force is also calculated for the configurational temperature.
///
/// @param		p1 pointer to the first atom
/// @param		p2 pointer to the second atom
/// @param		p3 pointer to the third atom
/// @param		p4 pointer to the forth atom
/// @param		dx12,dy12,dz12 position difference between atoms 1 and 2
/// @param		dx23,dy23,dz23 position difference between atoms 2 and 3
/// @param		dx34,dy34,dz34 position difference between atoms 2 and 3
/// @param		k dihedral constant
/// @param		equi equilibrium angle
/// @return		the value of the harmonic angle interaction energy

//================================================================================
real dihedralHarmonic (particleMD *p1, particleMD *p2, particleMD *p3, particleMD *p4, real dx12, real dy12, real dz12, real dx23, real dy23, real dz23, real dx34, real dy34, real dz34, real k, real equi)
//================================================================================
{
	real ang, c, k_eff;
	real fx1=0, fy1=0, fz1=0, fx4=0, fy4=0, fz4=0;
	real potE=0;
	real sign = 0.0;						  // to establish the quadrant of phi 
	int mode=0;                              //0==angularHarmonic 1==cosineHarmonic 2==cosineExpansion

	double c22,c23,c24,c33,c34,c44;
	double p,q,Q,t1,t2,t3,t4,t5,t6,t7;
	double beta1,beta2,beta3,beta4;

	// Calculate Rapaport's coefficients
	c22 = (dx12*dx12 + dy12*dy12 + dz12*dz12);
	c23 = (dx12*dx23 + dy12*dy23 + dz12*dz23);
	c24 = (dx12*dx34 + dy12*dy34 + dz12*dz34);
	c33 = (dx23*dx23 + dy23*dy23 + dz23*dz23);
	c34 = (dx23*dx34 + dy23*dy34 + dz23*dz34);
	c44 = (dx34*dx34 + dy34*dy34 + dz34*dz34);
	// p = -(u.v) where u = r12 x r23 and v = r23 x r34, normal vectors to planes including r12,r23 and r23,r34 
	// will be used to calculate angle between planes
	p=c24*c33-c23*c34;
	/// q = |u|^2 * |v|^2 , |u|^2 =  (c22*c33-c23*c23) and |v|^2 = (c33*c44-c34*c34)                                      
	q=(c22*c33-c23*c23)*(c33*c44-c34*c34);							

	// to establish the quadrant of phi
	sign = dx34*(dy12*dz23 - dy23*dz12) + dy23*(dx23*dz12 - dx12*dz23) + dz34*(dx12*dy23 - dx23*dy12);

	t1=p;
	t2=c22*c34-c23*c24;
	t3=c23*c23-c22*c33;
	// (r23 X r34)^2 = v^2
	t4=c33*c44-c34*c34;
	t5=c24*c34-c23*c44;
	t6=-t1;
	// (r12 X r23)^2 = u^2
	t7=c22*c33-c23*c23;
	beta2=c34/c33;
	beta3=c23/c33;
	beta1=-1.0-beta3;
	beta4=-1.0-beta2;

	if( sqrt(t7)>TOL && sqrt(t4)>TOL ) {
		// Trig
		Q = 1.0/sqrt(q);
		c = -1.0*p*Q;
		if( c>0.99999 && c<1.00001 ) ang=0.0;
		else if( c<-0.999999 && c>-1.00001 ) ang=M_PI;
		else ang = acos(c);

		if (sign<0) ang=0.0-ang;

		k_eff = k;

		// Angular-harmonic version of the force
		if(mode==0) {
			if(fabs(ang)<0.00001) k_eff*=1.0;
			else if(fabs(fabs(ang)-M_PI)<0.00001) k_eff*=0.0;
			else k_eff*=(ang-equi)/sin(ang);

			// compute the energy
			potE  = 0.5*k*(ang-equi)*(ang-equi);	
		}
		// Cosine-harmonic version of the force
		else if(mode==1) {
			k_eff *= -2.0*(c-cos(equi));

			// compute the energy
			potE  = k*(c-cos(equi))*(c-cos(equi));
		}
		// Cosine-expansion version of the force
		else if(mode==2) {
			printf("Error: mode=2 (cosine expansion) needs to be fixed.\n");
			exit(1);
			k_eff*=1.0;	

			// compute the energy
			potE  = k*(1.0-cos(ang-equi));              
		}
		// Calculate first and last force
		// first
		fx1 = -1.0*k_eff*Q*c33*( t1*dx12 + t2*dx23 + t3*dx34 )/t7;
		fy1 = -1.0*k_eff*Q*c33*( t1*dy12 + t2*dy23 + t3*dy34 )/t7;
		fz1 = -1.0*k_eff*Q*c33*( t1*dz12 + t2*dz23 + t3*dz34 )/t7;

		// last
		fx4 = -1.0*k_eff*Q*c33*( t4*dx12 + t5*dx23 + t6*dx34 )/t4;
		fy4 = -1.0*k_eff*Q*c33*( t4*dy12 + t5*dy23 + t6*dy34 )/t4;
		fz4 = -1.0*k_eff*Q*c33*( t4*dz12 + t5*dz23 + t6*dz34 )/t4;

		// Apply the force p1, p2, p3 and p4
		p1->ax += fx1;
		p1->ay += fy1;
		p1->az += fz1;
		p4->ax += fx4;
		p4->ay += fy4;
		p4->az += fz4;
		p2->ax += beta1*fx1+beta2*fx4;
		p2->ay += beta1*fy1+beta2*fy4;
		p2->az += beta1*fz1+beta2*fz4;
		p3->ax += beta3*fx1+beta4*fx4;
		p3->ay += beta3*fy1+beta4*fy4;
		p3->az += beta3*fz1+beta4*fz4;
	}
	return potE;
}

/// Compute the bendNematic (simple harmonic angle) interaction between
/// the bond vector of the polymer and the nematic background solvent.
/// This is a direct copy of bendHarmonic but modified to act between the monomers and nematic
/// The bond vector, the director and ks the strength of the BEND potential time scalar order
//   are passed as parameters.
///
/// @param		p1 pointer to the first atom
/// @param		p2 pointer to the second atom
/// @param		dx12,dy12,dz12 position difference between atoms 1 and 2
/// @param		dx23,dy23,dz23 ---> the director [nx,ny,nz] of the background nematic
/// @param		k bend spring constant
/// @param    S local scalar order parameter
/// @return		the value of the harmonic angle interaction energy

//================================================================================
real bendNematic (particleMD *p1, particleMD *p2, real dx12, real dy12, real dz12,
	        real dx23, real dy23, real dz23, real k, real S, real dt,
					struct spec *SP, struct cell *CL)
//================================================================================
{
	real r12r23, ir12, ir23, r12;
	real theta, c, isinc, dcx1, dcy1, dcz1;
	real ks;
	real fx1=0, fy1=0, fz1=0;				//Force
	real torque=0, f1=0, s=1;
	double m[] = {0, 0, 0}; 	//Effective magnetic field applied to local mesogens
	real potE=0;
	int bendStyle=0; //0==angularHarmonic 1==cosineHarmonic 2==cosineExpansion

	// calculate the distances squared and make sure these aren't too close to zero to explode
	ks=k*S;
	r12r23 = dx12*dx23 + dy12*dy23 + dz12*dz23;
	r12 = sqrt(dx12*dx12 + dy12*dy12 + dz12*dz12);
	ir23 = 1.0;		// The director should be unit vector

	// Because the interaction should be nematic in nature, need to flip director if anti-parallel
	if( r12r23<0 ) {
		dx23*=-1.0;
		dy23*=-1.0;
		dz23*=-1.0;
		r12r23*=-1.0;
	}
	// Calculate the forces
	if( r12>TOL ) {
		// Invert
		ir12 = 1.0/r12;
		// Effective magnetic field direction
		m[0] = dx12*ir12;
		m[1] = dy12*ir12;
		m[2] = dz12*ir12;
		// Trig
		c=r12r23*ir12*ir23;
		if( c>0.99999 && c<1.00001 ) theta=0.0;
		else if( c<-0.999999 && c>-1.00001 ) theta=0.0;
		else theta=acos(c);

		// compute the derivative of the cosine for the first *k
		// first
		dcx1=k*ir12*( dx23*ir23 - c*dx12*ir12 );
		dcy1=k*ir12*( dy23*ir23 - c*dy12*ir12 );
		dcz1=k*ir12*( dz23*ir23 - c*dz12*ir12 );

		if (bendStyle==0) { //angular harmonic
			// computes 1/sinc(theta)
			if(fabs(theta)<0.00001) isinc=1.0;
			else if(fabs(fabs(theta)-M_PI)<0.00001) isinc=0.0;
			else isinc=theta/sin(theta);

			//computes the forces
			fx1 = isinc*dcx1;
			fy1 = isinc*dcy1;
			fz1 = isinc*dcz1;

			// compute the energy
			potE  = 0.5*ks*theta*theta;
		}
		else if (bendStyle==1) { //cosine harmonic
			//computes the forces
			fx1 = -2.0*(c - 1.0)*dcx1;
			fy1 = -2.0*(c - 1.0)*dcy1;
			fz1 = -2.0*(c - 1.0)*dcz1;

			// compute the energy
			potE  = ks*(c-1.0)*(c-1.0);
		}
		else if (bendStyle==2) { // cosine expansion
			//computes the forces
			fx1 = dcx1;
			fy1 = dcy1;
			fz1 = dcz1;

			// compute the energy
			potE  = ks*(1.0-c);
		}
		// Apply the force p1 and p2
		p1->ax += fx1;
		p1->ay += fy1;
		p1->az += fz1;
		p2->ax -= fx1;
		p2->ay -= fy1;
		p2->az -= fz1;

		//compute the torque magnitude on the segment between the two monomers t = Fr sin(alpha)
		f1 = sqrt(fx1*fx1 + fy1*fy1 + fz1*fz1); //force magnitude
		if (f1 != 0.0){
			s = sqrt(1-pow( (dx12*fx1 + dy12*fy1 + dz12*fz1) / (f1 * r12), 2)); //sin(alpha)
			torque = f1*r12*0.5*s;

			m[0] = torque*m[0];
			m[1] = torque*m[1];
			m[2] = torque*m[2];

			// Apply the force to the nematic fluid as a local field
			magTorque_CL( CL,SP,(double)dt,m );
		}
	}

	return potE;
}





/// Compute the bendNematic (simple harmonic angle) interaction between
/// the bond vector of the polymer and the nematic background solvent (see Allen).
/// The two bond vectors and k the strength of the BEND potential are passed as
//  parameters. Divergence of the force is also calculated for the configurational temperature.
///
/// @param		p1 pointer to the first atom
/// @param		p2 pointer to the second atom
/// @param		dx12,dy12,dz12 position difference between atoms 1 and 2 from p2 to p1 i.e. dx12 = p1->wx - p2->wx
/// @param		nx,ny,nz the director of the background nematic solvent
/// @param		ks bend spring constant (times local scalar order parameter)
/// @return		the value of the harmonic angle interaction energy

// //================================================================================
// real bendNematic (particleMD *p1, particleMD *p2, real dx12, real dy12, real dz12,
// 				  real nx, real ny, real nz, real ks, real dt,struct spec *SP,struct cell *CL)
// //================================================================================
// {
// 	real r12n, ir12, irn;
// 	real s, c, signc, theta;
// 	double m[] = {0, 0, 0};
// 	real that[] = {0, 0, 0},mag=0;
// 	real f[] = {0.0, 0.0, 0.0};
// 	real potE=0;
//
// 	// calculate the distances squared
// 	r12n = dx12*nx + dy12*ny + dz12*nz;
// 	ir12 = sqrt(dx12*dx12 + dy12*dy12 + dz12*dz12);
// 	irn = sqrt(nx*nx + ny*ny + nz*nz);
// 	if( ir12>TOL && irn>TOL ) {
// 		ir12 = 1.0/ir12;
// 		irn = 1.0/irn;
// 		m[0] = dx12*ir12;
// 		m[1] = dy12*ir12;
// 		m[2] = dz12*ir12;
//
// 		// compute the angles
// 		c=r12n*ir12*irn;
// 		signc=(c > 0) ? 1 : ((c < 0) ? -1 : 0);
// 		if( c>0.99999 && c<1.000001 ) theta=0.0;
// 		else if( c<-0.999999 && c>-1.0000001 ) theta=0.0;
// 		else theta=acos(c);	//theta is ALWAYS positive
// 		s=sin(theta);
// 		// compute the theta-hat direction (sign should make up for signless theta)
// 		that[0] = nx-c*m[0];
// 		that[1] = ny-c*m[1];
// 		that[2] = nz-c*m[2];
// 		// Normalize and set sign
// 		mag=sqrt(that[0]*that[0]+that[1]*that[1]+that[2]*that[2]);
// 		if( mag>TOL ) {
// 			if(signc<0) {
// 				mag*=signc;
// 				theta=M_PI-theta;
// 				c=cos(theta);
// 				s=sin(theta);
// 			}
// 			that[0] /= mag;
// 			that[1] /= mag;
// 			that[2] /= mag;
//
// 			// compute the forces (should have a factor of 2 since ir12 is twice the distance BUT spread force across two particles)
// 			//Harmonic cosine potential
// 			f[0] = ks*s*that[0]*ir12;
// 			f[1] = ks*s*that[1]*ir12;
// 			f[2] = ks*s*that[2]*ir12;
// 			// Apply the forces to MD monomers
// 			p1->ax += f[0];
// 			p1->ay += f[1];
// 			p1->az += f[2];
// 			p2->ax -= f[0];
// 			p2->ay -= f[1];
// 			p2->az -= f[2];
// 			// Apply the force to the nematic fluid as a local field
// 			magTorque_CL( CL,SP,(double)dt,m );
// 			// compute the energy
// 			potE  = ks*(1-c);
// 		}
// 	}
//
// 	return potE;
// }

/// Apply the periodic boundary condition to all atoms. We simply call the
/// ApplyPBC function for each atom to make them all fit inside the simulation box.
///
/// @param		sim		a pointer to a simulation structure
/// @return 	void
/// @warning	This will only work if the atoms are at most one box length away
///				from a box boundary.

//================================================================================
void PeriodicBoundaries (simptr sim)
//================================================================================
{
	int	  	i, nAtom;
	particleMD	*atom;

	// local sim variables
	atom	= sim->atom.items;
	nAtom	= sim->atom.n;

	// apply periodic boundary condition
	for (i=0; i<nAtom; i++) {
		ApplyPBC (sim, &(atom[i].rx), &(atom[i].ry), &(atom[i].rz));
	}
}


/// Apply the periodic boundary conditions. Adjusts the values of positions
/// or distances pointed to by x, y and z according to the pbc condition
/// and the half box size set in the overall sim structure.
///
/// @param		sim	a pointer to a simulation structure
/// @param		px,py,pz pointers to the concerned position
/// @return 	void
/// @warning	Parameters are POINTERS to the positions, not the positions themselves.

//================================================================================
void ApplyPBC (simptr sim, real *px, real *py, real *pz)
//================================================================================
{
	// pbc along the x axis
	if (sim->pbcond & PBC_COND_x) {
		if 		(*px >= sim->box[x_])  *px -= sim->box[x_];
		else if	(*px <  0. ) *px += sim->box[x_];
	}

	// pbc along the y axis
	if (sim->pbcond & PBC_COND_y) {
		if	  	(*py >=  sim->box[y_]) *py -= sim->box[y_];
		else if	(*py < 0. ) *py += sim->box[y_];
	}

	// pbc along the z axis
	if (sim->pbcond & PBC_COND_z) {
		if	  	(*pz >=  sim->box[z_]) *pz -= sim->box[z_];
		else if	(*pz <  0. ) *pz += sim->box[z_];
	}
}
/// Apply the periodic boundary conditions to a difference. Adjusts the values of positions
/// or distances pointed to by x, y and z according to the pbc condition
/// and the half box size set in the overall sim structure.
///
/// @param		sim	a pointer to a simulation structure
/// @param		px,py,pz pointers to the concerned position
/// @return 	void
/// @warning	Parameters are POINTERS to the positions, not the positions themselves.

//================================================================================
void ApplyPBC_dr (simptr sim, real *dx, real *dy, real *dz)
//================================================================================
{
	real *box,*boxHalf;

	box = sim->box;
	boxHalf = sim->boxHalf;
	// pbc along the x axis
	if (sim->pbcond & PBC_COND_x) {
		if (*dx >= boxHalf[x_]) *dx -= box[x_];
		else if (*dx < -boxHalf[x_]) *dx += box[x_];
	}

	if (sim->pbcond & PBC_COND_y) {
		if (*dy >= boxHalf[y_]) *dy -= box[y_];
		else if	(*dy < -boxHalf[y_]) *dy += box[y_];
	}

	if (sim->pbcond & PBC_COND_z) {
		if (*dz >= boxHalf[z_]) *dz -= box[z_];
		else if	(*dz < -boxHalf[z_]) *dz += box[z_];
	}
}

/// Verifies wether we need to update the neighbor list. We update it when
/// the maximum displacement becomes larger than half the thickness of the
/// neighbor shell (typically a third of the neighborhood radius, rNebrShell
/// in the input file). The maximum displacement is saved in sim->drTotMax from
/// within the "VelocityVerletStep" function.
///
/// @param	sim		a pointer to a simulation structure
///7 @return void

//================================================================================
void CheckNebrList (simptr sim)
//================================================================================
{
	// check if we need to refresh the neighbor lists
	if (sim->drTotMax > 0.5*sim->rNebrShell) {
//		RefreshNebrList (sim);  //Note to tyler neighbor list disabled
		sim->drTotMax = 0;
	}
}


/// Updates the neighbor list by finding atom pairs in neighboring cells
/// that lie within the neighbor shell distance. The neighbor cells are
/// indexed with statically defined offsets. Neighbors are stored as POINTERS
/// to particle structures to avoid double indirection in the force
/// calculation (Steve Guillouzic, 2001). This is also what is done in more
/// advanced software (e.g. Espresso). Separate lists are kept for STD
/// (non-pbc) pairs and PBC pairs.
///
/// @param	sim		a pointer to a simulation structure
/// @return void

//================================================================================
void RefreshNebrList (simptr sim)
//================================================================================
{
	int			i, j, k, ***cell, *cellList;
	particleMD	*atom, *p1, *p2;
	list2STD	*nebrSTD;
	list2PBC    *nebrPBC;
	int			ncx, ncy, ncz, cx, cy, cz, nx, ny, nz;
	int			pbcIndex, pbcond;
	real		dx, dy, dz;
	real		r2, rNebr2;
	real		*pbc;
	static const int xoffset[] = {0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, -1, -1, -1};
	static const int yoffset[] = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	static const int zoffset[] = {0, 1, 1, 0, -1, -1, 0, 1, 1, 0, -1, -1, 0, 1};

	// local sim variables
	atom	 = sim->atom.items;
	nebrSTD  = &(sim->nebrSTD);
	nebrPBC  = &(sim->nebrPBC);
	cell	 = sim->cell;
	cellList = sim->cellList;
	rNebr2   = sim->rNebr2;
	ncx		 = sim->nCellAxis[x_];
	ncy		 = sim->nCellAxis[y_];
	ncz		 = sim->nCellAxis[z_];

	// make sure all atoms are inside the box
	PeriodicBoundaries (sim);

	// sort particles in cells
	RefreshCellList (sim);

	// loop over all cells
	nebrSTD->n=0;
	nebrPBC->n=0;

	for (cx=0; cx<ncx; cx++) {
		for (cy=0; cy<ncy; cy++) {
			for (cz=0; cz<ncz; cz++) {

				// loop over half of neighboring cells (14 in 3D)
				for (k=0; k<14; k++) {
					nx = cx + xoffset[k];
					ny = cy + yoffset[k];
					nz = cz + zoffset[k];

					// determine appropriate PBC
					pbcIndex = PBC_INDEX_NONE;
					pbcond   = 0;
					if (nx<0) {
						nx = ncx-1;
						pbcIndex -= PBC_INDEX_x;
						pbcond   |= PBC_COND_x;
					} else if (nx == ncx) {
						nx = 0;
						pbcIndex += PBC_INDEX_x;
						pbcond   |= PBC_COND_x;
					}
					if (ny<0) {
						ny = ncy-1;
						pbcIndex -= PBC_INDEX_y;
						pbcond   |= PBC_COND_y;
					} else if (ny == ncy) {
						ny = 0;
						pbcIndex += PBC_INDEX_y;
						pbcond   |= PBC_COND_y;
					}
					if (nz<0) {
						nz = ncz-1;
						pbcIndex -= PBC_INDEX_z;
						pbcond   |= PBC_COND_z;
					} else if (nz == ncz) {
						nz = 0;
						pbcIndex += PBC_INDEX_z;
						pbcond   |= PBC_COND_z;
					}

					// consider only aplicable periodic boundary condition
					if (pbcond && (pbcond & ~sim->pbcond)) continue;

					// loop over atoms in the current cell
					i = cell[cx][cy][cz];
					while (i >= 0) {

						// loop over atoms in the neighboring cell
						j = cell[nx][ny][nz];
						while (j >= 0) {

							// avoid double counting inside same cell
							if (k==0 && i <= j) {
								j = cellList[j];
								continue;
							}

							// calculate distance between atoms
							pbc = sim->pbc[pbcIndex];
							p1  = atom+i;
							p2  = atom+j;
							dx  = p2->rx - p1->rx + pbc[x_];
							dy  = p2->ry - p1->ry + pbc[y_];
							dz  = p2->rz - p1->rz + pbc[z_];
							r2  = dx*dx + dy*dy + dz*dz;

							// if we find a neighbor, register it
							if (r2 < rNebr2) {

								// STD or PBC pair
								if (pbcIndex == PBC_INDEX_NONE)  AddItem2STD (nebrSTD, p1, p2);
								else 							 AddItem2PBC (nebrPBC, p1, p2, pbc);
							}
							j = cellList[j];
						}
						i = cellList[i];
					}
				}
			}
		}
	}
}


/// Sort particles in cells according to their position. All particles MUST
/// lie within the simulation box for this to work. The list of atoms from
/// each cell is implemented using a single array in a linked list fashion
/// (See Knuth, Rapaport, or Allen & Tilsdeley).
///
/// @param	sim		a pointer to a simulation structure
/// @return void

//================================================================================
void RefreshCellList (simptr sim)
//================================================================================
{
	int			i, nAtom, ***cell, *cellList;
	int			cx, cy, cz, cxmax, cymax, czmax;
	particleMD	*atom;
	const real	*boxHalf, *cellInvWidth;

	// local sim variables
	atom		 = sim->atom.items;
	nAtom		 = sim->atom.n;
	cell		 = sim->cell;
	cellList	 = sim->cellList;
	cellInvWidth = sim->cellInvWidth;
	boxHalf	  	 = sim->boxHalf;

	// set maximum cell indices + 1
	cxmax = sim->nCellAxis[x_];
	cymax = sim->nCellAxis[y_];
	czmax = sim->nCellAxis[z_];

	// reset cell list
	for (i=0; i<sim->nCell; i++) {
		cell[0][0][i] = -1;
	}

	// sort atoms in appropriate cells
	for (i=0; i<nAtom; i++) {
// 		cx = (int) ((atom[i].rx+boxHalf[x_]) * cellInvWidth[x_]);
// 		cy = (int) ((atom[i].ry+boxHalf[y_]) * cellInvWidth[y_]);
// 		cz = (int) ((atom[i].rz+boxHalf[z_]) * cellInvWidth[z_]);
		cx = (int) ((atom[i].rx) * cellInvWidth[x_]);
		cy = (int) ((atom[i].ry) * cellInvWidth[y_]);
		cz = (int) ((atom[i].rz) * cellInvWidth[z_]);

		// this SHOULD never happen but it DOES because of finite precision arithmetic
		// (especially in single precision, e.g. (int) 102.9999 = 103). We just put
		// the offending particle in the last cell.
		if (cx == cxmax) cx--;
		if (cy == cymax) cy--;
		if (cz == czmax) cz--;

		// catch out of bounds atom
		if (cx<0 || cx>=cxmax || cy<0 || cy>=cymax || cz<0 || cz>=czmax ) {
			LOG ("\nAtom %d is outside simulation box\n", i);
			LOG ("box      (%f %f %f)\n",	sim->box[x_], sim->box[y_], sim->box[z_]);
			LOG ("boxHalf  (%f %f %f)\n",	boxHalf[x_], boxHalf[y_], boxHalf[z_]);
			LOG ("InvWidth (%f %f %f)\n\n", cellInvWidth[x_], cellInvWidth[y_], cellInvWidth[z_]);
			LOG ("pos      (%f %f %f)\n",	atom[i].rx, atom[i].ry, atom[i].rz);
			LOG ("cells    (%d %d %d)\n", 	cx, cy, cz);
			LOG ("cmax     (%d %d %d)\n\n", cxmax, cymax, czmax);
			LOG ("type  %d\n", 				atom[i].type);
			LOG ("group %#010X\n",			atom[i].group);
			LOG ("step  %d\n",				sim->step[count_]);
			error (ECELLSORT);
		}

		// add atom to the cell
		cellList[i] = cell[cx][cy][cz];
		cell[cx][cy][cz] = i;
	}
}
