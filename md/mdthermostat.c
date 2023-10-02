
//================================================================================
//
// name:   mdthermostat.c
// author: ftessier
// date:   2005-05-03 @ 11:04:36
//
// Thermostat functions for Molecular Dynamics
//
//================================================================================


#include <stdlib.h>
#include <math.h>
#include "mdtypes.h"
#include "mderror.h"
#include "mdutil.h"
#include "mdthermostat.h"


/// Attempt to keep the temperature constant by simply rescaling the particle
/// velocities. The particle velocities in the group specified by sim->groupTherm
/// are adjusted so that their average matches the kinetic E computed in the
/// function InitVelocities [inspired from Rapaport, p. 65]. Center of mass motion
///	is not handled, so this adjustment should NOT be applied to flowing particles, for example. May
/// be good to replace this with a more robust scheme like the Berendsen
/// thermostat [BerendsenPGDH-1984]. We don't worry about the small error due
/// to the constraints on some degrees of freedom; each particle should have
/// 3/2*kT kinetic energy on average.
///
/// @param		simvoid	a void cast of a pointer to a simulation structure
/// @return		void
/// @warning    Now that the DPD thermostat is implemented, we only use the
///				rescale thermostat for the wall particles, which are attached
///				to a fixed lattice.
/// @todo		Implement the idea of thermostat atom lists and a list of
///				thermostat at the simulation level.make

//================================================================================
void ThermostatRescale (void *simvoid)
//================================================================================
{
	int			i, nAtom, groupThermRescale;
	real		kT, kinETarget, kinETherm;
	real		vScale;
	particleMD	*atom;
	simptr		sim = (simptr) simvoid;

	// local sim variables
	atom	   			= sim->atom.items;
	nAtom	  			= sim->atom.n;
	groupThermRescale 	= sim->groupThermRescale[sim->phase];
	kinETherm  			= sim->kinETherm;
	kT		 			= sim->kT[sim->phase];

	// compute rescaling factor from the target kinetic energy
	kinETarget = 1.5 * kT;
    // sanity check on if kinETherm is non-zero
    // in reality, this shouldn't pose any problems as it will only be zero if a bitwise group check fails
    // this method also does a bitwise group check later, so this should not cause any issues.
    if (kinETherm > TOL) {
        vScale	   = sqrt (kinETarget / kinETherm);
    } else {
        vScale = 1;
    }

	// rescale velocities
	for (i=0; i<nAtom; i++) {
		if (atom[i].group & groupThermRescale) {
			atom[i].vx *= vScale;
			atom[i].vy *= vScale;
			atom[i].vz *= vScale;
		}
	}
}


/// Compute the forces needed to implement a DPD thermostat for a particle pair.
///
/// @param		p1,p2 pointers to two particles
/// @param		dx,dy,dz difference in position between the particles
/// @param		rCut2 square of the DPD interaction cutoff
/// @param		eta friction parameter
/// @param		sigma noise parameter
/// @param		dtSqrti inverse of the square root of the timestep dt
/// @param		groupThermDPD the group on which DPD should be applied
/// @return		void
/// @warning	Only consider fluid-fluid pairs or wall-wall pairs for thermostatting
///				when there is a flow, following a discussion with Kartunnen (2005-03-01).
///				If different types are considered, don't forget to add the reduced mass
///				in the value of sigma (Peters, Europhys. Lett. 66, 311)
/// @warning	It is not clear if the DPD interaction should be included in
///				the configurational temperature calculation (probably not).
/// @note		The SQRT_3 is there to provide a unit variance uniform distribution.
/// @note		The value of sigma should account for timestep artifacts via a timestep-dependent
/// 			correction factor (Peters, Europhys. Lett. 66, 311).

//================================================================================
void ThermostatDPD (particleMD *p1, particleMD *p2, real dx, real dy, real dz,
			   		real rCut2, real eta, real sigma, real dtSqrti, int groupThermDPD)
//================================================================================
{
	real 	r2, ri, fMag, fx, fy, fz;
	real 	dvx, dvy, dvz, dv2;

	// select only DPD pairs
	if (p1->group!=groupThermDPD || p2->group!=groupThermDPD) return;

	// calculate the distance squared
	r2 = dx*dx + dy*dy + dz*dz;

	// add the DPD thermostatting forces if we are inside the cutoff range
	if (r2 < rCut2) {
		ri = 1/sqrt(r2);

		// difference in velocities
		dvx = p2->vx - p1->vx;
		dvy = p2->vy - p1->vy;
		dvz = p2->vz - p1->vz;
		dv2 = dx*dvx + dy*dvy + dz*dvz;

		// dissipation force
		fMag = -eta * dv2/r2;
		fx = fMag * dx;
		fy = fMag * dy;
		fz = fMag * dz;

		// fluctuation force
		fMag = sigma * ri * SQRT_3 * (2.0*RandomReal()-1.0) * dtSqrti;
		fx += fMag * dx;
		fy += fMag * dy;
		fz += fMag * dz;

		// update particle accelerations
		p1->ax -= fx;
		p1->ay -= fy;
		p1->az -= fz;
		p2->ax += fx;
		p2->ay += fy;
		p2->az += fz;
/*
		#ifdef	TEMPERATURE_CONF
		// confgurational temperature
		p1->Tfx -= fx;
		p1->Tfy -= fy;
		p1->Tfz -= fz;
		p2->Tfx += fx;
		p2->Tfy += fy;
		p2->Tfz += fz;
		#endif
*/

	}
}

