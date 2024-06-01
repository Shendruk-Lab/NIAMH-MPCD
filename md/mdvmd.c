
//================================================================================
//
// name:   mdvmd.c
// author: tshen
// date:   2011-09-22 @ 13:59:36
//
// VMD functions to print out visualization files
//
//================================================================================


#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mdtypes.h"
#include "mderror.h"
#include "mdmemory.h"
#include "mdfiles.h"
#include "mdutil.h"
#include "mdsrd.h"
#include "mdvmd.h"


//================================================================================
void SetupVMD (simptr sim)
//================================================================================
{
	// Initializes all scenes (used to output visualization data).

	sceneptr s;

	// report
	LOG ("Setting up scenes\n");

	// include additioinal user histogram definitions
	#include "mdvmddef.h"
	// setup scene stream pointers and function pointers
	VMDSetupFiles (sim);
	VMDSetupFuncs (sim);
	// report
	s = sim->scenes;
	while (s) {
		if (s->file) {
			if (s->groupInc) LOG ("  %-25s<- included %#010X\n", s->label, s->groupInc);
			if (s->groupExc) LOG ("  %-25s-> excluded %#010X\n", s->label, s->groupExc);
		}
		s = s->next;
	}
}


//================================================================================
void VMDSetupFuncs (simptr sim)
//================================================================================
{
	// Links scenes to their handling functions. For now, we point everything to
	// the a default scene function, but eventually, we could do something more
	// complex, like with the histograms.

	sceneptr s;

	// local sim variables
	s = sim->scenes;

	// loop over all scenes
	while (s) {

		// scenefunc pointer
		s->scenefunc = &VMDFunction;

		// next scene
		s = s->next;
	}
}


//================================================================================
void VMDSetupFiles (simptr sim)
//================================================================================
{
	// Links scenes to their respective output file by saving the appropriate
	// stream into the stream member of the histogram structure. Data files MUST
	// be open and bear the same label as the scene for this setup to
	// work, otherwise the program will report an error and abort.

	sceneptr s = sim->scenes;

	while (s) {
		if (s->file) s->stream = GetSimStream (sim->files, s->label);
		s = s->next;
	}
}


//================================================================================
void VMDPrint (simptr sim, sceneptr s)
//================================================================================
{
	// Loop through the supplied scene list and check if they are scheduled
	// for data printing. If so, call the histfunc with PRINT if func pointer is
	// not null.

	int	phase;

	// local simulation variables
	phase = sim->phase;

	// check every histogram
	while (s) {
		if (s->active && s->stream) {
			if (sim->step[count_] == 0 || s->stepPrint[count_] == s->stepPrint[phase]) {
				s->stepPrint[count_] = 0;
				if (s->scenefunc) s->scenefunc ((void *) sim, s, ACTION_PRINT);
			}
		}
		s = s->next;
	}
}


//================================================================================
void VMDPrintObjects (simptr sim, FILE *stream)
//================================================================================
{
	// Write object properties to file. These object properties are used to
	// determine the visual characteristics of the particles in the scenic viewer.

	int 		i, nAtom, polyN;
	real		*box;

	// local simulation variables
	nAtom 	= sim->atom.n;
	box	= sim->box;
	polyN	= sim->polyN[0];

	// write object properties
	fprintf (stream, "unitcell %e %e %e\n",box[0],box[1],box[2]);
	fprintf (stream, "atom 0:%d radius 1 name polymer type %d\n",nAtom-1,OBJ_MONOMER);
	for( i=0;i<nAtom;i+=polyN) {
		fprintf( stream,"bond %d::%d\n", i,i+polyN-1);
	}
}


//================================================================================
int VMDSpaceDomain (sceneptr s, real s1, real s2, real s3)
//================================================================================
{
	// Check if the coordinates of the point lie within the scene's space
	// domain. Note that space domain has its own dimension and coordinate
	// system. Returns 0 if the point is outside the space domain, 1 otherwise.

	// points outside space domain (reject)
	if ( (s->sd >= 1) && (s1 < s->s1min || s1 > s->s1max) ) return 0;
	if ( (s->sd >= 2) && (s2 < s->s2min || s2 > s->s2max) ) return 0;
	if ( (s->sd >= 3) && (s3 < s->s3min || s3 > s->s3max) ) return 0;

	// points inside space domain (accept)
	return 1;
}


//================================================================================
sceneptr VMDNew (sceneptr *s)
//================================================================================
{
	// Adds a scene to the beginning of the linked list of histograms
	// pointed to by *s.

	sceneptr new;

	// allocate memory for the new simscene structure
	new = (sceneptr) mycalloc (1, sizeof(simhist));

	// insert the new scene at beginning of the list
	if (s) {
		new->next = *s;
		*s = new;
	}

	// return a pointer to the newly created scene
	return new;
}


//================================================================================
sceneptr VMDDuplicate (sceneptr s, char *newlabel)
//================================================================================
{
	// Duplicates the scene pointed to by s and inserts it AFTER s, giving
	// it the new label.

	sceneptr new;

	if (!s) return 0;

	// allocate memory for the new simscene structure
	new = (sceneptr) mycalloc (1, sizeof(simhist));

	// copy scene information and attribute new label
	memcpy   (new, s, sizeof(simscene));
	snprintf (new->label, STRMAX, "%s", newlabel);

	// insert new scene at beginning of the list
	new->next = s->next;
	s->next   = new;

	// return a pointer to the newly created scene
	return new;
}


//================================================================================
sceneptr VMDGet (sceneptr s, char *label)
//================================================================================
{
	// Finds the scene coresponding to the label in the linked list of scenes
	// pointed to by s. Returns a pointer to the scene, or generates a fatal
	// error if the requested scene is not found.

	int found=0;

	while (!found && s) {
		if (strcmp(label, s->label) == 0) found = 1;
		else s = s->next;
	}

	if (!found) {
		error (ESCENELIST);
	}

	return s;
}

//================================================================================
void VMDFunction (void *simvoid, sceneptr s, int action)
//================================================================================
{
	// General scene operation function. The requested action is specified
	// by the action flag. For the moment, only the PRINT action is supported.

	simptr		sim;
	int 		i, nAtom;
	particleMD	*atom, *p1;
	point		r, rs;
	real		*sc1=0, *sc2=0, *sc3=0;

	// local simulation variables
	sim	= (simptr) simvoid;
	atom	= sim->atom.items;
	nAtom 	= sim->atom.n;

	// execute requested action
	switch (action) {

		case ACTION_PRINT:

			// write objects before the first frame
			if (s->frame == 0) VMDPrintObjects (sim, s->stream);

			// increase frame count
			s->frame++;

			// print frame header
			fprintf (s->stream, "\ntimestep indexed\n");

			// initialize local point variables
			memset (&r,  0, sizeof(point));
			memset (&rs, 0, sizeof(point));
			r.coord  = CARTESIAN;
			rs.coord = s->scoord;

			// select space domain if applicable
			if (s->sdomain >= 0) {
				sc1 = &(rs.c1);
				sc2 = &(rs.c2);
				sc3 = &(rs.c3);
				CoordinateOrder (s->sdomain, &sc1, &sc2, &sc3);
			}

			// loop over all atoms
			for (i=0; i<nAtom; i++) {

				// particle pointer
				p1 = atom+i;

				// apply include-exclude rules
				if (!(p1->group & s->groupInc)) continue;
				if (  p1->group & s->groupExc) 	continue;


				// get particle coordinates and adjust for periodic boundaries
				r.c1 = p1->rx;
				r.c2 = p1->ry;
				r.c3 = p1->rz;
				ApplyPBC (sim, &r.c1, &r.c2, &r.c3);


				// filter points outside space domain if applicable
				if (s->sdomain >= 0) {
					CoordinateTransform (&r, &rs,0);
					if (!VMDSpaceDomain (s, *sc1, *sc2, *sc3)) continue;
				}

				// print the atom
				fprintf (s->stream, "%d %.6G %.6G %.6G", i, r.c1, r.c2, r.c3);
// 				// print link if there is one
// 				if (p1->next) {
// 					p1 = p1->next;
// 					// get particle coordinates and adjust for periodic boundaries
// 					r.c1 = p1->rx;
// 					r.c2 = p1->ry;
// 					r.c3 = p1->rz;
// 					ApplyPBC (sim, &r.c1, &r.c2, &r.c3);
// 					fprintf (s->stream, " 1 %.6G %.6G %.6G", r.c1, r.c2, r.c3);
// 				}
				fprintf (s->stream, "\n");
			}

			break;

		default:
			break;
	}
}


//================================================================================
void VMDFreeList (sceneptr s)
//================================================================================
{
	// Recursively frees all the memory allocated for the linked list of scenes.

	if (!s) return;
	VMDFreeList (s->next);

	// free histogram structure
	free (s);
}
