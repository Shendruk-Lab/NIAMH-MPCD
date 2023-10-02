
//================================================================================
//
// name:   mdhistogram.c
// author: ftessier
// date:   2005-05-03 @ 11:04:36
//
//================================================================================
///
/// @file mdhistogram.c
/// Functions to handle data histograms. This file contains an fairly advanced
/// interface to the manipulation of data histograms of different dimensions and
/// geometries. The list of data histograms is defined in a separate include file
/// for convenience. The basic principle is that once a histogram has been declared,
/// its histfunc function can be called with different macros to reset, print,
/// collect, and normalize the data. Each histogram also contains its own step
/// counters for collecting and printing data. Of course, it would be preferable
/// to keep the raw data from the simulation, but the prohibitive size of the
/// resulting data set makes this impractical. This is why histogram are used to
/// perform part of the collecting live, during the simulation.
///
//================================================================================


#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mdtypes.h"
#include "mderror.h"
#include "mdmemory.h"
#include "mdfiles.h"
#include "mdmeasure.h"
#include "mdutil.h"
#include "mdsrd.h"
#include "mdhistogram.h"
#include "../mpcd/headers/bc.h"

/// Initializes all data histograms.
///
/// @param		sim a pointer to a simulation structure
/// @return		void
/// @warning	The fields of the histogram must be initialized before the function
/// 			AllocateHistogram is called to allocate memory for the data arrays.
/// @warning	For clarity, the histogram definitions are put in an include file.
/// @warning	The number of bins 'ni' and the bin size values 'bini' are not independent:
/// 			if ni=0 then the bini is used to compute ni, otherwise the bini is
///				computed based on the histogram domain and the number of bins.

//================================================================================
void SetupHistograms (simptr sim)
//================================================================================
{
	histptr 	h;
// 	histptr 	h1;
// 	char		str[STRLEN];

	// report
	LOG ("Setting up histograms\n");

	// include additioinal user histogram definitions
	#include "mdhistogramdef.h"

	// setup histogram stream pointers and function pointers
	HistogramSetupFiles (sim);
	HistogramSetupFuncs (sim);

	// report
	h = sim->histograms;
	while (h) {
		if (h->active && h->file) {
			if (h->groupInc)	LOG ("  %-25s<- included %#010X\n", h->label, h->groupInc);
			if (h->groupExc)	LOG ("  %-25s-> excluded %#010X\n", h->label, h->groupExc);
		}
		h = h->next;
	}
}


/// Links histogram to their handling functions. Saves the function pointer in
/// the histfunc and func members of the histogram structure. This function is
/// needed mostly because function pointers change when we reload from a checkpoint
/// file.
///
/// @param		sim a pointer to a simulation structure
///	@return		void

//================================================================================
void HistogramSetupFuncs (simptr sim)
//================================================================================
{
	histptr h;

	// local sim variables
	h = sim->histograms;

	// loop over all histograms
	while (h) {

		// histfunc pointers
		h->histfunc = NULL;
		if (!strcmp(h->histfuncName,""));
		else if (!strcmp(h->histfuncName, "HistogramFunction"))			h->histfunc = &HistogramFunction;
		else error (EHIST);

		// func pointers
		h->func = NULL;
		if (!strcmp(h->funcName,""));
		else if (!strcmp(h->funcName, "HistogramFunction_density_rho_y"))	h->func = &HistogramFunction_density_rho_y;
		else if (!strcmp(h->funcName, "HistogramFunction_bond_orientation_rho_y"))	h->func = &HistogramFunction_bond_orientation_rho_y;
		else if (!strcmp(h->funcName, "HistogramFunction_bond_orientation"))	h->func = &HistogramFunction_bond_orientation;
		else if (!strcmp(h->funcName, "HistogramFunction_Tkin"))		h->func = &HistogramFunction_Tkin;
		else if (!strcmp(h->funcName, "HistogramFunction_Tconf"))		h->func = &HistogramFunction_Tconf;
		else if (!strcmp(h->funcName, "HistogramFunction_q"))			h->func = &HistogramFunction_q;
		else if (!strcmp(h->funcName, "HistogramFunction_fv"))			h->func = &HistogramFunction_fv;
		else if (!strcmp(h->funcName, "HistogramFunction_vx"))			h->func = &HistogramFunction_vx;
		else if (!strcmp(h->funcName, "HistogramFunction_dr_yz"))		h->func = &HistogramFunction_dr_yz;
		else if (!strcmp(h->funcName, "HistogramFunction_bond"))		h->func = &HistogramFunction_bond;
		else if (!strcmp(h->funcName, "HistogramFunction_fenex_N"))		h->func = &HistogramFunction_fenex_N;
		else if (!strcmp(h->funcName, "HistogramFunction_qeff_N"))		h->func = &HistogramFunction_qeff_N;
		else if (!strcmp(h->funcName, "HistogramFunction_hwall_N"))		h->func = &HistogramFunction_hwall_N;
		else if (!strcmp(h->funcName, "HistogramFunction_ellipsoid")) 	h->func = &HistogramFunction_ellipsoid;
		else if (!strcmp(h->funcName, "HistogramFunction_persist")) 	h->func = &HistogramFunction_persist;
		else error (EHIST);

		// next histogram
		h = h->next;
	}
}


/// Links histogram to their respective output file. Saves the appropriate
/// stream into the stream member of the histogram structure.
///
/// @param		sim a pointer to a simulation structure
/// @return		void
/// @warning	Data files MUST be open and bear the same label as the histogram
/// 			for this setup to work, otherwise the program will report an error
/// 			and abort.

//================================================================================
void HistogramSetupFiles (simptr sim)
//================================================================================
{
	histptr h;

	// local sim variables
	h = sim->histograms;

	// loop over all histograms
	while (h) {
		if (h->file) h->stream = GetSimStream (sim->files, h->label);
		h = h->next;
	}
}


/// Loop through the active histograms and check if they are scheduled for data
/// collection. If so, call the histfunc with ACTION_COLLECT if it is not NULL.
///
/// @param		sim a pointer to a simulation structure
///	@param		h a pointer to a histogram (list)

//================================================================================
void HistogramCollect (simptr sim, histptr h)
// void HistogramCollect (simptr sim, histptr h, bc WALL)
//================================================================================
{
	// local simulation variables
	int	phase = sim->phase;

	// check every histogram
	while (h) {
		if (h->active) {
			if (h->stepCollect[count_] == h->stepCollect[phase]) {
			h->stepCollect[count_] = 0;
//			transformCoord(sim,WALL);
			if (h->histfunc && h->n>0) h->histfunc ((void *) sim, h, ACTION_COLLECT);
//			untransformCoord(sim,WALL);
			}
		}
		h = h->next;
	}
}


/// Loop through active histograms list and check if they are scheduled for data
/// printing. If so, call the histfunc with ACTION_PRINT if it is not NULL.
///
/// @param		sim a pointer to a simulation structure
///	@param		h a pointer to a histogram (list)

//================================================================================
void HistogramPrint (simptr sim, histptr h)
//================================================================================
{
	int	phase;

	// local simulation variables
	phase = sim->phase;

	// check every histogram
	while (h) {
		if (h->active && h->stream) {
			if (h->stepPrint[count_] == h->stepPrint[phase]) {
				h->stepPrint[count_] = 0;
				if (h->histfunc) h->histfunc ((void *) sim, h, ACTION_PRINT);
			}
		}
		h = h->next;
	}
}


/// Add the value to an existing histogram h. First we calculate the bin
/// number that corresponds to the value, then we update the appropriate bin
/// and counts values. Note that if the value is outside the histogram
/// domain, k will take a negative value. In this case we don't do anything.
///
/// @param		h a pointer to a histogram structure
/// @param		h1,h2,h3 point at which to add the value
/// @param		value the value to add in the histogram
/// @return		void

//================================================================================
void HistogramAdd (histptr h, real h1, real h2, real h3, real value)
//================================================================================
{
	int k=0;

	// determine bin number
	k = HistogramIndexValue (h, h1, h2, h3);

	// add value and increment element count
	if (k >= 0) {
		h->bin[k]   += value;
		h->count[k] += 1.0;
		h->totCount += 1.0;
	}
}


/// Returns the flat index of an histogram bin associated with the supplied
/// position.
///
/// @param		h a pointer to a histogram struture
/// @param		h1,h2,h3 point for which we want the bin
/// @return		the flat index of the requested histogram bin
/// @warning	If the point lies outside the histogram domain we return -1,
/// 			so check the returned index before using it!
/// @note		This function properly handles the request in 1D, 2D and 3D,
/// 			according to the member 'd' of the histogram.

//================================================================================
int HistogramIndexValue (histptr h, real h1, real h2, real h3)
//================================================================================
{
	int k1=0, k2=0, k3=0;

	// disregard values outside histogram bounds
	if (				(h1 < h->min1 || h1 > h->max1) )	return -1;
	if ( (h->d >= 2) && (h2 < h->min2 || h2 > h->max2) ) 	return -1;
	if ( (h->d == 3) && (h3 < h->min3 || h3 > h->max3) ) 	return -1;

	// determine bin indices depending on dimension
	if (h->d >= 1) k1 = (int) ((h1 - h->min1)/(h->bin1));
	if (h->d >= 2) k2 = (int) ((h2 - h->min2)/(h->bin2));
	if (h->d == 3) k3 = (int) ((h3 - h->min3)/(h->bin3));

	// fix max condition and finite precision spills
	if (k1 == h->n1) k1 -= 1;
	if (k2 == h->n2) k2 -= 1;
	if (k3 == h->n3) k3 -= 1;

	// calculate and return flat index of bin
	if (h->d == 1) return (k1);
	if (h->d == 2) return (k1*h->n2 + k2);
	if (h->d == 3) return (k1*h->n2 + k2*h->n3 + k3);

	// assume failure if we get here!
	return -1;
}


/// Sets all the bins and counts of a histogram to 0. This function uses a memset
/// call to reset the bin values and bin counts, and also resets the total element
/// count and the measurement count.
///
/// @param		h a pointer to a histogram structure
/// @return		void

//================================================================================
void HistogramReset (histptr h)
//================================================================================
{
	// reset bin and count values using memset
	memset (h->bin,   0, h->n*sizeof(real));
	memset (h->count, 0, h->n*sizeof(real));

	// reset histogram counters
	h->totCount = 0;
	h->m = 0;
}


/// Normalize histogram according to the value of normMode in the histogram structure.
/// This function handles the normalization of histogram values, and is fairly involved
/// owing the the number of different possible normalization. New modes can easily be
/// added to this function by defining new corresponding NORM_ macros. By far the most
/// complicated case is the NORM_VOLUME, because the bin volume must be correctly
/// handled in curvilinear geometries.
///
/// @param		h a pointer to a histogram structure
/// @return		void
/// @warning	Make sure the normMode is correctly initialized in the histogram
/// 			structure.
/// @warning	This function modifies the values collected in the histogram.

//================================================================================
void HistogramNormalize (histptr h)
//================================================================================
{
	int		i, k;
	real	*b, *c, norm, rnorm;
	real 	r1, r2, power=2.0;

	// dereference histogram bin pointer
	b = h->bin;
	c = h->count;
	if (!b) return;

	switch (h->normMode) {

		case NORM_AVERAGE:

			// divide by number of elements in each bin
			for (i=0; i<h->n; i++) {
				norm = c[i];
				if (norm != 0) b[i] /= norm;
			}
			break;

		case NORM_MEASURE:

			// divide by number of measurements
			norm = h->m;
			if (norm != 0) for (i=0; i<h->n; i++) b[i] /= norm;
			break;

		case NORM_COUNT:

			// divide by total histogram element count
			norm = h->totCount;
			if (norm != 0) for (i=0; i<h->n; i++) b[i] /= norm;
			break;

		case NORM_USER:

			// divide by user supplied h->norm factor
			norm = h->norm;
			if (norm != 0) for (i=0; i<h->n; i++) b[i] /= norm;
			break;

		case NORM_INTEGRAL:

			// normalize histogram integral to 1 (bins all have same width)
			norm = h->totCount * HistogramGetBinSize (h);
			if (norm != 0) for (i=0; i<h->n; i++) b[i] /= norm;
			break;

		case NORM_VOLUME:

			// calculate constant part of bin volume
			norm = h->m * HistogramGetBinSize (h);

			// cartesian coordinates
			if (h->coord & CARTESIAN) {
				if (norm != 0) for (i=0; i<h->n; i++) b[i] /= norm;
			}

			// cylindrical and radial coordinates
			else {

				// radial coordinates
				if (h->domain == r_ || h->domain == rq_ || h->domain == rz_ || h->domain == rqz_) {
					if (h->coord & CYLINDRICAL) {
						power = 2.0;
						if (h->domain == r_  || h->domain == rz_)  norm /= h->bin1;
						if (h->domain == rq_ || h->domain == rqz_) norm /= 2 * h->bin1;
					}
					else if (h->coord & SPHERICAL) {
						norm *= 4.0/3.0;
						power = 3.0;
						if (h->domain == r_)  norm /= h->bin1;
						if (h->domain == rq_ || h->domain == rqp_) norm /= 2 * h->bin1;
						if (h->domain == rp_ || h->domain == rqp_) norm /= pi;
					}
					k=0;
					r1 = h->min1;
					r2 = r1 + h->bin1;
					for (i=0; i<h->n; i++) {
						rnorm = norm * (pow(r2,power)-pow(r1,power));
						if (rnorm != 0) b[i] /= rnorm;
						k++;
						if (k >= h->n2*h->n3) {
							r1 = r2;
							r2 = r1 + h->bin1;
							k  = 0;
						}
					}
				}

				// non-radial coordinates
				else if (norm != 0) for (i=0; i<h->n; i++) b[i] /= norm;
			}
			break;

		default:
			break;
	}
}


/// Calculates the orthogonal volume of histogram bins. This is useful for
/// example to normalize the integral of the histogram. The definition of
/// "volume" here of course depends on the dimensionality of the histogram.
///
/// @param		h a pointer to a histogram structure
/// @return		bin size (orthogonal multiplication of the bin widths)

//================================================================================
real HistogramGetBinSize (histptr h)
//================================================================================
{
	if	  	(h->d == 1) return (h->bin1);
	else if	(h->d == 2) return (h->bin1 * h->bin2);
	else if	(h->d == 3) return (h->bin1 * h->bin2 * h->bin3);

	return 0;
}


/// Check if the coordinates of the point lie within the histogram's space
/// domain. The space domain is a region in space that defines what elements are
/// to be considered for the histogram. For example, it is possible to consider
/// only the atoms that lie in a certain region, the space domain.
///
/// @param		h a pointer to a histogram strucure
/// @param 		s1,s2,s3 the point to analyze
/// @return		1 if the point is in the space domain, 0 otherwise.
/// @warning	The space domain has its own dimension and coordinate system, so
/// 			make sure to properly transform the point before calling this
/// 			function.

//================================================================================
int HistogramSpaceDomain (histptr h, real s1, real s2, real s3)
//================================================================================
{
	// points outside space domain (reject)
	if ( (h->sd >= 1) && (s1 < h->smin1 || s1 > h->smax1) ) return 0;
	if ( (h->sd >= 2) && (s2 < h->smin2 || s2 > h->smax2) ) return 0;
	if ( (h->sd >= 3) && (s3 < h->smin3 || s3 > h->smax3) ) return 0;

	// points inside space domain (accept)
	return 1;
}


/// Adds an histogram to the beginning of a linked list of histograms.
///
/// @param		*h the address of the histogram linked list
/// @return		the address of the new histogram
/// @warning	The new histogram is inserted at the BEGINNING of the linked list.

//================================================================================
histptr HistogramNew (histptr *h)
//================================================================================
{
	histptr new;

	// allocate memory for the new simhist structure
	new = (histptr) mycalloc (1, sizeof(simhist));

	// insert new histogram at beginning of the list
	if (h) {
		new->next = *h;
		*h = new;
	}

	// return a pointer to the newly created histogram
	return new;
}


/// Duplicates the histogram pointed to by h. The duplicate is inserted just AFTER h,
/// giving and takes on the new label.
///
/// @param		h a pointer to a histogram structure
/// @param		newlabel the (unique) label for the duplicate histogram
/// @return		a pointer to the new histogram

//================================================================================
histptr HistogramDuplicate (histptr h, char *newlabel)
//================================================================================
{

	histptr new;

	if (!h) return 0;

	// allocate memory for the new simfile structure
	new = (histptr) mycalloc (1, sizeof(simhist));

	// copy histogram information and attribute new label
	memcpy   (new, h, sizeof(simhist));
	snprintf (new->label, STRMAX, "%s", newlabel);

	// allocate new bins
	HistogramAllocateBins (new);

	// copy histogram data to new histogram
	memcpy (new->bin,   h->bin,   h->n*sizeof(real));
	memcpy (new->count, h->count, h->n*sizeof(real));

	// insert new histogram in list
	new->next = h->next;
	h->next   = new;

	// return a pointer to the newly created histogram
	return new;
}


//================================================================================
histptr HistogramGet (histptr h, char *label)
//================================================================================
{
	// Finds the histogram coresponding to the label in the linked of histograms
	// pointed to by h. Returns a pointer to the histogram, or generates a fatal
	// error if the requested histogram is not found.

	int found=0;

	while (!found && h) {
		if (strcmp(label, h->label) == 0) found = 1;
		else h = h->next;
	}

	if (!found) {
		error (EHISTLIST);
	}

	return h;
}


//================================================================================
void HistogramPrintData (simptr sim, histptr h, FILE *stream)
//================================================================================
{
	// Print histogram values in the specified stream.

	int i;

	// return if the histogram is empty
	if (h->totCount <= 0) return;

	// print frame header
	fprintf (stream, "# FRAME-BEGIN\n");
	fprintf (stream, "# frame=%d\n", h->frame++);
	fprintf (stream, "# tau=%.4g\n", sim->tNow);
	fprintf (stream, "# d=%d\n", h->d);
	fprintf (stream, "# coord=0x%X\n", h->coord);
	fprintf (stream, "# domain=0x%X\n", h->domain);

	// print histogram domain
	if (h->d >= 1) fprintf (stream, "# nx=%d xmin=%g xmax=%g xbin=%g\n", h->n1, h->min1, h->max1, h->bin1);
	if (h->d >= 2) fprintf (stream, "# ny=%d ymin=%g ymax=%g ybin=%g\n", h->n2, h->min2, h->max2, h->bin2);
	if (h->d == 3) fprintf (stream, "# nz=%d zmin=%g zmax=%g zbin=%g\n", h->n3, h->min3, h->max3, h->bin3);

	// print histogram values
	for (i=0; i<h->n; i++) fprintf (stream, "%g\n", h->bin[i]);

	// print frame trailer (blank line inserted to separate sets for xmgrace)
	fprintf (stream, "# FRAME-END\n#\n\n");
}


//================================================================================
void HistogramAllocateBins (histptr h)
//================================================================================
{
	// Generic function to allocate memory for histogram data. The parameters in
	// the histogram structure pointed to by h must be set BEFORE calling this
	// function. Note that n has precedence over binSize, i.e., binSize will
	// only be used if n<=0, otherwise n will be used to compute the appropriate
	// bin size.

	// check that the number of dimensions is ok
	if (h->d<1 || h->d>3) error (EHIST);

	// histogram sizes depending on histogram dimension
	if (h->d >= 1) {
		HistogramSetBinSize (&h->min1, &h->max1, &h->n1, &h->bin1);
		h->n = h->n1;
	}
	if (h->d >= 2) {
		HistogramSetBinSize (&h->min2, &h->max2, &h->n2, &h->bin2);
		h->n = h->n1 * h->n2;
	}
	if (h->d == 3) {
		HistogramSetBinSize (&h->min3, &h->max3, &h->n3, &h->bin3);
		h->n = h->n1 * h->n2 * h->n3;
	}

	// allocate space for bins and counts
	if (h->n > 0) {
		h->bin   = (real *) mycalloc ((size_t) h->n, sizeof(real));
		h->count = (real *) mycalloc ((size_t) h->n, sizeof(real));
	}

	// reset counters
	h->totCount = 0;
	h->m 	= 0;
}


//================================================================================
void HistogramSetBinSize (real *min, real *max, int *n, real *bin)
//================================================================================
{
	// Computes the bin size parameters depending on the histogram domain and
	// number of bins. If the number of bins is zero, the bin size "bin" is used
	// to compute the number of bins.

	if (*max < *min) error (EHIST);

	if (*n > 0) {
		*bin = (*max - *min) / *n;
	}
	else if (*bin > 0){
		*n   = (int) ceil( (*max-*min) / *bin );
		*max = *min + (*n)*(*bin);
	}
    else {
        *n = 1;
        *bin = (*max - *min);
    }
}


//================================================================================
void HistogramFreeList (histptr h)
//================================================================================
{
	// Recursively frees all the memory allocated for the linked list of
	// histograms, including that for the histogram data.

	if (!h) return;
	HistogramFreeList (h->next);

	// free bins and counts
	if (h->bin)   free (h->bin);
	if (h->count) free (h->count);

	// free histogram structure
	free (h);
}


/// General histogram operation function. The requested action is specified
/// by the action flag. The possible actions are defined as macros. This function
/// contains standard operations for RESET, NORMALIZE, PRINT and COLLECT action,
/// but they can all be overided by the private function (func) for each histogram
/// (which return 1 for override, 0 to carry on).
///
/// @param		simvoid a void cast of a pointer to a simulation strucutre
/// @param		h a pointer to a histogram structure
/// @param		action the action to perform
/// @return		void

//================================================================================
void HistogramFunction (void *simvoid, histptr h, int action)
//================================================================================
{
	simptr		sim;
	int 		i, nAtom, nPolymer, n;
	particleMD	*atom, *p1, *p2;
	point		r, rh, rs;
	real		*sc1=0, *sc2=0, *sc3=0;
	real		*hc1=0, *hc2=0, *hc3=0;
	real		dx,  dy,  dz, r2, rCut2;
	item2STD	*nebrSTD;
	item2PBC	*nebrPBC;
	itemPoly	*polymer;
	real		*pbc;

	// local simulation variables
	sim	  	 = (simptr) simvoid;
	atom	 = sim->atom.items;
	nAtom 	 = sim->atom.n;
	nebrSTD	 = sim->nebrSTD.items;
	nebrPBC	 = sim->nebrPBC.items;
	polymer  = sim->polymer.items;
	nPolymer = sim->polymer.n;
	rCut2	 = sim->rCut*sim->rCut;

	// execute requested action, giving priority to private implementations
	switch (action) {

		case ACTION_RESET:

			// pass RESET to histogram function
			if (h->func && h->func ((void *) sim, h, ACTION_RESET)) break;

			// perform default reset
			HistogramReset (h);
			break;

		case ACTION_NORMALIZE:

			// pass NORMALIZE to histogram function
			if (h->func && h->func ((void *) sim, h, ACTION_NORMALIZE)) break;

			// perform default normalize
			HistogramNormalize (h);
			break;

		case ACTION_PRINT:

			// pass PRINT to histogram function
			if (h->func && h->func ((void *) sim, h, ACTION_PRINT)) break;

			// perform default print
			h->histfunc ((void *) sim, h, ACTION_NORMALIZE);
			HistogramPrintData (sim, h, h->stream);
			h->histfunc ((void *) sim, h, ACTION_RESET);
			break;

		case ACTION_COLLECT:

			// pass COLLECT to histogram function
			if (h->func && h->func ((void *) sim, h, ACTION_COLLECT)) break;

			// initialize local point variables
			memset (&r,  0, sizeof(point));
			memset (&rh, 0, sizeof(point));
			memset (&rs, 0, sizeof(point));

			// set coordinate systems for r, rh and rs
			r.coord  = CARTESIAN;
			rh.coord = h->coord;
			rs.coord = h->scoord;

			// select histogram domain
			hc1 = &(rh.c1);
			hc2 = &(rh.c2);
			hc3 = &(rh.c3);
			CoordinateOrder (h->domain, &hc1, &hc2, &hc3);

			// select space domain if applicable
			if (h->sdomain >= 0) {
				sc1 = &(rs.c1);
				sc2 = &(rs.c2);
				sc3 = &(rs.c3);
				CoordinateOrder (h->sdomain, &sc1, &sc2, &sc3);
			}

			// set histogram pointers
			h->r	 = &r;
			h->rh	 = &rh;
			h->hc1	 = hc1;
			h->hc2   = hc2;
			h->hc3   = hc3;

			// loop
			switch (h->histLoop) {


				// loop over all atoms
				case LOOP_ATOM:

					// loop over all atoms
					for (i=0; i<nAtom; i++) {

						// particle pointer
						p1 = atom+i;

						// apply include-exclude rules
						if (!(p1->group & h->groupInc)) continue;
						if (  p1->group & h->groupExc) 	continue;

						// atom position
						r.c1 = p1->rx;
						r.c2 = p1->ry;
						r.c3 = p1->rz;

						// adjust position for periodic boundaries
						ApplyPBC (sim, &r.c1, &r.c2, &r.c3);

						// filter points outside space domain if applicable
						if (h->sdomain >= 0) {
							CoordinateTransform (&r, &rs, sim->box[y_]/2.0);
							if (!HistogramSpaceDomain (h, *sc1, *sc2, *sc3)) continue;
						}

						// particle pointer and default value to add
						h->p1    = p1;
						h->value = 1.0;

						// pass COMPUTE to histogram function
						if (h->func && h->func ((void *) sim, h, ACTION_COMPUTE)) continue;

						// add value to histogram
						CoordinateTransform (&r, &rh,sim->box[y_]/2.0);
						HistogramAdd (h, *hc1, *hc2, *hc3, h->value);
					}

					// increment measurement count
					h->m += 1;
					break;


				// loop over all neighbor pairs
				case LOOP_NEBR:

					// loop over all STD neighbor pairs
					n = sim->nebrSTD.n;
					for (i=0; i<n; i++) {

						// extract particle pointers
						p1 = nebrSTD[i].p1;
						p2 = nebrSTD[i].p2;

						// apply include-exclude rules
						if (!(p1->group & h->groupInc)) continue;
						if (!(p2->group & h->groupInc)) continue;
						if (  p1->group & h->groupExc) 	continue;
						if (  p2->group & h->groupExc) 	continue;

						// disregard if not in neighbor range
						dx = (p2->rx - p1->rx);
						dy = (p2->ry - p1->ry);
						dz = (p2->rz - p1->rz);
						r2 = dx*dx + dy*dy + dz*dz;
						if (r2 > rCut2) continue;

						// particle 1 position
						r.c1 = p1->rx;
						r.c2 = p1->ry;
						r.c3 = p1->rz;

						// adjust position for periodic boundaries
						ApplyPBC (sim, &r.c1, &r.c2, &r.c3);

						// filter points outside space domain if applicable
						if (h->sdomain >= 0) {
							CoordinateTransform (&r, &rs,sim->box[y_]/2.0);
							if (!HistogramSpaceDomain (h, *sc1, *sc2, *sc3)) continue;
						}

						// particle pointers and default value
						h->p1    = p1;
						h->p2    = p2;
						h->value = 1.0;

						// pass COMPUTE to histogram function
						if (h->func && h->func ((void *) sim, h, ACTION_COMPUTE)) continue;

						// add value to histogram
						CoordinateTransform (&r, &rh,sim->box[y_]/2.0);
						HistogramAdd (h, *hc1, *hc2, *hc3, h->value);
					}


					// loop over all PBC neighbor pairs
					n = sim->nebrPBC.n;
					for (i=0; i<n; i++) {

						// extract particle pointers
						p1  = nebrPBC[i].p1;
						p2  = nebrPBC[i].p2;
						pbc = nebrPBC[i].pbc;

						// apply include-exclude rules
						if (!(p1->group & h->groupInc)) continue;
						if (!(p2->group & h->groupInc)) continue;
						if (  p1->group & h->groupExc) 	continue;
						if (  p2->group & h->groupExc) 	continue;

						// disregard if not in neighbor range
						dx = (p2->rx - p1->rx) + pbc[x_];
						dy = (p2->ry - p1->ry) + pbc[y_];
						dz = (p2->rz - p1->rz) + pbc[z_];
						r2 = dx*dx + dy*dy + dz*dz;
						if (r2 > rCut2) continue;

						// particle 1 position
						r.c1 = p1->rx;
						r.c2 = p1->ry;
						r.c3 = p1->rz;

						// adjust position for periodic boundaries
						ApplyPBC (sim, &r.c1, &r.c2, &r.c3);

						// filter points outside space domain if applicable
						if (h->sdomain >= 0) {
							CoordinateTransform (&r, &rs, sim->box[y_]/2.0);
							if (!HistogramSpaceDomain (h, *sc1, *sc2, *sc3)) continue;
						}

						// particle pointers and default value
						h->p1    = p1;
						h->p2    = p2;
						h->value = 1.0;

						// pass COMPUTE to histogram function
						if (h->func && h->func ((void *) sim, h, ACTION_COMPUTE)) continue;

						// add value to histogram
						CoordinateTransform (&r, &rh, sim->box[y_]/2.0);
						HistogramAdd (h, *hc1, *hc2, *hc3, h->value);
					}

					// increment measurement count
					h->m += 1;
					break;


				// loop over all polymers
				case LOOP_POLYMER:

					// loop over all polymers
					for (i=0; i<nPolymer; i++) {

						// particle pointer
						p1 = polymer[i].p1;
						if (!p1) continue;

						// index of first monomer
						if (polymer[i].grafted) n=0;
						else n=1;

						// pointer to first monomer
						h->p2 = p1;

						// loop over all monomers
						while (p1) {

							// apply include-exclude rules
							if (!(p1->group & h->groupInc) || (p1->group & h->groupExc)) {
								p1 = p1->next;
								continue;
							}

							// particle position
							r.c1 = p1->rx;
							r.c2 = p1->ry;
							r.c3 = p1->rz;

							// adjust position for periodic boundaries
							ApplyPBC (sim, &r.c1, &r.c2, &r.c3);

							// filter points outside space domain if applicable
							if (h->sdomain >= 0) {
								CoordinateTransform (&r, &rs, sim->box[y_]/2.0);
								if (!HistogramSpaceDomain (h, *sc1, *sc2, *sc3)) {
									p1 = p1->next;
									continue;
								}
							}

							// particle pointer and monomer index
							h->p1    = p1;
							h->value = (real) n + 0.1;

							// pass COMPUTE to histogram function
							if (h->func && h->func ((void *) sim, h, ACTION_COMPUTE)) {
								p1 = p1->next;
								continue;
							}

							// add value to histogram
							CoordinateTransform (&r, &rh, sim->box[y_]/2.0);
							HistogramAdd (h, *hc1, *hc2, *hc3, h->value);

							// move on to next monomer
							p1 = p1->next;
							n++;
						}
					}
					break;


				default:
					break;
			}


		default:
			break;
	}
}


/// Calculates the local kinetic temperature in the system, based on local
/// relative velocities between neighboring atoms. See arxiv cond-mat/0502054,
/// for example. This is a func function which overrides only the actions that
/// it defines.
///
/// @param		simvoid a void cast of a pointer to a simulation structure
/// @param		h a pointer to a histogram structure
/// @param		action the action to perform
/// @return		the continuation flag (1 to break and 0 to continue)

//================================================================================
int HistogramFunction_Tkin (void *simvoid, histptr h, int action)
//================================================================================
{
	particleMD *p1, *p2;
	real	 M;
	real	 dvx, dvy, dvz, dv2;

	switch (action) {

		case ACTION_COMPUTE:
			// particle pointers
			p1 = h->p1;
			p2 = h->p2;
			// calculate reduced mass
			M = (p1->mass*p2->mass)/(p1->mass+p2->mass);
			// calculate dv and kT
			dvx = (p2->vx - p1->vx);
			dvy = (p2->vy - p1->vy);
			dvz = (p2->vz - p1->vz);
			dv2 = dvx*dvx + dvy*dvy + dvz*dvz;
			h->value = M*dv2/3.0;
			break;

		default:
			break;
	}

	// continue with default behavior
	return !DONE;
}


/// Calculates the local configurational temperature in the system (see Evans and
/// Delhomelle's papers for more information).
///
/// @param		simvoid a void cast of a pointer to a simulation structure
/// @param		h a pointer to a histogram structure
/// @param		action the action to perform
/// @return		the continuation flag (1 to break and 0 to continue)

//================================================================================
int HistogramFunction_Tconf (void *simvoid, histptr h, int action)
//================================================================================
{
	// Calculates the configurational temperature in the system

	int			i;
	histptr		hdivf;
	particleMD	*p;
	real		*bf2, *bdivf;
	real		f2, divf;
	real		*hc1, *hc2, *hc3;

	// local histogram variables
	hdivf 	= h->h1;
	hc1 	= h->hc1;
	hc2 	= h->hc2;
	hc3 	= h->hc3;

	switch (action) {

		case ACTION_RESET:
			// reset divf histograms
			HistogramReset (hdivf);
			break;

		case ACTION_NORMALIZE:
			// get f2 and divf bins
			bf2   = h->bin;
			bdivf = hdivf->bin;
			// normalize
			for (i=0; i<h->n; i++) {
				if (bdivf[i] != 0) bf2[i] /= bdivf[i];
			}
			// skip default action behavior
			return DONE;
			break;

		case ACTION_COMPUTE:
			// get f2 and divf values
			p    = h->p1;
			f2   = p->Tfx*p->Tfx + p->Tfy*p->Tfy + p->Tfz*p->Tfz;
			divf = p->Tdivf;
			// add values in histogram
			CoordinateTransform (h->r, h->rh, 0);
			HistogramAdd (h,     *hc1, *hc2, *hc3, f2);
			HistogramAdd (hdivf, *hc1, *hc2, *hc3, divf);
			// skip default action behavior
			return DONE;
			break;

		default:
			break;
	}

	// continue with default behavior
	return !DONE;
}


/// Returns the value 1 at the coordinates corresponding to the particle's
/// velocity. Useful for generating velocity distributions!
///
/// @param		simvoid a void cast of a pointer to a simulation structure
/// @param		h a pointer to a histogram structure
/// @param		action the action to perform
/// @return		the continuation flag (1 to break and 0 to continue)

//================================================================================
int HistogramFunction_fv (void *simvoid, histptr h, int action)
//================================================================================
{
	switch (action) {

		case ACTION_COMPUTE:
			h->r->c1 = h->p1->vx;
			h->r->c2 = h->p1->vy;
			h->r->c3 = h->p1->vz;
			h->value = 1.0;
			break;

		default:
			break;
	}

	// continue with default behavior
	return !DONE;
}


/// Sets the histogram value to the particle's velocity along the x axis
///
/// @param		simvoid a void cast of a pointer to a simulation structure
/// @param		h a pointer to a histogram structure
/// @param		action the action to perform
/// @return     the continuation flag (1 to break and 0 to continue)

//================================================================================
int HistogramFunction_vx (void *simvoid, histptr h, int action)
//================================================================================
{
	// Sets the value to the x-velocity of the particle.

	switch (action) {

		case ACTION_COMPUTE:
			h->value = h->p1->vx;
			break;

		default:
			break;
	}

	// continue with default behavior
	return !DONE;
}


//================================================================================
int HistogramFunction_q (void *simvoid, histptr h, int action)
//================================================================================
{
	// Sets the value to the charge of the particle.

	switch (action) {

		case ACTION_COMPUTE:
			h->value = (real) (h->p1->q);
			break;

		default:
			break;
	}

	// continue with default behavior
	return !DONE;
}


//================================================================================
int HistogramFunction_dr_yz (void *simvoid, histptr h, int action)
//================================================================================
{
	// Returns the *relative* radial displacement of the particle pointed to by
	// p, in the yz plane. Useful to check if the cylindrical capillary wall
	// extends or contracts significantly.

	real 		d0, dy, dz, dr;
	particleMD	*p1;

	switch (action) {

		case ACTION_COMPUTE:

			p1 = h->p1;
			dy = p1->ry - p1->y0;
			dz = p1->rz - p1->z0;
			d0 = sqrt (p1->y0*p1->y0 + p1->z0*p1->z0);
			dr = (dy*p1->y0 + dz*p1->z0) / (d0*d0);

			h->value = dr;
			break;

		default:
			break;
	}

	// continue with default action behavior
	return !DONE;
}



//================================================================================
int HistogramFunction_bond (void *simvoid, histptr h, int action)
//================================================================================
{
	// Sets the value to 1 at the coordinate of the the bond length if the
	// supplied particle is bonded to another one.

	particleMD *p1, *p2;

	switch (action) {

		case ACTION_COMPUTE:

			p1 = h->p1;
			if (p1->next) {
				p2 = p1->next;
				h->r->c1 = p2->wx - p1->wx;
				h->r->c2 = p2->wy - p1->wy;
				h->r->c3 = p2->wz - p1->wz;
				h->value = 1.0;
				return !DONE;
			}
			else {
				h->value = 0.0;
				return DONE;
			}
			break;

		default:
			break;
	}

	// continue with default behavior
	return !DONE;
}

//================================================================================
int HistogramFunction_bond_orientation (void *simvoid, histptr h, int action)
//================================================================================
{
	// Give the orientation of the bond with respect to flat plates along the xz plane

	real dx, dy, dz, dr_parallel;

	particleMD *p1, *p2;

	switch (action) {

		case ACTION_COMPUTE:

			p1 = h->p1;
			if (p1->next) {
				p2 = p1->next;
				h->r->c1 = (p2->rx + p1->rx)/2.0;
				h->r->c2 = (p2->ry + p1->ry)/2.0;
				h->r->c3 = (p2->rz + p1->rz)/2.0;
				dx = p2->wx - p1->wx;
				dy = p2->wy - p1->wy;
				dz = p2->wz - p1->wz;
				dr_parallel = sqrt(dx*dx+dz*dz);
				h->value = atan(dy/dr_parallel);
				return !DONE;
			}
			else {
				h->value = 0.0;
				return DONE;
			}
			break;

		default:
			break;
	}

	// continue with default behavior
	return !DONE;
}

//================================================================================
int HistogramFunction_bond_orientation_rho_y (void *simvoid, histptr h, int action)
//================================================================================
{
	// Give the orientation of the bond with respect to flat plates along the xz plane

	real dx, dy, dz, dr_parallel, x, y, z;
	simptr sim;
	sim	  	 = (simptr) simvoid;
	particleMD *p1, *p2;

	switch (action) {

		case ACTION_COMPUTE:

			p1 = h->p1;
			if (p1->next) {
				p2 = p1->next;
				x = (p2->rx + p1->rx)/2.0;
				y = (p2->ry + p1->ry)/2.0;
				z = (p2->rz + p1->rz)/2.0;
				if (x<-sim->boxHalf[x_]) x+=sim->box[x_];
				if (z<-sim->boxHalf[z_]) z+=sim->box[z_];
				h->r->c1 = sqrt(x*x+z*z);
				h->r->c2 = y;
				dx = p2->wx - p1->wx;
				dy = p2->wy - p1->wy;
				dz = p2->wz - p1->wz;
				dr_parallel = sqrt(dx*dx+dz*dz);
				h->value = atan(dy/dr_parallel);
				return !DONE;
			}
			else {
				h->value = 0.0;
				return DONE;
			}
			break;

		default:
			break;
	}

	// continue with default behavior
	return !DONE;
}

//================================================================================
int HistogramFunction_density_rho_y (void *simvoid, histptr h, int action)
//================================================================================
{
	// Give the orientation of the bond with respect to flat plates along the xz plane

	real x, y, z;
	simptr sim;
	sim	  	 = (simptr) simvoid;

	particleMD *p1;

	switch (action) {

		case ACTION_COMPUTE:

			p1 = h->p1;
			x = p1->rx;
			y = p1->ry;
			z = p1->rz;
			if (x<-sim->boxHalf[x_]) x+=sim->box[x_];
			if (z<-sim->boxHalf[z_]) z+=sim->box[z_];
			h->r->c1 = sqrt(x*x+z*z);
			h->r->c2 = y;
			h->value = 1;
			break;

		default:
			break;
	}

	// continue with default behavior
	return !DONE;
}


//================================================================================
int HistogramFunction_fenex_N (void *simvoid, histptr h, int action)
//================================================================================
{
	real		r2, r02, k;
	real		dx, dy, dz, f, fx;
	particleMD	*p1, *p2;
	simptr		sim;

	// local sim variables
	sim = (simptr) simvoid;
	r02 = sim->r0Fene;
	r02 = r02*r02;
	k   = sim->kFene;

	switch (action) {

		case ACTION_COMPUTE:

			p1 = h->p1;
			if (p1->next) {
				p2 = p1->next;
				h->r->c1  = h->value;
				dx = p2->wx - p1->wx;
				dy = p2->wy - p1->wy;
				dz = p2->wz - p1->wz;
				r2 = dx*dx + dy*dy + dz*dz;
				f  = k/(1-r2/r02);
				fx = f * dx;
				h->value = fx;
				return !DONE;
			}
			else {
				return DONE;
			}
			break;

		default:
			break;
	}

	// continue with default behavior
	return !DONE;
}


//================================================================================
int HistogramFunction_qeff_N (void *simvoid, histptr h, int action)
//================================================================================
{

	particleMD	*p1;


	switch (action) {

		case ACTION_COMPUTE:

			p1 = h->p1;
			if (p1) {
				h->r->c1  = h->value;
				h->value = p1->q;
				return !DONE;
			}
			else {
				return DONE;
			}
			break;

		default:
			break;
	}

	// continue with default behavior
	return !DONE;
}



//================================================================================
int HistogramFunction_hwall_N (void *simvoid, histptr h, int action)
//================================================================================
{
	real		r, capr;
	real		dy, dz;
	particleMD	*p1;
	simptr		sim;
	real halfHeight;

	// local sim variables
	sim  = (simptr) simvoid;
	capr = sim->caprInTrue;
	halfHeight = sim->box[y_]/2.0;

	switch (action) {

		case ACTION_COMPUTE:

			p1 = h->p1;
			if (p1) {
				if (sim->geometry==GEOM_CAPILLARY) {
					h->r->c1 = h->value;
					dy = p1->wy;
					dz = p1->wz;
					r  = sqrt(dy*dy + dz*dz);
					h->value = capr - r;
					return !DONE;
				}
				if (sim->geometry==GEOM_PLATES) {
					h->r->c1 = h->value;
					dy = fabs(p1->wy-halfHeight);
					h->value = fabs(halfHeight-dy);
					return !DONE;
				}
			}
			else {
				return DONE;
			}
			break;

		default:
			break;
	}

	// continue with default behavior
	return !DONE;
}


//================================================================================
int HistogramFunction_ellipsoid (void *simvoid, histptr h, int action)
//================================================================================
{
	real		x, y, z;
	real		x0, y0, z0, x1, y1, z1, yn, zn, rn;
	particleMD	*p1, *p2;
// 	real		capr;
// 	simptr		sim;

	// local sim variables
// 	sim  = (simptr) simvoid;
// 	capr = sim->caprInTrue;

	switch (action) {

		case ACTION_COMPUTE:

			p1 = h->p1;
			p2 = h->p2;
			if (p1 && p2) {

				// first monomer and current monomer positions
				x0 = p2->wx;
				y0 = p2->wy;
				z0 = p2->wz;
				x1 = p1->wx;
				y1 = p1->wy;
				z1 = p1->wz;

				// normal unit vector (capillary)
				rn = sqrt(y0*y0 + z0*z0);
				yn = -y0/rn;
				zn = -z0/rn;

				// rotated components
				x =  (x1-x0);
				y = -(y1-y0)*zn - (z1-z0)*yn;
				z =  (y1-y0)*yn + (z1-z0)*zn;

				// register for histogram adding
				h->r->c1 = x;
				h->r->c2 = y;
				h->r->c3 = z;
				h->value = 1.0;

				return !DONE;
			}
			else {
				return DONE;
			}
			break;

		default:
			break;
	}

	// continue with default behavior
	return !DONE;
}


//================================================================================
int HistogramFunction_persist (void *simvoid, histptr h, int action)
//================================================================================
{
	// Calculate the persistence length of a polymer by binning the bond
	// correlation vectors.

	real		r1, r2, s, cosQ;
	real		v1x, v1y, v1z;
	real		v2x, v2y, v2z;
	particleMD	*p1, *p2;

	// first initial vertex
	p1 = h->p1;

	// loop over chain
	while (p1 && p1->next) {

		// vector 1
		v1x = p1->next->wx - p1->wx;
		v1y = p1->next->wy - p1->wy;
		v1z = p1->next->wz - p1->wz;
		r1  = sqrt(v1x*v1x + v1y*v1y + v1z*v1z);

		// arc length
		s = 0;

		// second vertex
		p2 = p1;

		while (p2 && p2->next) {

			// vector 2
			v2x = p2->next->wx - p2->wx;
			v2y = p2->next->wy - p2->wy;
			v2z = p2->next->wz - p2->wz;
			r2  = sqrt(v2x*v2x + v2y*v2y + v2z*v2z);

			// dot product
			cosQ = (v1x*v2x + v1y*v2y + v1z*v2z) / (r1*r2);

			// add value in histogram
			HistogramAdd (h, s, 0, 0, cosQ);

			// next second vertex
			p2 = p2->next;
			s += r2;
		}

		// next initial vertex
		p1 = p1->next;
	}

	return 1;
}
void transformCoord(simptr sim,bc WALL) {
	int i;
	particleMD *atom,*p;
	int nAtom = sim->atom.n;

	atom = sim->atom.items;

	for (i=0; i<nAtom; i++) {
		p = atom+i;
		p->rx -= WALL.Q[0];
		p->ry -= WALL.Q[1];
		p->rz -= WALL.Q[2];
	}
}
void untransformCoord(simptr sim,bc WALL) {
	int i;
	particleMD *atom,*p;
	int nAtom = sim->atom.n;

	atom = sim->atom.items;

	for (i=0; i<nAtom; i++) {
		p = atom+i;
		p->rx += WALL.Q[0];
		p->ry += WALL.Q[1];
		p->rz += WALL.Q[2];
	}
}
