///
/// @file
/// @brief Contains the routines for reading files, including parsing input.
///
/// A variety of helper methods for reading files and parsing input. While there are some general methods, most methods
/// here are for either parsing legacy input files (ie, input.inp etc), or for parsing the all-in-one JSON input file.

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include<math.h>

# include "../headers/SRDclss.h"
# include "../headers/definitions.h"
# include "../headers/globals.h"
# include "../headers/pout.h"
# include "../headers/mtools.h"
# include "../headers/cJson.h"

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ******* ROUTINES FOR READING FILES ******* */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */

///
/// @brief If a read from file fileName fails, then write to fsynopsis and stop the program
///
/// Helper method to automate printing of read errors to terminal from legacy input files. Primarily for debugging
/// purposes.
///
/// @param flag Flag id. <0 if reached end of file early, 0 if failed to read, >0 if read successfully.
/// @param failure String showing the line that failed to read.
/// @param file
///
void checkRead( int flag,char failure[],char file[]) {
	if(flag<0) {
		printf( "\nError: Reached end of file before read %s from %s.\n",failure,file );
		exit( 1 );
	}
	else if(flag==0) {
		printf( "\nError: Failed to read %s from %s.\n",failure,file );
		exit( 1 );
	}
}

///
/// @brief LEGACY.
///
void readin( char fpath[],inputList *in,spec **SP,particleMPC **pSRD,cell ****CL,int *MD_mode ) {
/*
   By reading in the addresses of the variables as
   pointers this function sets the values to what
   is in the input file input.inp
*/
	FILE *finput;
	int i,j,MS,read;
	double MF;
	char STR[100],inSTR[100];

	strcpy( inSTR,fpath );
	strcat( inSTR,"input.inp" );
	finput = fopen( inSTR, "r" );
	if( !finput ) {					// file couldn't be opened
		printf( "Error:\tFile '%s' could not be opened.\n",inSTR );
		exit( 1 );
	}

	//Read the dimensionality of the system
	read=fscanf( finput,"%d %s",&DIM,STR );
	checkRead( read,"dimensionality",inSTR);
	//Read system dimensions
	for( i=0; i<_3D; i++ ) {
		read=fscanf( finput,"%i %s",&XYZ[i],STR );
		checkRead( read,"system size",inSTR);
	}
	for(i=0; i<_3D; i++ ) XYZ_P1[i] = XYZ[i]+1;
	//Read the thermal energy KBT
	read=fscanf( finput,"%lf %s",&(in->KBT),STR );
	checkRead( read,"thermal energy",inSTR);
	//Read if do initial transformation to rest frame
	read=fscanf( finput,"%d %s",&(in->RFRAME),STR );
	checkRead( read,"Galilean transform",inSTR);
	//Read how often transformation to rest frame
	read=fscanf( finput,"%d %s",&(in->zeroNetMom),STR );
	checkRead( read,"rest frame",inSTR);
	//Read if do random shift of cells
	read=fscanf( finput,"%d %s",&(in->GALINV),STR );
	checkRead( read,"random shift of cells",inSTR);
	//Read thermostat mode/technique
	read=fscanf( finput,"%d %s",&(in->TSTECH),STR );
	checkRead( read,"thermostat mode",inSTR);
	//Read collision operator technique/mode
	read=fscanf( finput,"%d %s",&(in->RTECH),STR );
	checkRead( read,"collision mode",inSTR);
	//Read liquid crystal (0 if not, 1 if LC)
	read=fscanf( finput,"%d %s",&(in->LC),STR );
	checkRead( read,"liquid crystal mode",inSTR);
	//Read the thermal relaxation time scale
	read=fscanf( finput,"%lf %s",&(in->TAU),STR );
	checkRead( read,"relaxation time",inSTR);
	//Read the rotation angle used
	read=fscanf( finput,"%lf %s",&(in->RA),STR );
	checkRead( read,"rotation angle",inSTR);
	//Read the Langevin thermostat friction coefficient
	read=fscanf( finput,"%lf %s",&(in->FRICCO),STR );
	checkRead( read,"friction coefficient",inSTR);
	//Read the liquid crystal mean-field potential
	read=fscanf( finput,"%lf %s",&(in->MFPOT),STR );
	checkRead( read,"LC mean-field potential",inSTR);
	//Read the constant external acceleration
	for( i=0; i<_3D; i++ ) {
		read=fscanf( finput,"%lf %s",&(in->GRAV[i]),STR );
		checkRead( read,"acceleration",inSTR);
	}
	//Read the constant external magnetic field
	for( i=0; i<_3D; i++ ) {
		read=fscanf( finput,"%lf %s",&(in->MAG[i]),STR );
		checkRead( read,"magnetic field",inSTR);
	}
	//Read the time (or number of loops) of the warmup
	read=fscanf( finput,"%d %s",&(in->warmupSteps),STR );
	checkRead( read,"warmup time",inSTR);
	//Read the total time (or number of loops) of the simulation
	read=fscanf( finput,"%d %s",&(in->simSteps),STR );
	checkRead( read,"total time",inSTR);
	//Read the time step of one iteration
	read=fscanf( finput,"%lf %s",&(in->dt),STR );
	checkRead( read,"time step",inSTR);
	//Read the random seed (0 if read from time)
	read=fscanf( finput,"%ld %s",&(in->seed),STR );
	checkRead( read,"random seed",inSTR);
	//Read the MD coupling mode
	read=fscanf(finput, "%d %s", MD_mode, STR );
	checkRead( read,"MD coupling",inSTR);
	//Read the number of MD steps per SRD step
	read=fscanf( finput,"%d %s",&(in->stepsMD),STR );
	checkRead( read,"MD time steps per SRD step",inSTR);
	//Read the species properties
	//Read the number of species
	read=fscanf( finput,"%d %s",&NSPECI,STR );
	checkRead( read,"number of species",inSTR);
	//Allocate the needed amount of memory for the species SP
	(*SP) = (spec*) calloc( NSPECI, sizeof( spec ) );
	for( i=0; i<NSPECI; i++ ) {
		//Read the species' mass
		read=fscanf( finput,"%lf %s",&MF,STR );
		checkRead( read,"mass",inSTR);
		(*SP+i)->MASS = MF;
		//Read the species' population
		read=fscanf( finput,"%i %s",&MS,STR );
		checkRead( read,"population",inSTR);
		(*SP+i)->POP = MS;
		//Read the species' position distribution function
		read=fscanf( finput,"%i %s",&MS,STR );
		checkRead( read,"positional distribution",inSTR);
		(*SP+i)->QDIST = MS;
		//Read the species' velocity distribution function
		read=fscanf( finput,"%i %s",&MS,STR );
		checkRead( read,"velocity distribution",inSTR);
		(*SP+i)->VDIST = MS;
		//Read the species' orienation distribution function
		read=fscanf( finput,"%i %s",&MS,STR );
		checkRead( read,"orientational distribution",inSTR);
		(*SP+i)->ODIST = MS;
		//Read the binary fluid interaction matrix for this species with all other species
		for( j=0; j<NSPECI; j++ ) {
			//Read the species' interaction matrix with other species
			read=fscanf( finput,"%lf %s",&MF,STR );
			checkRead( read,"interaction matrix",inSTR);
			(*SP+i)->M[j] = MF;
		}
		//Read the species' rotational friction coefficient
		read=fscanf( finput,"%lf %s",&MF,STR );
		checkRead( read,"rotational friction",inSTR);
		(*SP+i)->RFC = MF;
		//Read the species' tumbling parameter
		read=fscanf( finput,"%lf %s",&MF,STR );
		checkRead( read,"effective rod length",inSTR);
		(*SP+i)->LEN = MF;
		//Read the species' tumbling parameter
		read=fscanf( finput,"%lf %s",&MF,STR );
		checkRead( read,"tumbling parameter",inSTR);
		(*SP+i)->TUMBLE = MF;
		//Read the species' susceptibility to shear
		read=fscanf( finput,"%lf %s",&MF,STR );
		checkRead( read,"shear susceptibility",inSTR);
		(*SP+i)->CHIHI = MF;
		//Read the species' magnetic susceptibility anisotropy
		read=fscanf( finput,"%lf %s",&MF,STR );
		checkRead( read,"magnetic susceptibility",inSTR);
		(*SP+i)->CHIA = MF;
		//Read the species' activity
		read=fscanf( finput,"%lf %s",&MF,STR );
		checkRead( read,"activity",inSTR);
		(*SP+i)->ACT = MF;
		//Read the species' damping friction coefficient
		read=fscanf( finput,"%lf %s",&MF,STR );
		checkRead( read,"damping friction",inSTR);
		(*SP+i)->DAMP = MF;
	}
	//Total Number of particleMPCs
	GPOP = 0;
	for( i=0; i<NSPECI; i++ ) GPOP += (*SP+i)->POP;
	(*pSRD) = (particleMPC*) calloc( GPOP, sizeof( particleMPC ) );

	//Allocate memory for the cells
	//Allocate rows (x first)
	*CL = (cell***) calloc( XYZ_P1[0], sizeof( cell** ) );
	//For each x-element allocate the y columns
	for( i=0; i<XYZ_P1[0]; i++ ) {
		(*CL)[i] = (cell**) calloc( XYZ_P1[1], sizeof( cell* ) );
		//For each y-element allocate the z columns
		for( j=0; j<XYZ_P1[1]; j++ ) {
			(*CL)[i][j] = (cell*) calloc( XYZ_P1[2], sizeof( cell ) );
		}
	}

	fclose( finput );
}

///
/// @brief LEGACY.
///
void readpc( char fpath[],outputFlagsList *out ) {
/*
   By reading in the addresses of the variables as
   pointers this function sets the flags which tells
  the program what data to output
*/
	FILE *finput;
	char STR[100],inSTR[100];
	int read;

	strcpy( inSTR,fpath );
	strcat( inSTR,"printcom.inp" );
	finput = fopen( inSTR,"r" );
	if( !finput ) {				// file couldn't be opened
		printf( "Error:\tFile '%s' could not be opened.\n",inSTR );
		exit( 1 );
	}
	//Read debug mode
	read=fscanf( finput,"%d %s",&DBUG,STR );
	checkRead( read,"debug mode",inSTR);
	//Read how often the detailed species' trajectories are printed
	read=fscanf( finput,"%d %s",&(out->TRAJOUT),STR );
	checkRead( read,"detailed species' trajectories",inSTR);
	//Read which (number of) species whose detailed trajectories are printed
	read=fscanf( finput,"%d %s",&(out->printSP),STR );
	checkRead( read,"number of species whose detailed trajectories are printed",inSTR);
	//Read how often the coarse data is outputted
	read=fscanf( finput,"%d %s",&(out->COAROUT),STR );
	checkRead( read,"coarse data",inSTR);
	//Read how often the flow field is outputted
	read=fscanf( finput,"%d %s",&(out->FLOWOUT),STR );
	checkRead( read,"flow field",inSTR);
	//Read how often the total average MPCD velocity is outputted
	read=fscanf( finput,"%d %s",&(out->AVVELOUT),STR );
	checkRead( read,"total average MPCD velocity",inSTR);
	//Read how often the local director and scalar order parameter fields are outputted
	read=fscanf( finput,"%d %s",&(out->ORDEROUT),STR );
	checkRead( read,"director and scalar order parameter fields",inSTR);
	//Read how often the local order parameter tensor field is outputted
	read=fscanf( finput,"%d %s",&(out->QTENSOUT),STR );
	checkRead( read,"order parameter tensor field",inSTR);
	//Read how often the local order parameter tensor field in reciprical space is outputted
	read=fscanf( finput,"%d %s",&(out->QKOUT),STR );
	checkRead( read,"order parameter tensor field in reciprical space",inSTR);
	//Read how often the orientational energy field is outputted
	read=fscanf( finput,"%d %s",&(out->ENFIELDOUT),STR );
	checkRead( read,"orientational energy field",inSTR);
	//Read how often the colour/conc/phi field is outputted
	read=fscanf( finput,"%d %s",&(out->SPOUT),STR );
	checkRead( read,"colour/conc field",inSTR);
	//Read how often the pressure field is outputted
	read=fscanf( finput,"%d %s",&(out->PRESOUT),STR );
	checkRead( read,"pressure field",inSTR);
	//Read how often the orientational energy from neighbours is outputted
	read=fscanf( finput,"%d %s",&(out->ENNEIGHBOURS),STR );
	checkRead( read,"orientational energy from neighbours",inSTR);
	//Read how often the total average scalar order parameter is outputted
	read=fscanf( finput,"%d %s",&(out->AVSOUT),STR );
	checkRead( read,"average scalar order parameter",inSTR);
	//Read how often the standard deviation of the number per cell is outputted
	read=fscanf( finput,"%d %s",&(out->DENSOUT),STR );
	checkRead( read,"standard deviation of the number per cell",inSTR);
	//Read how often the total average enstrophy is outputted
	read=fscanf( finput,"%d %s",&(out->ENSTROPHYOUT),STR );
	checkRead( read,"total average enstrophy",inSTR);
	//Read how often the velocity distribution should be outputted
	read=fscanf( finput,"%d %s",&(out->HISTVELOUT),STR );
	checkRead( read,"velocity distribution",inSTR);
	//Read how often the speed distribution should be outputted
	read=fscanf( finput,"%d %s",&(out->HISTSPEEDOUT),STR );
	checkRead( read,"speed distribution",inSTR);
	//Read how often the vorticity distribution should be outputted
	read=fscanf( finput,"%d %s",&(out->HISTVORTOUT),STR );
	checkRead( read,"vorticity distribution",inSTR);
	//Read how often the enstrophy distribution should be outputted
	read=fscanf( finput,"%d %s",&(out->HISTENSTROUT),STR );
	checkRead( read,"enstrophy distribution",inSTR);
	//Read how often the director distribution should be outputted
	read=fscanf( finput,"%d %s",&(out->HISTDIROUT),STR );
	checkRead( read,"director distribution",inSTR);
	//Read how often the scalar order parameter distribution should be outputted
	read=fscanf( finput,"%d %s",&(out->HISTSOUT),STR );
	checkRead( read,"scalar order parameter distribution",inSTR);
	//Read how often the number per cell distribution should be outputted
	read=fscanf( finput,"%d %s",&(out->HISTNOUT),STR );
	checkRead( read,"number per cell distribution",inSTR);
	//Read how often the solid trajectory data is outputted
	read=fscanf( finput,"%d %s",&(out->SOLOUT),STR );
	checkRead( read,"solid trajectory",inSTR);
	//Read how often the topological charge field data is outputted
	read=fscanf( finput,"%d %s",&(out->TOPOOUT),STR );
	checkRead( read,"topological charge field",inSTR);
	//Read how often the defects positional data is outputted
	read=fscanf( finput,"%d %s",&(out->DEFECTOUT),STR );
	checkRead( read,"defect positions",inSTR);
	//Read how often the disclination tensor field data is outputted
	read=fscanf( finput,"%d %s",&(out->DISCLINOUT),STR );
	checkRead( read,"disclination tensor field",inSTR);
	//Read how often the energy data is outputted
	read=fscanf( finput,"%d %s",&(out->ENOUT),STR );
	checkRead( read,"energy",inSTR);
	//Read how often the velocity-velocity correlation data is outputted
	read=fscanf( finput,"%d %s",&(out->CVVOUT),STR );
	checkRead( read,"velocity-velocity correlation",inSTR);
	//Read how often the director-director correlation data is outputted
	read=fscanf( finput,"%d %s",&(out->CNNOUT),STR );
	checkRead( read,"director-director correlation",inSTR);
	//Read how often the vorticity-vorticity correlation data is outputted
	read=fscanf( finput,"%d %s",&(out->CWWOUT),STR );
	checkRead( read,"vorticity-vorticity correlation",inSTR);
	//Read how often the density-density correlation data is outputted
	read=fscanf( finput,"%d %s",&(out->CDDOUT),STR );
	checkRead( read,"density-density correlation",inSTR);
	//Read how often the order-order correlation data is outputted
	read=fscanf( finput,"%d %s",&(out->CSSOUT),STR );
	checkRead( read,"order-order correlation",inSTR);
	//Read how often the phase-phase (binary fluid) correlation data is outputted
	read=fscanf( finput,"%d %s",&(out->CPPOUT),STR );
	checkRead( read,"phase-phase (binary fluid) correlation",inSTR);
	//Read how often the energy spectrum data is outputted
	read=fscanf( finput,"%d %s",&(out->ENERGYSPECTOUT),STR );
	checkRead( read,"energy spectrum",inSTR);
	//Read how often the enstrophy spectrum data is outputted
	read=fscanf( finput,"%d %s",&(out->ENSTROPHYSPECTOUT),STR );
	checkRead( read,"enstrophy spectrum",inSTR);
	//Read how often the Binder cumulant is outputted
	read=fscanf( finput,"%d %s",&(out->BINDER),STR );
	checkRead( read,"Binder cumulant",inSTR);
	//Read the Binder cumulant bin size
	read=fscanf( finput,"%d %s",&(out->BINDERBIN),STR );
	checkRead( read,"Binder cumulant bin size",inSTR);
	//Read how often the swimmers' positions are outputted
	read=fscanf( finput,"%d %s",&(out->SWOUT),STR );
	checkRead( read,"swimmers' positions",inSTR);
	//Read how often the swimmers' orientations are outputted
	read=fscanf( finput,"%d %s",&(out->SWORIOUT),STR );
	checkRead( read,"swimmers' orientations",inSTR);
	//Read if the swimmers' run/tumble statistics are recorded
	read=fscanf( finput,"%d %s",&(out->RTOUT),STR );
	checkRead( read,"swimmers' run/tumble",inSTR);
	//Read if a synopsis printed
	read=fscanf( finput,"%d %s",&(out->SYNOUT),STR );
	checkRead( read,"synopsis",inSTR);
	//Read if checkpointing is on
	read=fscanf( finput,"%d %s",&(out->CHCKPNT),STR );
	checkRead( read,"checkpointing",inSTR);

	fclose( finput );
}

///
/// @brief LEGACY.
///
void bcin( FILE *fbc,bc *WALL,char fname[] ) {
/*
   By reading in the addresses of the variables as
   pointers this function sets the values to what
   is in the input file bc.inp
 */
	int x,read;
	double l;
	char LABEL[160];

	read=fscanf( fbc,"%s",LABEL );
	checkRead( read,"bc",fname);

	read=fscanf( fbc,"%d %s",&x,LABEL );
	checkRead( read,"bc",fname);
	WALL->COLL_TYPE = x;
	read=fscanf( fbc,"%d %s",&x,LABEL );
	checkRead( read,"bc",fname);
	WALL->PHANTOM = x;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->E = l;

	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->Q[0] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->Q[1] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->Q[2] = l;

	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->V[0] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->V[1] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->V[2] = l;

	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->O[0] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->O[1] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->O[2] = l;

	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->L[0] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->L[1] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->L[2] = l;

	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->G[0] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->G[1] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->G[2] = l;

	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->AINV[0] = l;
	if( fneq(WALL->AINV[0],0.0) ) WALL->A[0] = 1.0/WALL->AINV[0];
	else WALL->A[0] = 0.0;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->AINV[1] = l;
	if( fneq(WALL->AINV[1],0.0) ) WALL->A[1] = 1.0/WALL->AINV[1];
	else WALL->A[1] = 0.0;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->AINV[2] = l;
	if( fneq(WALL->AINV[2],0.0) ) WALL->A[2] = 1.0/WALL->AINV[2];
	else WALL->A[2] = 0.0;

	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->ROTSYMM[0] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->ROTSYMM[1] = l;
	read=fscanf( fbc,"%d %s",&x,LABEL );
	checkRead( read,"bc",fname);
	WALL->ABS = x;

	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->P[0] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->P[1] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->P[2] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->P[3] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->R = l;

	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->DN = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->DT = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->DVN = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->DVT = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->DVxyz[0] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->DVxyz[1] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->DVxyz[2] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->MVN = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->MVT = l;

	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->MUN = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->MUT = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->MUxyz[0] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->MUxyz[1] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->MUxyz[2] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->DUxyz[0] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->DUxyz[1] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->DUxyz[2] = l;

	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->KBT = l;
	read=fscanf( fbc,"%d %s",&x,LABEL );
	checkRead( read,"bc",fname);
	WALL->DSPLC = x;
	read=fscanf( fbc,"%d %s",&x,LABEL );
	checkRead( read,"bc",fname);
	WALL->INV = x;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->MASS = l;

	// Set the planar flag
	WALL->PLANAR = 0;							// sets default to non-planar
	if ( feq(WALL->B[0],0.0) ) {	// if not wavy, check if planar
		if( feq(WALL->P[0],1.0) && feq(WALL->P[1],1.0) && feq(WALL->P[2],1.0) ) {
			// Left or right wall
			if( fneq(WALL->A[0],0.0) && feq(WALL->A[1],0.0) && feq(WALL->A[2],0.0) ) WALL->PLANAR = 1;
			// Top or bottom wall
			else if( feq(WALL->A[0],0.0) && fneq(WALL->A[1],0.0) && feq(WALL->A[2],0.0) ) WALL->PLANAR = 1;
			// Far or near wall
			else if( feq(WALL->A[0],0.0) && feq(WALL->A[1],0.0) && fneq(WALL->A[2],0.0) ) WALL->PLANAR = 1;
		}
	}
}

///
/// @brief Determines if a boundary is a planar periodic boundary.
///
/// Performs a series of checks, in order, to verify if boundary condition is periodic. This is then marked
/// as a planar PBC in `XYZPBC`.
///
/// @param WALL Pointer to a particular boundary condition.
///
void setPBC( bc *WALL ) {
	int i;
	//Check if any axis is a planar PBC
	for( i=0; i<_3D; i++ ) {
		// Is the BC a plane?
		if( feq(WALL->P[i],1.0) ) {
			// Does it point in the correct direction?
			if( feq( fabs(WALL->AINV[i]),1.0 ) ) {
				// Does it keep the velocity pointed in the correct direction?
				if( WALL->MVN>0.0 && WALL->MVT>0.0 ) {
					//Does it modify the position?
					if( fneq(WALL->DN,0.0) || fneq(WALL->DT,0.0) ) {
						//Then it is a periodic boundary condition
						XYZPBC[i]=1;
					}
				}
			}
		}
	}
}

///
/// @brief LEGACY.
///
void readbc( char fpath[],bc **WALL ) {
/*
   By calling bcin this function sets each of
   the boundaries
*/
	FILE *fbc;
	int i,read;
	char STR[100],inSTR[100];

	strcpy( inSTR,fpath );
	strcat( inSTR,"bc.inp" );
	fbc = fopen(inSTR, "r");
	if( !fbc ) { // file couldn't be opened
		printf( "Error:\tFile '%s' could not be opened.\n",inSTR );
		exit( 1 );
	}

	read=fscanf( fbc,"%s",STR );
	checkRead( read,"bc",inSTR);
	//Read the number of boundary conditions being used
	read=fscanf( fbc,"%d %s",&NBC,STR );
	checkRead( read,"bc",inSTR);
	//Allocate the needed amount of memory for the boundary conditions WALLs
	(*WALL) = (bc*) calloc( NBC, sizeof( bc ) );
	//Read the parameters for each of the BCs
	for( i=0; i<NBC; i++ ) bcin( fbc,(*WALL+i),inSTR );
	//Determine if any BCs are periodic boundaries
	for( i=0; i<_3D; i++ ) XYZPBC[i]=0;
	for( i=0; i<NBC; i++ ) setPBC( (*WALL+i) );
	fclose( fbc );
}

///
/// @brief Reads a checkpoint to resume a simulation.
///
/// Reads the entirety of a checkpoint file and uses it to re-populate the simulation. The method iteratively goes
/// through all aspects of the simulation and reads them in, allocating memory as and when necessary.
/// This is used to resume an existing simulation. The only thing that is not checkpointed is the random number
/// generator state, which is re-seeded outside of this routine.
///
/// @param fpath Path to the directory where the checkpoint file is.
/// @param in Pointer to object containing input parameters (corresponding to legacy input.inp). Expected to be &in.
/// @param SP Pointer to the species list. Expected to be &SP.
/// @param pSRD Pointer to the particle list. Expected to be &pSRD.
/// @param CL Pointer to the cell array. Expected to be &CL.
/// @param MD_mode Pointer to the MD mode int flag. Expected to be &MD_mode.
/// @param WALL Pointer to the boundary condition list. Expected to be &WALL.
/// @param out Pointer to the object containing output flags (corresponding to legacy printcom.inp). Expected to be &out.
/// @param runtime Pointer to the current runtime of the simulation. Expected to be &runtime.
/// @param warmtime Pointer to the warmup timesteps performed by the simulation. Expected to be &warmtime.
/// @param theory Pointer to the object containing the kinetic theory parameters. Expected to be &theory.
/// @param AVVEL Pointer to the average velocity of the system. Expected to be &AVVEL.
/// @param AVS Pointer to the average scalar order parameter of the system. Expected to be &AVS.
/// @param avDIR Pointer to the average vector director of the system. Expected to be &avDIR.
/// @param S4 Pointer to the fourth moment of the scalar order parameter of the system. Expected to be &S4.
/// @param stdN Pointer to the standard deviation of the number of particles in the system. Expected to be &stdN.
/// @param KBTNOW Pointer to the current temperature of the system. Expected to be &KBTNOW.
/// @param AVV Pointer to the average velocity of the system. Expected to be &AVV.
/// @param AVNOW Pointer to the average velocity now of the system. Expected to be &AVNOW.
/// @param specS Pointer to the object containing the swimmer species hyperparameters. Expected to be &specS.
/// @param sw Pointer to the swimmer list. Expected to be &sw.
///
void readchckpnt(char fpath[], inputList *in, spec **SP, particleMPC **pSRD, cell ****CL, int *MD_mode, bc **WALL, outputFlagsList *out, int *runtime, int *warmtime, kinTheory *theory, double *AVVEL, double *AVS, double avDIR[_3D], double *S4, double *stdN, double *KBTNOW, double AVV[_3D], double AVNOW[_3D], specSwimmer *specS, swimmer **sw ) {
	FILE *finput;
	int i,j;
	char STR[100];

	strcpy( STR,fpath );
	strcat( STR,"checkpoint.dat" );
	finput = fopen( STR, "r" );
	if( !finput ) {					// file couldn't be opened
		printf( "Error:\tFile '%s' could not be opened.\n",STR );
		exit( 1 );
	}
	if(fscanf( finput,"%d",&(in->simSteps) ));		//Read time
	else printf("Warning: Failed to sim read time.\n");
	if(fscanf( finput,"%d %lf",&(in->warmupSteps),&(in->dt) ));		//Read time
	else printf("Warning: Failed to read time.\n");

	if(fscanf( finput,"%ld",&(in->seed) ));	//Read the random seed (0 if read from time)
	else printf("Warning: Failed to read random seed.\n");
	if(fscanf( finput,"%d %d %d %d %lf %lf",&DIM,&XYZ[0],&XYZ[1],&XYZ[2],&(in->KBT),KBTNOW ));
	else printf("Warning: Failed to read dimensionality, size or temperature.\n");
	for(i=0; i<_3D; i++ ) XYZ_P1[i] = XYZ[i]+1;
	if(fscanf( finput,"%d %d %d %d %d %d",&(in->RFRAME),&(in->zeroNetMom),&(in->GALINV),&(in->TSTECH),&(in->RTECH),&(in->LC) ));
	else printf("Warning: Failed to Galilean transform, rest frame, thermostat mode, collision mode or liquid crystal mode.\n");
	if(fscanf( finput,"%lf %lf %lf %lf",&(in->TAU),&(in->RA),&(in->FRICCO),&(in->MFPOT) ));		//Read the thermal relaxation time scale
	else printf("Warning: Failed to read relaxation time, rotation angle, friction coefficient or mean-field potential.\n");
	if(fscanf( finput,"%d %d %d",&(in->noHI),&(in->inCOMP),&(in->MULTIPHASE) ));		//Read no hydrodynamics, incompressibility and multi-phase
	else printf("Warning: Failed to read no hydrodynamics or incompressibility.\n");
	if(fscanf( finput,"%lf %lf %lf",&(in->GRAV[0]),&(in->GRAV[1]),&(in->GRAV[2]) ));	//Read the constant external acceleration
	else printf("Warning: Failed to read acceleration.\n");
	if(fscanf( finput,"%lf %lf %lf",&(in->MAG[0]),&(in->MAG[1]),&(in->MAG[2]) ));	//Read the constant external magnetic field
	else printf("Warning: Failed to read magnetic field.\n");


	if(fscanf(finput, "%d %d", MD_mode, &(in->stepsMD) ));	//Read the MD coupling mode
	else printf("Warning: Failed to read MD coupling.\n");
	if(fscanf( finput,"%d %d",&GPOP,&NSPECI ));	//Read the number of MPC particles
	else printf("Warning: Failed to read total number of particles or number of species.\n");

	if(fscanf( finput,"%d %d %lf %lf %d %d",runtime,warmtime,&(in->C),&(in->S),&(in->GRAV_FLAG),&(in->MAG_FLAG) ));//Read program variables
	else printf("Warning: Failed to read various program variables.\n");
	if(fscanf( finput,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&(theory->MFP), &(theory->VISC), &(theory->THERMD), &(theory->SDIFF), &(theory->SPEEDOFSOUND), &(theory->sumM), AVVEL, AVS, &avDIR[0], &avDIR[1], &avDIR[2], S4, stdN, &nDNST, &mDNST ));//Read program variables
	else printf("Warning: Failed to read various program variables.\n");
	if(fscanf( finput,"%lf %lf %lf %lf %lf %lf",&AVV[0], &AVV[1], &AVV[2], &AVNOW[0], &AVNOW[1], &AVNOW[2] ));//Read program variables
	else printf("Warning: Failed to read average velocities.\n");

	//Read output
	if(fscanf( finput,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",&DBUG, &(out->TRAJOUT), &(out->printSP), &(out->COAROUT), &(out->FLOWOUT), &(out->VELOUT), &(out->AVVELOUT), &(out->ORDEROUT), &(out->QTENSOUT), &(out->QKOUT), &(out->AVSOUT), &(out->SOLOUT), &(out->ENOUT), &(out->ENFIELDOUT), &(out->ENNEIGHBOURS), &(out->ENSTROPHYOUT), &(out->DENSOUT), &(out->CVVOUT), &(out->CNNOUT), &(out->CWWOUT), &(out->CDDOUT), &(out->CSSOUT), &(out->CPPOUT), &(out->BINDER), &(out->BINDERBIN), &(out->SYNOUT), &(out->CHCKPNT), &(out->CHCKPNTrcvr) ));
	else printf("Warning: Failed to read output.\n");
	if(fscanf( finput,"%d %d",&(out->SPOUT), &(out->PRESOUT) ));
	else printf("Warning: Failed to read output.\n");
	if(fscanf( finput,"%d %d %d %d %d %d %d",&(out->HISTVELOUT), &(out->HISTSPEEDOUT), &(out->HISTVORTOUT), &(out->HISTENSTROUT), &(out->HISTDIROUT), &(out->HISTSOUT), &(out->HISTNOUT) ));
	else printf("Warning: Failed to read histogram output.\n");
	if(fscanf( finput,"%d %d %d %d %d",&(out->ENERGYSPECTOUT), &(out->ENSTROPHYSPECTOUT), &(out->TOPOOUT), &(out->DEFECTOUT), &(out->DISCLINOUT) ));
	else printf("Warning: Failed to read histogram output.\n");
	if(fscanf( finput,"%d %d %d",&(out->SWOUT), &(out->SWORIOUT), &(out->RTOUT) ));
	else printf("Warning: Failed to read histogram output.\n");

	//Allocate the needed amount of memory for the species SP
	(*SP) = (spec*) calloc( NSPECI, sizeof( spec ) );
	for( i=0; i<NSPECI; i++ ) {
		if(fscanf( finput,"%lf %i %i %i %i %lf %lf %lf %lf %lf %lf %lf %lf %lf",&((*SP+i)->MASS), &((*SP+i)->POP), &((*SP+i)->QDIST), &((*SP+i)->VDIST), &((*SP+i)->ODIST), &((*SP+i)->RFC), &((*SP+i)->LEN), &((*SP+i)->TUMBLE), &((*SP+i)->CHIHI), &((*SP+i)->CHIA), &((*SP+i)->ACT), &((*SP+i)->SIGWIDTH), &((*SP+i)->SIGPOS), &((*SP+i)->DAMP) ));	//Read the species' mass
		else printf("Warning: Failed to read species %i.\n",i);
		for( j=0; j<NSPECI; j++ ) {
			//Read the species' interaction matrix with other species
			if(fscanf( finput,"%lf ",&((*SP+i)->M[j]) ));	//Read the species' interactions
			else printf("Warning: Failed to read species %d interaction with %d.\n",i,j);
		}
	}
	//Check total number of particleMPCs
	j = 0;
	for( i=0; i<NSPECI; i++ ) j += (*SP+i)->POP;
	if(GPOP!=j) {
		printf("Warning: GPOP does not match sum of species populations.\n");
		GPOP=j;
	}
	(*pSRD) = (particleMPC*) calloc( GPOP, sizeof( particleMPC ) );

	//Allocate memory for the cells
	//Allocate rows (x first)
	*CL = (cell***) calloc( XYZ_P1[0], sizeof( cell** ) );
	//For each x-element allocate the y columns
	for( i=0; i<XYZ_P1[0]; i++ ) {
		(*CL)[i] = (cell**) calloc( XYZ_P1[1], sizeof( cell* ) );
		//For each y-element allocate the z columns
		for( j=0; j<XYZ_P1[1]; j++ ) {
			(*CL)[i][j] = (cell*) calloc( XYZ_P1[2], sizeof( cell ) );
		}
	}
	//Read the BCs
	if(fscanf( finput,"%d",&NBC ));		//Read the number of BCs
	else printf("Warning: Failed to read number of BCs.\n");
	for( i=0; i<NBC; i++ ) {
		if(fscanf( finput,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&((*WALL+i)->COLL_TYPE), &((*WALL+i)->PHANTOM), &((*WALL+i)->E), &((*WALL+i)->Q[0]), &((*WALL+i)->Q[1]), &((*WALL+i)->Q[2]), &((*WALL+i)->V[0]), &((*WALL+i)->V[1]), &((*WALL+i)->V[2]), &((*WALL+i)->O[0]), &((*WALL+i)->O[1]), &((*WALL+i)->O[2]) ));
		else printf("Warning: Failed to read BC %d.\n",i);
		if(fscanf( finput,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &((*WALL+i)->L[0]), &((*WALL+i)->L[1]), &((*WALL+i)->L[2]), &((*WALL+i)->G[0]), &((*WALL+i)->G[1]), &((*WALL+i)->G[2]), &((*WALL+i)->A[0]), &((*WALL+i)->A[1]), &((*WALL+i)->A[2]), &((*WALL+i)->AINV[0]), &((*WALL+i)->AINV[1]), &((*WALL+i)->AINV[2]), &((*WALL+i)->P[0]),&((*WALL+i)->P[1]),&((*WALL+i)->P[2]),&((*WALL+i)->P[3]), &((*WALL+i)->R) ));
		else printf("Warning: Failed to read BC %d.\n",i);
		if(fscanf( finput,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &((*WALL+i)->DN), &((*WALL+i)->DT), &((*WALL+i)->DVN), &((*WALL+i)->DVT), &((*WALL+i)->DVxyz[0]), &((*WALL+i)->DVxyz[1]), &((*WALL+i)->DVxyz[2]), &((*WALL+i)->MVN), &((*WALL+i)->MVT), &((*WALL+i)->MUN), &((*WALL+i)->MUT), &((*WALL+i)->MUxyz[0]), &((*WALL+i)->MUxyz[1]), &((*WALL+i)->MUxyz[2]) ));
		else printf("Warning: Failed to read BC %d.\n",i);
		if(fscanf( finput,"%lf %lf %lf %lf %d %d %lf", &((*WALL+i)->DUxyz[0]), &((*WALL+i)->DUxyz[1]), &((*WALL+i)->DUxyz[2]), &((*WALL+i)->KBT), &((*WALL+i)->DSPLC), &((*WALL+i)->INV), &((*WALL+i)->MASS) ));
		else printf("Warning: Failed to read BC %d.\n",i);
		if(fscanf( finput,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &((*WALL+i)->W), &((*WALL+i)->VOL), &((*WALL+i)->Q_old[0]), &((*WALL+i)->Q_old[1]), &((*WALL+i)->Q_old[2]), &((*WALL+i)->O_old[0]), &((*WALL+i)->O_old[1]), &((*WALL+i)->O_old[2]), &((*WALL+i)->I[0][0]), &((*WALL+i)->I[0][1]), &((*WALL+i)->I[0][2]), &((*WALL+i)->I[1][0]), &((*WALL+i)->I[1][1]), &((*WALL+i)->I[1][2]), &((*WALL+i)->I[2][0]), &((*WALL+i)->I[2][1]), &((*WALL+i)->I[2][2]) ));
		else printf("Warning: Failed to read BC %d.\n",i);
		if(fscanf( finput,"%d %d %d %lf %lf", &((*WALL+i)->PLANAR), &((*WALL+i)->REORIENT), &((*WALL+i)->ABS), &((*WALL+i)->ROTSYMM[0]), &((*WALL+i)->ROTSYMM[1]) ));
		else printf("Warning: Failed to read BC %d.\n",i);
		if(fscanf( finput,"%lf %lf %lf %lf %lf %lf", &((*WALL+i)->dV[0]), &((*WALL+i)->dV[1]), &((*WALL+i)->dV[2]), &((*WALL+i)->dL[0]), &((*WALL+i)->dL[1]), &((*WALL+i)->dL[2]) ));
		else printf("Warning: Failed to read BC %d.\n",i);
	}

	//Read the MPCD particles
	for( i=0; i<GPOP; i++ ) {
		if(fscanf( finput,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&((*pSRD+i)->S_flag), &((*pSRD+i)->SPID), &((*pSRD+i)->q), &((*pSRD+i)->Q[0]), &((*pSRD+i)->Q[1]), &((*pSRD+i)->Q[2]), &((*pSRD+i)->V[0]), &((*pSRD+i)->V[1]), &((*pSRD+i)->V[2]), &((*pSRD+i)->U[0]), &((*pSRD+i)->U[1]), &((*pSRD+i)->U[2]), &((*pSRD+i)->T[0]), &((*pSRD+i)->T[1]), &((*pSRD+i)->T[2]) ));
		else printf("Warning: Failed to read MPCD particle %d.\n",i);
	}

	//Swimmers
	if(fscanf( finput,"%d %d %d %d %d %d %lf %lf %d %d",&NS, &(specS->TYPE), &(specS->QDIST), &(specS->ODIST), &(specS->headM), &(specS->middM), &(specS->iheadM), &(specS->imiddM), &(specS->HSPid), &(specS->MSPid) ));
	else printf("Warning: Failed to read swimmer-type variables.\n");
	if(fscanf( finput,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf", &(specS->FS), &(specS->TS), &(specS->DS), &(specS->sizeShrink), &(specS->springShrink), &(specS->fixDist), &(specS->k), &(specS->ro), &(specS->iro), &(specS->sig), &(specS->isig), &(specS->eps), &(specS->runTime), &(specS->tumbleTime), &(specS->shrinkTime), &(specS->MAGMOM) ));
	else printf("Warning: Failed to read swimmer-type variables.\n");

	//Allocate the memory for the swimmers
	(*sw) = (swimmer*) calloc( NS, sizeof( swimmer ) );

	for( i=0; i<NS; i++ ) {
		if(fscanf( finput,"%d %lf %lf %lf %d %d %lf %lf %lf %lf %lf\n", &((*sw+i)->RT), &((*sw+i)->n0[0]), &((*sw+i)->n0[1]), &((*sw+i)->n0[2]), &((*sw+i)->timeCNT), &((*sw+i)->timeRND), &((*sw+i)->ro), &((*sw+i)->iro), &((*sw+i)->sig), &((*sw+i)->isig), &((*sw+i)->k) ));
		else printf("Warning: Failed to read swimmer %d.\n",i);
		if(fscanf( finput,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &((*sw+i)->H.HorM), &((*sw+i)->H.Q[0]), &((*sw+i)->H.Q[1]), &((*sw+i)->H.Q[2]), &((*sw+i)->H.V[0]), &((*sw+i)->H.V[1]), &((*sw+i)->H.V[2]), &((*sw+i)->H.A[0]), &((*sw+i)->H.A[1]), &((*sw+i)->H.A[2]) ));
		else printf("Warning: Failed to read swimmer %d.\n",i);
		if(fscanf( finput,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &((*sw+i)->M.HorM), &((*sw+i)->M.Q[0]), &((*sw+i)->M.Q[1]), &((*sw+i)->M.Q[2]), &((*sw+i)->M.V[0]), &((*sw+i)->M.V[1]), &((*sw+i)->M.V[2]), &((*sw+i)->M.A[0]), &((*sw+i)->M.A[1]), &((*sw+i)->M.A[2]) ));
		else printf("Warning: Failed to read swimmer %d.\n",i);
	}

	fclose( finput );
}

///
/// @brief Read the arguments for the program.
///
/// Reads the arguments for the main program. All arguments are of form `-x [value]` where x is the argument and value
/// is the value of the argument. The arguments can be listed in any order.
///
/// @param argc The C argc variable.
/// @param argv The C argv variable.
/// @param ip The input file name, read from the arguments.
/// @param op The output directory, read from the arguments.
/// @param inMode The input mode, read from the arguments. Switches between legacy (.inp) and new (.json) input files.
///
void readarg( int argc, char* argv[], char ip[],char op[], int *inMode ) {
	int arg;
	int strln;

    // input/ output checkers
    int inProvided = 0;
    int outProvided = 0;

	// Default
	#ifdef DBG
		if( DBUG >= DBGINIT ) printf("Read arguments\n");
	#endif
	strcpy(ip,"mpcd/data/input.json");
	strcpy(op,"./"); // default output directory should be where you are
	for(arg=1; arg<argc; arg++) {
		// Check for a dash
		if(argv[arg][0]=='-') {
			switch (argv[arg][1]) {
				case 'i':
					*inMode = 0;
					arg++;
					sprintf(ip,"%s",argv[arg]);
                    inProvided = 1;
					break;
				case 'L': // legacy
					if (argv[arg][2] == 'i'){ // arg 'Li' for legacy input
						*inMode = 1;
						arg++;
						sprintf(ip,"%s",argv[arg]);
						// Make sure that the directory ends with a "/"
						strln = strlen(ip);
						if( ip[strln-1]!='/' ) strcat( ip,"/" );
                        inProvided = 1;
						break;
					}
				case 'o':
					arg++;
					sprintf(op,"%s",argv[arg]);
					// Make sure that the directory ends with a "/"
					strln = strlen(op);
					if( op[strln-1]!='/' ) strcat( op,"/" );
                    outProvided = 1;
					break;
				case 'v':
					printVersionSummary( );
				case 'h':
					printf("\nMulti-Particle Collision Dynamics\n");
					printf("Shendruk Lab, University of Edinburgh\n");
					printf("Originally by Tyler Shendruk\n");
					printf("for the Polymer Physics Research Group, University of Ottawa\n");
					printf("\nUsage:\n");
					printf("\t-i\t[path to JSON input file]\t\t\tRequired. See `sampleInputs` for examples.\n");
					printf("\t-o\t[path to output file directory]\t\t\tRequired.\n");
					printf("\t-Li\t(legacy) [path to input file directory]\t\tdefault=`mpcd/data/`\n");
					printf("\t-v\t(legacy) print version summary\n");
					printf("\t-h\tprint this help menu\n");
					exit(EXIT_SUCCESS);
					break;
				default:
					printf("Error:\n** Invalid arument: %s\n** Try -h\n",argv[arg]);
					exit(EXIT_FAILURE);
			}
		}
	}

    // do check to see if input and output were provided
    if (inProvided == 0){
        printf("Error: No input provided. See -h for help\n");
        exit(EXIT_FAILURE);
    }
    if (outProvided == 0){
        printf("Error: No output provided. See -h for help\n");
        exit(EXIT_FAILURE);
    }
}

///
/// @brief Checks if a given boundary condition contains the minimum primitives required by the input system.
///
/// The JSON input system allows for the user to ignore stating the value of certain JSON tags. For boundary conditions
/// however, there are some tags that are essential for a boundary to function. This method checks to ensure that these
/// tags are present in a given BC object.
///
/// The list of necessary tags is as follows:
/// - `"aInv"`
/// - `"P"`
/// - `"R"`
/// - `"DN"`
/// - `"MVN"`
/// - `"MVT"`
///
/// @param bc cJSON object containing the boundary condition to be checked.
/// @return 1 if valid, 0 if not.
///
int checkBC(cJSON *bc){
	char tagList[6][5] = {"aInv", "P", "R", "DN", "MVN", "MVT"};

	cJSON *temp = NULL;
	int i;
	for (i = 0; i < 6; i++){
		temp = cJSON_GetObjectItemCaseSensitive(bc, tagList[i]);
		if (temp == NULL) return 0;
	}

	return 1; // if you're here without returning then all succesful
}

///
/// @brief Main method for reading a JSON input file to populate the simulation parameters.
///
/// This (big) method performs all reading for the JSON input system. A summary of the method is as follows:
/// 1. Read in the JSON file and prepare it in cJSON format. Prepare helper linked lists.
/// 2. Perform "overrides", which are hacks to allow for shortcuts using the input system.
/// 3. Read in legacy input.inp parameters.
/// 4. Read in species hyper-parameters to populate `SP` array, and hence populate main particle list.
/// 5. Read in legacy printcom.inp parameters.
/// 6. Read in boundary conditions, corresponding to legacy bc.inp parameters.
/// 7. Read in legacy swimmers.inp parameters, and allocate memory for swimmers.
/// 8. Perform input verification, and clean-up before returning.
///
/// Many methods used in this method are defined in cJson.c.
///
/// @param fpath Path to the json input file.
/// @param in Pointer to the inputList struct to be populated. Expected to be `&in`.
/// @param SP Pointer to the particle species list to be populated. Expected to be `&SP`.
/// @param pSRD Pointer to the particle list. Expected to be `&pSRD`.
/// @param CL Pointer to the cell array. Expected to be `&CL`.
/// @param MDMode Pointer to the MD mode flag integer. Expected to be `&MDMode`.
/// @param out Pointer to the output flag structure. Expected to be `&out`.
/// @param WALL Pointer to the boundary condition list. Expected to be `&WALL`.
/// @param specS Pointer to the swimmer species list to be populated. Expected to be `&specS`.
/// @param sw Pointer to the swimmer list. Expected to be `&sw`.
/// @see cJson.c
///
void readJson( char fpath[], inputList *in, spec **SP, particleMPC **pSRD,
   cell ****CL, int *MDMode, outputFlagsList *out, bc **WALL,
   specSwimmer *specS, swimmer **sw){
	int i, j; // counting variables

	char* fileStr = NULL;
	if(getFileStr(fpath, &fileStr) != 0){ // read, return on error
		exit(EXIT_FAILURE);
	}

	#ifdef DBG // print json to user if in relevent debug mode
		if( DBUG > DBGTITLE ){
			printf("==== Read JSON =====\n%s", fileStr);
			printf("\n\n");
		}
	#endif

	//now can actually parse the json
	cJSON *jObj = cJSON_Parse(fileStr); // create the json object
	if (jObj == NULL) // error checking
	{
		const char *error_ptr = cJSON_GetErrorPtr();
		if (error_ptr != NULL)
		{
			fprintf(stderr, "Json read error. \nError before: %s\n", error_ptr);
		}
		exit(EXIT_FAILURE);
	}

	// set up input validation routines
	linkedList *jsonTagList = NULL;
	initLL(&jsonTagList);
	linkedList *arrayList = NULL;
	initLL(&arrayList);

	////////////////////////////////////////////////////////////////////////////
	// Perform parsing here.
	// This is done in the order declared in docs/InputGuide.md
	////////////////////////////////////////////////////////////////////////////

	// 0. Overrides ////////////////////////////////////////////////////////////
	// General override behaviour (ie, not BC and not species) should be set here

	// flags for overrides
	int domainWalls = 0; // Whether to add domain walls. 0 = off, 1 = PBC, 2 = solid
    float checkPointTimer = 0.0; // hours until MPCD will run a checkpoint

	// Get overrides
	domainWalls = getJObjInt(jObj, "domainWalls", 0, jsonTagList); // domainWalls
    checkPointTimer = getJObjDou(jObj, "checkpointTimerOut", 0.0, jsonTagList); // checkpointTimerOut

	// perform general overrides
	// NOTE: none here yet :)

	// 1. Old input.inp ////////////////////////////////////////////////////////
	// scroll up to void readin() to see better descriptions & definitions for these

	// dimensionality and domain bounds array
	cJSON *arrDomain = NULL;
	getCJsonArray(jObj, &arrDomain, "domain", jsonTagList, arrayList, 0);
	if(arrDomain != NULL) { // if the can be found in the json
		DIM = cJSON_GetArraySize(arrDomain);
		if(DIM > 3){ // check dimensionality is valid
			printf("Error: Dimensionality must be 3 or less.\n");
			exit(EXIT_FAILURE);
		}

		for (i = 0; i < _3D; i++) {
			if (i == 1 && DIM == 1) XYZ[i] = 1; // if 1D, set y to 1
			else if (i == 2 && (DIM == 2 || DIM == 1)) XYZ[i] = 1; // if 1D or 2D, set z to 1
			else XYZ[i] = cJSON_GetArrayItem(arrDomain, i)->valueint; // get the value
		}
	} else { // if array cannot be found then fallback to default
		DIM = 2;
		XYZ[0] = 30;
		XYZ[1] = 30;
		XYZ[2] = 1;
	}
	for(i=0; i<_3D; i++ ) XYZ_P1[i] = XYZ[i]+1; // add 1 to each dimension

	// first set of primitives
	in->KBT = getJObjDou(jObj, "kbt", 1, jsonTagList); // kbt
	in->dt = getJObjDou(jObj, "dt", 0.1, jsonTagList); // dt
	in->simSteps = getJObjInt(jObj, "simSteps", 2000, jsonTagList); // simSteps
	in->warmupSteps = getJObjInt(jObj, "warmUp", 0, jsonTagList); // warmupSteps
	in->RFRAME = getJObjInt(jObj, "rFrame", 1, jsonTagList); // RFRAME
	in->zeroNetMom = getJObjInt(jObj, "zeroNetMom", 0, jsonTagList); // zeroNetMom
	in->GALINV = getJObjInt(jObj, "galInv", 1, jsonTagList); // GALINV
	in->TSTECH = getJObjInt(jObj, "tsTech", 0, jsonTagList); // TSTECH
    const char* collOpTags[2] = {"rTech", "collOp"}; // possible tags for collision operator
	in->RTECH = getJObjIntMultiple(jObj, collOpTags, 2, 2, jsonTagList); // RTECH
	in->LC = getJObjInt(jObj, "lc", 0, jsonTagList); // LC
	in->TAU = getJObjDou(jObj, "tau", 0.5, jsonTagList); // TAU
	in->RA = getJObjDou(jObj, "rotAng", 1.570796, jsonTagList); // rotAng
	in->FRICCO = getJObjDou(jObj, "fricCoef", 1.0, jsonTagList); // fricCo
	in->MFPOT = getJObjDou(jObj, "mfpot", 10.0, jsonTagList); // mfpPot
	in->tolD = getJObjDou(jObj, "tolD", 0.01, jsonTagList); //defect tolerance
	in->noHI = getJObjInt(jObj, "noHI", 0, jsonTagList); // noHI
	in->inCOMP = getJObjInt(jObj, "incomp", 0, jsonTagList); // inCOMP
	in->MULTIPHASE = getJObjInt(jObj, "multiphase", 0, jsonTagList); // multiPhase

	// grav array
	cJSON *arrGrav = NULL;
	getCJsonArray(jObj, &arrGrav, "grav", jsonTagList, arrayList, 0);
	if (arrGrav != NULL) { // if grav has been found then....
		if (cJSON_GetArraySize(arrGrav) != _3D) { // check dimensionality is valid
			printf("Error: Grav must be a 3D array.\n");
			exit(EXIT_FAILURE);
		}

		for (i = 0; i < _3D; i++) { // get the value
			in->GRAV[i] = cJSON_GetArrayItem(arrGrav, i)->valuedouble;
		}
	} else { // if no grav specified then fallback
		in->GRAV[0] = 0;
		in->GRAV[1] = 0;
		in->GRAV[2] = 0;
	}

	// mag array
	cJSON *arrMag = NULL;
	getCJsonArray(jObj, &arrMag, "mag", jsonTagList, arrayList, 0);
	if (arrMag != NULL) { // if mag has been found then....
		if (cJSON_GetArraySize(arrMag) != _3D) { // check dimensionality is valid
			printf("Error: Mag must be a 3D array.\n");
			exit(EXIT_FAILURE);
		}

		for (i = 0; i < _3D; i++) { // get the value
			in->MAG[i] = cJSON_GetArrayItem(arrMag, i)->valuedouble;
		}
	} else { // if no grav specified then fallback
		in->MAG[0] = 0;
		in->MAG[1] = 0;
		in->MAG[2] = 0;
	}

	// second set of primitives
	in->seed = getJObjInt(jObj, "seed", 0, jsonTagList); // seed

	// Handle MD
	getJObjStr(jObj, "mdIn", "", &mdInputFile, jsonTagList);
    int mdCoupleMode = getJObjInt(jObj, "mdCoupleMode", 1, jsonTagList); // coupling mode for MD
	if (strcmp(mdInputFile, "") == 0){ // if no input file was found
		MDmode = 0;
	} else MDmode = mdCoupleMode; // otherwise set MD to correct coupling mode to enable it

	in->stepsMD = getJObjInt(jObj, "stepsMD", 20, jsonTagList); // stepsMD
    in->MFPLAYERH = getJObjInt(jObj, "mfpLayerH", 0, jsonTagList); // mfpLayerH

	// 2. Species //////////////////////////////////////////////////////////////
	// scroll up to void readin() to see better descriptions & definitions for these

	cJSON *arrSpecies = NULL;
	getCJsonArray(jObj, &arrSpecies, "species", jsonTagList, arrayList, 1);
	if(arrSpecies != NULL){ // if this can be found in the json
		NSPECI = cJSON_GetArraySize(arrSpecies); // get the number of species

		//Allocate the needed amount of memory for the species SP
		(*SP) = (spec*) calloc( NSPECI, sizeof( spec ) );
		for (i = 0; i < NSPECI; i++) { // loop through the species
			cJSON *objElem = cJSON_GetArrayItem(arrSpecies, i); // get the species object

			// now get first set of primitives
			(*SP+i)->MASS = getJObjDou(objElem, "mass", 1.0, jsonTagList); // mass

			// handle population related overrides
			double cellDens = getJObjDou(objElem, "dens", -1, jsonTagList);
			if (cellDens < 0){ // if cellDens is invalid
				(*SP+i)->POP = getJObjInt(objElem, "pop", XYZ[0]*XYZ[1]*XYZ[2]*20, jsonTagList); // pop
			} else { // otherwise, set population using per cell density
				(*SP+i)->POP = XYZ[0]*XYZ[1]*XYZ[2]*cellDens;
			}

			(*SP+i)->QDIST = getJObjInt(objElem, "qDist", 0, jsonTagList); // qDist
			(*SP+i)->VDIST = getJObjInt(objElem, "vDist", 0, jsonTagList); // vDist
			(*SP+i)->ODIST = getJObjInt(objElem, "oDist", 2, jsonTagList); // oDist

			//Read the binary fluid interaction matrix for this species with all other species
			cJSON *arrBFM = NULL;
			getCJsonArray(jObj, &arrBFM, "interMatr", jsonTagList, arrayList, 0);
			if (arrBFM != NULL) { // if grav has been found then....
				if (cJSON_GetArraySize(arrBFM) != NSPECI) { // check dimensionality is valid
					printf("Error: Interaction matrices must have columns of length equal to the number of species.\n");
					exit(EXIT_FAILURE);
				}

				for (j = 0; j < NSPECI; j++) { // get the value
					(*SP+i)->M[j] = cJSON_GetArrayItem(arrBFM, j)->valuedouble;
				}
			} else {
				for (j = 0; j < NSPECI; j++) { // get the value
					(*SP+i)->M[j] = 0;
				}
			}

			// get second set of primitives
			(*SP+i)->RFC = getJObjDou(objElem, "rfc", 0.01, jsonTagList); // rfCoef
			(*SP+i)->LEN = getJObjDou(objElem, "len", 0.007, jsonTagList); // len
			(*SP+i)->TUMBLE = getJObjDou(objElem, "tumble", 2.0, jsonTagList); // tumble
			(*SP+i)->CHIHI = getJObjDou(objElem, "shearSusc", 0.5, jsonTagList); // chiHi
			(*SP+i)->CHIA = getJObjDou(objElem, "magnSusc", 0.001, jsonTagList); // chiA
			(*SP+i)->ACT = getJObjDou(objElem, "act", 0.05, jsonTagList); // act
			(*SP+i)->SIGWIDTH = getJObjDou(objElem, "sigWidth", 1, jsonTagList); // sigWidth
			// error check, is SIGWIDTH 0?
			if ((*SP+i)->SIGWIDTH == 0) {
				printf("Error: SIGWIDTH cannot be 0.\n");
				exit(EXIT_FAILURE);
			}
			(*SP+i)->SIGPOS = getJObjDou(objElem, "sigPos", (*SP+i)->SIGWIDTH, jsonTagList); // sigPos
			(*SP+i)->MINACTRATIO = getJObjDou(objElem, "minActRatio", 0.0, jsonTagList); // minActRatio
			(*SP+i)->DAMP = getJObjDou(objElem, "damp", 0.0, jsonTagList); // damp
		}
	} else { // if nothing found in the JSON then fallback to the default
		// setting up a single species with default parameters
		//		note this is just copied from the above code w lines changed
		NSPECI = 1;

		(*SP) = (spec*) calloc( NSPECI, sizeof( spec ) );
		for (i = 0; i < NSPECI; i++) { // loop through the species
			// now get first set of primitives
			(*SP+i)->MASS = 1; // mass
			(*SP+i)->POP = XYZ[0]*XYZ[1]*XYZ[2]*20; // pop
			(*SP+i)->QDIST = 0; // qDist
			(*SP+i)->VDIST = 0; // vDist
			(*SP+i)->ODIST = 2; // oDist
			for (j = 0; j < NSPECI; j++) { // interaction matrix
				(*SP+i)->M[j] = 0;
			}
			(*SP+i)->RFC = 0.01; // rfCoef
			(*SP+i)->LEN = 0.007; // len
			(*SP+i)->TUMBLE = 2; // tumble
			(*SP+i)->CHIHI = 0.5; // chiHi
			(*SP+i)->CHIA = 0.001; // chiA
			(*SP+i)->ACT = 0.05; // act
			(*SP+i)->SIGWIDTH = 1; // sigwidth
			(*SP+i)->SIGPOS = 1; // sigpos
			(*SP+i)->MINACTRATIO = 0.0; // minActRatio
			(*SP+i)->DAMP = 0; // damp
		}
	}

	//Total Number of particleMPCs
	GPOP = 0;
	for( i=0; i<NSPECI; i++ ) GPOP += (*SP+i)->POP;
	(*pSRD) = (particleMPC*) calloc( GPOP, sizeof( particleMPC ) );

	//Allocate memory for the cells
	//Allocate rows (x first)
//	*CL = (cell***) malloc( XYZ_P1[0] * sizeof( cell** ) );
	*CL = calloc( XYZ_P1[0], sizeof( cell** ) );
	//For each x-element allocate the y columns
	for( i=0; i<XYZ_P1[0]; i++ ) {
//		(*CL)[i] = (cell**) malloc( XYZ_P1[1] * sizeof( cell* ) );
		(*CL)[i] = calloc( XYZ_P1[1], sizeof( cell* ) );
		//For each y-element allocate the z columns
		for( j=0; j<XYZ_P1[1]; j++ ) {
//			(*CL)[i][j] = (cell*) malloc( XYZ_P1[2] * sizeof( cell ) );
			(*CL)[i][j] = calloc( XYZ_P1[2], sizeof( cell ) );
		}
	}

	// 3. Printcom /////////////////////////////////////////////////////////////
	// scroll up to void readpc() to see better descriptions & definitions for these

	DBUG = getJObjInt(jObj, "debugOut", 3, jsonTagList); // dbug
	out->TRAJOUT = getJObjInt(jObj, "trajOut", 0, jsonTagList); // trajOut
	out->printSP = getJObjInt(jObj, "trajSpecOut", 0, jsonTagList); // printSP
	out->COAROUT = getJObjInt(jObj, "coarseOut", 0, jsonTagList); // coarOut
	out->FLOWOUT = getJObjInt(jObj, "flowOut", 0, jsonTagList); // flowOut
	out->VELOUT = getJObjInt(jObj, "velOut", 0, jsonTagList); // velOut
	out->AVVELOUT = getJObjInt(jObj, "avVelOut", 0, jsonTagList); // avVelOut
	out->ORDEROUT = getJObjInt(jObj, "dirSOut", 0, jsonTagList); // orderOut
	out->QTENSOUT = getJObjInt(jObj, "qTensOut", 0, jsonTagList); // qTensOut
	out->QKOUT = getJObjInt(jObj, "qkTensOut", 0, jsonTagList); // qKOut
	out->ENFIELDOUT = getJObjInt(jObj, "oriEnOut", 0, jsonTagList); // enFieldOut
	out->SPOUT = getJObjInt(jObj, "colourOut", 0, jsonTagList); // spOut
	out->PRESOUT = getJObjInt(jObj, "pressureOut", 0, jsonTagList); // presOut
	out->ENNEIGHBOURS = getJObjInt(jObj, "neighbourEnOut", 0, jsonTagList); // enNeighbours
	out->AVSOUT = getJObjInt(jObj, "avSOut", 0, jsonTagList); // avSOut
	out->DENSOUT = getJObjInt(jObj, "densSDOut", 0, jsonTagList); // densOut
	out->ENSTROPHYOUT = getJObjInt(jObj, "enstrophyOut", 0, jsonTagList); // enStrophOut
	out->HISTVELOUT = getJObjInt(jObj, "histVelOut", 0, jsonTagList); // histVelOut
	out->HISTSPEEDOUT = getJObjInt(jObj, "histSpeedOut", 0, jsonTagList); // histSpeedOut
	out->HISTVORTOUT = getJObjInt(jObj, "histVortOut", 0, jsonTagList); // histVortOut
	out->HISTENSTROUT = getJObjInt(jObj, "histEnsOut", 0, jsonTagList); // histEnstrophyOut
	out->HISTDIROUT = getJObjInt(jObj, "histDirOut", 0, jsonTagList); // histDirOut
	out->HISTSOUT = getJObjInt(jObj, "histSOut", 0, jsonTagList); // histSOut
	out->HISTNOUT = getJObjInt(jObj, "histNOut", 0, jsonTagList); // histNOut
	out->SOLOUT = getJObjInt(jObj, "solidTrajOut", 0, jsonTagList); // solOut
	out->TOPOOUT = getJObjInt(jObj, "topoFieldOut", 0, jsonTagList); // topoOut
	out->DEFECTOUT = getJObjInt(jObj, "defectsOut", 0, jsonTagList); // defectOut
	out->DISCLINOUT = getJObjInt(jObj, "disclinOut", 0, jsonTagList); // disclinationTensorFieldOut
	out->ENOUT = getJObjInt(jObj, "energyOut", 0, jsonTagList); // enOut
	out->CVVOUT = getJObjInt(jObj, "velCorrOut", 0, jsonTagList); // cvvOut
	out->CNNOUT = getJObjInt(jObj, "dirCorrOut", 0, jsonTagList); // cnnOut
	out->CWWOUT = getJObjInt(jObj, "vortCorrOut", 0, jsonTagList); // cwwOut
	out->CDDOUT = getJObjInt(jObj, "densCorrOut", 0, jsonTagList); // cddOut
	out->CSSOUT = getJObjInt(jObj, "orderCorrOut", 0, jsonTagList); // cssOut
	out->CPPOUT = getJObjInt(jObj, "phaseCorrOut", 0, jsonTagList); // cppOut
	out->ENERGYSPECTOUT = getJObjInt(jObj, "energySpecOut", 0, jsonTagList); // energySpectOut
	out->ENSTROPHYSPECTOUT = getJObjInt(jObj, "enstrophySpecOut", 0, jsonTagList); // enstrophySpectOut
	out->BINDER = getJObjInt(jObj, "binderOut", 0, jsonTagList); // binderOut
	out->BINDERBIN = getJObjInt(jObj, "binderBin", 0, jsonTagList); // binderBinOut
	out->SWOUT = getJObjInt(jObj, "swimQOut", 0, jsonTagList); // swOut
	out->SWORIOUT = getJObjInt(jObj, "swimOOut", 0, jsonTagList); // swOriOut
    const char* swimROutTags[2] = {"swimROut", "swimRTOut"}; // possible tags for collision operator
    out->RTOUT = getJObjIntMultiple(jObj, swimROutTags, 2, 0, jsonTagList); // RTECH
	out->SYNOUT = getJObjInt(jObj, "synopsisOut", 1, jsonTagList); // swSynOut
	out->CHCKPNT = getJObjInt(jObj, "checkpointOut", 0, jsonTagList); // chkpntOut

	// 3. Boundaries ///////////////////////////////////////////////////////////
	// scroll up to void bcin() to see better descriptions & definitions for these

	cJSON *arrBC = NULL;
	getCJsonArray(jObj, &arrBC, "BC", jsonTagList, arrayList, 1);
	if(arrBC != NULL){ // if this can be found in the json
		NBC = cJSON_GetArraySize(arrBC); // get the number of BCs

		//Allocate the needed amount of memory for the BCs
		(*WALL) = (bc*) calloc( NBC, sizeof( bc ) );
		for (i = 0; i < NBC; i++) { // loop through the BCs
			cJSON *objElem = cJSON_GetArrayItem(arrBC, i); // get the BC object
			bc *currWall = (*WALL + i); // get the pointer to the BC we want to write to

			// Do a check if the necessary parameters are present
			if (!checkBC(objElem)){
				printf("Error: BC %d is missing essential parameters\n", i);
				exit(EXIT_FAILURE);
			}

			// parse through first set of primitives
			currWall->COLL_TYPE = getJObjInt(objElem, "colType", 1, jsonTagList); // collType
			currWall->PHANTOM = getJObjInt(objElem, "phantom", 0, jsonTagList); // phantom
			currWall->E = getJObjDou(objElem, "E", -1.0, jsonTagList); // E

			// Q array
			cJSON *arrQ = NULL;
			getCJsonArray(objElem, &arrQ, "Q", jsonTagList, arrayList, 0);
			if (arrQ != NULL) { // if grav has been found then....
				if (cJSON_GetArraySize(arrQ) != _3D) { // check dimensionality is valid
					printf("Error: Q must be 3D.\n");
					exit(EXIT_FAILURE);
				}

				for (j = 0; j < _3D; j++) { // get the value
					currWall->Q[j] = cJSON_GetArrayItem(arrQ, j)->valuedouble;
				}
			} else {
				for (j = 0; j < _3D; j++) { // get the value
					currWall->Q[j] = 0;
				}
			}

			// V array
			cJSON *arrV = NULL;
			getCJsonArray(objElem, &arrV, "V", jsonTagList, arrayList, 0);
			if (arrV != NULL) { // if grav has been found then....
				if (cJSON_GetArraySize(arrV) != _3D) { // check dimensionality is valid
					printf("Error: V must be 3D.\n");
					exit(EXIT_FAILURE);
				}

				for (j = 0; j < _3D; j++) { // get the value
					currWall->V[j] = cJSON_GetArrayItem(arrV, j)->valuedouble;
				}
			} else {
				for (j = 0; j < _3D; j++) { // get the value
					currWall->V[j] = 0;
				}
			}

			// O array
			cJSON *arrO = NULL;
			getCJsonArray(objElem, &arrO, "O", jsonTagList, arrayList, 0);
			if (arrO != NULL) { // if grav has been found then....
				if (cJSON_GetArraySize(arrO) != _3D) { // check dimensionality is valid
					printf("Error: O must be 3D.\n");
					exit(EXIT_FAILURE);
				}

				for (j = 0; j < _3D; j++) { // get the value
					currWall->O[j] = cJSON_GetArrayItem(arrO, j)->valuedouble;
				}
			} else {
				for (j = 0; j < _3D; j++) { // get the value
					currWall->O[j] = 0;
				}
			}

			// L array
			cJSON *arrL = NULL;
			getCJsonArray(objElem, &arrL, "L", jsonTagList, arrayList, 0);
			if (arrL != NULL) { // if grav has been found then....
				if (cJSON_GetArraySize(arrL) != _3D) { // check dimensionality is valid
					printf("Error: L must be 3D.\n");
					exit(EXIT_FAILURE);
				}

				for (j = 0; j < _3D; j++) { // get the value
					currWall->L[j] = cJSON_GetArrayItem(arrL, j)->valuedouble;
				}
			} else {
				for (j = 0; j < _3D; j++) { // get the value
					currWall->L[j] = 0;
				}
			}

			// G array
			cJSON *arrG = NULL;
			getCJsonArray(objElem, &arrG, "G", jsonTagList, arrayList, 0);
			if (arrG != NULL) { // if grav has been found then....
				if (cJSON_GetArraySize(arrG) != _3D) { // check dimensionality is valid
					printf("Error: G must be 3D.\n");
					exit(EXIT_FAILURE);
				}

				for (j = 0; j < _3D; j++) { // get the value
					currWall->G[j] = cJSON_GetArrayItem(arrG, j)->valuedouble;
				}
			} else {
				for (j = 0; j < _3D; j++) { // get the value
					currWall->G[j] = 0;
				}
			}

			// aInv array - NECESSARY
			cJSON *arrAInv = NULL;
			getCJsonArray(objElem, &arrAInv, "aInv", jsonTagList, arrayList, 0);
			if (arrAInv != NULL) { // if grav has been found then....
				if (cJSON_GetArraySize(arrAInv) != _3D) { // check dimensionality is valid
					printf("Error: aInv must be 3D.\n");
					exit(EXIT_FAILURE);
				}

				for (j = 0; j < _3D; j++) { // get the value
					currWall->AINV[j] = cJSON_GetArrayItem(arrAInv, j)->valuedouble;
					if( fneq(currWall->AINV[j],0.0) ) currWall->A[j] = 1.0/currWall->AINV[j];
					else currWall->A[j] = 0.0;
				}
			} else {
				printf("Error: Could not find aInv in BC %d.\n", i);
				exit(EXIT_FAILURE);
			}

			// rotsymm array
			cJSON *arrRotSym = NULL;
			getCJsonArray(objElem, &arrRotSym, "rotSym", jsonTagList, arrayList, 0);
			if (arrRotSym != NULL) { // if grav has been found then....
				if (cJSON_GetArraySize(arrRotSym) != 2) { // check dimensionality is valid
					printf("Error: rotSym must have two components.\n");
					exit(EXIT_FAILURE);
				}

				for (j = 0; j < 2; j++) { // get the value
					currWall->ROTSYMM[j] = cJSON_GetArrayItem(arrRotSym, j)->valuedouble;
				}
			} else {
				for (j = 0; j < 2; j++) { // get the value
					currWall->ROTSYMM[j] = 4;
				}
			}

			currWall->ABS = getJObjInt(objElem, "abs", 0, jsonTagList); // abs

			// P array - NECESSARY
			cJSON *arrP = NULL;
			getCJsonArray(objElem, &arrP, "P", jsonTagList, arrayList, 0);
			if (arrP != NULL) { // if grav has been found then....
				if (cJSON_GetArraySize(arrP) != 4) { // check dimensionality is valid
					printf("Error: P must have 4 elements.\n");
					exit(EXIT_FAILURE);
				}

				for (j = 0; j < 4; j++) { // get the value
					currWall->P[j] = cJSON_GetArrayItem(arrP, j)->valuedouble;
				}
			} else {
				printf("Error: Could not find P in BC %d.\n", i);
				exit(EXIT_FAILURE);
			}

			// some more primitives
			currWall->R = getJObjDou(objElem, "R", 2, jsonTagList); // r - NECESSARY
			currWall->DN = getJObjDou(objElem, "DN", 1, jsonTagList); // dn - NECESSARY
			currWall->DT = getJObjDou(objElem, "DT", 0, jsonTagList); // dt
			currWall->DVN = getJObjDou(objElem, "DVN", 0, jsonTagList); // dvn
			currWall->DVT = getJObjDou(objElem, "DVT", 0, jsonTagList); // dvt

			// DVxyz array
			cJSON *arrDVxyz = NULL;
			getCJsonArray(objElem, &arrDVxyz, "DVxyz", jsonTagList, arrayList, 0);
			if (arrDVxyz != NULL) { // if grav has been found then....
				if (cJSON_GetArraySize(arrDVxyz) != _3D) { // check dimensionality is valid
					printf("Error: DVxyz must have two components.\n");
					exit(EXIT_FAILURE);
				}

				for (j = 0; j < _3D; j++) { // get the value
					currWall->DVxyz[j] = cJSON_GetArrayItem(arrDVxyz, j)->valuedouble;
				}
			} else {
				for (j = 0; j < _3D; j++) { // get the value
					currWall->DVxyz[j] = 0;
				}
			}

			currWall->MVN = getJObjDou(objElem, "MVN", 1, jsonTagList); // mvn - NECESSARY
			currWall->MVT = getJObjDou(objElem, "MVT", 1, jsonTagList); // mvt - NECESSARY
			currWall->MUN = getJObjDou(objElem, "MUN", 1, jsonTagList); // mun
			currWall->MUT = getJObjDou(objElem, "MUT", 1, jsonTagList); // mut

			// MUxyz array
			cJSON *arrMUxyz = NULL;
			getCJsonArray(objElem, &arrMUxyz, "MUxyz", jsonTagList, arrayList, 0);
			if (arrMUxyz != NULL) { // if grav has been found then....
				if (cJSON_GetArraySize(arrMUxyz) != _3D) { // check dimensionality is valid
					printf("Error: MUxyz must have two components.\n");
					exit(EXIT_FAILURE);
				}

				for (j = 0; j < _3D; j++) { // get the value
					currWall->MUxyz[j] = cJSON_GetArrayItem(arrMUxyz, j)->valuedouble;
				}
			} else {
				for (j = 0; j < _3D; j++) { // get the value
					currWall->MUxyz[j] = 1;
				}
			}

			// DUxyz array
			cJSON *arrDUxyz = NULL;
			getCJsonArray(objElem, &arrDUxyz, "DUxyz", jsonTagList, arrayList, 0);
			if (arrDUxyz != NULL) { // if grav has been found then....
				if (cJSON_GetArraySize(arrDUxyz) != _3D) { // check dimensionality is valid
					printf("Error: DUxyz must have two components.\n");
					exit(EXIT_FAILURE);
				}

				for (j = 0; j < _3D; j++) { // get the value
					currWall->DUxyz[j] = cJSON_GetArrayItem(arrDUxyz, j)->valuedouble;
				}
			} else {
				for (j = 0; j < _3D; j++) { // get the value
					currWall->DUxyz[j] = 0;
				}
			}

			currWall->KBT = getJObjDou(objElem, "kbt", 1, jsonTagList); // kbt
			currWall->DSPLC = getJObjInt(objElem, "dsplc", 0, jsonTagList); // dspc
			currWall->INV = getJObjInt(objElem, "inv", 0, jsonTagList); // inv
			currWall->MASS = getJObjDou(objElem, "mass", 1, jsonTagList); // mass

			// B array
			cJSON *arrB = NULL;
			getCJsonArray(objElem, &arrB, "wavy", jsonTagList, arrayList, 0);
			if (arrB != NULL) { // if wavewall parameters have been found then....
				if (cJSON_GetArraySize(arrB) != _3D) { // check dimensionality is valid
					printf("Error: B must be 3D.\n");
					exit(EXIT_FAILURE);
				}

				for (j = 0; j < _3D; j++) { // get the value
					currWall->B[j] = cJSON_GetArrayItem(arrB, j)->valuedouble;
				}
			} else for (j = 0; j < _3D; j++) { // get the value
					currWall->B[j] = 0.0;
			}

			// Handle BC overrides /////////////////////////////////////////////
			// anchoring
			if (getJObjInt(objElem, "homeotropic", 0, jsonTagList) == 1) { // homeotropic anchoring
				currWall->MUN = 1;
				currWall->MUT = 0;
			} else if (getJObjInt(objElem, "planar", 0, jsonTagList) == 1) { // planar anchoring
				currWall->MUN = 0;
				currWall->MUT = 1;
			}

			// Set the planar flag
			currWall->PLANAR = 0;							// sets default to non-planar
			if ( feq(currWall->B[0],0.0) ) {	// if not wavy, check if planar
				if( feq(currWall->P[0],1.0) && feq(currWall->P[1],1.0) && feq(currWall->P[2],1.0) ) {
					// Left or right wall
					if( fneq(currWall->A[0],0.0) && feq(currWall->A[1],0.0) && feq(currWall->A[2],0.0) ) currWall->PLANAR = 1;
					// Top or bottom wall
					else if( feq(currWall->A[0],0.0) && fneq(currWall->A[1],0.0) && feq(currWall->A[2],0.0) ) currWall->PLANAR = 1;
					// Far or near wall
					else if( feq(currWall->A[0],0.0) && feq(currWall->A[1],0.0) && fneq(currWall->A[2],0.0) ) currWall->PLANAR = 1;
				}
			}
		}
	} else { // otherwise default to periodic BCs about domain
		domainWalls = 1; // just trigger the domain walls override w PBCs
		NBC = 0;
	}

	// handle domainWalls override
	if (domainWalls == 1 || domainWalls == 2) {
		int oldBCNo = NBC; // number of BCs NOT including override created ones

		// realloc memory to store extra BCs
		NBC += DIM * 2; // creating 2 extra BCs for each dimension
		if (oldBCNo > 0) { // if wall already exists then realloc
			(*WALL) = (bc*) realloc(*WALL, NBC * sizeof(bc)); // realloc mem
		} else { // otherwise need to malloc
//			(*WALL) = (bc*) malloc(NBC * sizeof(bc)); // malloc mem
			(*WALL) = calloc(NBC, sizeof(bc)); // malloc mem
		}

		//set up PBCs on the xy plane based on the domain dimensions
		for (i = 0; i < 2 * DIM; i++) { // use i as the fundamental counter for setting these up
			bc *currWall = (*WALL + i + oldBCNo); // get the pointer to the BC we want to write to

			// handle the things that change for each wall first
			for (j = 0; j < _3D; j++) { // set A, Ainv to default vals
				currWall->AINV[j] = 0;
				currWall->A[j] = 0.0;
			}
			switch (i) { // set up walls based on index
				case 0: // left wall
					currWall->AINV[0] = 1; // aInv array - NECESSARY
					currWall->A[0] = 1.0/currWall->AINV[0];
					currWall->R = 0; // r
					currWall->DN = XYZ[0]; // dn
					break;

				case 1: // right wall
					currWall->AINV[0] = -1; // aInv array - NECESSARY
					currWall->A[0] = 1.0/currWall->AINV[0];
					currWall->R = -XYZ[0]; // r
					currWall->DN = XYZ[0]; // dn
					break;

				case 2: // bottom wall
					currWall->AINV[1] = 1; // aInv array - NECESSARY
					currWall->A[1] = 1.0/currWall->AINV[1];
					currWall->R = 0; // r
					currWall->DN = XYZ[1]; // dn
					break;

				case 3: // top wall
					currWall->AINV[1] = -1; // aInv array - NECESSARY
					currWall->A[1] = 1.0/currWall->AINV[1];
					currWall->R = -XYZ[1]; // r
					currWall->DN = XYZ[1]; // dn
					break;

				case 4: // near wall
					currWall->AINV[2] = 1; // aInv array - NECESSARY
					currWall->A[2] = 1.0/currWall->AINV[2];
					currWall->R = 0; // r
					currWall->DN = XYZ[2]; // dn
					break;

				case 5: // far wall
					currWall->AINV[2] = -1; // aInv array - NECESSARY
					currWall->A[2] = 1.0/currWall->AINV[2];
					currWall->R = -XYZ[2]; // r
					currWall->DN = XYZ[2]; // dn
					break;
			}

			// check if we're using PBCs or solid BCs
			if (domainWalls == 1) { // if PBCs, set appropriate flags
				currWall->PHANTOM = 0;
				currWall->MVN = 1;
				currWall->MVT = 1;
			} else { // otherwise, set flags for solid walls
				currWall->PHANTOM = 1;
				currWall->DN = 0; // override the value of DN from earlier, needs to be 0 for solid
				currWall->MVN = -1.0;
				currWall->MVT = -1.0;
			}

			// set all the default values
			currWall->COLL_TYPE = 1; // collType
			currWall->E = -1.0; // E
			for (j = 0; j < _3D; j++) { // Q array
				currWall->Q[j] = 0;
			}
			for (j = 0; j < _3D; j++) { // V array
				currWall->V[j] = 0;
			}
			for (j = 0; j < _3D; j++) { // O array
				currWall->O[j] = 0;
			}
			for (j = 0; j < _3D; j++) { // L array
				currWall->L[j] = 0;
			}
			for (j = 0; j < _3D; j++) { // G array
				currWall->G[j] = 0;
			}
			for (j = 0; j < 2; j++) { // rotsymm array
				currWall->ROTSYMM[j] = 4;
			}
			currWall->ABS = 0; // abs
			for (j = 0; j < 4; j++) { // P array - NECESSARY
				currWall->P[j] = 1;
			}
			currWall->DT = 0; // dt
			currWall->DVN = 0; // dvn
			currWall->DVT = 0; // dvt
			for (j = 0; j < _3D; j++) { // DVxyz array
				currWall->DVxyz[j] = 0;
			}
			currWall->MUN = 1; // mun
			currWall->MUT = 1; // mut
			for (j = 0; j < _3D; j++) { // MUxyz array
				currWall->MUxyz[j] = 1;
			}
			for (j = 0; j < _3D; j++) { // DUxyz array
				currWall->DUxyz[j] = 0;
			}
			currWall->KBT = 1; // kbt
			currWall->DSPLC = 0; // dspc
			currWall->INV = 0; // inv
			currWall->MASS = 1; // mass

			// Set the planar flag
			currWall->PLANAR = 0;							// sets default to non-planar
			if ( feq(currWall->B[0],0.0) ) {	// if not wavy, check if planar
				if( feq(currWall->P[0],1.0) && feq(currWall->P[1],1.0) && feq(currWall->P[2],1.0) ) {
					// Left or right wall
					if( fneq(currWall->A[0],0.0) && feq(currWall->A[1],0.0) && feq(currWall->A[2],0.0) ) currWall->PLANAR = 1;
					// Top or bottom wall
					else if( feq(currWall->A[0],0.0) && fneq(currWall->A[1],0.0) && feq(currWall->A[2],0.0) ) currWall->PLANAR = 1;
					// Far or near wall
					else if( feq(currWall->A[0],0.0) && feq(currWall->A[1],0.0) && fneq(currWall->A[2],0.0) ) currWall->PLANAR = 1;
				}
			}
		}
	}
	//Determine if any BCs are periodic boundaries
	for( i=0; i<_3D; i++ ) XYZPBC[i]=0;
	for( i=0; i<NBC; i++ ) setPBC( (*WALL+i) );

    // handle checkpoint timer override
    if (checkPointTimer != 0.0) {
        out->CHCKPNT = 1; // just set this as a flag to enable behaviour
        out->CHCKPNTTIMER = checkPointTimer;
    }

	// 4. Swimmers /////////////////////////////////////////////////////////////
	// look at void readswimmers() in swimmers.c to see better descriptions & definitions for these

	specS->TYPE = getJObjInt(jObj, "typeSwim", 2, jsonTagList); // type
	NS = getJObjInt(jObj, "nSwim", 0, jsonTagList); // ns
	specS->QDIST = getJObjInt(jObj, "qDistSwim", 0, jsonTagList); // qDisp
	specS->ODIST = getJObjInt(jObj, "oDistSwim", 0, jsonTagList); // oDisp
	specS->headM = getJObjInt(jObj, "headMSwim", 20, jsonTagList); // headMass
	specS->middM = getJObjInt(jObj, "midMSwim", 20, jsonTagList); // middleMass
	specS->HSPid = getJObjInt(jObj, "hspIdSwim", 1, jsonTagList); // hspId
	specS->MSPid = getJObjInt(jObj, "mspIdSwim", 1, jsonTagList); // mspId
	specS->FS = getJObjDou(jObj, "fsSwim", 20, jsonTagList); // fs
	specS->DS = getJObjDou(jObj, "dsSwim", 1, jsonTagList); // ds
	specS->TS = getJObjDou(jObj, "tsSwim", 0, jsonTagList); // ts
	specS->sizeShrink = getJObjDou(jObj, "sizeShrinkSwim", 0.1, jsonTagList); // sizeShrink
	specS->springShrink = getJObjDou(jObj, "springShrinkSwim", 0.1, jsonTagList); // springShrink
	specS->k = getJObjDou(jObj, "kSwim", 30, jsonTagList); // k
	specS->ro = getJObjDou(jObj, "roSwim", 4, jsonTagList); // ro
	specS->sig = getJObjDou(jObj, "sigSwim", 4, jsonTagList); // sig
	specS->eps = getJObjDou(jObj, "epsSwim", 1, jsonTagList); // eps
	specS->runTime = getJObjDou(jObj, "runTSwim", 0, jsonTagList); // runTime
	specS->tumbleTime = getJObjDou(jObj, "tumTSwim", 0, jsonTagList); // tumbleTime
	specS->shrinkTime = getJObjDou(jObj, "shrTSwim", 2, jsonTagList); // shrinkTime
	specS->MAGMOM = getJObjDou(jObj, "magMomSwim", 1, jsonTagList); // magMom
	specS->fixDist = getJObjDou(jObj, "fixDistSwim", 0, jsonTagList); // fixDist

	//Allocate the memory for the swimmers
	(*sw) = (swimmer*) calloc( NS, sizeof( swimmer ) );

	////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////

	// input verification step
   	verifyJson(jObj, jsonTagList, arrayList);

	// clear memory
	free(fileStr);
   	cJSON_Delete(jObj); // free the json object
	freeLL(jsonTagList); // free the linked lists
	freeLL(arrayList);

	return;
}
