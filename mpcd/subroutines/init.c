///
/// @file
///
/// @brief Set of tools to initialize the system
///
/// Set of tools to initialize the system
///

# include<math.h>
# include<time.h>
# include<string.h>
# include<stdlib.h>

# include "../headers/definitions.h"
# include "../headers/SRDclss.h"
# include "../headers/globals.h"
# include "../headers/rand.h"
# include "../headers/pout.h"
# include "../headers/mtools.h"
# include "../headers/ctools.h"
# include "../headers/pout.h"
# include "../headers/bc.h"
# include "../headers/mpc.h"
# include "../headers/therm.h"
# include "../headers/lc.h"
# include "../headers/swimmers.h"

# include "../../md/mdtypes.h"
# include "../../md/mdsrd.h"
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* *************** OPEN FILES *************** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */

///
/// @brief Function to open a file for reading and writing.
///
/// This function is the basic tool to open files for reading and writing.
/// If the file exists, its content is deleted.
/// It is employed by many other methods, for example opencoarse().
///
/// @param fout Return pointer to the file being opened.
/// @param dir Path to the directory of the file being opened.
/// @param filestring Name of the file being opened without its extension.
/// @param fileextension Extension of the file being opened.
/// @see opencoarse()
///
void openBasic( FILE **fout,char dir[],char filestring[],char fileextension[] ) {
	char filename[200];

	strcpy( filename,dir );
	strcat( filename,filestring );
	strcat( filename,fileextension );
	*fout = fopen(filename,"w+" );
	if( !*fout ){ // file couldn't be opened
		printf( "Error: File '%s'could not be opened.\n",filename );
		exit( 1 );
	}
	outheader( *fout,0 );
}

///
/// @brief Function to open the checkpoint file for writing.
///
/// This function opens the checkpoint.dat file for writing.
///
/// @param fout Return pointer to the file being opened.
/// @param dir Path to the directory of the checkpoint.dat file.
///
void openCheckpoint( FILE **fout,char dir[] ) {
	char filename[200];
	char filechckpoint[]="checkpoint";
	char fileextension[]=".dat";

	strcpy( filename,dir );
	strcat( filename,filechckpoint );
	strcat( filename,fileextension );
	*fout = fopen(filename,"w" );
	if( !*fout ){ // file couldn't be opened
		printf( "Error: File '%s'could not be opened.\n",filename );
		exit( 1 );
	}
}

///
/// @brief Function that opens the output file for the i-th species for reading and writing.
///
/// This functions opens up the output file for the i-th species for reading and writing.
/// In addition, it sets up the header for the file and formarts it.
///
/// @param i Index specifying the species associated to the file being opened.
/// @param fdetail Array of return pointers to the list of files associated to all species.
/// @param dir Path to the directory of the checkpoint.dat file.
/// @param fileprefix Name of the file being opened.
/// @param filesuffix Suffix specifying the species associated to the file being opened, updated to "i" within the function.
/// @param fileextension Extension of the file being opened.
///
void opendetails( int i,FILE *fdetail[],char dir[],char fileprefix[],char filesuffix[],char fileextension[] ) {
	char filename[200];
	strcpy( filename,dir );
	strcat( filename,fileprefix );
	sprintf( filesuffix,"%i",i );
	strcat( filename,filesuffix );
	strcat( filename,fileextension );
	fdetail[i] = fopen( filename,"w+" );
	if( !fdetail[i] ){ // file couldn't be opened
		printf( "Error: File '%s'could not be opened.\n",filename );
		exit( 1 );
	}
	outheader( fdetail[i],i );
	fprintf( fdetail[i],"SPECIES: %i\n",i );
	coordheader( fdetail[i] );
}

///
/// @brief Function that initializes the coarse grained output file.
///
/// This function initializes the coarse grained output file.
/// It opens it up for writing and reading while formatting it with its header.
///
/// @param f Return pointer to the coarse grained output file being opened.
/// @param dir Path to the directory of the coarse grained output file.
/// @param fname Name of the coarse grained output file.
/// @param ext Extension of the coarse grained output file.
///
void opencoarse( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	coarseheader( *f );
}

///
/// @brief Function that initializes the global average velocity MPCD output file.
///
/// This function initializes the global average velocity MPCD output file.
/// It opens it up for writing and reading while formatting it with its header.
///
/// @param f Return pointer to the global average velocity MPCD output file being opened.
/// @param dir Path to the directory of the global average velocity MPCD output file.
/// @param fname Name of the coarse grained output file.
/// @param ext Extension of the coarse grained output file.
///
void openavvel( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	// 	avvelheader( *favvel );
	avvelWithGradVelheader( *f );
}

///
/// @brief Function that initializes the director output file.
///
/// This function initializes the director output file.
/// It opens it up for writing and reading while formatting it with its header.
///
/// @param f Return pointer to the director output file being opened.
/// @param dir Path to the directory of the director output file.
/// @param fname Name of the director output file.
/// @param ext Extension of the director output file.
///
void openorder( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	orderheader( *f );
}

///
/// @brief Function that initializes the tensor order parameter output file.
///
/// This function initializes the tensor order parameter output file.
/// It opens it up for writing and reading while formatting it with its header.
///
/// @param f Return pointer to the tensor order parameter output file being opened.
/// @param dir Path to the directory of the tensor order parameter output file.
/// @param fname Name of the tensor order parameter output file.
/// @param ext Extension of the tensor order parameter output file.
///
void openorderQ( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	orderQheader( *f );
}

///
/// @brief Function that initializes the tensor order parameter (in reciprocal space) output file.
///
/// This function initializes the tensor order parameter (in reciprocal space) output file.
/// It opens it up for writing and reading while formatting it with its header.
///
/// @param f Return pointer to the tensor order parameter (in reciprocal space) output file being opened.
/// @param dir Path to the directory of the tensor order parameter (in reciprocal space) output file.
/// @param fname Name of the tensor order parameter (in reciprocal space) output file.
/// @param ext Extension of the tensor order parameter (in reciprocal space) output file.
///
void openorderQK( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	orderQKheader( *f );
}

///
/// @brief Function that initializes the mean scalar order parameter output file.
///
/// This function initializes the mean scalar order parameter output file.
/// It opens it up for writing and reading while formatting it with its header.
///
/// @param f Return pointer to the mean scalar order parameter output file being opened.
/// @param dir Path to the directory of the mean scalar order parameter output file.
/// @param fname Name of the mean scalar order parameter output file.
/// @param ext Extension of the mean scalar order parameter output file.
///
void openavs( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	avsheader( *f );
}

///
/// @brief Function that initializes the density variations output file.
///
/// This function initializes the density variations output file.
/// It opens it up for writing and reading while formatting it with its header.
///
/// @param f Return pointer to the density variations output file being opened.
/// @param dir Path to the directory of the density variations output file.
/// @param fname Name of the density variations output file.
/// @param ext Extension of the density variations output file.
///
void opendensSTD( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	densheader( *f );
}

///
/// @brief Function that initializes the velocity distribution output file.
///
/// This function initializes the velocity distribution output file.
/// It opens it up for writing and reading while formatting it with its header.
///
/// @param f Return pointer to the velocity distribution output file being opened.
/// @param dir Path to the directory of the velocity distribution output file.
/// @param fname Name of the velocity distribution output file.
/// @param ext Extension of the velocity distribution output file.
///
void openhistVel( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	histVelheader( *f );
}

///
/// @brief Function that initializes the speed distribution output file.
///
/// This function initializes the speed distribution output file.
/// It opens it up for writing and reading while formatting it with its header.
///
/// @param f Return pointer to the speed distribution output file being opened.
/// @param dir Path to the directory of the speed distribution output file.
/// @param fname Name of the speed distribution output file.
/// @param ext Extension of the speed distribution output file.
///
void openhistSpeed( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	histSpeedheader( *f );
}

///
/// @brief Function that initializes the vorticity distribution output file.
///
/// This function initializes the vorticity distribution output file.
/// It opens it up for writing and reading while formatting it with its header.
///
/// @param f Return pointer to the vorticity distribution output file being opened.
/// @param dir Path to the directory of the vorticity distribution output file.
/// @param fname Name of the vorticity distribution output file.
/// @param ext Extension of the vorticity distribution output file.
///
void openhistVort( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	histVortheader( *f );
}

///
/// @brief Function that initializes the enstrophy distribution output file.
///
/// This function initializes the enstrophy distribution output file.
/// It opens it up for writing and reading while formatting it with its header.
///
/// @param f Return pointer to the enstrophy distribution output file being opened.
/// @param dir Path to the directory of the enstrophy distribution output file.
/// @param fname Name of the enstrophy distribution output file.
/// @param ext Extension of the enstrophy distribution output file.
///
void openhistEnstrophy( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	histEnstrheader( *f );
}

///
/// @brief Function that initializes the director distribution output file.
///
/// This function initializes the director distribution output file.
/// It opens it up for writing and reading while formatting it with its header.
///
/// @param f Return pointer to the director distribution output file being opened.
/// @param dir Path to the directory of the director distribution output file.
/// @param fname Name of the director distribution output file.
/// @param ext Extension of the director distribution output file.
///
void openhistDir( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	histDirheader( *f );
}

///
/// @brief Function that initializes the scalar order parameter distribution output file.
///
/// This function initializes the scalar order parameter distribution output file.
/// It opens it up for writing and reading while formatting it with its header.
///
/// @param f Return pointer to the scalar order parameter distribution output file being opened.
/// @param dir Path to the directory of the scalar order parameter distribution output file.
/// @param fname Name of the scalar order parameter distribution output file.
/// @param ext Extension of the scalar order parameter distribution output file.
///
void openhistS( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	histSheader( *f );
}

///
/// @brief Function that initializes the density distribution output file.
///
/// This function initializes the density distribution output file.
/// It opens it up for writing and reading while formatting it with its header.
///
/// @param f Return pointer to the density distribution output file being opened.
/// @param dir Path to the directory of the density distribution output file.
/// @param fname Name of the density distribution output file.
/// @param ext Extension of the density distribution output file.
///
void openhistDens( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	histNheader( *f );
}

///
/// @brief Function that initializes the mean enstrophy output file.
///
/// This function initializes the mean enstrophy output file.
/// It opens it up for writing and reading while formatting it with its header.
///
/// @param f Return pointer to the mean enstrophy output file being opened.
/// @param dir Path to the directory of the mean enstrophy output file.
/// @param fname Name of the mean enstrophy output file.
/// @param ext Extension of the mean enstrophy output file.
///
void openavenstrophy( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	avenstrophyheader( *f );
}

///
/// @brief Function that initializes the flow field output file.
///
/// This function initializes the flow field output file.
/// It opens it up for writing and reading while formatting it with its header.
///
/// @param f Return pointer to the flow field output file being opened.
/// @param dir Path to the directory of the flow field output file.
/// @param fname Name of the flow field output file.
/// @param ext Extension of the flow field output file.
///
void openflow( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	flowheader( *f );
}

///
/// @brief Function that initializes the velocity field output file.
///
/// This function initializes the velocity field output file.
/// It opens it up for writing and reading while formatting it with its header.
///
/// @param f Return pointer to the flow field output file being opened.
/// @param dir Path to the directory of the flow field output file.
/// @param fname Name of the flow field output file.
/// @param ext Extension of the flow field output file.
///
void openvel( FILE **f,char dir[],char fname[],char ext[] ) {
    openBasic( f,dir,fname,ext );
    flowheader( *f );
}

///
/// @brief Function that initializes the energy output file.
///
/// This function initializes the energy output file.
/// It opens it up for writing and reading while formatting it with its header.
///
/// @param f Return pointer to the energy output file being opened.
/// @param dir Path to the directory of the energy output file.
/// @param fname Name of the energy output file.
/// @param ext Extension of the energy output file.
///
void openenergy( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	energyheader( *f );
}

///
/// @brief Function that initializes the energy field output file.
///
/// This function initializes the energy field output file.
/// It opens it up for writing and reading while formatting it with its header.
///
/// @param f Return pointer to the energy field output file being opened.
/// @param dir Path to the directory of the energy field output file.
/// @param fname Name of the energy field output file.
/// @param ext Extension of the energy field output file.
///
void openenergyfield( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	energyfieldheader( *f );
}

///
/// @brief Function that initializes the energy from neighbours output file.
///
/// This function initializes the energy from neighbours output file.
/// It opens it up for writing and reading while formatting it with its header.
///
/// @param f Return pointer to the energy from neighbours output file being opened.
/// @param dir Path to the directory of the energy from neighbours output file.
/// @param fname Name of the energy from neighbours output file.
/// @param ext Extension of the energy from neighbours output file.
///
void openenergyneighbours( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	energyneighboursheader( *f );
}

///
/// @brief Function that initializes the synopsis output file.
///
/// This function initializes the synopsis output file.
/// It opens it up for writing and reading while formatting it with its header.
/// when called for the first time, as specified by the `firsttime` parameter, it truncates the file to
/// zero length. When called for the second time, it is opened in apending mode.
///
/// @param fsynopsis Return pointer to the synopsis output file being opened.
/// @param dir Path to the directory of the synopsis output file.
/// @param firsttime Integer specifying if it is the first time opening the synopsis file.
///
void opensynopsis( FILE **fsynopsis,char dir[],int firsttime ) {
	char filename[200];
	char filesynopsis[]="synopsis";
	char fileextension[]=".dat";

	strcpy( filename,dir );
	strcat( filename,filesynopsis );
	strcat( filename,fileextension );
	if( firsttime ) *fsynopsis = fopen( filename,"w+" );
	else *fsynopsis = fopen( filename,"a" );
	if( !*fsynopsis ){ // file couldn't be opened
		printf("Error: File '%s'could not be opened.\n",filename );
		exit( 1 );
	}
	if( firsttime ) outheader( *fsynopsis,0 );
}

///
/// @brief Function that initializes the solids' trajectories (or boundary condition (BC) motion) output files.
///
/// This function initializes the solids' trajectories (or boundary condition (BC) motion) output files.
/// (one file for each BC).
///
/// @param bc Integer specifying the boundary condition whose motion is being outputed.
/// @param fsolids Array of return pointers to the list of files associated to all BCs.
/// @param dir Path to the directory of the BC motion output file.
/// @param filesolids Name of the BC motion output file.
/// @param filesuffix Suffix specifying the BC associated to the file being opened, updated to "bc" within the function.
/// @param fileextension Extension of the BC motion output file.
///
void opentraj( int bc,FILE *fsolids[],char dir[],char filesolids[],char filesuffix[],char fileextension[] ){
	char filename[200];

	strcpy( filename,dir );
	strcat( filename,filesolids );
	sprintf (filesuffix,"%i",bc );
	strcat( filename,filesuffix );
	strcat( filename,fileextension );
	fsolids[bc] = fopen( filename,"w+" );
	if( !fsolids[bc] ){ // file couldn't be opened
		printf( "Error: File '%s'could not be opened.\n",filename );
		exit( 1 );
	}
	outheader( fsolids[bc],bc );
	solidsheader( fsolids[bc] );
}

///
/// @brief Function thet opens the MPCD particles positions file for reading.
///
/// This function opens the MPCD particles positions files for reading. The files (one for each species) must exist.
///
/// @param i Index of the relevant species.
/// @param fin Array of return pointers to the list of files associated to all species.
/// @param dir Path to the directory of the MPCD particles positions file.
/// @param fileprefix Name of the directory of the MPCD particles positions file.
/// @param filesuffix Suffix specifying the species associated to the file being opened, updated to "i" within the function.
/// @param fileextension Extension of the MPCD particles position file.
///
void openplace( int i,FILE *fin[],char dir[],char fileprefix[],char filesuffix[16],char fileextension[] ){
	char filename[200];

	strcpy( filename,dir );
	strcat( filename,fileprefix );
	sprintf( filesuffix,"%i",i );
	strcat( filename,filesuffix );
	strcat( filename,fileextension );
	fin[i] = fopen( filename,"r" );
	if( !fin[i] ) { // file couldn't be opened
		printf( "Error: File '%s'could not be opened.\n",filename );
		exit( 1 );
	}
}

///
/// @brief Function that initializes a correlation output file.
///
/// This function initializes a correlation output file.
/// It opens it up for writing and reading while formatting it with its header.
///
/// @param f Return pointer to the correlation output file being opened.
/// @param dir Path to the directory of the correlation output file.
/// @param fname Name of the correlation output file.
/// @param ext Extension of the correlation output file.
///
void opencorr( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	corrheader( *f );
}

///
/// @brief Function that initializes the energy spectrum output file.
///
/// This function initializes the energy spectrum output file.
/// It opens it up for writing and reading while formatting it with its header.
///
/// @param f Return pointer to the energy spectrum output file being opened.
/// @param dir Path to the directory of the energy spectrum output file.
/// @param fname Name of the energy spectrum output file.
/// @param ext Extension of the energy spectrum output file.
///
void openenergyspect( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	energyspectheader( *f );
}

///
/// @brief Function that initializes the enstrophy spectrum output file.
///
/// This function initializes the enstrophy spectrum output file.
/// It opens it up for writing and reading while formatting it with its header.
///
/// @param f Return pointer to the enstrophy spectrum output file being opened.
/// @param dir Path to the directory of the enstrophy spectrum output file.
/// @param fname Name of the enstrophy spectrum output file.
/// @param ext Extension of the enstrophy spectrum output file.
///
void openenstrophyspect( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	enstrophyspectheader( *f );
}

///
/// @brief Function that initializes the topological charge field output file.
///
/// This function initializes the topological charge field output file (only for 2D systems).
/// It opens it up for writing and reading while formatting it with its header.
///
/// @param f Return pointer to the topological charge field output file being opened.
/// @param dir Path to the directory of the topological charge field output file.
/// @param fname Name of the topological charge field output file.
/// @param ext Extension of the topological charge field output file.
///
void opentopo( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	if(DIM==_3D) printf("Warning: Topological charge field is only outputted for 2D!\n");
	else topoheader( *f );
}

///
/// @brief Function that initializes the defect tracker output file.
///
/// This function initializes the defect tracker output file (only for 2D systems).
/// It opens it up for writing and reading while formatting it with its header.
///
/// @param f Return pointer to the defect tracker output file being opened.
/// @param dir Path to the directory of the defect tracker output file.
/// @param fname Name of the defect tracker output file.
/// @param ext Extension of the defect tracker output file.
///
void opendefect( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	if(DIM==_3D) printf("Warning: Defects are only outputted for 2D!\n");
	else defectheader( *f );
}

///
/// @brief Function that initializes the disclination tensor output file.
///
/// This function initializes the disclination tensor output file.
/// It opens it up for writing and reading while formatting it with its header.
///
/// @param f Return pointer to the disclination tensor output file being opened.
/// @param dir Path to the directory of the disclination tensor output file.
/// @param fname Name of the disclination tensor output file.
/// @param ext Extension of the disclination tensor output file.
///
void opendisclin( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	disclinTensorheader( *f );
}

///
/// @brief Function that initializes the phi/color/species-type field output file.
///
/// This function initializes the  phi/color/species-type field output file.
/// It opens it up for writing and reading while formatting it with its header.
///
/// @param f Return pointer to the  phi/color/species-type field output file being opened.
/// @param dir Path to the directory of the  phi/color/species-type field output file.
/// @param fname Name of the  phi/color/species-type field output file.
/// @param ext Extension of the  phi/color/species-type field output file.
///
void openmultiphase( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	multiphaseheader( *f );
}

///
/// @brief Function that initializes the pressure field output file.
///
/// This function initializes the pressure field output file.
/// It opens it up for writing and reading while formatting it with its header.
///
/// @param f Return pointer to the pressure field output file being opened.
/// @param dir Path to the directory of the pressure field output file.
/// @param fname Name of the pressure field output file.
/// @param ext Extension of the pressure field output file.
///
void openpressure( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	pressureheader( *f );
}

///
/// @brief Function that initializes the Binder cumulant output file.
///
/// This function initializes the Binder cumulant output file.
/// It opens it up for writing and reading while formatting it with its header.
///
/// @param f Return pointer to the Binder cumulant output file being opened.
/// @param dir Path to the directory of the Binder cumulant output file.
/// @param fname Name of the Binder cumulant output file.
/// @param ext Extension of the Binder cumulant output file.
/// @param binSize The size of bins.
///
void openbinder( FILE **f,char dir[],char fname[],char ext[],int binSize ) {
	openBasic( f,dir,fname,ext );
	binderheader( *f,binSize );
}

///
/// @brief Function that initializes the swimmer output file.
///
/// This function initializes the swimmer output file.
/// It opens it up for writing and reading while formatting it with its header.
///
/// @param f Return pointer to the swimmer output file being opened.
/// @param dir Path to the directory of the swimmer output file.
/// @param fname Name of the swimmer output file.
/// @param ext Extension of the swimmer output file.
///
void openswimmer( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	swimmerheader( *f );
}

///
/// @brief Function that initializes the swimmer orientation output file.
///
/// This function initializes the swimmer orientation output file.
/// It opens it up for writing and reading while formatting it with its header.
///
/// @param f Return pointer to the swimmer orientation output file being opened.
/// @param dir Path to the directory of the swimmer orientation output file.
/// @param fname Name of the swimmer orientation output file.
/// @param ext Extension of the swimmer orientation output file.
///
void openswimmerOri( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	swimmeroriheader( *f );
}

///
/// @brief Function that initializes the swimmer run and tumble output file.
///
/// This function initializes the swimmer run and tumble output file.
/// It opens it up for writing and reading while formatting it with its header.
///
/// @param f Return pointer to the swimmer run and tumble output file being opened.
/// @param dir Path to the directory of the swimmer run and tumble output file.
/// @param fname Name of the swimmer run and tumble output file.
/// @param ext Extension of the swimmer run and tumble output file.
///
void openruntumble( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	runtumbleheader( *f );
}
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ***************** THEORY ***************** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */

///
/// @brief Function that calculates several coefficients characterizing the system.
///
/// This function calculates: the mean free path, the speed of sound, the dynamic viscocity,
/// the thermal diffusion coefficient and the self-diffusion coefficient.
/// It assumes that the cell size is a=1. The calculation follows Noguchi & Gompper,
/// Transport coefficients of off-lattice mesoscale-hydrodynamics simulation techniques, PRE 78, 016706 (2008).
/// If a synopsis file is requested, this function also prints the above information to the file.
///
/// @param MFP Return pointer to the mean free path.
/// @param VISC Return pointer to the viscocity.
/// @param THERMD Return pointer to the thermal diffusion coefficient.
/// @param SDIFF Return pointer to the self-diffusion coefficient.
/// @param SPEEDOFSOUND Return pointer to the speed of sound.
/// @param RA Rotation angle for SRD.
/// @param FRICCO Friction coefficient for Langevin thermostat.
/// @param KBT Temperature (a third of thermal energy).
/// @param dt MPCD time-step value.
/// @param sumM Sum of all masses.
/// @param RTECH Rotation technique.
/// @param SYNOUT Integer specifying if the synopsis file should be outputted (1 for yes, 0 for no).
/// @param fsynopsis Pointer to the synopsis file.
///
void theory_trans( double *MFP,double *VISC,double *THERMD,double *SDIFF,double *SPEEDOFSOUND,double RA,double FRICCO,double KBT,double dt,double sumM,int RTECH,int SYNOUT,FILE *fsynopsis ) {
	double a=1.0;							//MPCD cell size
	double A,B,CM,SM;								//Correlation factors from Table 1 of Nguchi & Gompper
	double VISCKIN,VISCCOL;		//Kinetic and collisional parts of DYNAMIC viscosity
	double THERMDKIN,THERMDCOL;		//Kinetic and collisional parts of thermal diffusion coefficient
	double inv_nDNST;		//Inverse of the number density
	double M;						//Just a value that comes up often
	double inv_DIM;		//Inverse of the dimensionality
	double avMASS = sumM / (double)GPOP;	//Average mass of a particle

	inv_DIM=1.0/(float)DIM;
	inv_nDNST = 1.0/nDNST;
	M = nDNST / (nDNST -1.0 + smrtPow( e,-1.0*nDNST ) );
	//Calculate mean free path
	*MFP = dt * sqrt( KBT/avMASS );
	//Calculate the speed of sound
	*SPEEDOFSOUND = sqrt( 2.0*KBT/avMASS );

	//Calculate the correlation factors
	if(RTECH==ORTHAXIS) {
		A=2.0*inv_DIM*(1.0-cos(RA));
		if(DIM==_2D) B=1.0-cos(2.0*RA);
		else B=(2.0/5.0)*(2.0-cos(RA)-cos(2.0*RA));
	}
	else if(RTECH==ARBAXIS ){
		A=2.0*inv_DIM*(1.0-sin(RA)/RA);
		if(DIM==_2D) B=1.0-0.5*sin(2.0*RA)/RA;
		else B=(2.0/5.0)*(2.0-sin(RA)/RA-0.5*sin(RA)/RA);
	}
	else if(RTECH==MPCAT || RTECH==RAT ) {
		A=1.0;
		B=1.0;
	}
	else if(RTECH==LANG || RTECH==RLANG) {
		A=(FRICCO*dt/avMASS)/(1.0+0.5*FRICCO*dt/avMASS);
		B=(2.0*FRICCO*dt/avMASS)/(1.0+0.5*FRICCO*dt/avMASS)/(1.0+0.5*FRICCO*dt/avMASS);
	}
	else {
		// All other versions do not have theoretically known transport coefficients
		A=1.0;
		B=1.0;
	}
	//Calculate the viscosity
	// From https://journals.aps.org/pre/abstract/10.1103/PhysRevE.78.016706
	VISCKIN = 0.0;
	VISCCOL = 0.0;
	//Kinetic part of viscosity
	if(RTECH==ORTHAXIS || RTECH==ARBAXIS || RTECH==MPCAT || RTECH==LANG) {
		//All of the versions without angular-momentum conservation have the same form
		CM=B/M;
		VISCKIN=nDNST*KBT*dt*smrtPow(a,-DIM)*( 1.0/CM-0.5 );
	}
	else if(RTECH==RAT || RTECH==RLANG) {
		//All of the versions with angular-momentum conservation have the same form
		CM=B*(1.0-exp(-nDNST)*(1.0+nDNST)) + (A+inv_DIM*B)*nDNST*exp(-nDNST)/(float)(DIM+2.0);
		CM+=0.5*(A*DIM-0.5*B*(3.0*DIM+2.0))*(1.0-exp(-nDNST)*(1.0+nDNST+0.5*nDNST*nDNST))/nDNST;
		VISCKIN=nDNST*KBT*dt*smrtPow(a,-DIM)*( 1.0/CM-0.5 );
	}
	else {
		//Approximate with only the scaling result
		VISCKIN=nDNST*KBT*dt*smrtPow(a,-DIM);
	}

	//Collisional part of viscosity
	if(RTECH==ORTHAXIS || RTECH==ARBAXIS || RTECH==MPCAT || RTECH==LANG) {
		//All of the versions without angular-momentum conservation have the same form
		VISCCOL=(A*nDNST/M)*avMASS/(12.0*dt*smrtPow(a,DIM-2));
	}
	else if(RTECH==RAT || RTECH==RLANG) {
		//All of the versions with angular-momentum conservation have the same form
		VISCCOL=(A*avMASS)/(24.0*dt*smrtPow(a,DIM-2))*( nDNST-7.0/5.0+exp(-nDNST)*(7.0/5.0+2.0*nDNST/5.0+(inv_DIM-0.3)*nDNST*nDNST) );
	}
	else {
		//Approximate with only the scaling result
		VISCCOL=avMASS/(dt*smrtPow(a,DIM-2));
	}
	//Total viscosity
	*VISC=VISCKIN+VISCCOL;

	//Calculate the self diffusion coefficient
	if(RTECH==ORTHAXIS || RTECH==ARBAXIS || RTECH==MPCAT || RTECH==LANG) {
		//All of the versions without angular-momentum conservation have the same form
		SM=A/M;
		*SDIFF = (KBT*dt/avMASS)*(SM-0.5);
	}
	else if(RTECH==RAT || RTECH==RLANG) {
		//All of the versions with angular-momentum conservation have the same form
		SM=0.5*exp(-nDNST)*( 0.5*inv_DIM*nDNST*(DIM-1.0)*(DIM-2.0) + DIM - 1.0 + (DIM+1.0)*inv_nDNST );
		SM+=1.0-0.5*(DIM+1.0)*inv_nDNST;
		SM*=A;
		*SDIFF = (KBT*dt/avMASS)*(SM-0.5);
	}
	else {
		//Approximate with only the scaling result
		*SDIFF = (KBT*dt/avMASS);
	}

	if(RTECH==ORTHAXIS || RTECH==ARBAXIS ) {
		//Calculate the thermal diffusion coefficient
		if( DIM ==  _2D ) {
			THERMDKIN = 2. / (1.-cos(RA));
			THERMDKIN -= 1.;
			THERMDKIN += ( 4.*inv_nDNST )*( 1.-0.25 / ( sin( 0.5*RA )*sin( 0.5*RA ) ) );
			THERMDKIN *= 0.5 * KBT * dt;
		}
		else if( DIM == _3D ) {
			THERMDKIN = 0.5 * (2.+cos(RA)) / (1.-cos(RA));
			THERMDKIN += (3.*inv_nDNST) * (0.8-0.25 / (sin( 0.5*RA )*sin( 0.5*RA )));
			THERMDKIN *= KBT*dt;
		}
		else {
			printf( "DIM must be 2 or 3 to calculated thermal diffusivity.\n" );
			THERMDKIN = 0.;
		}
		THERMDCOL = (1./(12.*(DIM+2.)*dt)) * ((nDNST)-0.555555*nDNST*nDNST) * (1.-cos(RA));
		if( THERMDCOL <= 0. ){
			THERMDCOL = 1. + smrtPow( e,-nDNST ) * ( log( nDNST )-1. ) ;
			THERMDCOL -= inv_nDNST;
			THERMDCOL -= smrtPow( inv_nDNST,2. );
			THERMDCOL -= 2. * smrtPow( inv_nDNST,3. );
			THERMDCOL *= (inv_nDNST) / (3.*(DIM+2.)*dt);
			THERMDCOL *= (1.-cos(RA));
		}
	}
	else {
		THERMDKIN = 0.5;
		THERMDCOL = 0.5;
	}
	*THERMD = THERMDKIN+THERMDCOL;
	if( SYNOUT == OUT ) {
		fprintf( fsynopsis,"Mean Free Path: %lf\n",*MFP );
		fprintf( fsynopsis,"Dynamic Viscosity: %lf\n",*VISC );
		fprintf( fsynopsis,"\tkinetic contribution: %lf\n\tcollisional contribution: %lf\n",VISCKIN,VISCCOL );
		fprintf( fsynopsis,"Kinematic Viscosity: %lf\n",(*VISC)*inv_nDNST/avMASS );
		fprintf( fsynopsis,"\tkinetic contribution: %lf\n\tcollisional contribution: %lf\n",VISCKIN*inv_nDNST/avMASS,VISCCOL*inv_nDNST/avMASS );
		fprintf( fsynopsis,"Self Diffusion Coefficient: %lf\n",*SDIFF );
		fprintf( fsynopsis,"Schmidt number: %lf\n",(*VISC)/(*SDIFF)/mDNST );
		fprintf( fsynopsis,"Speed of sound: %lf\n",*SPEEDOFSOUND );
		fprintf( fsynopsis,"Thermal Diffusion Coefficient: %lf\n",*THERMD );
	}
}

///
/// @brief Function that calculates the fluid number density.
///
/// This function calculates the number density of the fluid, either per unit area (2D) or volume (3D).
/// The volume occupied by boundary conditions, such as colloids, is excluded from the calculation.
///
/// @param WALL Array of the system's boundary conditions.
/// @return The fluid number density.
///
double ndensity( bc WALL[] ) {
/*
   Calculates the number density of the fluid.
   Either per unit volume or area
*/
	int i;
	double V;	//The control volume.
	double D;	//The number density
	V = (double)( XYZ[0] * XYZ[1] * XYZ[2] );
	for( i=0; i<NBC; i++ ) V -= WALL[i].VOL;
	D = GPOP / V;
	return D;
}

///
/// @brief Function that calculates the fluid number density.
///
/// This function calculates the number density of the fluid, either per unit area (2D) or volume (3D).
/// The volume occupied by boundary conditions, such as colloids, is excluded from the calculation.
///
/// @param WALL Array of the system's boundary conditions.
/// @param MASS The system total mass.
/// @return The fluid mass density.
///
double mdensity( bc WALL[],double MASS ) {
	int i;
	double V;	//The control volume.
	double D;	//The mass density
	V = (double)( XYZ[0] * XYZ[1] * XYZ[2] );
	for( i=0; i<NBC; i++ ) V -= WALL[i].VOL;
	D = MASS / V;
	return D;
}

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ***************** ZEROING **************** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */


///
/// @brief Function that zeros everything in all the particles.
///
/// This functios zeros all parameters for all MPCD particles.
///
/// @param pp Return pointer to the first MPCD particle in the array.
///
void zeroparticles( particleMPC *pp ) {
	int i;

	for( i=0; i<GPOP; i++ ) {
        zerovec( (pp+i)->Q, _3D );
        zerovec( (pp+i)->V, _3D );
        zerovec( (pp+i)->U, _3D );
        zerovec( (pp+i)->T, _3D );
	}
}

///
/// @brief Function that zeros counters.
///
/// This function zeros the current un-thermostated temperature, the current average flow velocity
/// and the average of the scalar order parameter.
///
/// @param KBTNOW Return pointer to the current un-thermostated temperature.
/// @param AVNOW Return pointer to the current average flow velocity.
/// @param AVS Return pointer to the average scalar order parameter.
///
void zerocnt( double *KBTNOW,double AVNOW[],double *AVS ) {
	*KBTNOW = 0.;
	*AVS = 0.;
    zerovec( AVNOW, _3D );
}

///
/// @brief Function that zeros a vector histogram.
///
/// This function zeros a vector histogram.
///
/// @param HIST Vector histogram being zeroed.
///
void zeroHISTVEC( int HIST[_3D][BINS] ) {
    int i;
    for (i=0; i < _3D; i++) {
        zerovec( (double *) HIST[i], BINS);
    }
}

///
/// @brief Function that zeros a scalar histogram.
///
/// This function zeros a scalar histogram.
///
/// @param HIST Scalar histogram being zeroed.
///
void zeroHISTSCALAR( int HIST[BINS] ) {
    zerovec( (double *) HIST, BINS);
}

///
/// @brief Function that zeros all the contents of a cell list.
///
/// This functions zeros the entire content of the receiving cell list.
/// This includes: population, mass, scalar order parameter, center of mass, velocity of center
/// of mass, flow velocity, director, velocity gradient tensor, moment of inertia, streaming and
/// collisional parts of the stress tensor.
/// Pointers to the first SRD, MD and swimmer particles are set to NULL.
///
/// @param CL Return pointer to the cell list being zeroed.
///
void zerocell( cell ***CL ) {
	int i,j,k,l;
	for( i=0; i<XYZ_P1[0]; i++ ) for( j=0; j<XYZ_P1[1]; j++ ) for( k=0; k<XYZ_P1[2]; k++ ) {
		CL[i][j][k].POP = 0;
		CL[i][j][k].POPSRD = 0;
		CL[i][j][k].POPSW = 0;
		CL[i][j][k].POPMD = 0;
		CL[i][j][k].MASS = 0.0;
		CL[i][j][k].S = 0.0;

        zerovec( CL[i][j][k].CM, _3D );
        zerovec( CL[i][j][k].VCM, _3D );
        zerovec( CL[i][j][k].FLOW, _3D );
        zerovec( CL[i][j][k].DIR, _3D );

        for (l=0; l<_3D; l++) { // zero matrices
            zerovec( CL[i][j][k].E[l], _3D );
            zerovec( CL[i][j][k].I[l], _3D );
            zerovec( CL[i][j][k].Ps[l], _3D );
            zerovec( CL[i][j][k].Pc[l], _3D );
        }

		//The list doesn't have anyone in it yet so it doesn't point anywhere
		CL[i][j][k].pp = NULL;
		CL[i][j][k].MDpp = NULL;
		CL[i][j][k].sp = NULL;
	}
}

///
/// @brief Function that zeros the collisional pressure term.
///
/// This function zeros the collisional pressure term.
///
/// @param CL Return pointer to the cell list whose collisional pressure term is being zeroed.
///
void zeroPressureColl( cell *CL ) {
    int i;
    for (i=0; i<DIM; i++) {
        zerovec( CL->Pc[i], DIM );
    }
}

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************** INITIALIZING ************** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */

///
/// @brief Function that initializes program variables and physical parameters.
///
/// This functions initializes many program variables and some physical parameters. These include setting
/// the seed for the random engine, time clock, sum of all masses, average speed (assuming an isotropic distribution),
/// cosine and sine of the rotation angle for SRD, external magnetic field, number density and mass density.
/// It also sets the physical parameters of the boundary conditions: their volume and moment of inertia.
/// In addition, it zeros the runtimes, warmtimes, average director, active velocities,
/// everything the cell lists and in the particles.
///
/// @param seed Return pointer to seed.
/// @param to CPU time.
/// @param co Clock time.
/// @param runtime Return pointer to runtime.
/// @param warmtime Return pointer to warmtime.
/// @param sumM Return pointer to sum of all masses.
/// @param AV Return pointer to average velocity.
/// @param avDIR Return pointer to average director.
/// @param SP Array of all species.
/// @param C Return pointer to the cosine of the rotation angle for SRD.
/// @param S Return pointer to the sine of the rotation angle for SRD.
/// @param RA Rotation angle for SRD.
/// @param AVVEL Return pointer to the average speed.
/// @param KBT Temperature (a third of thermal energy).
/// @param WALL Return pointer to array of boundary conditions.
/// @param MAG Return pointer to constant external magnetic field.
/// @param CL Return pointer to the array of all cell lists.
/// @param pp Return pointer to the array of all MPCD particles.
///
void initvar( unsigned long *seed,time_t *to,clock_t *co,int *runtime,int *warmtime,double *sumM,double AV[_3D],double avDIR[_3D],spec SP[],double *C,double *S,double RA,double *AVVEL,double KBT,bc WALL[],double MAG[_3D],cell ***CL,particleMPC *pp ) {
	int i;

	*seed = RandomSeedSRD (*seed);			//note to tyler

	//init_genrand( time( NULL ) );	//Initialize random number generator.
	*to = time( NULL );
	*co = clock( );
	*runtime = 0;
	*warmtime = 0;
	for( i=0; i<_3D; i++ ) AV[i] = 0.;
	for( i=0; i<_3D; i++ ) avDIR[i] = 0.;
	*sumM = 0.;
	for( i=0; i<NSPECI; i++ ) *sumM += (SP[i].POP) * (SP[i].MASS);
	*C = cos( RA );
	*S = sin( RA );
	*AVVEL = sqrt( DIM*KBT * GPOP / (*sumM) );	//If isotropic velocity dist.
	//Calculate physical parameters of BCs
	for( i=0; i<NBC; i++ ) {		//Set the velocity of the walls to zero
		dim_vol( &WALL[i],XYZ,DIM );
		mominert( &WALL[i],XYZ,DIM );
	}

	nDNST = ndensity( WALL );
	mDNST = mdensity( WALL,*sumM );
	for( i=0; i<DIM; i++ ) MAG[i]/=(nDNST);//Want torque per unit volume so divide field by number density
	//Zero everything in the cell lists
	zerocell( CL );
	//Zero everything in the particles
	zeroparticles( pp );
}

///
/// @brief Function that initiates the position of a MPCD particle.
///
/// This function initiates the position for a MPCD particle, either randomly
/// or reading it from a file.
///
/// @param Q Return pointer to the position of the MPCD particle.
/// @param PL Integer that is 0 if position is determined randomly, 1 if read from a file.
/// @param fin File from which to read the particle position.
///
void place( double Q[],int PL,FILE *fin ) {
	int d;

	if( PL == PRF ) for( d=0; d<DIM; d++ ) Q[d] = genrand_real( ) * XYZ[d];
	else if( PL == READ ) for( d=0; d<DIM; d++ ) {
		if(fscanf( fin, "%lf",&Q[d] ));
		else printf("Warning: Failed to read place input file.\n");
	}
	else{
		printf( "Error: Particle placement distribution unacceptable.\n" );
		exit( 1 );
	}
}

///
/// @brief Function that randomly places an MPCD particle.
///
/// This function randomly places an MPCD particle within one SRD cell of its current position.
///
/// @param p Return pointer to MPCD particle whose position is being altered.
///
void replace( particleMPC *p ) {
	int i;
	for( i=0; i<DIM; i++ ) p->Q[i] = (double)XYZ[i] * genrand_real( );
}

///
/// @brief Function that sets the initial velocity of an MPCD particle.
///
/// This subroutines sets the initial velocity of an MPCD particle.
/// It can do so either: i) randomly with a uniform distribution, ii) randomly with a gaussian
/// distribution, iii) randomly with a spherically symmetric and separable gaussian distribution,
/// iv) with all particles having the average speed (but travelling along each axis in either direction)
/// or v)reading it from a file.  For all random methods, the average speed is fixed at sqrt(DIM*KBT/M).
///
/// @param V Return pointer to an MPCD particle's velocity.
/// @param KBT Temperature (a third of thermal energy).
/// @param PL Integer indicating assignment method: 0 if uniformly random, 1 if read from file, 2 for spherically symmetric gaussian,
///           3 for average speed in each axis (random direction) and 4 for Gaussian distribution.
/// @param MASS Mass of the particle.
/// @param fin File from which the velocity is read.
///
void push( double V[],double KBT,int PL,double MASS,FILE *fin ) {
	int d;
	double normalize;
	if( PL == RANDVEL ) for( d=0; d<DIM; d++ ) V[d] = sqrt( KBT/MASS ) * (2. * genrand_real() - 1.);
							//The particleMPCs have an average speed of sqrt(DIM*KBT/M) BUT each compenent makes this up so v^2=vx^2+vy^2+vz^2 where vx=vy=vz so must divide by sqrt(3) to get component thus no 3.
							//Also for uniform must consider from zero to 2*sqrt(KBT/M).
							//To get negative we make the interval twice as big and subtract off the original max.
	else if(PL == READ) for( d=0; d<DIM; d++ ) {
		if(fscanf( fin, "%lf",&V[d] ));
		else printf("Warning: Failed to read push input file.\n");
	}
	else if( PL == HOMOAXIS ) for( d=0; d<DIM; d++ ) {
		V[d] = sqrt( DIM*KBT/MASS );			//Each particleMPC has the average energy
		V[d] *= genrand_pmOne();			//But could travel in either direction
	}
	else if( PL == GAUSS ) for( d=0; d<DIM; d++ ) {
		V[d] = genrand_gaussMB( KBT,MASS );	//Velocity distribution is Gaussian
	}
	else if( PL == HOMO ) {
		for( d=0; d<DIM; d++ ) 	V[d] = genrand_gaussMB( KBT,MASS );	//Velocity distribution is Gaussian (spherically symmetric and separable)
		norm( V,DIM );
		normalize = sqrt((double) DIM *KBT/MASS);
		for( d=0; d<DIM; d++ ) 	V[d] *= normalize;
	}
	else{
		printf( "Error: Particle velocity distribution unacceptable.\n" );
		exit( 1 );
	}
}

///
/// @brief Function that sets the initial orientation of an MPCD particle.
///
/// This subroutines sets the initial orientation. It can do so by:
/// i) set it randomly, ii) aligned to one of the axes, iii) in a 45 degree configuration,
/// iv) randomly parallel to one of the planes formef by the axes, v) pointing parallel to
/// the ray starting at the origin and pointing towards the top right corner ot vi) forming a
/// defect pair (+1/2 and -1/2) configuration.
///
/// @param U Return pointer to the MPCD particle's orientation.
/// @param Q Position of the MPCD particle.
/// @param PL Integer determining the desired configuration. 0 for random orientation, 1, 2 and 3 for
///           allignment parallel to the x, y and z axes, respectively, 4 for a 45 degree configuration,
///           5,6 and 7 for random alignment in the XY, XZ, YZ axes, respectively, 8 for pointing parallel
///           to line from origin to top right corner and 9 for the defect pair configuration.
///
void orient( double U[],double Q[],int PL ) {
	int d;
	double noise;
	if( PL==RANDORIENT ) {
		genrand_coneNP( U,pi,DIM );
// 		if(DIM==_2D) U[0]*=genrand_pmOne();
	}
	else if( PL==ALIGNX ) {
		//There is a problem sometimes if it is PERFECTLY aligned so give a bit of noise (+/- 1% max)
		noise = 1.0 + 0.02*( genrand_real()-0.5 );
		U[0]=0.9999*noise;
		U[1]=0.01414178*noise;
		if( DIM>=_3D ) U[2]=0.000008*noise;
	}
	else if( PL==ALIGNY ) {
		//For some reason the collision operation isn't happy if every nematogen is exactly aligned along y
		//so give a bit of noise (+/- 1% max)
		noise = 1.0 + 0.02*( genrand_real()-0.5 );
		U[1]=0.9999*noise;
		U[0]=0.01414178*noise;
		if( DIM>=_3D ) U[2]=0.000008*noise;
	}
	else if( PL==ALIGNZ ) {
		//There is a problem sometimes if it is PERFECTLY aligned so give a bit of noise (+/- 1% max)
		noise = 1.0 + 0.02*( genrand_real()-0.5 );
		U[2]=0.9999*noise;
		U[0]=0.007071*noise;
		U[1]=0.007071*noise;
	}
	else if( PL==ALIGN45 ) {
		if( DIM==_2D ) for( d=0; d<DIM; d++ ) U[d]=1./sqrt(2.);
		else if( DIM==_3D ) for( d=0; d<DIM; d++ ) U[d]=1./sqrt(3.);
		else if( DIM==_1D ) {
			printf( "Error: 1D liquid crystal.\n" );
			exit( 1 );
		}
	}
	else if( PL==PLANEZ ) {
		//Just 2D as far as the orientation is concerned
		genrand_coneNP( U,pi,2 );
	}
	else if( PL==PLANEY ) {
		//Just 2D but switch the ordering of the components to put in XZ-plane
		genrand_coneNP( U,pi,2 );
		U[2]=U[1];
		U[1]=0.0;
	}
	else if( PL==PLANEX ) {
		//Just 2D but switch the ordering of the components to put in YZ-plane
		genrand_coneNP( U,pi,2 );
		U[2]=U[0];
		U[0]=0.0;
	}
	else if( PL==ALIGNTR ) {
		// Orientation points from the origin to the top right corner (furthest point)
		// of the system size, with unit size

		double L2 = 0.0;
		for( d=0; d<DIM; d++ ) L2 += XYZ[d]*XYZ[d];
		L2 = sqrt(L2);

		for( d=0; d<DIM; d++ ) U[d] = XYZ[d]/L2;
	}
	else if( PL==ALIGNDEFECTPAIR ) {
		// Orientation points as if around two oppositely charged defects
		// of the system size, with unit size
		double QP[_2D],QM[_2D];		//Position of the Plus and Minus defects
		double UP[_2D],UM[_2D];		//Orientations of the defects
		double kP=0.5,kM=-0.5;		//Charges of the defects
		double phiP,phiM;			//Director angles at the point Q
		double wP,wM,wY;			//Weightings of U for each defect at the point Q and the walls(linear)
		double r,R;

		//Linear weighting slope
		R=0.5*sqrt( 2.25*((double)XYZ[0])*((double)XYZ[0]) + ((double)XYZ[1])*((double)XYZ[1]) );
		// Do the +1/2 defect on the left side
		QP[0] = ((double)XYZ[0])*0.25;
		QP[1] = ((double)XYZ[1])*0.5;
		phiP = atan2(QP[1]-Q[1],QP[0]-Q[0]);
		UP[0] = cos(kP*phiP);
		UP[1] = sin(kP*phiP);
		r=sqrt( pow(QP[0]-Q[0],2.0)+pow(QP[1]-Q[1],2.0) );
		wP=1.0-r/R;
		// Do the -1/2 defect on the right side
		QM[0] = ((double)XYZ[0])*0.75;
		QM[1] = ((double)XYZ[1])*0.5;
		phiM = atan2(-QM[1]+Q[1],-QM[0]+Q[0]);
		UM[0] = cos(kM*phiM);
		UM[1] = sin(kM*phiM);
		r=sqrt( pow(QM[0]-Q[0],2.0)+pow(QM[1]-Q[1],2.0) );
		wM=1.0-r/R;
		if(Q[1]<=0.5*((double)XYZ[1])) wY=1.0-2.0*Q[1]/(double)XYZ[1];
		else wY=1.0-2.0*((double)XYZ[1]-Q[1])/(double)XYZ[1];

		// Set the director field with a linear weighting from each
		U[0] = wP*UP[0] + wM*UM[0] + wY;
		U[1] = wP*UP[1] + wM*UM[1];
		if( DIM>_2D ) U[2] = 0.0;
		norm( U,DIM );
	}

	else{
		printf( "Error: Particle orientation distribution unacceptable.\n" );
		exit( 1 );
	}
}

///
/// @brief Function that checks if an MPCD particle is within an BC.
///
/// This subroutine checks if a particle is within an obstacle (or BC).
/// If that is the case, then the particle's position is shifted to a neighboring SRD cell.
///
/// @param i Index of the particle whose position is being checked.
/// @param pp Return pointer to the first particle of the particle array.
/// @param SP Array of all species.
/// @param WALL Array of all boundary conditions.
/// @return If the position of the particle had to be shifted, the function returns i-1, so
///         the test can be checked again. If its not being shifted it returns i, so we can progress
///         to the next particle.
///
int checkplaceMPC( int i,particleMPC *pp,spec SP[],bc WALL[] ) {
	double shift[_3D];
	int j,k;

	for( j=0; j<NBC; j++ ) {
		//Zero the shift vector (just in case)
		for( k=0; k<_3D; k++ ) shift[k] = 0.;
		shiftBC( shift,&WALL[j],(pp+i) );
		rotateBC( &WALL[j],(pp+i),0 );
        WALL[j].W = calcW( WALL[j],*(pp+i) );
		rotatebackBC( &WALL[j],(pp+i),0 );
		shiftbackBC( shift,&WALL[j] );
		//If W<=0 then the particleMPC is inside an obstacle and must be replaced
		if( WALL[j].W <= TOL ) {
			replace( (pp+i) );
			//push( (pp+i)->V,KBT,SP[(pp+i)->SPID].VDIST, SP[(pp+i)->SPID].M,NULL );
			i--;
			return i;
		}
	}
	return i;
}

///
/// @brief Function that checks if an MPCD particle is within an BC.
///
/// This subroutine checks if a particle is within an obstacle (or BC).
/// If that is the case, then the particle's position is shifted.  The process
/// is repeated until the particle is no longer within an obstacle.
///
/// @param pp Return pointer to the MPCD particle whose position is being checked.
/// @param WALL Array of all boundary conditions (obstacles).
///
void replacePos_WithCheck( particleMPC *pp,bc WALL[] ) {
	double shift[_3D];
	int j,k;
	int flag=1;

	while( flag ) {
		replace( pp );
		//push( pp->V,KBT,SP[pp->SPID].VDIST, SP[pp->SPID].M,NULL );
		flag=0;
		for( j=0; j<NBC; j++ ) {
			//Zero the shift vector (just in case)
			for( k=0; k<_3D; k++ ) shift[k] = 0.;
			shiftBC( shift,&WALL[j],pp );
			rotateBC( &WALL[j],pp,0 );
			WALL[j].W = calcW( WALL[j],*pp );
			rotatebackBC( &WALL[j],pp,0 );
			shiftbackBC( shift,&WALL[j] );
			//If W<=0 then the particleMPC is inside an obstacle and must be replaced
			if( WALL[j].W <= TOL ) flag=1;
		}
	}
}

///
/// @brief  Function that checks if an MPCD particles is within an obstacle.
///
/// This function checks if an MPCD is within an obstacle. If they are,
/// it shifts their position until they are not. The position of the MPCD particles are checked
/// and shifted using checkplaceMPC.  If MPC in MD mode is being run, in addition we check that
/// the MPCD particle is not too close (considering periodic boundary conditions) to the MD particles
/// (as set by 1.25*rCut). If they are, the MPCD particle is shifted until its not.
///
/// @param IN Index of the MPCD particle being checked.
/// @param pp Return pointer to the first particle of the MPCD particle array.
/// @param SP Array of all species.
/// @param WALL Array of boundary conditions (obstacles).
/// @param simMD A pointer to the MD simulation.
/// @param KBT Temperature.
/// @param MD_mode Integer describing the MD simulation mode.
/// @see checkplaceMPC()
/// @return If the particle had to beshifted, it returns IN-1, so the check can be performed again.
///         If its not, it returns IN, so the next particle can be checked.
///
int checkplace( int IN,particleMPC *pp,spec SP[],bc WALL[],simptr simMD,double KBT,int MD_mode ) {
	particleMD *atom;
	double d[_3D];
	int i=IN,j;

	i=checkplaceMPC( i,pp,SP,WALL );

	if(MD_mode == MPCinMD ) {
		atom = simMD->atom.items;
		for( j=0; j<simMD->atom.n; j++ ) {
			d[0] = (pp+i)->Q[0] - (atom+j)->rx;
			d[1] = (pp+i)->Q[1] - (atom+j)->ry;
			d[2] = (pp+i)->Q[2] - (atom+j)->rz;

			//Owen made hard code it in instead of call the routine I made. I'm mopey
			if (simMD->pbcond & PBC_COND_x) {
				if (d[0] >= 0.5*XYZ[0]) d[0] -= XYZ[0];
				else if	(d[0] < -0.5*XYZ[0]) d[0] += XYZ[0];
			}
			if (simMD->pbcond & PBC_COND_y) {
				if (d[1] >= 0.5*XYZ[1]) d[1] -= XYZ[1];
				else if	(d[1] < -0.5*XYZ[1]) d[1] += XYZ[1];
			}
			if (simMD->pbcond & PBC_COND_z) {
				if (d[2] >= 0.5*XYZ[2]) d[2] -= XYZ[2];
				else if	(d[2] < -0.5*XYZ[2]) d[2] += XYZ[2];
			}
			if( sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]) <= (simMD->rCut)*1.25 ) {
				replace( (pp+i) );
				//push( (pp+i)->V,KBT,SP[(pp+i)->SPID].VDIST, SP[(pp+i)->SPID].M,NULL );
				i--;
				return i;
			}
		}
	}
	return i;
}

///
/// @brief Initialization of the MPCD particle's position, velocity and orientation.
///
/// This subroutine initializes the MPCD particle's position, velocity and orientations by calling
/// place(), push() and orient(). Positions and velocities can be read from a file.
/// After having their position set, they are all checked for conflicts with obstacles using checkplace().
/// Finally, it flags all particles as streaming.
///
/// @param dir Directory of the position files.
/// @param SP Array of all species.
/// @param pp Return pointer of the fist MPCD particle in the MPCD particle array.
/// @param KBT Temperature.
/// @param AVVEL Return pointer to average velocity.
/// @param WALL Array of all boundary conditions (obstacles).
/// @param simMD A pointer to the MD simulation.
/// @param MD_mode Integer describing the MD simulation mode.
/// @param LC Integer determining if the species is isotropic (i.e., not a liquid crystal). 0 means isotropic.
/// @see place()
/// @see push()
/// @see return()
/// @see checkplace()
///
void setcoord(char dir[], spec SP[], particleMPC *pp, double KBT, double AVVEL[], bc WALL[], simptr simMD, int MD_mode, int LC ) {
	int i,j,k=0;
	FILE *fin[NSPECI];
	char fileprefix[] = "placeSP";
	char filesuffix[16] = "0000";
	char fileextension[] = ".inp";

	for( i=0; i<NSPECI; i++ ) for( j=0; j<SP[i].POP; j++ ) {
		(pp+k)->SPID = i;
		k++;
	}
	for( j=0; j<DIM; j++ ) AVVEL[j] = 0.0;
	//Open files incase the particleMPCs' coordinates are being read in
	//Initialize the detailed output files
	for( i=0; i<NSPECI; i++ ) if( SP[i].POP > 0 && SP[i].QDIST == READ ) openplace( i,fin,dir,fileprefix,filesuffix,fileextension );
	//Set  particleMPC coordinates
	for( i=0; i<GPOP; i++ ) {
		//zero (dimensions not used will stay zeroed)
		for( j=0; j<_3D; j++ ) {
			(pp+i)->Q[j] = 0.0;
			(pp+i)->V[j] = 0.0;
			(pp+i)->U[j] = 0.0;
		}
		//Set particleMPC position, velocity and orientation
		place( (pp+i)->Q,SP[(pp+i)->SPID].QDIST, fin[(pp+i)->SPID] );

		push( (pp+i)->V,KBT,SP[(pp+i)->SPID].VDIST, SP[(pp+i)->SPID].MASS,fin[(pp+i)->SPID] );
		//Shift first mode of the velocity dist by the average velocity (push() centres about zero)
		for( j=0; j<DIM; j++ ) AVVEL[j] += (pp+i)->V[j];

		if( LC>ISOF ) orient( (pp+i)->U,(pp+i)->Q,SP[(pp+i)->SPID].ODIST );
	}
	for( j=0; j<DIM; j++ ) AVVEL[j] /= (double)GPOP;
	//Check particleMPC coordinates
	for( i=0; i<GPOP; i++ ) i = checkplace(i, pp, SP, WALL, simMD, KBT, MD_mode );
	//They will all need to stream
	for( i=0; i<GPOP; i++ ) (pp+i)->S_flag = STREAM;

	//Close the input files
	for( i=0; i<NSPECI; i++ ) if( SP[i].POP > 0 && SP[i].QDIST == READ ) fclose( fin[i] );
}

/// @brief Routine that checks for odd input.
///
/// This routine goes over various inputs, checking they are sound. In particular, it checks
/// the thermostat, the number of species, the number of boundary conditions, dimensionality,
/// the sensibility of the BCs, a proper values for flags setting up: i) liquid crystal (LC parameter)
/// ii) abscence of hydrodynamic interactions (noHI parameter), iii) incompressibility (inCOMP parameter),
/// iv) multiple phases (MULTIPHASE).  It also checks for a non-zero value of friction for nematogens, for a
/// non-zero mean field potential, if modelling a liquid crystal, and that the parameters for swimmers
/// are sensible.
///
/// @param fsynopsis Return point to the synopsis file.
/// @param SYNOUT Integer indicating if a synopsis file is needed.
/// @param in List of inputs.
/// @param SP Pointer to the first species of the array of all species.
/// @param WALL Pointer to the first boundary condition (BC) of the array of all BCs.
/// @param SS Species of swimmer.
///
void checkSim( FILE *fsynopsis,int SYNOUT,inputList in,spec *SP,bc *WALL,specSwimmer SS ) {
	int i,j;

	//Check thermostat
	if( in.TSTECH==BEREND && in.TAU<0.5*in.dt ) {
		#ifdef DBG
			if( DBUG >= DBGWARN ) printf( "Error: TAU = %lf. For Berendsen thermostat TAU must be close to 1*dt.\n",in.TAU );
		#endif
		if(SYNOUT == OUT) fprintf(fsynopsis,"Warning: TAU = %lf. For Berendsen thermostat TAU must be close to 1*dt.\n",in.TAU );
		exit(1);
	}
	if( in.TSTECH == HEYES && (in.TAU<0.05 || in.TAU>0.3) ) {
		#ifdef DBG
			if( DBUG >= DBGWARN ) printf( "Error: TAU = %lf. For Heyes thermostat must have 0.05<TAU/dt<0.3.\n",in.TAU );
		#endif
		if(SYNOUT == OUT) fprintf(fsynopsis,"Warning: TAU = %lf. For Heyes thermostat must have 0.05<TAU/dt<0.3.\n",in.TAU );
		exit(1);
	}
	if( ( in.TSTECH==VSC || in.TSTECH==BEREND ) && ( in.RTECH==MPCAT || in.RTECH==RAT ) ) {
		#ifdef DBG
			if( DBUG >= DBGWARN ) printf( "Error: RTECH = %d but TSTECH = %d. For Andersen MPCD, the thermostat must be turned off or apply only to the centre of mass.\n",in.RTECH,in.TSTECH );
		#endif
		if(SYNOUT == OUT) fprintf(fsynopsis,"Warning: RTECH = %d but TSTECH = %d. For Andersen MPCD, the thermostat must be turned off or apply only to the centre of mass.\n",in.RTECH,in.TSTECH );
		exit(1);
	}
	if( in.TSTECH==HEYES && ( in.RTECH==MPCAT || in.RTECH==RAT ) ) {
		#ifdef DBG
			if( DBUG >= DBGWARN ) printf( "Warning: RTECH = %d but TSTECH = %d. For Andersen MPCD, the including the Heyes thermostat is redundant.\n",in.RTECH,in.TSTECH );
		#endif
		if(SYNOUT == OUT) fprintf(fsynopsis,"Warning: RTECH = %d but TSTECH = %d. For Andersen MPCD, the including the Heyes thermostat is redundant.\n",in.RTECH,in.TSTECH );
	}
	if( in.RTECH==VICSEK || in.RTECH==CHATE ) for( i=0; i<NSPECI; i++ ) if( SP[i].ACT>1.0 || SP[i].ACT<0.0 ) {
		#ifdef DBG
			if( DBUG >= DBGWARN ) printf( "Error: RTECH = %d (Vicsek or Chate algorithm) but ACT=%lf (the noise) is not within [0,1].\n",in.RTECH,SP[i].ACT );
		#endif
		if(SYNOUT == OUT) fprintf(fsynopsis,"Error: RTECH = %d (Vicsek or Chate algorithm) but ACT=%lf (the noise) is not within [0,1].\n",in.RTECH,SP[i].ACT );
		exit(1);
	}
	if( in.TSTECH==MAXV && in.RTECH!=MPCAT ) {
		#ifdef DBG
			if( DBUG >= DBGWARN ) printf( "Warning: RTECH = %d but TSTECH = %d. This thermostat is a maximum velocity controller designed to work with Andersen MPCD.\n",in.RTECH,in.TSTECH );
		#endif
		if(SYNOUT == OUT) fprintf(fsynopsis,"Warning: RTECH = %d but TSTECH = %d. This thermostat is a maximum velocity controller designed to work with Andersen MPCD.\n",in.RTECH,in.TSTECH );
	}
	// Check number of species
	if( NSPECI > MAXSPECI ) {
		printf("Error: The number of species cannot exceed %d unless definitions.h is altered to increase MAXSPECI.\n", MAXSPECI);
		if(SYNOUT == OUT) fprintf(fsynopsis,"Error: The number of species cannot exceed %d unless definitions.h is altered to increase MAXSPECI.\n", MAXSPECI);
		exit(1);
	}
	// Check number of BC
	if( NBC > MAXBC ) {
		printf("Error: The number of BCs cannot exceed %d unless definitions.h is altered to increase MAXBC.\n", MAXBC);
		if(SYNOUT == OUT) fprintf(fsynopsis,"Error: The number of BCs cannot exceed %d unless definitions.h is altered to increase MAXBC.\n", MAXBC);
		exit(1);
	}
	// Check dimensionality
	if( DIM == _2D ) if ( XYZ[2] != 1 ) {
		printf("Error: In 2D, the simulation must be in the z-plane. Set DZ=1\n");
		if(SYNOUT == OUT) fprintf(fsynopsis,"Error: In 2D, the simulation must be in the z-plane. Set DZ=1\n");
		exit(1);
	}
	if( DIM > _3D ) {
		printf("Error: Hyper-dimensional solvent simulations not supported.\n");
		if(SYNOUT == OUT) fprintf(fsynopsis,"Error: In 2D, the simulation must be in the z-plane. Set DZ=1\n");
		exit(1);
	}
	// Initialize BC
	//Check that BCs make sense
	for( i=0; i<NBC; i++ ) {
		if( WALL[i].DSPLC==0 ) for( j=0; j<DIM; j++ ) {
			if( fneq(WALL[i].V[j],0.0) ) {
				printf( "Warning:\tBC %d immobile; setting velocity to zero.\n",i );
				if(SYNOUT == OUT) fprintf(fsynopsis,"Warning:\tBC %d immobile; setting velocity to zero.\n",i );
				WALL[i].V[j]=0.0;
			}
			if( fneq(WALL[i].L[j],0.0) ) {
				printf( "Warning:\tBC %d immobile; setting angular velocity to zero.\n",i );
				if(SYNOUT == OUT) fprintf(fsynopsis,"Warning:\tBC %d immobile; setting angular velocity to zero.\n",i );
				WALL[i].L[j]=0.0;
			}
			if( fneq(WALL[i].G[j],0.0) ) {
				printf( "Warning:\tBC %d immobile; setting acceleration to zero.\n",i );
				if(SYNOUT == OUT) fprintf(fsynopsis,"Warning:\tBC %d immobile; setting acceleration to zero.\n",i );
				WALL[i].V[j]=0.0;
			}
		}
		if( WALL[i].ABS==0 ) {
			// The user should likely use the abs() operator for BCs that aren't a plane or don't have even powers
			for( j=0; j<DIM; j++ ) {
				if( fneq(WALL[i].P[j],1.0) ) if( fneq( fmod(WALL[i].P[j],2.0),0.0 ) ) {
					printf( "Warning:\tBC %d does not use abs() but likely should when P[%d]=%lf is not even or unity.\n",i,j,WALL[i].P[j] );
					printf( "\tIf the user means to create a closed obstacle or particle this is likely to generate errors.\n" );
					printf( "\t***However***, this is left to the user's discretion.\n" );
					if(SYNOUT == OUT) {
						fprintf(fsynopsis,"Warning:\tBC %d does not use abs() but likely should when P[%d]=%lf is not even or unity.\n",i,j,WALL[i].P[j] );
						fprintf( fsynopsis,"\tIf the user means to create a closed obstacle or particle this is likely to generate errors.\n" );
						fprintf( fsynopsis,"\t***However***, this is left to the user's discretion.\n" );
					}
				}
			}
		}
	}
	if( !( in.LC==ISOF || in.LC==LCL || in.LC==LCG) ){
		printf( "Error: Unrecognized value of LC=%d.\n",in.LC );
		exit( 1 );
	}
	if( in.LC==LCG ) for( i=0; i<NSPECI; i++ ) if( SP[i].ODIST==RANDORIENT ) {
		printf( "Warning: Using global S (LC=%d) but initiated in isotropic phase. Simulation may not reach nematic phase.\n",in.LC );
		if(SYNOUT == OUT) fprintf( fsynopsis,"Warning: Using global S (LC=%d) but initiated in isotropic phase. Simulation may not reach nematic phase.\n",in.LC );
	}
	if( !( in.noHI==HION || in.noHI==HIOFF) ){
		printf( "Error: Unrecognized value of noHI=%d.\n",in.noHI );
		exit( 1 );
	}
	if( !( in.inCOMP==INCOMPOFF || in.inCOMP==INCOMPSWAP || in.inCOMP==INCOMPVIRIAL || in.inCOMP==INCOMPSUB) ){
		printf( "Error: Unrecognized value of inCOMP=%d.\n",in.inCOMP );
		exit( 1 );
	}
	if( !( in.MULTIPHASE==MPHOFF || in.MULTIPHASE==MPHPOINT || in.MULTIPHASE==MPHSURF ) ){
		printf( "Error: Unrecognized value of MULTIPHASE=%d.\n",in.MULTIPHASE );
		exit( 1 );
	}
	if( ( in.MULTIPHASE==MPHOFF && NSPECI>1 ) ){
		printf( "Warning: MULTIPHASE=%d (off) but more than one species present (NSPECI=%d).\n",in.MULTIPHASE,NSPECI );
	}
	//Check that nematogens have non-zero friction
	if( in.LC>ISOF ) for( i=0; i<NSPECI; i++ ) if( feq(SP[i].RFC,0.0) ) {
		printf( "Warning:\tSpecies %d has zero rotational friction coefficient\n",i );
		printf( "\t\tThis is acceptable iff no magnetic field is applied\n");
		if(SYNOUT == OUT) fprintf( fsynopsis,"Warning:\tSpecies %d has zero rotational friction coefficient\n\t\tThis is acceptable iff no magnetic field is applied\n",i );
	}
	//Check that if running as a LC that the LC mean field potential MFPOT is greater than zero
	if( in.LC>ISOF && in.MFPOT<=0.0 ) {
		#ifdef DBG
			if( DBUG >= DBGWARN ) printf( "Error: Running as nematic liquid crystal (LC=%d) but mean field = %lf. Must be greater than 0 or run as isotropic\n",in.LC,in.MFPOT );
		#endif
		if(SYNOUT == OUT) fprintf(fsynopsis,"Error: Running as nematic liquid crystal (LC=%d) but mean field = %lf. Must be greater than 0 or run as isotropic\n",in.LC,in.MFPOT );
		exit(1);
	}
// 	if( in.LC==LCG || in.LC==LCL ) for( i=0; i<NSPECI; i++ ) if( SP[i].RFC*SP[i].TUMBLE>1.0 ) {
// 		printf( "Warning: Tumbling parameter (%lf) times rotational friction coefficient (%lf) greater than 1. Simulation likely to crash.\n",SP[i].TUMBLE,SP[i].RFC );
// 		if(SYNOUT == OUT) fprintf( fsynopsis,"Warning: Tumbling parameter (%lf) times rotational diffusion coefficient (%lf) greater than 1. Simulation likely to crash.\n",SP[i].TUMBLE,SP[i].RFC );
// 	}
	if( in.LC!=ISOF ) for( i=0; i<NSPECI; i++ ) if( SP[i].ODIST==ALIGNZ && DIM<_3D ) {
		#ifdef DBG
			if( DBUG >= DBGWARN ) printf( "Error: \tSpecies %d orientations initialized along z-axis but not 3D\n",i );
		#endif
		if(SYNOUT == OUT) fprintf(fsynopsis,"Error: \tSpecies %d orientations initialized along z-axis but not 3D\n",i );
		exit(1);
	}
	if( in.LC!=ISOF ) if( in.RTECH==CHATE || in.RTECH==CHATE_MPCAT || in.RTECH==CHATE_LANG || in.RTECH==DIPOLE_VCM ) {
		printf( "Error: RTECH is a Chate-style active nematic but LC should be %d (LC=%d and RTECH=%d).\n",ISOF,in.LC,in.RTECH );
		if(SYNOUT == OUT) fprintf( fsynopsis,"Error: RTECH is a Chate-style active nematic but LC should be %d (LC=%d and RTECH=%d).\n",ISOF,in.LC,in.RTECH );
		exit( 1 );
	}
	if( (in.RTECH==DIPOLE_DIR_SUM || in.RTECH==DIPOLE_DIR_AV) && in.LC==ISOF ) {
		printf( "Error: RTECH is a dipole force along the director (RTECH=%d) but LC is for an isotropic fluid (LC=%d).\n",in.RTECH,in.LC );
		if(SYNOUT == OUT) fprintf( fsynopsis, "Error: RTECH is a dipole force along the director (RTECH=%d) but LC is for an isotropic fluid (LC=%d).\n",in.RTECH,in.LC );
		exit( 1 );
	}
	// Check binary fluid
	// if( in.binaryDELTA > 1.0 || in.binaryDELTA < 0.0 ) {
	// 	printf("Warning: The binary fluid fractional energy difference (binaryDELTA) should be between zero and unity.\n");
	// 	if(SYNOUT == OUT) fprintf(fsynopsis,"Warning: The binary fluid fractional energy difference (binaryDELTA) should be between zero and unity.\n");
	// 	exit(1);
	// }
	// Check swimmers
	if( SS.sizeShrink>1.0 ) {
		printf("Warning: The swimmers grow (smaller rotational diffusion) during tumble phase rather than shrink (sizeShrink=%lf).\n",SS.sizeShrink);
		if(SYNOUT == OUT) fprintf(fsynopsis,"Warning: The swimmers grow (smaller rotational diffusion) during tumble phase rather than shrink (sizeShrink=%lf).\n",SS.sizeShrink);
	}
	if( feq(SS.DS,0.0) ) {
		printf("Warning: The swimmers' dipole strength is zero so they do not geneate dipolar flow.\n");
		if(SYNOUT == OUT) fprintf(fsynopsis,"Warning: The swimmers' dipole strength is zero so they do not geneate dipolar flow.\n");
	}
	if( !feq(fabs(SS.TS),0.0) && DIM<_3D ) {
		printf("Warning: The swimmers' rotlet dipole strength is non-zero but the simulation is 2D so no rotlet possible.\n");
		if(SYNOUT == OUT) fprintf(fsynopsis,"Warning: The swimmers' rotlet dipole strength is non-zero but the simulation is 2D so no rotlet possible.\n");
	}
	if( SS.runTime<0.0 || SS.tumbleTime<0.0 ) {
		printf("Error: The run or tumble time must be greater than zero (or zero for no run/tumble dynamics).\n");
		if(SYNOUT == OUT) fprintf(fsynopsis,"Error: The run or tumble time must be greater than zero (or zero for no run/tumble dynamics).\n");
		exit( 1 );
	}
}

///
/// @brief Function that initializes output files as requested by the input file.
///
/// This subroutines initializes the output files that are requested by the input file. It does so
/// by checking for the corresponding flags and the open methods, e.g., openflow().
///
/// @param op Path to the output directory.
/// @param outFlag Structute listing output flags with values read (obtaind) from the input file.
/// @param outFile Structure listing outpout files.
/// @param in Structure containing the input lists, read from the input file.
/// @param SP Pointer of the first species in the array of all species.
/// @param WALL Array of all boundary conditions (obstacles).
/// @see openflow()
///
void initOutput( char op[],outputFlagsList *outFlag,outputFilesList *outFile,inputList in,spec *SP, bc WALL[] ) {

	int i;
	char filecoarse[]="coarsegrain";
	char fileavvel[]="avVel";
	char fileorder[]="directorfield";
	char fileorderQ[]="ordertensor";
	char fileorderQK[]="recipOrder";
	char fileavs[]="avS";
	char filedensSTD[]="densSTD";
	char fileenstrophy[]="avEnstrophy";
	char fileflow[]="flowfield";
	char filevel[]="velfield";
	char filesolids[]="solidtraj";
	char filetopo[]="topochargefield";
	char filedefect[]="defects";
	char filedisclin[]="disclinTensorfield";
	char fileprefix[]="detailedSP";
	char filehistVel[]="distVel";
	char filehistVort[]="distVort";
	char filehistDens[]="distDens";
	char filehistS[]="distS";
	char filehistDir[]="distDir";
	char filehistSpeed[]="distSpeed";
	char filehistEnstr[]="distEnstrophy";
	char fileenergy[]="energy";
	char fileenergyfield[]="enfield";
	char fileenneighbours[]="enneighbours";
	char filecorrVV[]="corrVelVel";
	char filecorrNN[]="corrDirDir";
	char filecorrDD[]="corrDensDens";
	char filecorrSS[]="corrOrderOrder";
	char filecorrPP[]="corrPhiPhi";
	char filecorrWW[]="corrVortVort";
	char fileenergyspect[]="energySpectrum";
	char fileenstrophyspect[]="enstrophySpectrum";
	char fileBinder[]="binderCumulant";
	char fileswimmer[]="swimmers";
	char fileswimmerori[]="swimmersOri";
	char fileruntumble[]="runtumble";
	char filemultiphase[]="multiphase";
	char filepressure[]="pressure";

	char filesuffix[]="0000";
	char fileextension[]=".dat";

	#ifdef DBG
		if( DBUG >= DBGINIT ) printf("Initialization\n");
		if( DBUG >= DBGINIT ) printf("Initialize Output Files\n");
	#endif
	//Don't bother with LC stuff if its not being used
	if(in.LC==ISOF) {
		outFlag->ORDEROUT=0;
		outFlag->QTENSOUT=0;
		outFlag->QKOUT=0;
		outFlag->CNNOUT=0;
		outFlag->CSSOUT=0;
	}
	//Initialize the detailed output files
	if( (outFlag->TRAJOUT)>=OUT ) for(i=0;i<NSPECI;i++) if(SP[i].POP>=1) opendetails( i,outFile->fdetail,op,fileprefix,filesuffix,fileextension );
	//Initialize the course grained output file
	if( (outFlag->COAROUT)>=OUT ) opencoarse( &(outFile->fcoarse),op,filecoarse,fileextension );
	//d
	if( (outFlag->AVVELOUT)>=OUT ) openavvel( &(outFile->favvel),op,fileavvel,fileextension );
	//Initialize the director output file
	if( (outFlag->ORDEROUT)>=OUT ) openorder( &(outFile->forder),op,fileorder,fileextension );
	//Initialize the tensor order parameter output file
	if( (outFlag->QTENSOUT)>=OUT ) openorderQ( &(outFile->forderQ),op,fileorderQ,fileextension );
	if( (outFlag->QKOUT)>=OUT ) openorderQK( &(outFile->forderQK),op,fileorderQK,fileextension );
	//Initialize the scalar order parameter output file
	if( (outFlag->AVSOUT)>=OUT ) openavs( &(outFile->favs),op,fileavs,fileextension );
	//Initialize the scalar order parameter output file
	if( (outFlag->DENSOUT)>=OUT ) opendensSTD( &(outFile->fdensSTD),op,filedensSTD,fileextension );
	//Initialize the enstrophy output file
	if( (outFlag->ENSTROPHYOUT)>=OUT ) openavenstrophy( &(outFile->fenstrophy),op,fileenstrophy,fileextension );
	//Initialize the flow and velocity field output files
	if( (outFlag->FLOWOUT)>=OUT ) openflow( &(outFile->fflow),op,fileflow,fileextension );
	if( (outFlag->VELOUT)>=OUT ) openvel( &(outFile->fvel),op,filevel,fileextension );
	//Initialize the distribution output files
	if( (outFlag->HISTVELOUT)>=OUT ) openhistVel( &(outFile->fhistVel),op,filehistVel,fileextension );
	if( (outFlag->HISTSPEEDOUT)>=OUT ) openhistSpeed( &(outFile->fhistSpeed),op,filehistSpeed,fileextension );
	if( (outFlag->HISTVORTOUT)>=OUT ) openhistVort( &(outFile->fhistVort),op,filehistVort,fileextension );
	if( (outFlag->HISTENSTROUT)>=OUT ) openhistEnstrophy( &(outFile->fhistEnstr),op,filehistEnstr,fileextension );
	if( (outFlag->HISTDIROUT)>=OUT ) openhistDir( &(outFile->fhistDir),op,filehistDir,fileextension );
	if( (outFlag->HISTSOUT)>=OUT ) openhistS( &(outFile->fhistS),op,filehistS,fileextension );
	if( (outFlag->HISTNOUT)>=OUT ) openhistDens( &(outFile->fhistDens),op,filehistDens,fileextension );
	//Initialize the energy output files
	if( (outFlag->ENOUT)>=OUT ) openenergy( &(outFile->fenergy),op,fileenergy,fileextension );
	//Initialize the energy field output files
	if( (outFlag->ENFIELDOUT)>=OUT ) openenergyfield( &(outFile->fenergyfield),op,fileenergyfield,fileextension );
	//Initialize the energy from neighbours output files
	if( (outFlag->ENNEIGHBOURS)>=OUT ) openenergyneighbours( &(outFile->fenneighbours),op,fileenneighbours,fileextension );
	//Initialize the spatial correlation output files
	if( (outFlag->CVVOUT)>=OUT ) opencorr( &(outFile->fcorrVV),op,filecorrVV,fileextension );
	if( (outFlag->CWWOUT)>=OUT ) opencorr( &(outFile->fcorrWW),op,filecorrWW,fileextension );
	if( (outFlag->CDDOUT)>=OUT ) opencorr( &(outFile->fcorrDD),op,filecorrDD,fileextension );
	if( (outFlag->CPPOUT)>=OUT ) opencorr( &(outFile->fcorrPP),op,filecorrPP,fileextension );
	if( (outFlag->CNNOUT)>=OUT ) opencorr( &(outFile->fcorrNN),op,filecorrNN,fileextension );
	if( (outFlag->CSSOUT)>=OUT ) opencorr( &(outFile->fcorrSS),op,filecorrSS,fileextension );
	//Initialize the spectrum output files
	if( (outFlag->ENERGYSPECTOUT)>=OUT ) openenergyspect( &(outFile->fenergyspect),op,fileenergyspect,fileextension );
	if( (outFlag->ENSTROPHYSPECTOUT)>=OUT ) openenstrophyspect( &(outFile->fenstrophyspect),op,fileenstrophyspect,fileextension );
	//Initialize the Binder cumulant output files
	if( (outFlag->BINDER)>=OUT ) openbinder( &(outFile->fbinder),op,fileBinder,fileextension,outFlag->BINDERBIN );
	//Initialize swimmer output files
	if( (outFlag->SWOUT)>=OUT ) openswimmer( &(outFile->fswimmers),op,fileswimmer,fileextension );
	//Initialize swimmer output files
	if( (outFlag->SWORIOUT)>=OUT ) openswimmerOri( &(outFile->fswimmersOri),op,fileswimmerori,fileextension );
	//Initialize swimmer run/tumble output files
	if( (outFlag->RTOUT)>=OUT ) openruntumble( &(outFile->fruntumble),op,fileruntumble,fileextension );
	//Initialize the synopsis output files
	if( (outFlag->SYNOUT)>=OUT ) opensynopsis( &(outFile->fsynopsis),op,1 );
	//Initialize the solids' trajectories (or BC motion) output files
	if( (outFlag->SOLOUT)>=OUT ) for ( i=0; i<NBC; i++ ) if ( WALL[i].DSPLC ) opentraj(i,outFile->fsolids,op,filesolids,filesuffix,fileextension);
	//Initialize the topological charge output file
	if( (outFlag->TOPOOUT)>=OUT ) opentopo( &(outFile->ftopo),op,filetopo,fileextension );
	//Initialize the defect trajectories output file
	if( (outFlag->DEFECTOUT)>=OUT ) opendefect( &(outFile->fdefects),op,filedefect,fileextension );
	//Initialize the disclination tensor output file
	if( (outFlag->DISCLINOUT)>=OUT ) opendisclin( &(outFile->fdisclination),op,filedisclin,fileextension );
	//Initialize the phi/color/species-type field output file
	if( (outFlag->SPOUT)>=OUT ) openmultiphase( &(outFile->fmultiphase),op,filemultiphase,fileextension );
	//Initialize the pressure field output file
	if( (outFlag->PRESOUT)>=OUT ) openpressure( &(outFile->fpressure),op,filepressure,fileextension );

	if( (outFlag->SYNOUT)==OUT) {
		fprintf(outFile->fsynopsis,"Output:\n" );
		fprintf(outFile->fsynopsis,"\tTrajectories:\t%d\n",outFlag->TRAJOUT);
		fprintf(outFile->fsynopsis,"\tCoarse-Grained Flow:\t%d\n",outFlag->COAROUT);
		fprintf(outFile->fsynopsis,"\tGlobal Average velocity:\t%d\n",outFlag->AVVELOUT);
		fprintf(outFile->fsynopsis,"\tFlow:\t\t%d\n",outFlag->FLOWOUT);
		fprintf(outFile->fsynopsis,"\tPrint distributions:\n");
		fprintf(outFile->fsynopsis,"\t\tVel: %d\n\t\tSpeed: %d\n\t\tVorticity: %d\n\t\tEnstrophy: %d\n\t\tDirector: %d\n\t\tScalar order parameter: %d\n\t\tDensity: %d\n",outFlag->HISTVELOUT,outFlag->HISTSPEEDOUT,outFlag->HISTVORTOUT,outFlag->HISTENSTROUT,outFlag->HISTDIROUT,outFlag->HISTSOUT,outFlag->HISTNOUT);
		fprintf(outFile->fsynopsis,"\tGlobal average scalar order parameter:\t%d\n",outFlag->AVSOUT);
		fprintf(outFile->fsynopsis,"\tAverage enstrophy:\t%d\n",outFlag->ENSTROPHYOUT);
		fprintf(outFile->fsynopsis,"\tStandard deviation of density:\t%d\n",outFlag->DENSOUT);
		fprintf(outFile->fsynopsis,"\tDirector and scalar order parameter:\t%d\n",outFlag->ORDEROUT);
		fprintf(outFile->fsynopsis,"\tOrder parameter tensor:\t%d\n",outFlag->QTENSOUT);
		fprintf(outFile->fsynopsis,"\tOrder parameter tensor (reciprocal space):\t%d\n",outFlag->QKOUT);
		fprintf(outFile->fsynopsis,"\tEnergy:\t\t%d\n",outFlag->ENOUT);
		fprintf(outFile->fsynopsis,"\tEnergy field:\t\t%d\n",outFlag->ENFIELDOUT);
		fprintf(outFile->fsynopsis,"\tEnergy neighbours:\t\t%d\n",outFlag->ENNEIGHBOURS);
		if(DIM==_3D) {
			fprintf(outFile->fsynopsis,"\tTopological charge field:\t\t%d\tNot outputted in 3D!\n",outFlag->TOPOOUT);
			fprintf(outFile->fsynopsis,"\tDefect positions:\t\t%d\tNot outputted in 3D!\n",outFlag->DEFECTOUT);
		}
		else {
			fprintf(outFile->fsynopsis,"\tTopological charge field:\t\t%d\n",outFlag->TOPOOUT);
			fprintf(outFile->fsynopsis,"\tDefect positions:\t\t%d\n",outFlag->DEFECTOUT);
		}
		fprintf(outFile->fsynopsis,"\tDisclination tensor field:\t\t%d\n",outFlag->DISCLINOUT);
		fprintf(outFile->fsynopsis,"\tPhi/colour/species-type field:\t\t%d\n",outFlag->SPOUT);
		fprintf(outFile->fsynopsis,"\tPressure field:\t\t%d\n",outFlag->PRESOUT);
		fprintf(outFile->fsynopsis,"\tVelocity-velocity correlation:\t\t%d\n",outFlag->CVVOUT);
		fprintf(outFile->fsynopsis,"\tDirector-director correlation:\t\t%d\n",outFlag->CNNOUT);
		fprintf(outFile->fsynopsis,"\tVorticity-vorticity correlation:\t\t%d\n",outFlag->CWWOUT);
		fprintf(outFile->fsynopsis,"\tDensity-density correlation:\t\t%d\n",outFlag->CDDOUT);
		fprintf(outFile->fsynopsis,"\tOrder-order correlation:\t\t%d\n",outFlag->CSSOUT);
		fprintf(outFile->fsynopsis,"\tPhase-phase (binary fluid) correlation:\t\t%d\n",outFlag->CPPOUT);
		fprintf(outFile->fsynopsis,"\tBinder cumulant:\t\t%d --- bin size:\t\t%d\n",outFlag->BINDER,outFlag->BINDERBIN);
		fprintf(outFile->fsynopsis,"\tSolid BC:\t%d\n",outFlag->SOLOUT);
		fprintf(outFile->fsynopsis,"\tSwimmers:\t%d\n",outFlag->SWOUT);
		fprintf(outFile->fsynopsis,"Files initialized.\n" );
	}
}

///
/// @brief Function that initializes the simulation.
///
/// This subroutine initializes the simulation by calling all specific initializers, such as setcoord(), initvar(),
/// zerocnt(), etc..
///
/// @param CL Return pointer to array of all cell lists.
/// @param SRDparticles Return pointer to array of MPCD particles.
/// @param SP Array of all subspecies.
/// @param WALL Array of all boundary conditions (obstacles).
/// @param simMD A pointer to the MD simulation.
/// @param specS A pointer to the first element in the array of all swimming species.
/// @param swimmers A pointer to the first element in the array of all swimmers.
/// @param argc Number of terminal arguments (i.e., path to input and output).
/// @param argv Pointers to strings of terminal arguments (i.e., path to input and output).
/// @param in List of inputs.
/// @param to CPU time.
/// @param co Wall time.
/// @param runtime Return pointer to runtime.
/// @param warmtime Return pointer to warmtime.
/// @param AVVEL Return pointer to average velocity.
/// @param theory Return pointer to structure containing theoretical system parameters.
/// @param KBTNOW Return pointer to current temperature.
/// @param AVS Return pointer to average scalar order parameter.
/// @param S4 Return pointer to fourth moment of the scalar order parameter.
/// @param stdN Return pointer to density fluctuations.
/// @param AVNOW Return pointer to current average of flow velocity.
/// @param AVV Return pointer to past average flow velocity.
/// @param avDIR Return pointer to average director.
/// @param outFlags List of output flags.
/// @param MD_mode Integer specifying MD mode.
/// @param fsynopsis Synopsis file.
/// @param ip Path to input directory.
/// @see setcoord()
/// @see initvar()
/// @see zerocnt()
///
void initializeSIM(cell ***CL, particleMPC *SRDparticles, spec SP[], bc WALL[], simptr simMD, specSwimmer *specS, swimmer *swimmers, int argc, char* argv[], inputList *in, time_t *to, clock_t *co, int *runtime, int *warmtime, double *AVVEL, kinTheory *theory, double *KBTNOW, double *AVS, double *S4, double *stdN, double AVNOW[_3D], double AVV[_3D], double avDIR[_3D], outputFlagsList outFlags, int MD_mode, FILE *fsynopsis, char ip[] ) {
	int i,j;
	#ifdef DBG
		if( DBUG >= DBGINIT ) printf("\tInitialize Parameters\n");
	#endif
	initvar( &(in->seed),to,co,runtime,warmtime,&(theory->sumM),AVNOW,avDIR,SP,&(in->C),&(in->S),in->RA,AVVEL,in->KBT,WALL,in->MAG,CL,SRDparticles );

	maxXYZ=(int) sqrt( (double)(XYZ[0]*XYZ[0]+XYZ[1]*XYZ[1]+XYZ[2]*XYZ[2]) );

	in->GRAV_FLAG = 0;
	if( fneq(in->GRAV[0],0.0) || fneq(in->GRAV[1],0.0) || fneq(in->GRAV[2],0.0) ) in->GRAV_FLAG = 1;
	for( i=0; i<NBC; i++ ) if( fneq((WALL+i)->G[0],0.0) || fneq((WALL+i)->G[1],0.0) || fneq((WALL+i)->G[2],0.0) ) in->GRAV_FLAG = 1;
	// Flag whether or not to do BC re-orientations
	for( i=0; i<NBC; i++ ) {
		(WALL+i)->REORIENT = 1;
		// If the colloid cannot move and doesn't have an initial non-zero orientation, then never bother with rotations
		if( (WALL+i)->DSPLC==0 && feq((WALL+i)->O[0],0.0) && feq((WALL+i)->O[1],0.0) && feq((WALL+i)->O[2],0.0) ) (WALL+i)->REORIENT = 0;
		// If the colloid has circular symmetry, then never bother with rotations
		else if( DIM==_2D && feq((WALL+i)->P[0],2.0) && feq((WALL+i)->P[1],2.0) && feq((WALL+i)->A[0],(WALL+i)->A[1]) ) (WALL+i)->REORIENT = 0;
		// If the colloid has spherical symmetry, then never bother with rotations
		else if( feq((WALL+i)->P[0],2.0) && feq((WALL+i)->P[1],2.0) && feq((WALL+i)->P[2],2.0) && feq((WALL+i)->A[0],(WALL+i)->A[1]) && feq((WALL+i)->A[1],(WALL+i)->A[2]) ) (WALL+i)->REORIENT = 0;
	}
	in->MAG_FLAG = 0;
	if( fneq(in->MAG[0],0.0) || fneq(in->MAG[1],0.0) || fneq(in->MAG[2],0.0) ) in->MAG_FLAG = 1;

	if( in->TSTECH==MAXV ) for( j=0; j<DIM; j++ ) {
		AVV[j]=in->GRAV[j];
		in->GRAV[j]=0.0;
	}
	else for( j=0; j<DIM; j++ ) AVV[j]=0.0;

	#ifdef DBG
		if( DBUG >= DBGINIT ) printf( "\tInitialize System\n" );
		if( DBUG >= DBGINIT ) printf( "\tInitialize MPCD\n" );
	#endif
	//Intialize positions, velocity and angular velocity
	#ifdef DBG
		if( DBUG >= DBGINIT ) printf( "\tPlace MPCD particle\n" );
	#endif
	setcoord(ip, SP, SRDparticles, in->KBT, AVV, WALL, simMD, MD_mode, in->LC );
	//Intialize positions and orientations of swimmers
	if( NS>0 ) {
		#ifdef DBG
			if( DBUG >= DBGINIT ) printf( "\tPlace swimmers\n" );
		#endif
		setswimmers( specS,swimmers,WALL,in->stepsMD,in->dt );
	}

	if(outFlags.SYNOUT == OUT) fprintf(fsynopsis,"\nMPCD particles placed.\n" );
	//Calculate the theoretical properties of the SRD gas
	theory_trans( &(theory->MFP),&(theory->VISC),&(theory->THERMD),&(theory->SDIFF),&(theory->SPEEDOFSOUND),in->RA,in->FRICCO,in->KBT,in->dt,theory->sumM,in->RTECH,outFlags.SYNOUT,fsynopsis );
		//Zero counters
	zerocnt( KBTNOW,AVNOW,AVS );
	if( in->TSTECH==MAXV ) for( j=0; j<DIM; j++ ) AVNOW[j]=AVV[j];
	//Bin particles
	#ifdef DBG
		if( DBUG >= DBGINIT ) printf( "\tBin MPCD particles\n" );
	#endif
	binin( SRDparticles,CL );
	bininSwimmers( *specS,swimmers,CL );
	if(outFlags.SYNOUT == OUT) fprintf(fsynopsis,"\nMPCD particles binned for first time.\n" );
	if( MD_mode ) bininMD(simMD, CL );
	localPROP( CL,SP,*specS,in->RTECH,in->LC );
	*S4=0.;
	*stdN=0.;
	if( outFlags.AVSOUT>=OUT ) {
		*AVS = avOrderParam( SRDparticles,in->LC,avDIR );
		*S4 = avS4( SRDparticles,in->LC,avDIR );
		*stdN = stdNum( CL,GPOP,XYZ,XYZ_P1 );
		}
	/* ****************************************** */
	/* ******** GALILEAN TRANSFORMATION ********* */
	/* ****************************************** */
	*KBTNOW = TEMP( SRDparticles,SP,WALL,AVV );
	#ifdef DBG
		if( DBUG >= DBGINIT ) printf( "Actual System Temperature before moving to rest frame: %lf\n",*KBTNOW );
	#endif
	if( in->RFRAME==1 && GPOP<=2 ) {
		in->RFRAME=0;
		printf( "Simulation of 2 particles or less in the rest frame is boring\nDon't do Galilean Transformation\n" );
	}
	#ifdef DBG
		if( DBUG >= DBGINIT ) {
			printf( "Inputted System Temperature (units of KB): %lf\n",in->KBT );
			if( in->RFRAME==1) printf( "Galilean Transformation to System Rest Frame\n" );
		}
	#endif
	//Do the Galilean transformation of the system to its rest frame i.e. remove system's net momentum
	if( in->RFRAME ) galileantrans(SRDparticles, WALL, simMD, SP, in->KBT, AVV, GPOP, NBC, MD_mode, DIM );
	//Now that the initial shift is done we use RFRAME to signal when it should happen periodically (with zeroNetMom)
	//But don't want to do it if accelerating, duh
	if( !(in->zeroNetMom) ) in->RFRAME = 0;
	else if( in->RFRAME ) if( fneq(in->GRAV[0],0.0) || fneq(in->GRAV[1],0.0) || fneq(in->GRAV[2],0.0) ) in->RFRAME = 0;
	if(outFlags.SYNOUT == OUT) fprintf(fsynopsis,"\nGalilean transformation completed.\n" );
	/* ****************************************** */
	/* *********** TEMPERATURE SCALING ********** */
	/* ****************************************** */
	*KBTNOW = TEMP( SRDparticles,SP,WALL,AVV );
	#ifdef DBG
		if( DBUG >= DBGINIT ) {
			printf( "System Temperature before scaling: %lf\n",*KBTNOW );
			printf( "Apply thermostat %d.\n",in->TSTECH );
		}
	#endif
	//Scale temperature TSTECH=VSC here so that temperature moved instantaneously to initialized value
	scaleT( in->KBT,*KBTNOW,in->dt,in->TAU,AVV,AVNOW,VSC,SP,in->LC,WALL,SRDparticles,CL );
	*KBTNOW = TEMP( SRDparticles,SP,WALL,AVV );
	#ifdef DBG
		if( DBUG >= DBGINIT ) printf( "System Temperature after moving to rest frame and scaling: %lf\n",*KBTNOW );
	#endif
	if(outFlags.SYNOUT == OUT) fprintf(fsynopsis,"Temperature rescaled to input value.\n" );
	avVel( CL,AVNOW );
}

///
/// @brief Function that initializes the simulation out of a checkpoint.
///
/// This function initializes the recovery of a simulation out of a checkpoint. An MD simulation cannot be recovered.
///
/// @param CL Return pointer to array of all cells list.
/// @param SRDparticles Return pointer to first element in array of all MPCD particles.
/// @param SP Array of all particle subspecies.
/// @param specS Array of all species of swimmers.
/// @param RTECH Integer specifying rotation technique.
/// @param LC Integer specifying type of Liquid Crystal.
/// @param MD_mode Integer specifying type of MD simulation.
/// @param SYNOUT Integer specifying if a synopsis file is requires.
/// @param fsynopsis Synopsis file.
///
void initializeRecovery(cell ***CL, particleMPC *SRDparticles, spec SP[], specSwimmer specS, int RTECH, int LC, int MD_mode, int SYNOUT, FILE *fsynopsis ) {
	//int i;
	if(SYNOUT == OUT) fprintf(fsynopsis,"\nSimulation recovered from checkpoint.\n" );
	// MD isn't checkpointed so can't recover MD simulation
	if(MD_mode != noMD ) {
		printf("Error: Cannot recover MD simulation.\n");
		exit(EXIT_FAILURE);
	}
	//maxXYZ=0;
	//for( i=0;i<_3D;i++ ) if( XYZ[i]>maxXYZ ) maxXYZ=XYZ[i];
	maxXYZ=(int) sqrt( (double)(XYZ[0]*XYZ[0]+XYZ[1]*XYZ[1]+XYZ[2]*XYZ[2]) );
	//Bin particles
	#ifdef DBG
		if( DBUG >= DBGINIT ) printf( "\tBin MPCD particles\n" );
	#endif
	binin( SRDparticles,CL );
	if(SYNOUT == OUT) fprintf(fsynopsis,"\nMPCD particles binned for first time.\n" );
	localPROP( CL,SP,specS,RTECH,LC );
}
