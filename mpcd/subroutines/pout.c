///
/// @file
///
/// @brief Prints data output files for NIAMH-MPCD.
///
/// A collection of functions for constructing and printing different raw data outputs to .dat files. The types of files produced must be specified in input.json.
///

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "../headers/definitions.h"
# include "../headers/globals.h"
# include "../headers/SRDclss.h"
# include "../headers/therm.h"
# include "../headers/lc.h"
# include "../headers/mtools.h"
# include "../headers/ctools.h"
# include "../headers/mpc.h"
# include "../headers/init.h"
# include "../headers/swimmers.h"

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* **************** PRINTING **************** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */

///
/// @brief This function prints a descriptive header for all data output files.
///
/// The header is printed to all data output files, crediting the creator.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param SP This is the subpopulation of species.
///
void outheader( FILE *fout,int SP ) {
	fprintf( fout," **********************************************\n" );
	fprintf( fout," **********************************************\n" );
	fprintf( fout," **********************************************\n" );
	fprintf( fout," ***************** NIAMH-MPCD *****************\n" );
	fprintf( fout," **********************************************\n" );
	fprintf( fout," **************** Shendruk Lab ****************\n" );
	fprintf( fout," ********** University of Edinburgh ***********\n" );
	fprintf( fout," ************ Code begun Oct. 2008 ************\n" );
	fprintf( fout," **********************************************\n" );
	fprintf( fout," **********************************************\n" );
	fprintf( fout," **********************************************\n\n" );
}

///
/// @brief Prints column headers for data output files.
///
/// Column headers are produced for the relevant .dat files to display raw data in a table format.
/// Time, t, is the first column header.
/// QX, QY, and QZ are indices of spatial positions in Cartesian co-ordinates.
/// VX, VY, VZ, and UX, UY, UZ are velocity values.
/// |V| is the speed.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
void coordheader( FILE *fout ) {
	fprintf( fout,"t\t QX\t\t QY\t\t QZ\t\tVX\t\tVY\t\tVZ\t\t|V|\t\tUX\t\tUY\t\tUZ\n" );
}

///
/// @brief Prints column headers for flow field data output files.
///
/// Column headers are produced for flow field .dat files to display raw data in a table format.
/// Time, t, is the first column header.
/// QX, QY, and QZ are indices of spatial positions in Cartesian co-ordinates.
/// VcmX, VcmY, and VcmZ are centre of mass velocity values in the Cartesian directions.
/// POP is the cell population
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
void coarseheader( FILE *fout ) {
	int n;
	fprintf( fout,"t\t\t\t\tQX\t\tQY\t\tQZ\tVcmX\t\t\tVcmY\t\t\tVcmZ\t\t\t\tPOP" );
	for( n=0; n<NSPECI; n++ ) fprintf( fout,"\t\tSP%d",n );
	fprintf( fout,"\n" );
}

///
/// @brief Prints column headers for director field data output files.
///
/// Column headers are produced for director field .dat files to display raw data in a table format.
/// Time, t, is the first column header.
/// QX, QY, and QZ are indices of spatial positions in Cartesian co-ordinates.
/// NX, NY, and NZ are director values in the Cartesian directions.
/// S is the scalar order parameter.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
void orderheader( FILE *fout ) {
	fprintf( fout,"t\t\tQX\t\tQY\t\tQZ\t\tNX\t\tNY\t\tNZ\t\tS\n" );
}

///
/// @brief Prints column headers for tensor order parameter data output files.
///
/// Column headers are produced for tensor order parameter .dat files to display raw data in a table format.
/// X, Y, and Z are indices of spatial positions in Cartesian co-ordinates.
/// QXX, QXY, QXZ, QYX, QYZ, QZX, QZY, and QZZ are components of the order parameter tensor.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
void orderQheader( FILE *fout ) {
	fprintf( fout,"X\tY\tZ\tQXX\tQXY\tQXZ\tQYX\tQYY\tQYZ\tQZX\tQZY\tQZZ\n" );
}

///
/// @brief Prints column headers for tensor order parameter data output files.
///
/// Column headers are produced for reciprocal space of tensor order parameter .dat files to display raw data in a table format.
/// Time, t, is the first column header.
/// K123_X, K123_Y, and K123_Z are the Fourier transformed wave vectors in Cartesian space.
/// |QXX|2 to |QZZ|2 are squared components of the order parameter tensor as they may be complex numbers.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
void orderQKheader( FILE *fout ) {
	fprintf( fout,"t\tK123_X\ttK123_Y\ttK123_Z\t|QXX|2\t|QXY|2\t|QXZ|2\t|QYX|2\t|QYY|2\t|QYZ|2\t|QZX|2\t|QZY|2\t|QZZ|2\n" );
}

///
/// @brief Prints column headers for average cell velocity data output files.
///
/// Column headers are produced for average velocity .dat files to display raw data in a table format.
/// Time, t, is the first column header.
/// VcmX, VcmY, and VcmZ are centre of mass velocities in Cartesian co-ordinates.
/// KBT is thermal energy and is only considered if `COLL_TYPE` is set to thermal collisions in input.json.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
void avvelheader( FILE *fout ) {
	fprintf( fout,"t\t VcmX\t\tVcmY\t\tVcmZ\t\tKBT\n" );
}

///
/// @brief Prints column headers for average cell velocity and velocity gradient tensor data output files.
///
/// Column headers are produced for average velocity and velocity gradient .dat files to display raw data in a table format.
/// Time, t, is the first column header.
/// VcmX, VcmY, and VcmZ are centre of mass velocities in Cartesian co-ordinates.
/// KBT is thermal energy and is only considered if `COLL_TYPE` is set to thermal collisions in input.json
/// dVXX, dVYX, and dVZX are derivatives of x, y, and z velocities with respect to x.
/// dVXY, dVYY, and dVZY are derivatives of x, y, and z velocities with respect to y.
/// dVXZ, dVYZ, and dVZZ are derivatives of x, y, and z velocities with respect to z.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
void avvelWithGradVelheader( FILE *fout ) {
	fprintf( fout,"t\t VcmX\t\tVcmY\t\tVcmZ\t\tKBT\t\tdVXX\t\tdVXY\t\tdVXZ\t\tdVYX\t\tdVYY\t\tdVYZ\t\tdVZX\t\tdVZY\t\tdVZZ\t\n" );
}

///
/// @brief Prints column headers for global average orientation output file.
///
/// Column headers are produced for average orientation .dat files to display raw data in a table format.
/// Time, t, is the first column header.
/// NcmX, NcmY, and NcmZ are centre of mass orinetations in Cartesian co-ordinates.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
void avOriheader( FILE *fout ) {
	fprintf( fout,"t\t NcmX\t\tNcmY\t\tNcmZ\t\n" );
}

///
/// @brief Prints column headers for radial correlation data output files.
///
/// Column headers are produced for radial correlation .dat files to display raw data in a table format.
/// Time, t, is the first column header.
/// dr is the radial separation that the correlation function is taken over.
/// C is the value of the correlation function at separation dr.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
void corrheader( FILE *fout ) {
	fprintf( fout,"t\tdr\t C\n" );
}

///
/// @brief Prints column headers for energy spectra data output files.
///
/// Column headers are produced for energy spectra .dat files to display raw data in a table format.
/// Time, t, is the first column header.
/// k is the wave number.
/// E is the corresponding energy value.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
void energyspectheader( FILE *fout ) {
	fprintf( fout,"t\t k\t\t E\n" );
}

///
/// @brief Prints column headers for enstrophy spectra data output files.
///
/// Column headers are produced for enstrophy spectra .dat files to display raw data in a table format.
/// Time, t, is the first column header.
/// k is the wave number.
/// Omega is the corresponding enstropy value.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
void enstrophyspectheader( FILE *fout ) {
	fprintf( fout,"t\t k\t\t Omega\n" );
}

///
/// @brief Prints column headers for topological charge data output files.
///
/// Column headers are produced for topological charge .dat files to display raw data in a table format.
/// Time, t, is the first column header.
/// QX, QY, and QZ are indices of spatial positions in Cartesian co-ordinates.
/// Charge is the resultant topological charge.
/// angle is the orientational angle of a topological defect.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
void topoheader( FILE *fout ) {
	fprintf( fout,"t\t QX\t\t QY\t\t QZ\t\t charge\t\t angle\n" );
}

///
/// @brief Prints column headers for 2D defect data output files.
///
/// Column headers are produced for 2D defect tracking .dat files to display raw data in a table format.
/// Time, t, is the first column header.
/// numDefects relates to the number of defects on a line at time t.
/// QX, QY, and QZ are indices of spatial positions in Cartesian co-ordinates.
/// Charge is the resultant topological charge.
/// angle is the orientational angle of a topological defect.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
void defectheader( FILE *fout ) {
	fprintf( fout,"t\t numDefects\t \n QX\t\t QY\t\t charge\t\t angle\n" );
}

///
/// @brief Prints column headers for 3D defect disclination data output files.
///
/// Column headers are produced for 3D defect tracking .dat files to display raw data in a table format.
/// X, Y, and Z are indices of spatial positions in Cartesian co-ordinates.
/// DXX to DZZ are disclination tensor components.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
void disclinTensorheader( FILE *fout ) {
	fprintf( fout,"X\tY\tZ\tDXX\tDXY\tDXZ\tDYX\tDYY\tDYZ\tDZX\tDZY\tDZZ\n" );
}

///
/// @brief Prints column headers for multiphase data output files.
///
/// Column headers are produced for multiphase .dat files to display raw data in a table format.
/// Time, t, is the first column header.
/// QX, QY, and QZ are indices of spatial positions in Cartesian co-ordinates.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
void multiphaseheader( FILE *fout ) {
	int i;
	fprintf( fout,"t\t\tQX\t\tQY\t\tQZ" );
	for( i=0; i<NSPECI; i++ ) fprintf( fout,"\t\tN_%d",i );
	fprintf( fout,"\n" );
}

///
/// @brief Prints column headers for pressure tensor data output files.
///
/// Column headers are produced for pressure tensor .dat files to display raw data in a table format.
/// Time, t, is the first column header.
/// QX, QY, and QZ are indices of spatial positions in Cartesian co-ordinates.
/// Pxx to Pzz are components of the pressure tensor P.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
void pressureheader( FILE *fout ) {
	fprintf( fout,"t\t\tQX\t\tQY\t\tQZ\tPxx\tPxy\tPxz\tPyx\tPyy\tPyz\tPzx\tPzy\tPzz\n" );
}

///
/// @brief Prints column headers for multiphase pressure data output files.
///
/// Column headers are produced.
/// Time, t, is the first column header.
/// The next column denote each species whose pressure will be outputted below.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
void sppressureheader( FILE *fout ) {
	fprintf( fout,"t\t" );
	for (int q=0;q<NSPECI;q++) {
		fprintf( fout,"\tavP%d",q);
		fprintf( fout,"\tstdP%d",q);
	}
	fprintf( fout,"\n");
}

///
/// @brief Prints column headers for binning data output files for histogram use.
///
/// Column headers are produced for histogram bins in .dat files in order to display raw data in a table format.
/// Bin size is the size of each histogram bin.
/// Time, t, is the first column header.
/// BinderCumulant is the number of counts for the given bin.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param binSize This is the size of bins for the binder.
///
void binderheader( FILE *fout,int binSize ) {
	fprintf( fout,"Bin Size:\t%d\n",binSize );
	fprintf( fout,"t\tBinderCumulant\n" );
}

///
/// @brief Prints column headers for average scalar order parameter data output files.
///
/// Column headers are produced for average scalar order parameter .dat files to display raw data in a table format.
/// Time, t, is the first column header.
/// S and S4 are the scalar order parameter and fourth moment of the scalar order paremeter respectively.
/// nX, nY, and nZ are the director values in each of the Cartesian directions.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
void avsheader( FILE *fout ) {
	fprintf( fout,"t\t\t S\t\t S4\t\t nX\t\t nY\t\t nZ\n" );
}

///
/// @brief Prints column headers for density data output files.
///
/// Column headers are produced for density .dat files to display raw data in a table format.
/// Time, t, is the first column header.
/// densSTD is the standard deviation of the density.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
void densheader( FILE *fout ) {
	fprintf( fout,"t\t\t densSTD\n" );
}

///
/// @brief Prints column headers for average enstrophy data output files.
///
/// Column headers are produced for average enstrophy .dat files to display raw data in a table format.
/// Time, t, is the first column header.
/// Enstrophy is the average enstrophy.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
void avenstrophyheader( FILE *fout ) {
	fprintf( fout,"t\t\t enstrophy\n" );
}

///
/// @brief Prints column headers for velocity data output files.
///
/// Column headers are produced for velocity .dat files to display raw data in a table format.
/// Time, t, is the first column header.
/// QX, QY, and QZ are indices of spatial positions in Cartesian co-ordinates.
/// VcmX, VcmY, and VcmZ are centre of mass velocities in Cartesian co-ordinates.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
void flowheader( FILE *fout ) {
	fprintf( fout,"   t\t   QX\t   QY\t   QZ\tVcmX\t\tVcmY\t\tVcmZ\n" );
}

///
/// @brief Prints column headers for angular velocity and orientation data output files.
///
/// Column headers are produced for angular velocity .dat files to display raw data in a table format.
/// Time, t, is the first column header.
/// QX, QY, and QZ are indices of spatial positions in Cartesian co-ordinates.
/// VX, VY, and VZ are velocities in Cartesian co-ordinates.
/// OX, OY, and OZ are orientation values in Cartesian co-ordinates.
/// LX, LY, and LZ are angular velocities in Cartesian co-ordinates.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
void solidsheader( FILE *fout ) {
	fprintf( fout,"t \tQX\t\tQY\t\tQZ\t\tVX\t\tVY\t\tVZ\t\tOX\t\tOY\t\tOZ\t\tLX\t\tLY\t\tLZ\n" );
}

///
/// @brief Prints column headers for velocity probability distribution histogram data output files.
///
/// Column headers are produced for velocity probability distribution histogram .dat files to display raw data in a table format.
/// Time, t, is the first column header.
/// V is the speed.
/// PX, PY, and PZ are the probability distributions in Cartesian co-ordinates.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
void histVelheader( FILE *fout ) {
	fprintf( fout,"t\t V\t\t\tPX\t\tPY\t\tPZ\n" );
}

///
/// @brief Prints column headers for vorticity probability distribution histogram data output files.
///
/// Column headers are produced for vorticity probability distribution histogram .dat files to display raw data in a table format.
/// Time, t, is the first column header.
/// W is the vorticity.
/// PX, PY, and PZ are the probability distributions in Cartesian co-ordinates.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
void histVortheader( FILE *fout ) {
	fprintf( fout,"t\t w\t\t\tPX\t\tPY\t\tPZ\n" );
}

///
/// @brief Prints column headers for director probability distribution histogram data output files.
///
/// Column headers are produced for director probability distribution histogram .dat files to display raw data in a table format.
/// Time, t, is the first column header.
/// n is the average director value.
/// PX, PY, and PZ are the probability distributions in Cartesian co-ordinates.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
void histDirheader( FILE *fout ) {
	fprintf( fout,"t\t n\t\t\tPX\t\tPY\t\tPZ\n" );
}

///
/// @brief Prints column headers for speed probability distribution histogram data output files.
///
/// Column headers are produced for speed probability distribution histogram .dat files to display raw data in a table format.
/// Time, t, is the first column header.
/// |V| is the speed.
/// P is the probability distribution.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
void histSpeedheader( FILE *fout ) {
	fprintf( fout,"t\t |V|\t\tP\n" );
}

///
/// @brief Prints column headers for enstrophy probability distribution histogram data output files.
///
/// Column headers are produced for enstrophy probability distribution histogram .dat files to display raw data in a table format.
/// Time, t, is the first column header.
/// |w| is the enstrophy.
/// P is the probability distribution.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
void histEnstrheader( FILE *fout ) {
	fprintf( fout,"t\t |w|\t\tP\n" );
}

///
/// @brief Prints column headers for director probability distribution histogram data output files.
///
/// Column headers are produced for director probability distribution histogram .dat files to display raw data in a table format.
/// Time, t, is the first column header.
/// stdN is the standard deviation of the director orientation.
/// P is the probability distribution.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
void histNheader( FILE *fout ) {
	fprintf( fout,"t\t stdN\t\tP\n" );
}

///
/// @brief Prints column headers for scalar order parameter probability distribution histogram data output files.
///
/// Column headers are produced for scalar order parameter probability distribution histogram .dat files to display raw data in a table format.
/// Time, t, is the first column header.
/// S is the scalar order parameter.
/// P is the probability distribution.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
void histSheader( FILE *fout ) {
	fprintf( fout,"t\t S\t\tP\n" );
}

///
/// @brief Prints column headers for energy contribution data output files.
///
/// Column headers are produced for energy .dat files to display raw data in a table format.
/// Time, t, is the first column header.
/// MPC_kin is the kinetic contributions to the energy, with BC_kin the kinetic contribution at the boundary.
/// MPC_nem is the nematic contributions to the energy.
/// BC_rot is the rotational energy contribution at the boundary.
/// Total is the total energy.
/// KBT is the thermal energy contribution and is only considered if `COLL_TYPE` is set to thermal collisions in input.json.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
/// @brief Prints column headers for energy contribution data output files.
///
/// Column headers are produced for energy .dat files to display raw data in a table format.
/// Time, t, is the first column header.
/// MPCD_kin is the kinetic contributions to the energy, with BC_kin the kinetic contribution at the boundary.
/// MPCD_nem is the nematic contributions to the energy.
/// BC_rot is the rotational energy contribution at the boundary.
/// Total is the total energy.
/// KBT is the thermal energy contribution and is only considered if `COLL_TYPE` is set to thermal collisions in input.json.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
void energyheader( FILE *fout ) {
	fprintf( fout,"t\t\tMPCD_kin\t\tMPCD_nem\t\tBC_kin\t\tBC_rot\t\tTotal\t\tKBT\n" );
}

///
/// @brief Prints column headers for energy field data output files.
///
/// Column headers are produced for energy field .dat files to display raw data in a table format.
/// Time, t, is the first column header.
/// QX, QY, and QZ are indices of spatial positions in Cartesian co-ordinates.
/// MPCD_kin is the kinetic contributions to the energy.
/// MPCD_nem is the nematic contributions to the energy.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
void energyfieldheader( FILE *fout ) {
	fprintf( fout,"QX\tQY\tQZ\tMPCD_kin\tMPCD_nem\n" );
}

///
/// @brief Prints column headers for swimmer data output files.
///
/// Column headers are produced for swimmer .dat files to display raw data in a table format.
/// Time, t, is the first column header.
/// Swimmer consists of a head (H) and a middle (M).
/// HX to HZ and MX to MZ are head and middle positions in Cartesian co-ordinates respectively.
/// HVX, HVY, and HVZ are head velocities in Cartesian co-ordinates.
/// MVX, MVY, and MVZ are middle velocities in Cartesian co-ordinates.
/// RTphase is the run-tumble phase of the swimmer.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
void swimmerheader( FILE *fout ) {
	fprintf( fout,"t\t\tHX\tHY\tHZ\tHVX\tHVY\tHVZ\tMX\tMY\tMZ\tMVX\tMVY\tMVZ\tRTphase\n" );
}

///
/// @brief Prints column headers for swimmer orientation data output files.
///
/// Column headers are produced for swimmer orientation .dat files to display raw data in a table format.
/// Time, t, is the first column header.
/// nX, nY, and nZ are the swimmer orientation in Cartesian co-ordinates.
/// RTphase is the run-tumble phase of the swimmer.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
void swimmeroriheader( FILE *fout ) {
	fprintf( fout,"t\t\tnX\tnY\tnZ\tRTphase\n" );
}

///
/// @brief Prints column headers for centre of masses of multiphase fluids.
///
/// Column headers are produced.
/// Time, t, is the first column header.
/// X, Y, and Z are the centres of masses' Cartesian co-ordinates.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
void comheader( FILE *fout ) {
	fprintf( fout,"t\t\t" );
	for (int k=0;k<NSPECI;k++) fprintf( fout,"X%d\tY%d\tZ%d\tvarXX%d\tvarXY%d\tvarXZ%d\tvarYX%d\tvarYY%d\tvarYZ%d\tvarZX%d\tvarZY%d\tvarZZ%d\t",k,k,k,k,k,k,k,k,k,k,k,k );
	fprintf( fout,"\n" );
}

///
/// @brief Prints column headers for run-tumble data output files.
///
/// Column headers are produced for run-tumble .dat files to display raw data in a table format.
/// The first column, RTphase, is the run-tumble phase of the swimmer.
/// dt_cnt is the timestep count.
/// dAng is the change in angle of the swimmer.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
void runtumbleheader( FILE *fout ) {
	fprintf( fout,"RTphase\t\tdt_cnt\tdAng\n" );
}

///
/// @brief Prints column headers for neighbouring cell nematic energy data output files.
///
/// Column headers are produced for neighbouring cell nematic energy .dat files to display raw data in a table format.
/// Time, t, is the first column header.
/// MPCD_nem is the nematic contributions to the energy.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
///
void energyneighboursheader( FILE *fout ) {
	fprintf( fout,"t\ttMPCD_nem\n" );
}

///
/// @brief Prints Cartesian co-ordinates of MPCD particles to data output files.
///
/// This function prints the Cartesian co-ordinates of MPC particles to output files that require it and are requested by input.json.
///
/// @param fout This is a pointer to the output .dat file name to be produced, covering all species specified.
/// @param pr This is an index value for the subpopulation SP.
/// @param T Timestep, the output rate of which is specified in input.json.
/// @param p List of MPCD particle index numbers.
/// @param SP Subpopulation of species.
/// @see outputResults()
///
void coordout( FILE *fout[MAXSPECI],int pr,double T,particleMPC p[],spec SP[] ) {
	int i,j;
	double v;
	for( i=0; i<NSPECI; i++ ) {
		//Do this for the species that are to be printed out and only for the species with population and mass greater than zero
		if( i<pr && SP[i].POP!=0 && fneq(SP[i].MASS,0.0) ) {
			//fprintf(fout[i],"\nTIME STEP: %i\n",T);//Print to file
			for( j=0; j<GPOP; j++ ) if( (p+j)->SPID == i ) {
				fprintf( fout[i],"%12.5e\t",T );
				fprintf( fout[i],"%12.5e\t%12.5e\t%12.5e\t",p[j].Q[0],p[j].Q[1],p[j].Q[2] );
				fprintf( fout[i],"%12.5e\t%12.5e\t%12.5e\t",p[j].V[0],p[j].V[1],p[j].V[2] );
				v = sqrt( p[j].V[0]*p[j].V[0]+p[j].V[1]*p[j].V[1]+p[j].V[2]*p[j].V[2] );
				fprintf( fout[i],"%12.5e\t",v );
				fprintf( fout[i],"%12.5e\t%12.5e\t%12.5e\n",p[j].U[0],p[j].U[1],p[j].U[2] );
			}
			fprintf( fout[i],"\n" );
			#ifdef FFLSH
				fflush(fout[i]);
			#endif
		}
	}
}

///
/// @brief Prints coarse grained data to data output files.
///
/// This function prints the coarse grained cell and velocity data to output files that require it and are requested by input.json.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param t Timestep, the output rate which is specified in input.json.
/// @param CL This is a pointer to the co-ordinates and cell of each particle in the MPCD list.
/// @see cellout()
/// @see outputResults()
///
void coarseout( FILE *fout,double t,cell ***CL ) {
	int i,j,k,n;
	for( i=0; i<XYZ[0]; i++ ) for( j=0; j<XYZ[1]; j++ ) for( k=0; k<XYZ[2]; k++ ) {
		fprintf( fout,"%.2f\t",t );
		fprintf( fout,"%5d\t%5d\t%5d\t",i,j,k );
		if( CL[i][j][k].POP == 0 ) {
			fprintf( fout, "%12.5e\t%12.5e\t%12.5e\t%5i",0.,0.,0.,0 );
			for( n=0; n<NSPECI; n++ ) fprintf( fout, "\t%5i",0 );
			fprintf( fout,"\n" );
		}
		else {
			fprintf( fout, "%12.5e\t%12.5e\t%12.5e\t%5i",CL[i][j][k].VCM[0],CL[i][j][k].VCM[1],CL[i][j][k].VCM[2],CL[i][j][k].POP );
			for( n=0; n<NSPECI; n++ ) fprintf( fout, "\t%5i",CL[i][j][k].SP[n] );
			fprintf( fout,"\n" );
		}
	}
	#ifdef FFLSH
		fflush(fout);
	#endif
}

///
/// @brief Produces co-ordinates of MPCD and MD particles.
///
/// This function produces the co-ordinates, resident cell, and resident cell population of particles in the list of MPCD and MD particles including swimmers.
///
/// @param CL This is a pointer to the co-ordinates and cell of each particle.
/// @see outputResults()
///
void cellout( cell ***CL ) {
	int i,j,k,l;
	particleMPC *pMPC;	//Temporary pointer to MPCD particles
	particleMD *pMD;		//Temporary pointer to MD particles
	smono *pSW;					//Temporary pointer to swimmer monomers

	for( i=0; i<XYZ_P1[0]; i++ ) for( j=0; j<XYZ_P1[1]; j++ ) for( k=0; k<XYZ_P1[2]; k++ ) {
		l=0;
		// MPCD
		pMPC = CL[i][j][k].pp;
		pMD = CL[i][j][k].MDpp;
		pSW = CL[i][j][k].sp;

		if( pMPC != NULL || pMD != NULL || pSW != NULL ) printf( "In Cell [%d,%d,%d]:\n",i,j,k );
		while( pMPC != NULL ) {
			l++;
			printf( "\tMPCD Particle %d:",l );
			printf( "\t\t\tQ=(%lf,%lf,%lf)\n",pMPC->Q[0],pMPC->Q[1],pMPC->Q[2] );
			printf( "\t\t\tV=(%lf,%lf,%lf) \n",pMPC->V[0],pMPC->V[1],pMPC->V[2] );
			//Increment link in list
			pMPC = pMPC->next;
		}
		// MD
		while( pMD != NULL ) {
			l++;
			printf( "\tMD Particle %d:",l );
			printf( "\t\t\tQ=(%lf,%lf,%lf)\n",pMD->rx,pMD->ry,pMD->rz );
			printf( "\t\t\tV=(%lf,%lf,%lf) \n",pMD->vx,pMD->vy,pMD->vz );
			//Increment link in list
			pMD = pMD->nextSRD;
		}
		while( pSW != NULL ) {
			l++;
			printf( "\tSwimmer monomer %d:",l );
			printf( "\t\t\tQ=(%lf,%lf,%lf)\n",pSW->Q[0],pSW->Q[1],pSW->Q[2] );
			printf( "\t\t\tV=(%lf,%lf,%lf) \n",pSW->V[0],pSW->V[1],pSW->V[2] );
			//Increment link in list
			pSW = pSW->next;
		}
		printf( "\tP=%d\n",CL[i][j][k].POP );
	}
}

///
/// @brief Outputs entire list of MPCD particles and MD partciles.
///
/// This function outputs the list of all particles co-ordinates and cell information as an array of lists.
///
/// @param CL This is a pointer to the co-ordinates and cell of each particle.
/// @param XYZ_p1 This is three-dimensional list of particle positions.
/// @see cellout()
///
void listout( cell ***CL,int XYZ_p1[_3D] ) {
	int a,b,c,d;
	particleMPC *pMPC;
	particleMD *pMD;
	smono *pSW;

	printf( "Local properties:\n" );
	for( a=0; a < XYZ_p1[0]; a++ ) for(b=0; b < XYZ_p1[1]; b++ ) for(c=0; c < XYZ_p1[2]; c++ ) {
		d=0;
		printf( "\tCell [%d,%d,%d]:\n",a,b,c );
		if( CL[a][b][c].pp != NULL ) {
			pMPC = CL[a][b][c].pp;
			while( pMPC != NULL ) {
				d++;
				printf( "\t\tMPCD Particle %d:\n",d );
				printf( "\t\t\tQ=(%lf,%lf,%lf)\n",pMPC->Q[0],pMPC->Q[1],pMPC->Q[2] );
				printf( "\t\t\tV=(%lf,%lf,%lf)\n",pMPC->V[0],pMPC->V[1],pMPC->V[2] );
				pMPC = pMPC->next;
			}
		}
		d = 0;
		if( CL[a][b][c].MDpp != NULL ) {
			pMD = CL[a][b][c].MDpp;
			while( pMD != NULL ) {
				d++;
				printf( "\t\tMD Particle %d:\n",d );
				printf( "\t\t\tQ=(%lf,%lf,%lf)\n",pMD->rx,pMD->ry,pMD->rz );
				printf( "\t\t\tV=(%lf,%lf,%lf)\n",pMD->vx,pMD->vy,pMD->vz );
				pMD = pMD->nextSRD;
			}
		}
		d = 0;
		if( CL[a][b][c].sp != NULL ) {
			pSW = CL[a][b][c].sp;
			while( pSW != NULL ) {
				d++;
				printf( "\t\tSwimmer Particle %d:\n",d );
				printf( "\t\t\tQ=(%lf,%lf,%lf)\n",pSW->Q[0],pSW->Q[1],pSW->Q[2] );
				printf( "\t\t\tV=(%lf,%lf,%lf)\n",pSW->V[0],pSW->V[1],pSW->V[2] );
				pSW = pSW->next;
			}
		}
		if( CL[a][b][c].pp != NULL || CL[a][b][c].MDpp != NULL || CL[a][b][c].sp != NULL ) {
			printf( "\t\tPopulation: %d\n",CL[a][b][c].POP );
			printf( "\t\tMass: %lf\n",CL[a][b][c].MASS );
// 			printf( "\t\tThermal Energy: %lf\n",CL[a][b][c].KBT );
			printf( "\t\tCentre of Mass Velocity: (%lf,%lf,%lf)\n",CL[a][b][c].VCM[0],CL[a][b][c].VCM[1],CL[a][b][c].VCM[2] );
		}
	}
}

///
/// @brief Outputs position and velocity of MPCD particles to the terminal.
///
/// This function prints MPCD particle respective positions and velocities to terminal.
///
/// @param p This is an index for each MPCD particle.
///
void pcoord( particleMPC p ) {
	printf( "\tQ=(%6.12e,%6.12e,%6.12e)\n",p.Q[0],p.Q[1],p.Q[2] );
	printf( "\tV=(%6.12e,%6.12e,%6.12e)\n",p.V[0],p.V[1],p.V[2] );
	printf( "\tU=(%6.12e,%6.12e,%6.12e)\n",p.U[0],p.U[1],p.U[2] );
}

///
/// @brief Outputs boundary information to the terminal.
///
/// This function prints boundary positions, velocity, and angular velocity to terminal.
///
/// @param WALL This is a pointer obtaining information on boundary conditions.
///
void bccoord( bc WALL ) {
	printf( "\tQ=(%6.12e,%6.12e,%6.12e)\n",WALL.Q[0],WALL.Q[1],WALL.Q[2] );
	printf( "\tV=(%6.12e,%6.12e,%6.12e)\n",WALL.V[0],WALL.V[1],WALL.V[2] );
	printf( "\tO=(%6.12e,%6.12e,%6.12e)\n",WALL.O[0],WALL.O[1],WALL.O[2] );
	printf( "\tL=(%e,%e,%e)\n",WALL.L[0],WALL.L[1],WALL.L[2] );
}

///
/// @brief Outputs position and velocity of MD partciles to the terminal.
///
/// This function prints MD particle respective positions and velocities to terminal.
///
/// @param p This is an index for each MD particle.
///
void mdcoord( particleMD p ) {
	printf( "\tQ=(%lf,%lf,%lf)\n",p.rx,p.ry,p.rz );
	printf( "\tV=(%lf,%lf,%lf)\n",p.vx,p.vy,p.vz );
}

///
/// @brief Outputs position and velocity of swimmers to the terminal.
///
/// This function prints swimmers respective positions and velocities to terminal.
///
/// @param sw This is an index for each swimmer.
///
void swcoord( swimmer sw ) {
	printf( "\tH Q=(%lf,%lf,%lf) ",sw.H.Q[0],sw.H.Q[1],sw.H.Q[2] );
	printf( "\tV=(%lf,%lf,%lf)\n",sw.H.V[0],sw.H.V[1],sw.H.V[2] );
	printf( "\tM Q=(%lf,%lf,%lf) ",sw.M.Q[0],sw.M.Q[1],sw.M.Q[2] );
	printf( "\tV=(%lf,%lf,%lf)\n",sw.M.V[0],sw.M.V[1],sw.M.V[2] );
}

///
/// @brief Generic function to print vectors.
///
/// This is a generic function that prints vectors.
///
/// @param VEC This is the vector that will be printed.
/// @param dimension This is the dimensionality of the vector that will be printed.
///
void pvec( double VEC[],int dimension ) {
	int i;
	printf( " (" );
	for( i=0; i<(dimension-1); i++ ) printf( "%lf,",VEC[i] );
	printf( "%lf)\n",VEC[dimension-1] );
}

///
/// @brief Generic function to print 3D tensors.
///
/// This is a generic function that prints three-dimensional tensors.
///
/// @param TENS This is the 3D tensor that will be printed.
///
void ptens3D( double TENS[][_3D] ) {
	int i,j;
	int dimension=3;
	printf( " [ " );
	for( j=0; j<(dimension); j++ ) {
		printf( " (" );
		for( i=0; i<(dimension-1); i++ ) printf( "%lf,",TENS[j][i] );
		printf( "%lf)\n",TENS[j][dimension-1] );
	}
	printf( " ]\n" );
}

///
/// @brief Generic function to print 2D tensors.
///
/// This is a generic function that prints two-dimensional tensors.
///
/// @param TENS This is the 2D tensor that will be printed.
///
void ptens2D( double TENS[][_2D] ) {
	int i,j;
	int dimension=2;
	printf( " [ " );
	for( j=0; j<(dimension); j++ ) {
		printf( " (" );
		for( i=0; i<(dimension-1); i++ ) printf( "%lf,",TENS[j][i] );
		printf( "%lf)\n",TENS[j][dimension-1] );
	}
	printf( " ]\n" );
}

///
/// @brief Generic function to print 2D or 3D tensors.
///
/// This is a generic function that prints two-dimensional and three-dimensional tensors.
///
/// @param TENS This is pointer to the 2D or 3D tensor that will be printed.
/// @param dimension This is the dimensionality of the tensor.
///
void ptens( double **TENS,int dimension ) {
	int i,j;
	double T2[_2D][_2D],T3[_3D][_3D];
	if(dimension==_3D) {
		for( i=0; i<dimension; i++ ) for( j=0; j<dimension; j++ ) T3[i][j]=TENS[i][j];
		ptens3D( T3 );
	}
	else if(dimension==_2D) {
		for( i=0; i<dimension; i++ ) for( j=0; j<dimension; j++ ) T2[i][j]=TENS[i][j];
		ptens2D( T2 );
	}
	else printf( "Warning: can only print tensors of order 2 or 3\n" );
}

///
/// @brief Outputs position and velocity of MPCD particles to the terminal.
///
/// This function prints MPCD particle respective positions and velocities to terminal.
///
/// @param POS This is the three-dimensional co-ordinates of MPCD particles.
/// @param VEL This is the three-dimensional velocities of MPCD particles.
/// @param ANG This is the three-dimensional angular velocities of MPCD particles.
/// @param dimension This is the dimensionality.
///
void pvcoord( double POS[_3D],double VEL[_3D],double ANG[_3D],int dimension ) {
	printf( "\tQ=" );
	pvec( POS,dimension );
	printf( "\tV=" );
	pvec( VEL,dimension );
	printf( "\tW=" );
	pvec( ANG,dimension );

}

///
/// @brief Outputs all MPCD particles information to the terminal.
///
/// This function prints all MPCD particle information to terminal.
///
/// @param p This is an index for each MPCD particle.
///
void pall( particleMPC p[] ) {
	int i;
	for( i=0; i<GPOP; i++ ) {
		printf( "Particle %i\t",i );
		pcoord( p[i] );
	}
}

///
/// @brief Outputs input data and parameters.
///
/// This function prints relevent input information to the terminal.
///
/// @param in This is a pointer that fetches information from input.json.
/// @param AVVEL This is the average speed in any direction.
/// @param SP This is a pointer that fetches species information such as population.
/// @param theorySP This is a pointer that fetches theoretical information for each species.
/// @param theoryGl This is the theoretical information for the global system.
///
void listinput( inputList in,double AVVEL,spec SP[],kinTheory theorySP[],kinTheory theoryGl ) {
	int i,n;
	#ifdef DBG
		if( DBUG >= DBGINIT ){
			printf( "\tDebug Mode: %i\n",DBUG );
			printf( "\tDimensions: %i\n",DIM );
			printf( "\tRotation Technique: %i\n",in.RTECH );
			printf( "\tLiquid Crystal: %i\n",in.LC );
			printf( "\tRotation angle: %lf\n",in.RA );
			printf( "\tSystem Size: [%i,%i,%i]\n",XYZ[0],XYZ[1],XYZ[2] );
			printf( "\tAccessible volume: %lf\n",VOL );
			printf( "\tNo hydodynamic interaction: %i\n",in.noHI);
			printf( "\tIncompressibility correction: %i\n",in.inCOMP);
			printf( "\tMultiphase interactions: %i\n",in.MULTIPHASE);
			printf( "\tConstant external acceleration:" );
			pvec( in.GRAV,DIM );
			printf( "\tConstant external magnetic field:" );
			pvec( in.MAG,DIM );
			printf( "\tAverage speed in any direction: %lf\n",AVVEL );
			printf( "\tNumber of Species: %d\n",NSPECI );
			printf( "\tWarmup Iterations: %d\n",in.warmupSteps );
			printf( "\tSimulation Iterations: %d\n",in.simSteps );
			printf( "\tTime Step: %lf\n",in.dt );
			printf( "\n\tThe number of boundaries: %i\n",NBC );
			printf( "\tInputted System Temperature (units of KB): %lf\n",in.KBT );
			printf( "\tRemove system's net momentum (1=yes, 0=no): %i\n",in.RFRAME );
			printf( "\tThermostat Method: %i\n",in.TSTECH );
			printf( "\tThermal Relaxation Scale: %lf\n",in.TAU );
			printf( "\tSystem-wide properties:\n" );
			printf( "\t\tSystem population: %i\n",GPOP );
			printf( "\t\tSystem total mass: %lf\n",GMASS );
			printf( "\t\tParticle number density: %lf\n",GnDNST );
			printf( "\t\tMass density: %lf\n",GmDNST );
			printf( "\t\tMean Free Path: %lf\n",theoryGl.MFP );
			printf( "\t\tKinematic Viscosity: %lf\n",theoryGl.VISC );
			printf( "\t\tSelf Diffusion Coefficient: %lf\n",theoryGl.SDIFF );
			printf( "\t\tSchmidt number: %lf\n",theoryGl.VISC/theoryGl.SDIFF/GmDNST );
			printf( "\t\tSpeed of sound: %lf\n",theoryGl.SPEEDOFSOUND );
			printf( "\t\tThermal Diffusion Coefficient: %lf\n",theoryGl.THERMD );
			printf( "\tSpecies properties:\n" );
			for( i=0; i<NSPECI; i++ ) {
				printf( "\t\tSpecies ID: %i\n",i );
				printf( "\t\tPopulation: %i\n",SP[i].POP );
				printf( "\t\tAccessible volume: %lf\n",SP[i].VOL );
				printf( "\t\tParticle mass: %lf\n",SP[i].MASS );
				if(in.LC) {
					printf( "\t\tLiquid crystal properties:\n" );
					printf( "\t\t\tRotational friction coefficient: %lf\n",SP[i].RFC );
					printf( "\t\t\tTumbling parameter: %lf\n",SP[i].TUMBLE );
					printf( "\t\t\tHydrodynamic susceptibility: %lf\n",SP[i].CHIHI );
					printf( "\t\t\tMagnetic susceptibility: %lf\n",SP[i].CHIA );
					printf( "\t\t\tLength: %lf\n",SP[i].LEN );
					printf( "\t\t\tActivity: %lf\n",SP[i].ACT );
					printf( "\t\t\t\tActive dipole sigmoid width: %lf\n",SP[i].SIGWIDTH );
					printf( "\t\t\t\tActive dipole sigmoid position: %lf\n",SP[i].SIGPOS );
					printf( "\t\t\tMean field potential: %lf\n",SP[i].MFPOT );
				}
				printf( "\t\tInteraction matrix: " );
				pvec( SP[i].M,NSPECI );
				printf( "\t\tSpecies total mass: %lf\n",SP[i].DAMP );
				printf( "\t\tDamping/friction coefficient: %lf\n",theorySP[i].sumM );
				printf( "\t\tThe following assume the species are perfectly separated from each other\n" );
				printf( "\t\t\tParticle number density: %lf\n",SP[i].nDNST );
				printf( "\t\t\tMass density: %lf\n",SP[i].mDNST );
				printf( "\t\t\tMean Free Path: %lf\n",theorySP[i].MFP );
				printf( "\t\t\tKinematic Viscosity: %lf\n",theorySP[i].VISC );
				printf( "\t\t\tSelf Diffusion Coefficient: %lf\n",theorySP[i].SDIFF );
				printf( "\t\t\tSchmidt number: %lf\n",theorySP[i].VISC/theorySP[i].SDIFF/SP[i].mDNST );
				printf( "\t\t\tSpeed of sound: %lf\n",theorySP[i].SPEEDOFSOUND );
				printf( "\t\t\tThermal Diffusion Coefficient: %lf\n",theorySP[i].THERMD );
			}
		}
	#endif
	n = 0;
	for( i=0; i<NSPECI; i++ ) n += SP[i].POP;
	if( n != GPOP ){
		printf( "Error:\tGPOP does not equal sum of species' pop.\n" );
		exit( 1 );
	}
	if( DIM != _3D && XYZ[2] != 1 ){
		printf( "Error:\tZ component contradicts stated number of dimensions (%i)\n",DIM );
		exit( 1 );
	}
}

///
/// @brief Parameters and input data ouput to synopsis data file.
///
/// This function prints relevent information to the synopsis data file. All data is given in these <a href="https://journals.aps.org/pre/abstract/10.1103/PhysRevE.78.016706">units</a>.
///
/// @param in This is a pointer that fetches information from input.json.
/// @param SP This is a pointer that fetches species information such as population.
/// @param WALL This is a pointer obtaining information on boundary conditions.
/// @param SS This is a pointer that fetches swimmer specifications such as initial conditions, type, and run-tumble conditions.
/// @param out This is a flag that determines if data should be output or not.
/// @param theorySP This is a pointer that fetches theoretical information for each species.
/// @param theoryGl This is the theoretical information for the global system.
/// @param fsynopsis This is a pointer to the synopsis.dat output file.
///
void stateinput( inputList in,spec SP[],bc WALL[],specSwimmer SS,outputFlagsList out,kinTheory theorySP[],kinTheory theoryGl,FILE *fsynopsis ) {
	int i,j;

	if( out.SYNOUT == OUT ) {
		fprintf( fsynopsis,"\nBasic Units:\n" );
		fprintf( fsynopsis,"Length: a = 1, MPCD cell size\n" );
		fprintf( fsynopsis,"Mass: m = 1, MPCD particle mass\n" );
		fprintf( fsynopsis,"Energy: kT = 1, thermal energy\n" );
		fprintf( fsynopsis,"Derived Units:\n" );
		fprintf( fsynopsis,"Time: tau = a * sqrt(m/kT) = 1\n" );
		fprintf( fsynopsis,"Density units: 1/a^d = 1/a^%i\n",DIM );
		fprintf( fsynopsis,"Diffusion const: a * sqrt(kT/m) = a^2/tau\n" );
		fprintf( fsynopsis,"Stress: kT/a^d = kT/a^%i\n",DIM );
		fprintf( fsynopsis,"Dynamic viscosity: sqrt(m*kT)/a^(d-1) = kT*tau/a^d = kT*tau/a^%i\n",DIM );
		fprintf( fsynopsis,"Kinetic viscosity: kT*tau/m\n" );

		fprintf( fsynopsis,"\nUser defined variables:\n" );
		fprintf( fsynopsis,"Dimensionality: %i\n",DIM );
		fprintf( fsynopsis,"System dimensions: (%i,%i,%i)\n",XYZ[0],XYZ[1],XYZ[2] );
		fprintf( fsynopsis,"Volume of the control volum: %lf\n",VOL );
		fprintf( fsynopsis,"Rotation technique: %i\n",in.RTECH );
		fprintf( fsynopsis,"Nematic Liquid Crystal: ");
		if(in.LC) fprintf( fsynopsis,"YES\n" );
		else fprintf( fsynopsis,"NO\n" );
		fprintf( fsynopsis,"Thermal energy: %lf\n",in.KBT );
		fprintf( fsynopsis,"Remove system's net momentum (1=yes, 0=no): %i\n",in.RFRAME );
		fprintf( fsynopsis,"Randomly shift the MPCD cells for Galilean invariance (1=yes, 0=no): %i\n",in.GALINV );
		fprintf( fsynopsis,"Thermostating method: %i\n",in.TSTECH );
		fprintf( fsynopsis,"Thermal relaxation time scale: %lf\n",in.TAU );
		fprintf( fsynopsis,"Rotation angle: %lf\n",in.RA );
		fprintf( fsynopsis,"Langevin friction coefficient: %lf\n",in.FRICCO );
		fprintf( fsynopsis,"Hydodynamic interactions: ");
		if(in.noHI) fprintf( fsynopsis,"OFF\n" );
		else fprintf( fsynopsis,"On\n" );
		fprintf( fsynopsis,"Incompressibility correction: ");
		if(in.inCOMP) fprintf( fsynopsis,"OFF\n" );
		else fprintf( fsynopsis,"On\n" );
		fprintf( fsynopsis,"Multiphase interactions mode: %d\n",in.MULTIPHASE );
		fprintf( fsynopsis,"External acceleration: (%lf,%lf,%lf)\n",in.GRAV[0],in.GRAV[1],in.GRAV[2] );
		fprintf( fsynopsis,"External magnetic field: (%lf,%lf,%lf)\n",in.MAG[0],in.MAG[1],in.MAG[2] );
		fprintf( fsynopsis,"Total simulation iterations: %d\n",in.simSteps );
		fprintf( fsynopsis,"Warmup iterations: %d\n",in.warmupSteps );
		fprintf( fsynopsis,"Time step: %lf\n",in.dt );
		fprintf( fsynopsis,"Random seed: %ld\n",in.seed );

		fprintf( fsynopsis,"\nSpecies variables:\n" );
		fprintf( fsynopsis,"Number of species: %d\n",NSPECI );
		for( i=0; i<NSPECI; i++ ) {
			fprintf( fsynopsis,"Species: %i\n",i );
			fprintf( fsynopsis,"\tMass: %lf\n\tPopulation: %i\n",SP[i].MASS,SP[i].POP );
			fprintf( fsynopsis,"\tVolume accessible: %lf\n\tParticle Number Density: %lf\n\tMass Density: %lf\n",SP[i].VOL,SP[i].nDNST,SP[i].mDNST );
			fprintf( fsynopsis,"\tRotational Friction Coefficient: %lf\n",SP[i].RFC);
			fprintf( fsynopsis,"\tEffective rod-length to couple MPCD torque to BC force: %lf\n",SP[i].LEN);
			fprintf( fsynopsis,"\tTumbling parameter: %lf\n",SP[i].TUMBLE);
			fprintf( fsynopsis,"\tMagnetic Sysceptibility: %lf\n\tShear Sysceptibility: %lf\n",SP[i].CHIA,SP[i].CHIHI );
			fprintf( fsynopsis,"\tActivity: %lf\n\tActivity-sigmoid width: %lf\n\tActivity-sigmoid position: %lf\n\tMinimum proportion for activity %lf\n",SP[i].ACT,SP[i].SIGWIDTH,SP[i].SIGPOS,SP[i].MINACTRATIO );
			fprintf( fsynopsis,"\tMean field potential: %lf\n",SP[i].MFPOT );
			fprintf( fsynopsis,"\tDamping friction: %lf\n",SP[i].DAMP );
			fprintf( fsynopsis,"\tPos. dist.: %i\n\tVel. dist.: %i\n\tOri. dist.: %i\n",SP[i].QDIST,SP[i].VDIST,SP[i].ODIST );
			fprintf( fsynopsis,"\tInteraction matrix: [" );
			for( j=0; j<NSPECI-1; j++ ) fprintf( fsynopsis,"%lf, ",SP[i].M[j] );
			fprintf( fsynopsis,"%lf]\n",SP[i].M[NSPECI-1] );
			fprintf( fsynopsis,"\tAccessible volume: %lf\n\tNumber density: %lf\n\tMass density: %lf\n",SP[i].VOL,SP[i].nDNST,SP[i].mDNST );
		}
		fprintf( fsynopsis,"\nBC variables:\n" );
		fprintf( fsynopsis,"Number of BCs: %i\n",NBC );
		fprintf( fsynopsis,"Periodic BC axes: [%i,%i,%i]\n",XYZPBC[0],XYZPBC[1],XYZPBC[2] );
		for( i=0; i<NBC; i++ ) {
			fprintf( fsynopsis,"BC %i:\n",i );
			fprintf( fsynopsis,"\tCoefficient of Restitution: %lf\n",WALL[i].E );
			fprintf( fsynopsis,"\tCentre: (%lf,%lf,%lf)\n",WALL[i].Q[0],WALL[i].Q[1],WALL[i].Q[2] );
			fprintf( fsynopsis,"\tInitial Velocity: (%lf,%lf,%lf)\n",WALL[i].V[0],WALL[i].V[1],WALL[i].V[2] );
			fprintf( fsynopsis,"\tInitial Orientation: (%lf,%lf,%lf)\n",WALL[i].O[0],WALL[i].O[1],WALL[i].O[2] );
			fprintf( fsynopsis,"\tInitial Angular Velocity: (%lf,%lf,%lf)\n",WALL[i].L[0],WALL[i].L[1],WALL[i].L[2] );
			fprintf( fsynopsis,"\tNormal: (%lf,%lf,%lf)\n",WALL[i].A[0],WALL[i].A[1],WALL[i].A[2] );
			fprintf( fsynopsis,"\tPowers: (%lf,%lf,%lf)\n",WALL[i].P[0],WALL[i].P[1],WALL[i].P[2] );
			fprintf( fsynopsis,"\tRadius (with power %lf): %lf\n",WALL[i].P[3],WALL[i].R );
			fprintf( fsynopsis,"\tAbsolute value of terms (0=false; 1=true):  %d\n",WALL[i].ABS );
			fprintf( fsynopsis,"\t%lf-fold and %lf-fold rotational symmetries\n",WALL[i].ROTSYMM[0],WALL[i].ROTSYMM[1] );
			fprintf( fsynopsis,"\tSurface rotations (0=false; 1=true): %d\n",WALL[i].REORIENT );
			fprintf( fsynopsis,"\tTransformation POS: \t+%lf norm, +%lf tang\n",WALL[i].DN,WALL[i].DT );
			fprintf( fsynopsis,"\tTransformation VEL: \t+%lf norm, +%lf tang\n",WALL[i].DVN,WALL[i].DVT );
			fprintf( fsynopsis,"\t\t\t\t\t*%lf norm, *%lf tang\n",WALL[i].MVN,WALL[i].MVT );
			fprintf( fsynopsis,"\tTransformation ORI: \t+%lf norm, +%lf tang\n",WALL[i].MUN,WALL[i].MUT );
			fprintf( fsynopsis,"\tTransformation ORI: \t%lf x, %lf y, %lf z\n",WALL[i].MUxyz[0],WALL[i].MUxyz[1],WALL[i].MUxyz[2] );
			fprintf( fsynopsis,"\tMove flag: %i\n",WALL[i].DSPLC );
			fprintf( fsynopsis,"\tBC mass: %lf\n",WALL[i].MASS );
			fprintf( fsynopsis,"\tBC volume: %lf\n",WALL[i].VOL );
			fprintf( fsynopsis,"\tBC inertia tensor:\t[%lf, %lf %lf]\n",WALL[i].I[0][0],WALL[i].I[1][0],WALL[i].I[2][0] );
			fprintf( fsynopsis,"\t\t\t\t\t\t\t\t[%lf, %lf %lf]\n",WALL[i].I[0][1],WALL[i].I[1][1],WALL[i].I[2][1] );
			fprintf( fsynopsis,"\t\t\t\t\t\t\t\t[%lf, %lf %lf]\n",WALL[i].I[0][2],WALL[i].I[1][2],WALL[i].I[2][2] );
			fprintf( fsynopsis,"\tInteraction matrix with SRD particles: [" );
			for( j=0; j<NSPECI-1; j++ ) fprintf( fsynopsis,"%i, ",WALL[i].INTER[j] );
			fprintf( fsynopsis,"%i]\n",WALL[i].INTER[NSPECI-1] );
			fprintf( fsynopsis,"\tInteraction matrix with MD particles: %i\n",WALL[i].INTER[MAXSPECI+0] );
			fprintf( fsynopsis,"\tInteraction matrix with swimmer particles: %i\n",WALL[i].INTER[MAXSPECI+1] );
		}
		fprintf( fsynopsis,"\nSwimmer variables:\n" );
		fprintf( fsynopsis,"\tType: %d\n",SS.TYPE );
		fprintf( fsynopsis,"\tNumber of swimmers: %d\n",NS );
		fprintf( fsynopsis,"\tInitialize position: %d\n",SS.QDIST );
		fprintf( fsynopsis,"\tInitialize orientation: %d\n",SS.ODIST );
		fprintf( fsynopsis,"\tHead Mass: %d\n",SS.headM );
		fprintf( fsynopsis,"\tHead Species id: %d\n",SS.HSPid );
		fprintf( fsynopsis,"\tMiddle Mass: %d\n",SS.middM );
		fprintf( fsynopsis,"\tMiddle Species id: %d\n",SS.MSPid );
		fprintf( fsynopsis,"\tSwimming propulsion force: %lf\n",SS.FS );
		fprintf( fsynopsis,"\tSwimming dipole strength: %lf",SS.DS );
		if(SS.DS>0.0) fprintf( fsynopsis," --- Pusher\n" );
		else fprintf( fsynopsis," --- Puller\n" );
		fprintf( fsynopsis,"\tTumbling shrink size: %lf",SS.sizeShrink );
		fprintf( fsynopsis,"\tTumbling shrink spring const: %lf",SS.springShrink );
		fprintf( fsynopsis,"\tSpring const: %lf\n",SS.k );
		fprintf( fsynopsis,"\tSpring separation: %lf\n",SS.ro );
		fprintf( fsynopsis,"\tLJ sigma: %lf\n",SS.sig );
		fprintf( fsynopsis,"\tLJ energy: %lf\n",SS.eps );
		if(SS.dep==0){fprintf( fsynopsis," --- No short range attractive potential.\n" ); }
		if(SS.dep==1) {fprintf( fsynopsis," --- AO potential turned on.\n" );fprintf( fsynopsis,"\tAO potential depth: %lf\n",SS.depth );fprintf( fsynopsis,"\tAO potential range: %lf\n",SS.range );}
		if(SS.dep==2) {fprintf( fsynopsis," --- Square attractive potential turned on.\n" );fprintf( fsynopsis,"\t Square potential depth: %lf\n",SS.depth );fprintf( fsynopsis,"\t Square potential range: %lf\n",SS.range ); }
		fprintf( fsynopsis,"\tFixed distance from bottom wall (if type==DUMBBELL_NEARWALL): %lf\n",SS.fixDist );
		fprintf( fsynopsis,"\tAverage run time: %lf\n",SS.runTime );
		fprintf( fsynopsis,"\tAverage tumble time: %lf\n",SS.tumbleTime );
		fprintf( fsynopsis,"\tMagnetic moment strength: %lf\n",SS.MAGMOM );

		fprintf( fsynopsis,"\nOutput variables:\n" );
		fprintf( fsynopsis,"Debug mode: %i\n",DBUG );
		fprintf( fsynopsis,"Print detailed trajectories every %i time steps\n",out.TRAJOUT );
		fprintf( fsynopsis,"Number of species with detailed output: %i\n",out.printSP );
		fprintf( fsynopsis,"Print coarse data every %i time steps\n",out.COAROUT );
		fprintf( fsynopsis,"Print flow data: %i\n",out.FLOWOUT );
        fprintf( fsynopsis,"Print velocity data: %i\n",out.VELOUT );
		fprintf( fsynopsis,"Print flow around first swimmer data: %i\n",out.SWFLOWOUT );
		fprintf( fsynopsis,"Print averaged flow data: %i\n",out.AVVELOUT );
		fprintf( fsynopsis,"Print averaged orientation data: %i\n",out.AVORIOUT );
		fprintf( fsynopsis,"Print energy data: %i\n",out.ENOUT );
		fprintf( fsynopsis,"Print director and scalar order parameter fields: %i\n",out.ORDEROUT );
		fprintf( fsynopsis,"Print tensor order parameter data: %i\n",out.QTENSOUT );
		fprintf( fsynopsis,"Print reciprocal space order parameter data: %i\n",out.QKOUT );
		fprintf( fsynopsis,"Print averaged order parameter data: %i\n",out.AVSOUT );
		fprintf( fsynopsis,"Print averaged center of mass of multiphase species: %i\n",out.AVCOMOUT );
		fprintf( fsynopsis,"Print standard deviation of density data: %i\n",out.DENSOUT );
		fprintf( fsynopsis,"Print averaged enstrophy data: %i\n",out.ENSTROPHYOUT );
		fprintf( fsynopsis,"Velocity-velocity correlation:\t\t%d\n",out.CVVOUT);
		fprintf( fsynopsis,"Director-director correlation:\t\t%d\n",out.CNNOUT);
		fprintf( fsynopsis,"Vorticity-vorticity correlation:\t\t%d\n",out.CWWOUT);
		fprintf( fsynopsis,"Density-density correlation:\t\t%d\n",out.CDDOUT);
		fprintf( fsynopsis,"Order-order correlation:\t\t%d\n",out.CSSOUT);
		fprintf( fsynopsis,"Phase-phase (binary fluid) correlation:\t\t%d\n",out.CPPOUT);
		fprintf( fsynopsis,"Binder cumulant:\t\t%d --- bin size:\t\t%d\n",out.BINDER,out.BINDERBIN);
		fprintf( fsynopsis,"How often  solids' trajectories outputted: %i\n",out.SOLOUT );
		fprintf( fsynopsis,"Print distributions:\n" );
		fprintf( fsynopsis,"\tVel: %d\n\tSpeed: %d\n\tVorticity: %d\n\tEnstrophy: %d\n\tDirector: %d\n\tScalar order parameter: %d\n\tDensity: %d\n",out.HISTVELOUT,out.HISTSPEEDOUT,out.HISTVORTOUT,out.HISTENSTROUT,out.HISTDIROUT,out.HISTSOUT,out.HISTNOUT );
		fprintf( fsynopsis,"Synopsis of Simulation: %i\n",out.SYNOUT );
	}
}

///
/// @brief Outputs data for velocity probability distributions.
///
/// This function prints velocity probability distributions for a given number of bins to produce a histogram.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param vel This is the velocity in three dimensions for the corresponding bin it resides in.
/// @param minRange This is the minimum velocity.
/// @param maxRange This is the maximum velocity.
/// @param t This is the time.
/// @see outputResults()
///
void histVelout( FILE *fout,int vel[_3D][BINS],double minRange,double maxRange,double t ) {
	int i;
	double dv = (maxRange-minRange)/((float)BINS - 1.0);
	fprintf( fout,"%12.5e\n",t );
	for( i=0; i<BINS; i++ ) fprintf( fout,"\t%12.5e\t%7i\t%7i\t%7i\n",minRange+i*dv,vel[0][i],vel[1][i],vel[2][i] );
	fprintf( fout,"\n" );
	#ifdef FFLSH
		fflush(fout);
	#endif
}

///
/// @brief Outputs data for speed probability distributions.
///
/// This function prints speed probability distributions for a given number of bins to produce a histogram.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param speed This is the speed for the corresponding bin it resides in.
/// @param minRange This is the minimum speed.
/// @param maxRange This is the maximum speed.
/// @param t This is the time.
/// @see outputResults()
///
void histSpeedout( FILE *fout,int speed[BINS],double minRange,double maxRange,double t ) {
	int i;
	double dv = (maxRange-minRange)/((float)BINS - 1.0);

	fprintf( fout,"%12.5e\n",t );
	for( i=0; i<BINS; i++ ) fprintf( fout,"\t%12.5e\t%7i\n",minRange+i*dv,speed[i] );

	fprintf( fout,"\n" );
	#ifdef FFLSH
		fflush(fout);
	#endif
}

///
/// @brief Outputs data for vorticity probability distributions.
///
/// This function prints vorticity probability distributions for a given number of bins to produce a histogram.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param vort This is the vorticity in three dimensions for the corresponding bin it resides in.
/// @param minRange This is the minimum vorticity.
/// @param maxRange This is the maximum vorticity.
/// @param t This is the time.
/// @see outputResults()
///
void histVortout( FILE *fout,int vort[_3D][BINS],double minRange,double maxRange,double t ) {
	int i;
	double dw = (maxRange-minRange)/((float)BINS - 1.0);

	fprintf( fout,"%12.5e\n",t );
	for( i=0; i<BINS; i++ ) fprintf( fout,"\t%12.5e\t%7i\t%7i\t%7i\n",minRange+i*dw,vort[0][i],vort[1][i],vort[2][i] );

	fprintf( fout,"\n" );
	#ifdef FFLSH
		fflush(fout);
	#endif
}

///
/// @brief Outputs data for enstrophy probability distributions.
///
/// This function prints enstrophy probability distributions for a given number of bins to produce a histogram.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param enstrophy This is the enstrophy for the corresponding bin it resides in.
/// @param minRange This is the minimum enstrophy.
/// @param maxRange This is the maximum enstrophy.
/// @param t This is the time.
/// @see outputResults()
///
void histEnstrout( FILE *fout,int enstrophy[BINS],double minRange,double maxRange,double t ) {
	int i;
	double dw2 = (maxRange-minRange)/((float)BINS - 1.0);

	fprintf( fout,"%12.5e\n",t );
	for( i=0; i<BINS; i++ ) fprintf( fout,"\t%12.5e\t%7i\n",minRange+i*dw2,enstrophy[i] );

	fprintf( fout,"\n" );
	#ifdef FFLSH
		fflush(fout);
	#endif
}

///
/// @brief Outputs data for director orientation probability distributions.
///
/// This function prints director orientation probability distributions for a given number of bins to produce a histogram.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param dir This is the director orientation in three dimensions for the corresponding bin it resides in.
/// @param minRange This is the minimum enstrophy.
/// @param maxRange This is the maximum enstrophy.
/// @param t This is the time.
/// @see outputResults()
///
void histDirout( FILE *fout,int dir[_3D][BINS],double minRange,double maxRange,double t ) {
	int i;
	double dn = (maxRange-minRange)/((float)BINS - 1.0);

	fprintf( fout,"%12.5e\n",t );
	for( i=0; i<BINS; i++ ) fprintf( fout,"\t%12.5e\t%7i\t%7i\t%7i\n",minRange+i*dn,dir[0][i],dir[1][i],dir[2][i] );

	fprintf( fout,"\n" );
	#ifdef FFLSH
		fflush(fout);
	#endif
}

///
/// @brief Outputs data for scalar order parameter probability distributions.
///
/// This function prints scalar order parameter probability distributions for a given number of bins to produce a histogram.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param S This is the scalar order parameter for the corresponding bin it resides in.
/// @param minRange This is the minimum scalar order parameter.
/// @param maxRange This is the maximum scalar order parameter.
/// @param t This is the time.
/// @see outputResults()
///
void histSout( FILE *fout,int S[BINS],double minRange,double maxRange,double t ) {
	int i;
	double dS = (maxRange-minRange)/((float)BINS - 1.0);

	fprintf( fout,"%12.5e\n",t );
	for( i=0; i<BINS; i++ ) fprintf( fout,"\t%12.5e\t%7i\n",minRange+i*dS,S[i] );

	fprintf( fout,"\n" );
	#ifdef FFLSH
		fflush(fout);
	#endif
}

///
/// @brief Outputs data for particle density probability distributions.
///
/// This function prints particle density probability distributions for a given number of bins to produce a histogram.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param dens This is the particle density for the corresponding bin it resides in.
/// @param minRange This is the minimum particle density.
/// @param maxRange This is the maximum particle density.
/// @param t This is the time.
/// @see outputResults()
///
void histNout( FILE *fout,int dens[BINS],double minRange,double maxRange,double t ) {
	int i;
	double dp = (maxRange-minRange)/((float)BINS - 1.0);

	fprintf( fout,"%12.5e\n",t );
	for( i=0; i<BINS; i++ ) fprintf( fout,"\t%12.5e\t%7i\n",minRange+i*dp,dens[i] );

	fprintf( fout,"\n" );
	#ifdef FFLSH
		fflush(fout);
	#endif
}

///
/// @brief Outputs energy data.
///
/// This function collates and prints nematic, kinetic and boundary energy contributions.
/// The thermal energy and rotational energy at the boundary are also output.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param pp This is a pointer to MPCD particle indices.
/// @param pSP This is a pointer to subpopulation indices.
/// @param WALL This is a pointer to boundary position information.
/// @param t This is the time.
/// @param KBT This is a pointer that fetches thermal energy and is only considered if `COLL_TYPE` is set to thermal collisions in input.json.
/// @param wmf This is the mean-field potential.
/// @see outputResults()
///
void enout( FILE *fout,particleMPC *pp,spec *pSP,bc WALL[],double t,double KBT,double wmf ) {
	int i,j,k;
	double MPC_K=0.0,BC_K=0.0,BC_R=0.0,TE=0.0,E=0.0;

	for( i=0; i<GPOP; i++ ) {
		E = 0.0;
		for( j=0; j<DIM; j++ ) E += (pp+i)->V[j] * (pp+i)->V[j];
		E *= 0.5 * (pSP+(pp+i)->SPID)->MASS;
		MPC_K += E;
		TE += E;
	}
	for( i=0; i<NBC; i++ ) if( WALL[i].DSPLC ) {
		E = 0.0;
		//Kinetic
		for( j=0; j<DIM; j++ ) E += WALL[i].V[j]*WALL[i].V[j];
		E *= 0.5 * WALL[i].MASS;
		BC_K += E;
		TE += E;
		//Rotational
		E = 0.0;
		for( j=0; j<_3D; j++ ) for( k=0; k<_3D; k++ ) E += WALL[i].L[j]*WALL[i].I[j][k]*WALL[i].L[j];
		E *= 0.5;
		BC_R += E;
		TE += E;
	}
	//Nematic potential energy
	TE += wmf;
// 	KBT = TEMP( pp,pSP,WALL );

	fprintf( fout,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",t,MPC_K,wmf,BC_K,BC_R,TE,KBT );
	#ifdef FFLSH
		fflush(fout);
	#endif
}

///
/// @brief Outputs orientation interaction energy and Kinetic energy within a cell.
///
/// This function calculates and prints the orientation interaction energy and kinetic energy within a cell.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param CL This is a pointer to the co-ordinates and cell of each particle.
/// @param SP This is a pointer to species subpopulation indices.
/// @param LC This is a flag that states if the system is a liquid crystal.
/// @see outputResults
///
void enfieldout( FILE *fout,cell ***CL,spec *SP,int LC ) {
	int a,b,c,d,id;
	double enK,wmf,S,un,DIR[_3D],u[_3D],m;
	double invdim=1./((double)DIM);
	particleMPC *tmpc;	//Temporary pointer to MPCD particles

	for( a=0; a<XYZ[0]; a++ ) for( b=0; b<XYZ[1]; b++ ) for( c=0; c<XYZ[2]; c++ ) {
		wmf=0.;
		enK=0.;
		if( CL[a][b][c].POPSRD > 1 ) {
			S = CL[a][b][c].S;
			for( d=0; d<DIM; d++ ) DIR[d] = CL[a][b][c].DIR[d];
			tmpc = CL[a][b][c].pp;
			while( tmpc != NULL ) {
				id = tmpc->SPID;
				//Nematic energy
				if( LC ) {
					for( d=0; d<DIM; d++ ) u[d] = tmpc->U[d];
					un = dotprod( u,DIR,DIM );
					wmf += ( S*un*un  + (1.-S)*invdim )*( (SP+id)->MFPOT );
				}
				//Kinetic energy
				m = (SP+id)->MASS;
				for( d=0; d<DIM; d++ ) u[d] = tmpc->V[d];
				un = dotprod( u,u,DIM );
				enK += 0.5*m*un;

				tmpc = tmpc->next;
			}
		}
		fprintf( fout, "%5i\t%5i\t%5i\t%e\t%e\n",a,b,c,enK,wmf );
	}
}

///
/// @brief Outputs average cell orientation interaction energy with neighbouring cells.
///
/// This function calculates the total director of all the cells under consideration and computes the energy based on the local director and tensor order parameter.
/// The value computed is the average cell orientation interaction energy with neighbouring cells and is printed out.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param t This is time.
/// @param CL This is a pointer to the co-ordinates and cell of each particle.
/// @param SP This is a pointer to species subpopulation indices.
/// @param LC This is a flag that states if the system is a liquid crystal.
/// @see outputResults()
///
void enneighboursout( FILE *fout,double t,cell ***CL,spec *SP,int LC ) {
	int a,b,c,d,id;
	double avMFPOT,wmf,un,sumWMF;
	double local_DIR[DIM],nnn_DIR[DIM],local_S,nnn_S;
	double **Q,eigval[DIM];
	//double invDIM=1.0/((double)DIM);
	particleMPC *tmpc;	//Temporary pointer to MPCD particles

	sumWMF=0.;
	//Allocate memory for tensor order parameter Q
	Q = calloc ( DIM, sizeof( *Q ) );
	for( a=0; a<DIM; a++ ) Q[a] = calloc ( DIM, sizeof( *Q[a] ) );
	for( a=0; a<DIM; a++ ) for( b=0; b<DIM; b++ ) Q[a][b] = 0.0;

	if( LC ) for( a=0; a<XYZ[0]; a++ ) for( b=0; b<XYZ[1]; b++ ) for( c=0; c<XYZ[2]; c++ ) if( CL[a][b][c].POPSRD>1 ) {
		//Local values
		for( d=0; d<DIM; d++ ) local_DIR[d]=CL[a][b][c].DIR[d];
		local_S=CL[a][b][c].S;
		avMFPOT=0.;
		if( CL[a][b][c].POPSRD > 1 ) {
			tmpc = CL[a][b][c].pp;
			while( tmpc != NULL ) {
				id = tmpc->SPID;
				avMFPOT+=(SP+id)->MFPOT;
				tmpc = tmpc->next;
			}
			avMFPOT/=(float)CL[a][b][c].POPSRD;
		}
		//Next-nearest values
		//Calculate the tensor order parameter from the cell and its neighbous
		tensOrderParamNNN( CL,Q,LC,a,b,c );
		// From the tensor order parameter find eigenvalues and vectors --- Q is written over as normalized eigenvectors
		solveEigensystem( Q,DIM,eigval );
		//The scalar order parameter is the largest eigenvalue which is given first, ie eigval[0]
		if(DIM==_3D) nnn_S = -1.*(eigval[1]+eigval[2]);
		else nnn_S=eigval[0];

		if( nnn_S<1./(1.-DIM) ){
			printf("Warning: Global scalar order parameter < 1/(1-DIM)\n");
			printf("Eigenvalues=");
			pvec(eigval,DIM);
			printf("Eigenvectors=");
			for( d=0; d<DIM; d++ ) pvec(Q[d],DIM);
		}
		// The director is the eigenvector corresponding to the largest eigenvalue
		for( d=0; d<DIM; d++ ) nnn_DIR[d] = Q[0][d];

		//We have the local director and order parameter and the next-nearest equivalents
		//From these we can see the energy
		un = dotprod( local_DIR,nnn_DIR,DIM );
		wmf = local_S*un*un;
		//wmf += (1.-local_S)*invDIM;	//Don't include the constant (wrt u.n) term
		wmf*=avMFPOT;
		sumWMF+=wmf;
	}
	fprintf( fout, "%12.5e\t%12.5e\n",t,sumWMF );

	for( d=0; d<DIM; d++ ) free( Q[d] );
	free( Q );
}

///
/// @brief Outputs coarse grained average velocity to file.
///
/// This function simply prints the coarse grained average velocity to the output data file.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param t This is time.
/// @param vel This is the velocity in three dimensions.
/// @param KBT This is a pointer that fetches thermal energy and is only considered if `COLL_TYPE` is set to thermal collisions in input.json.
/// @see outputResults()
///
void avvelout( FILE *fout,double t,double vel[_3D],double KBT ) {
	fprintf( fout, "%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\n",t,vel[0],vel[1],vel[2],KBT );
	#ifdef FFLSH
		fflush(fout);
	#endif
}

///
/// @brief Outputs average orientation to file.
///
/// This function simply prints the average orientation to the output data file.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param t This is time.
/// @param ori This is the velocity in three dimensions.
/// @see outputResults()
///
void avoriout( FILE *fout,double t,double ori[_3D] ) {
	fprintf( fout, "%12.5e\t%12.5e\t%12.5e\t%12.5e\n",t,ori[0],ori[1],ori[2] );
	#ifdef FFLSH
		fflush(fout);
	#endif
}

///
/// @brief Outputs coarse grained average velocity and velocity gradient tensor to file.
///
/// This function simply prints the coarse grained average velocity and velocity gradient tensor to the output data file.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param t This is time.
/// @param vel This is the velocity in three dimensions.
/// @param KBT This is a pointer that fetches thermal energy and is only considered if `COLL_TYPE` is set to thermal collisions in input.json.
/// @param gradVel This is velocity gradient tensor.
/// @see outputResults()
///
void avveloutWithGradVel( FILE *fout,double t,double vel[_3D],double KBT,double gradVel[_3D][_3D] ) {
	fprintf( fout, "%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t",t,vel[0],vel[1],vel[2],KBT );
	fprintf( fout, "%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\n",gradVel[0][0],gradVel[0][1],gradVel[0][2],gradVel[1][0],gradVel[1][1],gradVel[1][2],gradVel[2][0],gradVel[2][1],gradVel[2][2] );
	#ifdef FFLSH
		fflush(fout);
	#endif
}

///
/// @brief Outputs coarse grained scalar order parameter to file.
///
/// This function simply prints the coarse grained scalar order parameter to the output data file.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param t This is time.
/// @param S This is a pointer to the scalar order parameter.
/// @param S4 This is a pointer to the fourth moment of the scalar order parameter.
/// @param DIR This is a pointer that fetches the director orientation in three dimensions.
/// @see outputResults()
///
void avsout( FILE *fout,double t,double S,double S4,double DIR[] ) {
	fprintf( fout, "%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\n",t,S,S4,DIR[0],DIR[1],DIR[2] );
	#ifdef FFLSH
		fflush(fout);
	#endif
}

///
/// @brief Outputs centres of masses of multiphase fluid to file.
///
/// This function simply prints the centres of masses of multiphase fluid to file..
///
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param t This is time.
/// @param CoM This is a pointer to the center of mass of a species from a multiphase fluid.
/// @param stdCoM This is a pointer to the std tensor of the center of mass for a given species from a multiphase fluid.
/// @see outputResults()
///
void avcomout( FILE *fout,double t,double CoM[MAXSPECI][_3D],double stdCoM[MAXSPECI][_3D][_3D]) {
	int i;

	fprintf( fout, "%12.5e",t );
	for (i=0;i<NSPECI;i++) fprintf( fout, "\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e",CoM[i][0],CoM[i][1],CoM[i][2],stdCoM[i][0][0],stdCoM[i][0][1],stdCoM[i][0][2],stdCoM[i][1][0],stdCoM[i][1][1],stdCoM[i][1][2],stdCoM[i][2][0],stdCoM[i][2][1],stdCoM[i][2][2] );
	fprintf( fout, "\n" );
	#ifdef FFLSH
		fflush(fout);
	#endif
}

///
/// @brief Outputs coarse grained standard deviations of density to file.
///
/// This function simply prints the coarse grained standard deviations of density to the output data file.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param t This is time.
/// @param stdN This is a pointer to the standard deviations of density.
/// @see outputResults()
///
void densSTDout( FILE *fout,double t,double stdN ) {
	fprintf( fout, "%12.5e\t%12.5e\n",t,stdN );
	#ifdef FFLSH
		fflush(fout);
	#endif
}

///
/// @brief Outputs coarse grained average enstrophy to file.
///
/// This function simply prints the coarse grained average enstrophy to the output data file.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param t This is time.
/// @param E This is a pointer to the average enstrophy.
/// @see outputResults()
///
void avenstrophyout( FILE *fout,double t,double E ) {
	fprintf( fout, "%12.5e\t%12.5e\n",t,E );
	#ifdef FFLSH
		fflush(fout);
	#endif
}

///
/// @brief Outputs binder cumulants for each given timestep.
///
/// This function simply prints the binder cumulants to the output data file.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param t This is time.
/// @param UL This is a pointer to the binder cumulant.
/// @see outputResults()
///
void binderout( FILE *fout,double t,double UL ) {
	fprintf( fout, "%12.5e\t%12.5e\n",t,UL );
	#ifdef FFLSH
		fflush(fout);
	#endif
}

///
/// @brief Outputs flow field data calculated from cell velocities.
///
/// This function computes the centre of mass velocity in the `x`,`y`, and `z` directions for each cell, as well as the average velocity.
/// The sum of all centre of mass velocities are computed into an average in the `x`, `y`, and `z` directions.
/// The centre of mass velocities and averages are printed to the output file.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param CL This is a pointer to the co-ordinates and cell of each particle.
/// @param interval is the time interval used for normalisation.
/// @param t This is the time step.
/// @see outputResults()
///
void flowout( FILE *fout,cell ***CL,int interval, double t) {
	int h=0,i=0,j=0,k=0;
	double av[_3D];
    zerovec(av, _3D);
	// for( i=0; i<_3D; i++ ) av[i] = 0.0;
	double dint = (double)interval;

	for( i=0; i<_3D; i++ ) av[i] = 0.0;			//In some compilers, this routine might write nonsense in the z-component in 2D otherwise
	for( i=0; i<XYZ[0]; i++ ) for( j=0; j<XYZ[1]; j++ ) for( k=0; k<XYZ[2]; k++ ) {
	// for( i=0; i<XYZ_P1[0]; i++ ) for( j=0; j<XYZ_P1[1]; j++ ) for( k=0; k<XYZ_P1[2]; k++ ) {
		for( h=0; h<DIM; h++ ) av[h] = CL[i][j][k].FLOW[h]/dint;		//Normalize the sum to get the average
		fprintf( fout,"%12.5e\t", t); // print time
		fprintf( fout, "%5d\t%5d\t%5d\t",i,j,k );
		fprintf( fout, "%12.5e\t%12.5e\t%12.5e\n",av[0],av[1],av[2] );
		// for( h=0; h<DIM; h++ ) CL[i][j][k].FLOW[h] = 0.0;		//Reset sum
		for( h=0; h<_3D; h++ ) CL[i][j][k].FLOW[h] = 0.0;		//Reset sum; no harm in doing it in 3D. Might store nonsense in the z-component in 2D otherwise
	}
	#ifdef FFLSH
		fflush(fout);
	#endif
}

///
/// @brief Outputs velocity field data calculated from cell velocities.
///
/// This method is very similar to flowout(), but instead of printing the time averaged flow velocity, it prints the
/// explicit current velocity of the cell.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param CL This is a pointer to the co-ordinates and cell of each particle.
/// @param t This is the time step.
/// @see flowout()
/// @see outputResults()
///
void velout( FILE *fout,cell ***CL, double t) {
/*
    Prints an instantaneous vcm field of the cells
*/
    int h,i,j,k;
    double vel[_3D];

	for( i=0; i<_3D; i++ ) vel[i] = 0.0;			//In some compilers, this routine might write nonsense in the z-component in 2D otherwise
    for( i=0; i<XYZ[0]; i++ ) for( j=0; j<XYZ[1]; j++ ) for( k=0; k<XYZ[2]; k++ ) {
        for( h=0; h<DIM; h++ ) vel[h] = CL[i][j][k].VCM[h];
        fprintf( fout,"%12.5e\t", t); // print time
        fprintf( fout, "%5d\t%5d\t%5d\t",i,j,k );
        fprintf( fout, "%12.5e\t%12.5e\t%12.5e\n",vel[0],vel[1],vel[2] );
    }
#ifdef FFLSH
    fflush(fout);
#endif
}

///
/// @brief Outputs flow field around the first swimmer, with a 'hijack' method described in sumSWFLOW().
///
/// This function computes the average of the centre-of-mass velocity of located at each grid position relative to the first swimmer.
/// Depending on swimmer orientation, cells will contain more or less information, especially at the corners: a weight is given to counter this effect.
/// The centre of mass velocities and averages are printed to the output file.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param CL This is a pointer to the co-ordinates and cell of each particle.
/// @param interval This is the time interval used for normalisation.
/// @param t This is the time step.
/// @see outputResults()
/// @see sumSWFLOW()
///
void swflowout( FILE *fout,cell ***CL,int interval, double t) {
	int h,i,j,k;
	double av[_3D];
	for( i=0; i<_3D; i++ ) av[i] = 0.0;

	for( i=0; i<XYZ[0]; i++ ) for( j=0; j<XYZ[1]; j++ ) for( k=0; k<XYZ[2]; k++ ) {
		if (CL[i][j][k].SWFLOW[3]>0.0) for( h=0; h<DIM; h++ ) av[h] = CL[i][j][k].SWFLOW[h]/CL[i][j][k].SWFLOW[3];		//Normalize the sum to get the average, with weights contained in SWFLOW[3]
		else for( h=0; h<DIM; h++ ) av[h] = 0.0;
		fprintf( fout,"%12.5e\t", t); // print time
		fprintf( fout, "%5d\t%5d\t%5d\t",i,j,k );
		fprintf( fout, "%12.5e\t%12.5e\t%12.5e\n",av[0],av[1],av[2] );
		for( h=0; h<7; h++ ) CL[i][j][k].SWFLOW[h] = 0.0;		//Reset sum
	}
	#ifdef FFLSH
		fflush(fout);
	#endif
}

///
/// @brief Prints solid object data.
///
/// This function simply prints the solid object data to the output data file.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param WALL This is a pointer to all of the walls (boundary conditions).
/// @param t This is the time step.
/// @see outputResults()
///
void solidout( FILE *fout,bc WALL,double t ) {
	fprintf( fout,"%12.5e\t",t );
	fprintf( fout,"%12.5e\t%12.5e\t%12.5e\t",WALL.Q[0],WALL.Q[1],WALL.Q[2] );
	fprintf( fout,"%12.5e\t%12.5e\t%12.5e\t",WALL.V[0],WALL.V[1],WALL.V[2] );
	fprintf( fout,"%12.5e\t%12.5e\t%12.5e\t",WALL.O[0],WALL.O[1],WALL.O[2] );
	fprintf( fout,"%12.5e\t%12.5e\t%12.5e\n",WALL.L[0],WALL.L[1],WALL.L[2] );
	#ifdef FFLSH
		fflush(fout);
	#endif
}

///
/// @brief Outputs topological charge and defect angle data to file for 2D defects only.
///
/// This function only works for 2D defects and serves two functions. If `TOPOUT` is flagged in input.json, the topological charge and angle of defect are calculated by `topoChargeLocal` and `topoAngleLocal` and the results are printed to a .dat file.
/// If `DEFECTOUT` is flagged in input.json, the charge and angle are are grouped into clusters of non-zero values to identify defects. Each defect is given an identity, `defID`, with the topological charge and angle averaged for the whole defect.
///
/// @param ftopo This is a pointer to the topological charge .dat file to be produced.
/// @param TOPOOUT This is a flag that determines whether the topological charge and angle should be printed, specified in input.json.
/// @param fdefect This is a pointer to the topological defect .dat file to be produced.
/// @param DEFECTOUT This is a flag that determines whether the defect topological charge and angle should be printed, specified in input.json.
/// @param t This is the time step.
/// @param CL This is a pointer to the co-ordinates and cell of each particle.
/// @param tolD This is a cutoff that acts as the tolerance of the defect tracker.
/// @see topoChargeLocal()
/// @see topoAngleLocal()
/// @see outputResults()
///
void topoChargeAndDefectsOut( FILE *ftopo,int TOPOOUT,FILE *fdefect,int DEFECTOUT,double t,cell ***CL,double tolD){
	//FIXME:
	int i,j,k,cntD;
	double m,cmx,cmy,avC,avAx,avAy,avA;

	double topoC[XYZ[0]][XYZ[1]]; //init topo charge array
	for( i=0; i<XYZ[0]; i++ ) for( j=0; j<XYZ[1]; j++ ) topoC[i][j] = 0.0;
	double topoAngle[XYZ[0]][XYZ[1]]; //init topo angle array
	for( i=0; i<XYZ[0]; i++ ) for( j=0; j<XYZ[1]; j++ ) topoAngle[i][j] = 0.0;

	//loop through non-CB boundary cells and calculate topo charge
	for( i=1; i<XYZ[0]-1; i++ ) for( j=1; j<XYZ[1]-1; j++ ) topoC[i][j] = topoChargeLocal(CL, i, j, 0);
	//loop through non-CB boundary cells and calculate topo angle
	for( i=2; i<XYZ[0]-2; i++ ) for( j=2; j<XYZ[1]-2; j++ ){
		//FIXME: Too lazy to handle derivatives properly at the boundaries, so we just ignoring another layer there instead. Oopsies. Same goes for the loop above.
		//Tyler: I think this is reasonable since we don't want to assume PBCs
		if( fabs(topoC[i][j])>TOL ) topoAngle[i][j] = topoAngleLocal(CL, i, j, 0, topoC[i][j]);
	}
	if( TOPOOUT ) {
		//Output
		for( i=0; i<XYZ[0]; i++ ) for( j=0; j<XYZ[1]; j++ ) for( k=0; k<XYZ[2]; k++ ) {
			fprintf( ftopo,"%.2f\t%5d\t%5d\t%5d\t",t,i,j,k );
			if( CL[i][j][k].POPSRD == 0 ) fprintf( ftopo, "%06.3f\t%12.5e\n", 0.0, 0.0);
			else fprintf( ftopo, "%06.3f\t%12.5e\n",topoC[i][j], topoAngle[i][j]);
		}
		#ifdef FFLSH
			fflush(ftopo);
		#endif
	}
	if( DEFECTOUT ) {
		// Group topological charges into "clusters" of nearby non-zero values
		int defID[XYZ[0]][XYZ[1]]; //ID value for each "cluster"
		for( j=0; j<XYZ[1]; j++ ) for( i=0; i<XYZ[0]; i++ ) defID[i][j]=0; //ID is zero everywhere there is NO "cluster"
		cntD=0;	// Counts the number of "clusters" (blurry defects)
		for( j=0; j<XYZ[1]; j++ ) for( i=0; i<XYZ[0]; i++ ) if( fabs(fabs(topoC[i][j])-0.5)<tolD ) {
			// First MPCD cell
			if( i==0 && j==0 ) {
				//Simply set to new cluster
				cntD+=1;
				defID[i][j]=cntD;
			}
			// First row
			else if( j==0 ) {
				//Check backwards (that has an ID already and the same topological charge [product positive])
				if( defID[i-1][j]>0 && topoC[i-1][j]*topoC[i][j]>0.0 ) defID[i][j]=defID[i-1][j];
				else { //New cluster
					cntD+=1;
					defID[i][j]=cntD;
				}
			}
			// First column
			else if( i==0 ) {
				//Check directly below
				if( defID[i][j-1]>0 && topoC[i][j]*topoC[i][j]>0.0 ) defID[i][j]=defID[i][j-1];
				//Check below and forward
				else if( defID[i+1][j-1]>0 && topoC[i][j]*topoC[i][j]>0.0 ) defID[i][j]=defID[i+1][j-1];
				else { //New cluster
					cntD+=1;
					defID[i][j]=cntD;
				}
			}
			// Bulk
			else {
				//Check below and back
				if( defID[i-1][j-1]>0 && topoC[i-1][j-1]*topoC[i][j]>0.0 ) defID[i][j]=defID[i-1][j-1];
				//Check directly below
				else if( defID[i][j-1]>0 && topoC[i][j-1]*topoC[i][j]>0.0 ) defID[i][j]=defID[i][j-1];
				//Check below and forward
				else if( defID[i+1][j-1]>0 && topoC[i+1][j-1]*topoC[i][j]>0.0 ) defID[i][j]=defID[i+1][j-1];
				//Check directly back
				else if( defID[i-1][j]>0 && topoC[i-1][j]*topoC[i][j]>0.0 ) defID[i][j]=defID[i-1][j];
				else { //New cluster
					cntD+=1;
					defID[i][j]=cntD;
				}
			}
		}
		//Average each cluster's values and output their positions
		fprintf( fdefect,"%.2f\t%d\n",t,cntD );
		for( k=1; k<=cntD; k++ ) {
            cmx=0.0; // CoM pos
			cmy=0.0;
			avC=0.0; // charge
            avAx=0.0; // x component of vector angle (for average)
            avAy=0.0; // y component of vector angle (for average)
            avA=0.0; // reset average angle (to dump)
			m=0.0; // count

			for( i=0; i<XYZ[0]; i++ ) for( j=0; j<XYZ[1]; j++ ) if( defID[i][j]==k ) {
				m+=1.0;
				cmx+=(double)i + 0.5;
				cmy+=(double)j + 0.5;
				avC+=topoC[i][j];

				//Wrapping defect orientation on proper period prior to average
				double locAngle;
				if( topoC[i][j]>0.0 ){
					// +1/2 defects have a 2 pi symmetry so can be handled reasonably
					locAngle = fmod(topoAngle[i][j], 2.0*M_PI);
					locAngle = (locAngle < 0.0) ? (2.0*M_PI + locAngle) : locAngle;
				} else {
					// -1/2 defects need additional considerations due to 3-fold symmetry
					locAngle = fmod(topoAngle[i][j], 2.0*M_PI/3.0);
					locAngle = 3.0*((locAngle < 0.0) ? (2.0*M_PI/3.0 + locAngle) : locAngle);
				}
                avAx += cos(locAngle);
                avAy += sin(locAngle);
			}
			if( m>TOL ){
				cmx/=m;
				cmy/=m;
				avC/=m;
			}

			//Scaling back to proper period
			avA = atan2(-avAy, -avAx) + M_PI;
			avA = (avC>0.0) ? avA : avA/3.0;

			fprintf( fdefect,"%12.5e\t%12.5e\t%06.3f\t%12.5e\n",cmx,cmy,avC,avA );
		}
		fprintf( fdefect,"\n" );
		#ifdef FFLSH
			fflush(fdefect);
		#endif
	}
}

///
/// @brief Outputs disclination tensor data to file.
///
/// This function calculates and prints three-dimensional disclination tensor data to file based on methods by <a href="https://pubs.rsc.org/en/content/articlelanding/2022/SM/D1SM01584B">C Schimming and J Vinals</a>.
/// The disclination tensor is calculated from the tensor ordere parameter `Q`, which in turn is calculated by 'tensOrderParam'. 
///
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param t This is the time step.
/// @param CL This is a pointer to the co-ordinates and cell of each particle.
/// @param LC This is a flag that states if the system is a liquid crystal.
/// @see tensOrderParam()
/// @see outputResults()
///
void disclinationTensorOut( FILE *fout,double t,cell ***CL,int LC ) {
	printf( "Warning:\tdisclinationTensorOut() not yet implemented .\n" );

	int i,j,k;
	double **Q,**D;

	//Allocate memory for tensor order parameter Q and disclination tensor D
	Q = calloc ( _3D, sizeof( *Q ) );
	D = calloc ( _3D, sizeof( *D ) );
	for( i=0; i<_3D; i++ ) {
		Q[i] = calloc ( _3D, sizeof( *Q[i] ) );
		D[i] = calloc ( _3D, sizeof( *D[i] ) );
	}
	//Zero
	for( i=0; i<_3D; i++ ) for( j=0; j<_3D; j++ ) {
		Q[i][j] = 0.0;
		D[i][j] = 0.0;
	}

	for( i=0; i<XYZ[0]; i++ ) for( j=0; j<XYZ[1]; j++ ) for( k=0; k<XYZ[2]; k++ ) {
		// Find the local tensor order parameter
		tensOrderParam( &CL[i][j][k],Q,LC );
		//Calculate disclination tensor from Q

		// Output the disclination tensor
		// fprintf( fout,"%5d\t%5d\t%5d\t",i,j,k );
		// if( CL[i][j][k].POP == 0 ) fprintf( fout, "%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\n" ,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 );
		// else fprintf( fout, "%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\n", D[0][0],D[0][1],D[0][2],D[1][0],D[1][1],D[1][2],D[2][0],D[2][1],D[2][2] );
	}

	for( i=0; i<_3D; i++ ) {
		free( Q[i] );
		free( D[i] );
	}
	free( Q );
	free( D );
	#ifdef FFLSH
		fflush(fout);
	#endif
}

///
/// @brief Outputs scalar order parameter and director field data to file.
///
/// This function simply prints the scalar order parameter and director field to an output data file.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param t This is the time step.
/// @param CL This is a pointer to the co-ordinates and cell of each particle.
/// @param LC This is a flag that states if the system is a liquid crystal.
/// @see outputResults()
///
void orderout( FILE *fout,double t,cell ***CL,int LC ) {
	int i,j,k;
	for( i=0; i<XYZ[0]; i++ ) for( j=0; j<XYZ[1]; j++ ) for( k=0; k<XYZ[2]; k++ ) {
		//Output
		fprintf( fout,"%.2f\t",t );
		fprintf( fout,"%5d\t%5d\t%5d\t",i,j,k );
		if( CL[i][j][k].POPSRD == 0 ) fprintf( fout, "%12.5e\t%12.5e\t%12.5e\t%12.5e\n",0.0,0.0,0.0,0.0 );
		else fprintf( fout, "%12.5e\t%12.5e\t%12.5e\t%12.5e\n",CL[i][j][k].DIR[0],CL[i][j][k].DIR[1],CL[i][j][k].DIR[2],CL[i][j][k].S );
	}
	#ifdef FFLSH
		fflush(fout);
	#endif
}

///
/// @brief Outputs multiphase data to file.
///
/// This function simply prints multiphase data to an output data file.
/// The species type, colour, and multiphase phi parameter are printed.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param t This is the time step.
/// @param CL This is a pointer to the co-ordinates and cell of each particle.
/// @see outputResults()
///
void multiphaseout( FILE *fout,double t,cell ***CL ) {
	int i,j,k,n;
	for( i=0; i<XYZ[0]; i++ ) for( j=0; j<XYZ[1]; j++ ) for( k=0; k<XYZ[2]; k++ ) {
		//Output
		fprintf( fout,"%.2f\t",t );
		fprintf( fout,"%5d\t%5d\t%5d",i,j,k );
		for( n=0; n<NSPECI; n++ ) fprintf( fout, "\t\t%d",CL[i][j][k].SP[n] );
		fprintf( fout, "\n" );
	}
	#ifdef FFLSH
		fflush(fout);
	#endif
}

///
/// @brief Outputs pressure field data to file.
///
/// This function simply prints pressure field data to an output data file.
/// The streaming pressure `Ps` and collisional pressure `Pc` are the contributing factors computed.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param t This is the time step.
/// @param CL This is a pointer to the co-ordinates and cell of each particle.
/// @see outputResults()
///
void pressureout( FILE *fout,double t,cell ***CL ) {
	int i,j,k;

	for( i=0; i<XYZ[0]; i++ ) for( j=0; j<XYZ[1]; j++ ) for( k=0; k<XYZ[2]; k++ ) {
		//Output
		fprintf( fout,"%.2f\t",t );
		fprintf( fout,"%5d\t%5d\t%5d\t",i,j,k );
		// for( l=0; l<DIM; l++ ) for( m=0; m<DIM; m++ ) printf( "%lf\n",CL[i][j][k].Ps[l][m] );
		if( CL[i][j][k].POP == 0 ) fprintf( fout, "%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\n",0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 );
		else {
			//Print the pressure
			fprintf( fout, "%12.5e\t%12.5e\t%12.5e\t",CL[i][j][k].Ps[0][0]+CL[i][j][k].Pc[0][0],CL[i][j][k].Ps[0][1]+CL[i][j][k].Pc[0][1],CL[i][j][k].Ps[0][2]+CL[i][j][k].Pc[0][2] );
			fprintf( fout, "%12.5e\t%12.5e\t%12.5e\t",CL[i][j][k].Ps[1][0]+CL[i][j][k].Pc[1][0],CL[i][j][k].Ps[1][1]+CL[i][j][k].Pc[1][1],CL[i][j][k].Ps[1][2]+CL[i][j][k].Pc[1][2] );
			fprintf( fout, "%12.5e\t%12.5e\t%12.5e\n",CL[i][j][k].Ps[2][0]+CL[i][j][k].Pc[2][0],CL[i][j][k].Ps[2][1]+CL[i][j][k].Pc[2][1],CL[i][j][k].Ps[2][2]+CL[i][j][k].Pc[2][2] );
		}
	}
	#ifdef FFLSH
		fflush(fout);
	#endif
}

///
/// @brief Outputs species-specific pressure field data to file.
///
/// This function simply prints species-specific pressure field data to an output data file.
/// The streaming pressure `Ps` and collisional pressure `Pc` are the contributing factors computed.
/// The species that is the most present has its pressure ouputted.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param t This is the time step.
/// @param CL This is a pointer to the co-ordinates and cell of each particle.
/// @see outputResults()
///
void sppressureout( FILE *fout,double t,cell ***CL ) {
	int i,j,k,pop,maxp,maxid;
	double avP[NSPECI],stdP[NSPECI],cnt[NSPECI],cutoff=0.7,P;


	fprintf( fout,"%.2f",t );
	for( i=0; i<XYZ[0]; i++ ) for( j=0; j<XYZ[1]; j++ ) for( k=0; k<XYZ[2]; k++ ) {

		if( CL[i][j][k].POP > 0 ){
			pop=CL[i][j][k].POP;
			maxp=0;

			for (int q=0;q<NSPECI;q++) {
				if (CL[i][j][k].SP[q]>maxp){
					maxp=CL[i][j][k].SP[q];
					maxid=q;
				}
			}	

			if (cutoff<=(float)maxp/(float)pop){
				P=0;
				cnt[maxid]++;
				for (int q=0;q<DIM;q++) P+=CL[i][j][k].Ps[q][q]+CL[i][j][k].Pc[q][q];
				P/=(float)DIM;
				avP[maxid]+=P;
				stdP[maxid]+=P*P;
			}
		}
	}

	for (i=0;i<NSPECI;i++){
		avP[i]/=cnt[i];
		stdP[i]/=cnt[i];
		avP[i]=fabs(avP[i]);
		stdP[i]=sqrt(fabs(stdP[i]-avP[i]*avP[i]));
		fprintf( fout,"\t%12.5e\t%12.5e",avP[i],stdP[i]);
	}
	fprintf( fout,"\n");
	#ifdef FFLSH
		fflush(fout);
	#endif
}

///
/// @brief Outputs tensor order parameter `Q` data to file.
///
/// This function obtains the tensor order parameter from `tensOrderParam`.
/// The tensor order parameter value is then printed to an output data file.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param t This is the time step.
/// @param CL This is a pointer to the co-ordinates and cell of each particle.
/// @param LC This is a flag that states if the system is a liquid crystal.
/// @see outputResults()
///
void orderQout( FILE *fout,double t,cell ***CL,int LC ) {
	int i,j,k;
	double **Q;

	//Allocate memory for tensor order parameter Q
	Q = calloc ( _3D, sizeof( *Q ) );
	for( i=0; i<_3D; i++ ) Q[i] = calloc ( _3D, sizeof( *Q[i] ) );
	for( i=0; i<_3D; i++ ) for( j=0; j<_3D; j++ ) Q[i][j] = 0.0;

	for( i=0; i<XYZ[0]; i++ ) for( j=0; j<XYZ[1]; j++ ) for( k=0; k<XYZ[2]; k++ ) {
		//// Find the tensor order parameter based on self, neighbours and next-nearest neighbours (NNN)
		//tensOrderParamNNN( CL,Q,LC,i,j,k );
		// Find the local tensor order parameter
		tensOrderParam( &CL[i][j][k],Q,LC );
		//Output
		fprintf( fout,"%5d\t%5d\t%5d\t",i,j,k );
		if( CL[i][j][k].POPSRD == 0 ) fprintf( fout, "%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\n" ,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 );
		else fprintf( fout, "%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\n", Q[0][0],Q[0][1],Q[0][2],Q[1][0],Q[1][1],Q[1][2],Q[2][0],Q[2][1],Q[2][2] );
	}

	for( i=0; i<_3D; i++ ) free( Q[i] );
	free( Q );
	#ifdef FFLSH
		fflush(fout);
	#endif
}

///
/// @brief Outputs Fourier transform of tensor order parameter in space to data file.
///
/// This function calculates the real and imaginary parts of the tensor order parameter as a function of the Fourier transformed spatial wave number.
/// The wave number and squared modulus of the tensor order paramter are printed to file.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param t This is the time step.
/// @param pMPC This a pointer to MPCD particles.
/// @param CL This is a pointer to the co-ordinates and cell of each particle.
/// @param LC This is a flag that states if the system is a liquid crystal.
/// @see findRotationMatrix()
/// @see dotprodMatVec()
/// @see dotprodMatMat()
/// @see outputResults()
///
void orderQKout( FILE *fout,double t,particleMPC pMPC[],cell ***CL,int LC ) {
	int i,j,k,a,b,n;
	double S,ReQ[_3D][_3D],ImQ[_3D][_3D],temp_ReQ[_3D][_3D],temp_ImQ[_3D][_3D],modQ2[_3D][_3D];
	double rotMat[_3D][_3D],rotMatTranspose[_3D][_3D];
	double invL[_3D],DIR[_3D],xhat[_3D],zhat[_3D];
	double waveNum[_3D],K123[_3D],Kprime[_3D],k3;
	double U[_3D],pos[_3D],kr,ckr,skr;
	double c1=1./((double)DIM-1.);
	double c2=VOL/((double)GPOP);
	double fDIM=(double)DIM;
	c2*=c2;

	for( i=0; i<_3D; i++ ) {
		invL[i]=2.*pi/((double)XYZ[i]);
		U[i]=0.;
		pos[i]=0.;
		DIR[i]=0.;
		xhat[i]=0.;
		zhat[i]=0.;
	}
	xhat[0]=1.;
	zhat[2]=1.;
	for( a=0; a<_3D; a++ ) for( b=0; b<_3D; b++ ) {
		ReQ[a][b] = 0.0;
		ImQ[a][b] = 0.0;
		modQ2[a][b] = 0.0;
	}

	for( i=0; i<XYZ[0]; i++ ) {
		waveNum[0] = invL[0]*i;
		for( j=0; j<XYZ[1]; j++ ) {
			waveNum[1] = invL[1]*j;
			for( k=0; k<XYZ[2]; k++ ) {
				waveNum[2] = invL[2]*k;
				//Zero
				for( a=0; a<DIM; a++ ) for( b=0; b<DIM; b++ ) {
					ReQ[a][b] = 0.0;
					ImQ[a][b] = 0.0;
				}
				for( a=0; a<DIM; a++ ) DIR[a] = CL[i][j][k].DIR[a];
				//Calculate order paramter
				for( n=0; n<GPOP; n++ ) {
					for( a=0; a<DIM; a++ ) {
						U[a] = (pMPC+n)->U[a];
						pos[a] = (pMPC+n)->Q[a];
					}
					kr=dotprod( waveNum, pos,DIM );
					ckr=cos(kr);
					skr=sin(kr);
					for( a=0; a<DIM; a++ ) for( b=0; b<DIM; b++ ) {
						S = c1 * fDIM * U[a] * U[b];
						if( a==b ) S-=c1;
						ReQ[a][b] += S*ckr;
						ImQ[a][b] += S*skr;
					}
				}
				//Do the rotation to 13-frame
				//Rotate the director to be parallel to the z-direction
				findRotationMatrix( rotMat,DIR,zhat );
				dotprodMatVec( rotMat,waveNum,Kprime,DIM );
				//Find the projection of the new K onto the xy-plane
				k3=Kprime[2];
				Kprime[2]=0.;
				//Find the Q tensor in this frame
				for( a=0; a<_3D; a++ ) for( b=0; b<_3D; b++ ) rotMatTranspose[a][b] = rotMat[b][a];
				dotprodMatMat( rotMat,ReQ,temp_ReQ,DIM );
				dotprodMatMat( temp_ReQ,rotMatTranspose,ReQ,DIM );
				dotprodMatMat( rotMat,ImQ,temp_ImQ,DIM );
				dotprodMatMat( temp_ImQ,rotMatTranspose,ImQ,DIM );
				//Rotate away the y component
				//Find the rotation matrix that rotates kappa onto x-hat
				findRotationMatrix( rotMat,Kprime,xhat );
				dotprodMatVec( rotMat,Kprime,K123,DIM );
				K123[2]=k3;
				//Find the Q tensor in this frame
				for( a=0; a<_3D; a++ ) for( b=0; b<_3D; b++ ) rotMatTranspose[a][b] = rotMat[b][a];
				dotprodMatMat( rotMat,ReQ,temp_ReQ,DIM );
				dotprodMatMat( temp_ReQ,rotMatTranspose,ReQ,DIM );
				dotprodMatMat( rotMat,ImQ,temp_ImQ,DIM );
				dotprodMatMat( temp_ImQ,rotMatTranspose,ImQ,DIM );
				//Find the modulus squared
				for( a=0; a<DIM; a++ ) for( b=0; b<DIM; b++ ) modQ2[a][b] = c2*( ReQ[a][b]*ReQ[a][b] + ImQ[a][b]*ImQ[a][b] );
				//Output
				fprintf( fout,"%12.5e\t%12.5e\t%12.5e\t",K123[0],K123[1],K123[2] );
				fprintf( fout, "%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\n", modQ2[0][0],modQ2[0][1],modQ2[0][2],modQ2[1][0],modQ2[1][1],modQ2[1][2],modQ2[2][0],modQ2[2][1],modQ2[2][2] );
			}
		}
	}
	#ifdef FFLSH
		fflush(fout);
	#endif
}

///
/// @brief Outputs autocorrelation function data to file.
///
/// This function simply prints the autocorrelation function data to a data output file.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param corr This is the autocorrelation data.
/// @param t This is the time step.
/// @see outputResults()
///
void corrout( FILE *fout,double corr[],double t ) {
	int i;

	//Output
	fprintf( fout,"%12.5e\n",t );
	for( i=0; i<maxXYZ; i++ ) {
		// fprintf( fout,"%12.5e\t%12.5e\n",(double)i,corr[i] );
		fprintf( fout,"\t%d\t%12.5e\n",i,corr[i] );
	}
	fprintf( fout,"\n" );
	#ifdef FFLSH
		fflush(fout);
	#endif
}

///
/// @brief Outputs spectra data to file.
///
/// This function simply prints the energy or enstrophy spectrum data to a data output file.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param spect This is the spectrum data.
/// @param t This is the time step.
/// @see outputResults()
///
void spectout( FILE *fout,double spect[],double t ) {
	int i;
	double k,pi2;

	pi2=2.0*pi;
	//Output
	fprintf( fout,"%12.5e\n",t );
	// for( i=1; i<maxXYZ; i++ ) {
	// 	k=pi2/((double)i);
	// 	fprintf( fout,"\t%12.5e\t%12.5e\n", k,spect[i] );
	// }
	for( i=maxXYZ-1; i>0; i-- ) {
		k=pi2/((double)i);
		fprintf( fout,"\t%12.5e\t%12.5e\n", k,spect[i] );
	}
	fprintf( fout,"\n" );
	#ifdef FFLSH
		fflush(fout);
	#endif
}

///
/// @brief A checkpointing function to allow MPCD to re-initialise a simulation from these parameters.
///
/// This function outputs entire simulation data so that it can be used as a set of re-initialisation parameter for another simulation.
/// Total time and input parameter, as well as the all the output parameters are output.
/// Information on swimmers, MPCD particle species, MD particles and boundaries are also output.
///
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param in This is the list of inputs from input.json.
/// @param SP This is the species-wide information about MPCD particles.
/// @param pSRD This is a list of information for all MPCD particles.
/// @param MD_mode This is a flag to determine if MD mode is on.
/// @param WALL This is a pointer to boundary position information.
/// @param outFlag This is a flag for .dat files to be output.
/// @param runtime This is the length of time the simulation runs.
/// @param warmtime This is the length of warm up time of the simulation.
/// @param AVVEL This is a pointer to the average speed.
/// @param AVS This is is a pointer to the average scalar order parameter.
/// @param avDIR This is the average director orientation in three-dimensions.
/// @param S4 This is a pointer to the fourth moment of the scalar order parameter.
/// @param stdN This is the standard deviation of the density.
/// @param KBTNOW This is a pointer to the current un-thermostated temperature.
/// @param AVV This is a pointer to the past average flow velocities.
/// @param AVNOW This is a pointer to the current average flow velocities.
/// @param theorySP These are theoretical values for each species based off input.json.
/// @param theoryGl These are the global theoretical values based off input.json.
/// @param specS This is the swimmer species.
/// @param sw This is a pointer to the list of swimmers.
///
void checkpoint(FILE *fout, inputList in, spec *SP, particleMPC *pSRD, int MD_mode, bc *WALL, outputFlagsList outFlag, int runtime, int warmtime, double AVVEL, double AVS, double avDIR[_3D], double S4, double stdN, double KBTNOW, double AVV[_3D], double AVNOW[_3D], kinTheory theorySP[], kinTheory theoryGl, specSwimmer specS, swimmer *sw ) {
	int i,j;

	fprintf( fout,"%d\n",in.simSteps );					//total time (or number of iterations)
	fprintf( fout,"%d %lf\n",in.warmupSteps,in.dt );	//Warmup time iterations, and time step

	fprintf( fout,"%ld\n",in.seed );					//Random seed (0 if read from time)
	fprintf( fout,"%d %d %d %d %lf %lf\n",DIM,XYZ[0],XYZ[1],XYZ[2],in.KBT,KBTNOW );
	fprintf( fout,"%d %d %d %d %d %d\n",in.RFRAME,in.zeroNetMom,in.GALINV,in.TSTECH,in.RTECH,in.LC );
	fprintf( fout,"%lf %lf %lf\n",in.TAU,in.RA,in.FRICCO );
	fprintf( fout,"%d %d %d %d\n",in.noHI,in.inCOMP,in.MULTIPHASE,in.MFPLAYERH );
	fprintf( fout,"%lf %lf %lf\n",in.GRAV[0],in.GRAV[1],in.GRAV[2] );	//Acceleration (external force)
	fprintf( fout,"%lf %lf %lf\n",in.MAG[0],in.MAG[1],in.MAG[2] );		//External magnetic field
	fprintf(fout, "%d %d\n", MD_mode, in.stepsMD );						//MD coupling mode and number of MD steps per SRD step
	fprintf( fout,"%d %d\n",GPOP,NSPECI);								//Total number of particles and number of species

	fprintf( fout,"%d %d %lf %lf %d %d %lf\n",runtime,warmtime,in.C,in.S,in.GRAV_FLAG,in.MAG_FLAG,in.tolD );
	fprintf( fout,"%lf %lf %lf %lf %lf %lf %lf %lf\n", AVVEL, AVS, avDIR[0], avDIR[1], avDIR[2], S4, stdN, VOL );
	fprintf( fout,"%lf %lf %lf %lf %lf %lf\n",AVV[0], AVV[1], AVV[2], AVNOW[0], AVNOW[1], AVNOW[2] );
	fprintf( fout,"%lf %lf %lf %lf %lf %lf\n", theoryGl.MFP, theoryGl.VISC, theoryGl.THERMD, theoryGl.SDIFF, theoryGl.SPEEDOFSOUND, theoryGl.sumM );
	fprintf( fout,"%d %d\n",in.chckpntIn,MDmode );		//MD and Checkpointing strings don't currently get read
	// fprintf( fout,"%d %s\n",in.chckpntIn,in.chckpntInputFile );		//Checkpointing
	// fprintf( fout,"%d %s\n",MDmode,mdInputFile );					//MD

	//Output variables
	fprintf( fout,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",DBUG,outFlag.TRAJOUT,outFlag.printSP,outFlag.COAROUT,outFlag.FLOWOUT,outFlag.VELOUT,outFlag.AVVELOUT,outFlag.ORDEROUT,outFlag.QTENSOUT,outFlag.QKOUT,outFlag.AVSOUT,outFlag.SOLOUT,outFlag.ENOUT,outFlag.ENFIELDOUT,outFlag.ENNEIGHBOURS,outFlag.ENSTROPHYOUT,outFlag.DENSOUT,outFlag.AVCOMOUT,outFlag.CVVOUT,outFlag.CNNOUT,outFlag.CWWOUT,outFlag.CDDOUT,outFlag.CSSOUT,outFlag.CPPOUT,outFlag.BINDER,outFlag.BINDERBIN,outFlag.SYNOUT,outFlag.CHCKPNT,outFlag.CHCKPNTrcvr );
	fprintf( fout,"%d %d\n",outFlag.SPOUT,outFlag.PRESOUT );
	fprintf( fout,"%d %d %d %d %d %d %d\n",outFlag.HISTVELOUT,outFlag.HISTSPEEDOUT,outFlag.HISTVORTOUT,outFlag.HISTENSTROUT,outFlag.HISTDIROUT,outFlag.HISTSOUT,outFlag.HISTNOUT );
	fprintf( fout,"%d %d %d %d %d\n",outFlag.ENERGYSPECTOUT,outFlag.ENSTROPHYSPECTOUT,outFlag.TOPOOUT,outFlag.DEFECTOUT,outFlag.DISCLINOUT );
	fprintf( fout,"%d %d %d\n",outFlag.SWOUT,outFlag.SWORIOUT,outFlag.RTOUT );

	//Species of MPCD particles
	for( i=0; i<NSPECI; i++ ) {
		fprintf( fout,"%lf %i %i %i %i ",(SP+i)->MASS,(SP+i)->POP,(SP+i)->QDIST,(SP+i)->VDIST,(SP+i)->ODIST );
		fprintf( fout,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ",(SP+i)->RFC, (SP+i)->LEN, (SP+i)->TUMBLE, (SP+i)->CHIHI, (SP+i)->CHIA, (SP+i)->ACT, (SP+i)->MFPOT, (SP+i)->SIGWIDTH, (SP+i)->SIGPOS, (SP+i)->DAMP );
		fprintf( fout,"%lf %lf %lf %lf %lf\n",(SP+i)->VOL,(SP+i)->nDNST,(SP+i)->mDNST,(SP+i)->MINACTRATIO,(SP+i)->BS );
		for( j=0; j<NSPECI; j++ ) fprintf( fout,"%lf ",(SP+i)->M[j] );			//Binary fluid control parameters
		fprintf( fout,"\n" );
		fprintf( fout,"%lf %lf %lf %lf %lf %lf\n", theorySP[i].MFP, theorySP[i].VISC, theorySP[i].THERMD, theorySP[i].SDIFF, theorySP[i].SPEEDOFSOUND, theorySP[i].sumM );
	}
	fprintf( fout,"%lf %lf %lf %d\n",GnDNST,GmDNST,GMASS,maxXYZ );

	//BCs
	fprintf( fout,"%d\n",NBC );
	for( i=0; i<NBC; i++ ) {
		fprintf( fout,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",WALL[i].COLL_TYPE, WALL[i].PHANTOM, WALL[i].E, WALL[i].Q[0], WALL[i].Q[1], WALL[i].Q[2], WALL[i].V[0], WALL[i].V[1], WALL[i].V[2], WALL[i].O[0], WALL[i].O[1], WALL[i].O[2] );
		fprintf( fout,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", WALL[i].L[0], WALL[i].L[1], WALL[i].L[2], WALL[i].G[0], WALL[i].G[1], WALL[i].G[2], WALL[i].A[0], WALL[i].A[1], WALL[i].A[2], WALL[i].AINV[0], WALL[i].AINV[1], WALL[i].AINV[2],WALL[i].P[0],WALL[i].P[1],WALL[i].P[2],WALL[i].P[3], WALL[i].R, WALL[i].B[0],WALL[i].B[1],WALL[i].B[2] );
		fprintf( fout,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", WALL[i].DN, WALL[i].DT, WALL[i].DVN, WALL[i].DVT, WALL[i].DVxyz[0], WALL[i].DVxyz[1], WALL[i].DVxyz[2], WALL[i].MVN, WALL[i].MVT, WALL[i].MUN, WALL[i].MUT, WALL[i].MUxyz[0], WALL[i].MUxyz[1], WALL[i].MUxyz[2] );
		fprintf( fout,"%lf %lf %lf %lf %d %d %lf\n", WALL[i].DUxyz[0], WALL[i].DUxyz[1], WALL[i].DUxyz[2], WALL[i].KBT, WALL[i].DSPLC, WALL[i].INV, WALL[i].MASS );
		fprintf( fout,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", WALL[i].W, WALL[i].VOL, WALL[i].Q_old[0], WALL[i].Q_old[1], WALL[i].Q_old[2], WALL[i].O_old[0], WALL[i].O_old[1], WALL[i].O_old[2], WALL[i].I[0][0], WALL[i].I[0][1], WALL[i].I[0][2], WALL[i].I[1][0], WALL[i].I[1][1], WALL[i].I[1][2], WALL[i].I[2][0], WALL[i].I[2][1], WALL[i].I[2][2] );
		fprintf( fout,"%d %d %d %lf %lf\n", WALL[i].PLANAR, WALL[i].REORIENT, WALL[i].ABS, WALL[i].ROTSYMM[0], WALL[i].ROTSYMM[1] );
		fprintf( fout,"%lf %lf %lf %lf %lf %lf\n", WALL[i].dV[0], WALL[i].dV[1], WALL[i].dV[2], WALL[i].dL[0], WALL[i].dL[1], WALL[i].dL[2] );
		for( j=0; j<MAXSPECI+2; j++ ) fprintf( fout,"%i ", WALL[i].INTER[j] );		//BC particle interaction flags
		fprintf( fout,"\n" );
	}

	//MPCD particles
	for( i=0; i<GPOP; i++ ) {
		fprintf( fout,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",pSRD[i].S_flag, pSRD[i].SPID, pSRD[i].q, pSRD[i].Q[0], pSRD[i].Q[1], pSRD[i].Q[2], pSRD[i].V[0], pSRD[i].V[1], pSRD[i].V[2], pSRD[i].U[0], pSRD[i].U[1], pSRD[i].U[2], pSRD[i].T[0], pSRD[i].T[1], pSRD[i].T[2] );
	}

	//Swimmers
	fprintf( fout,"%d %d %d %d %d %d %lf %lf %d %d\n", NS,specS.TYPE, specS.QDIST, specS.ODIST, specS.headM, specS.middM, specS.iheadM, specS.imiddM, specS.HSPid, specS.MSPid );
	fprintf( fout,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf %lf %lf %d %lf\n", specS.FS, specS.TS, specS.DS, specS.sizeShrink, specS.springShrink, specS.fixDist, specS.k, specS.ro, specS.iro, specS.sig, specS.isig, specS.eps, specS.dep, specS.range, specS.depth, specS.runTime, specS.tumbleTime, specS.shrinkTime, specS.MAGMOM );
	for( i=0; i<NS; i++ ) {
		fprintf( fout,"%d %lf %lf %lf %d %d %lf %lf %lf %lf %lf\n",(sw+i)->RT,(sw+i)->n0[0],(sw+i)->n0[1],(sw+i)->n0[2],(sw+i)->timeCNT,(sw+i)->timeRND,(sw+i)->ro,(sw+i)->iro,(sw+i)->sig,(sw+i)->isig,(sw+i)->k );
		fprintf( fout,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",(sw+i)->H.HorM,(sw+i)->H.Q[0],(sw+i)->H.Q[1],(sw+i)->H.Q[2],(sw+i)->H.V[0],(sw+i)->H.V[1],(sw+i)->H.V[2],(sw+i)->H.A[0],(sw+i)->H.A[1],(sw+i)->H.A[2] );
		fprintf( fout,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",(sw+i)->M.HorM,(sw+i)->M.Q[0],(sw+i)->M.Q[1],(sw+i)->M.Q[2],(sw+i)->M.V[0],(sw+i)->M.V[1],(sw+i)->M.V[2],(sw+i)->M.A[0],(sw+i)->M.A[1],(sw+i)->M.A[2] );
	}

	fflush(fout); // force flush
}

///
/// @brief This function runs a checkpoint operation.
///
/// This function runs a checkpointing operation that clears up code in `mpcd.c` during checkpointing.
/// Checkpoints may either be time-based or not and is based off `checkpoint()`.
///
/// @param op This is the path to the output file.
/// @param lastCheckpoint This is the last checkpoint time.
/// @param fout This is a pointer to the output .dat file name to be produced.
/// @param in This is the list of inputs from input.json.
/// @param SP This is the species-wide information about MPCD particles.
/// @param pSRD This is a list of information for all MPCD particles.
/// @param MD_mode This is a flag to determine if MD mode is on.
/// @param WALL This is a pointer to boundary position information.
/// @param outFlag This is a flag for .dat files to be output.
/// @param runtime This is the length of time the simulation runs.
/// @param warmtime This is the length of warm up time of the simulation.
/// @param AVVEL This is a pointer to the average speed.
/// @param AVS This is is a pointer to the average scalar order parameter.
/// @param avDIR This is the average director orientation in three-dimensions.
/// @param S4 This is a pointer to the fourth moment of the scalar order parameter.
/// @param stdN This is the standard deviation of the density.
/// @param CoM These are the coordinates of the centres of masses of each species.
/// @param KBTNOW This is a pointer to the current un-thermostated temperature.
/// @param AVV This is a pointer to the past average flow velocities.
/// @param AVNOW This is a pointer to the current average flow velocities.
/// @param theory These are theoretical values based off input.json.
/// @param specS This is the swimmer species.
/// @param sw This is a pointer to the list of swimmers.
/// @see checkpoint()
/// @see openCheckpoint()
///
void runCheckpoint(char op[STRLN], time_t *lastCheckpoint, FILE *fout, inputList in, spec *SP, particleMPC *pSRD, int MD_mode, bc *WALL, outputFlagsList outFlag, int runtime, int warmtime, double AVVEL, double AVS, double avDIR[_3D], double S4, double stdN, double KBTNOW, double AVV[_3D], double AVNOW[_3D], kinTheory theorySP[], kinTheory theoryGl, specSwimmer specS, swimmer *sw ) {
    // if time-based checkpointing has been enabled, see if a checkpoint needs to be made
    // otherwise return early
    if (outFlag.CHCKPNTTIMER != 0.0) {
        time_t currTime = time(NULL);
        if (currTime - *lastCheckpoint >= outFlag.CHCKPNTTIMER*60*60) {
            // if time diff is more than the set checkpointing time
            #ifdef DBG
            	if( DBUG >= DBGRUN ) printf( "\nTimer based checkpoint triggered." );
            #endif
            *lastCheckpoint = currTime;
        } else {
            return; // early return, no checkpoint needed
        }
    }
    #ifdef DBG
    if( DBUG >= DBGRUN ) printf( "\nCheckpointing.\n" );
    #endif
    // normal checkpoint
    openCheckpoint( &(fout),op );
    checkpoint(fout, in, SP, pSRD, MD_mode, WALL, outFlag, runtime, warmtime, AVVEL, AVS, avDIR, S4, stdN, KBTNOW, AVV, AVNOW, theorySP, theoryGl, specS, sw);
    fclose( fout );
}

///
/// @brief This function outputs all results except histograms.
///
/// This function calls relevent functions to calculate values and output the results into data files. This includes all correlation functions and averaged data requested in the input.
///
/// @param CL This is a pointer to the co-ordinates and cell of each particle.
/// @param SRDparticles This is a pointer to the array of particles.
/// @param SP This is the species-wide information about MPCD particles.
/// @param WALL This is a pointer to boundary position information.
/// @param simMD This is a pointer to the MD simulation.
/// @param SS This is a pointer to the swimmer species.
/// @param swimmers This is the list of swimmers.
/// @param AVNOW This is a pointer to the current average flow velocities.
/// @param AVV This is a pointer to the past average flow velocities.
/// @param avDIR This is the average director orientation in three-dimensions.
/// @param runtime This is the length of time the simulation runs.
/// @param in This is the list of inputs from input.json.
/// @param AVVEL This is a pointer to the average speed.
/// @param KBTNOW This is a pointer to the current un-thermostated temperature.
/// @param AVS This is is a pointer to the average scalar order parameter.
/// @param S4 This is a pointer to the fourth moment of the scalar order parameter.
/// @param stdN This is the standard deviation of the density.
/// @param MD_mode This is a flag to determine if MD mode is on.
/// @param outFlag This is a flag for .dat files to be output.
/// @param outFiles This is the list of output files.
/// @see solidout()
/// @see bin()
/// @see binSwimmers()
/// @see binMD()
/// @see localPROP()
/// @see avVel()
/// @see localVelGrad()
/// @see galileantrans()
/// @see zeroExtraDims()
/// @see avOrderParam()
/// @see avS4()
/// @see avsout()
/// @see densSTDout()
/// @see binderCumulant()
/// @see binderout()
/// @see avveloutWithGradVel()
/// @see avEnstrophy()
/// @see avenstrophyout()
/// @see enout()
/// @see enfieldout()
/// @see enneighboursout()
/// @see swimout()
/// @see swimoriout()
/// @see corrout()
/// @see velvelCorr()
/// @see dirdirCorr()
/// @see orderorderCorr()
/// @see normCorr()
/// @see FTspectrum()
/// @see spectout()
/// @see coordout()
/// @see flowout()
/// @see coarseout()
/// @see orderout()
/// @see orderQout()
/// @see disclinationTensorOut()
/// @see multiphaseout()
/// @see pressureout()
/// @see orderQKout()
///
void outputResults( cell ***CL,particleMPC *SRDparticles,spec SP[],bc WALL[],simptr simMD,specSwimmer SS, swimmer swimmers[],double AVNOW[_3D],double AVV[_3D],double avDIR[_3D], int runtime, inputList in, double AVVEL, double KBTNOW,double *AVS,double *S4,double *stdN,double *CoM[MAXSPECI][_3D],int MD_mode,outputFlagsList outFlag,outputFilesList outFiles ) {
	int a,b,c,i,j,k;
	double time_now = runtime*in.dt;					//Simulation time
	double wmf;
	double corr[maxXYZ],spect[maxXYZ];				//Correlation functions and energy spectra
	double UL;																//Binder cumulant
	double avGradVel[_3D][_3D];	
	double com[MAXSPECI][_3D];					//Center of mass vector for each species
	double sigcom[MAXSPECI][_3D][_3D];	        // Center of mass std tensor for each species							//Velocity gradient
	double AVORI[_3D];
	/* ****************************************** */
	/* ************** BC trajectory ************* */
	/* ****************************************** */
	for( i=0; i<NBC; i++ ) if( (WALL+i)->DSPLC ) {
		// Write values
		if( outFlag.SOLOUT>=OUT && runtime%outFlag.SOLOUT==0 ) solidout(outFiles.fsolids[i],WALL[i],time_now);
	}
	/* ****************************************** */
	/* ************** BIN and CALC ************** */
	/* ****************************************** */
	#ifdef DBG
		if( DBUG >= DBGTITLE ) printf( "Bin Particles.\n" );
	#endif
	// Bin SRD particles
	bin( CL,SP,WALL,in.KBT,in.LC,0 );
	// Bin swimmer monomers
	binSwimmers( CL,0 );
	// Bin MD particles
	if( MD_mode ) binMD(CL );
	//Calculate the local properties of each cell (VCM,in.KBT,POPulation,Mass)
	localPROP( CL,SP,SS,in.RTECH,in.LC );
	avVel( CL,AVNOW );
	avOri( SRDparticles,AVORI );
	//Calculate velocity gradient
	if( (outFlag.AVVELOUT>=OUT && runtime%outFlag.AVVELOUT==0) || (outFlag.ENSTROPHYOUT>=OUT && runtime%outFlag.ENSTROPHYOUT==0) || (outFlag.HISTVORTOUT>=OUT && runtime%outFlag.HISTVORTOUT==0)  || (outFlag.HISTENSTROUT>=OUT && runtime%outFlag.HISTENSTROUT==0) || (outFlag.CWWOUT>=OUT && runtime%outFlag.CWWOUT==0) ) {
		//Velocity gradient
		localVelGrad( CL );
	}
	/* ****************************************** */
	/* ************ ZERO NET MOMENTUM *********** */
	/* ****************************************** */
	if( in.RFRAME && runtime%in.zeroNetMom==0 ) {
		#ifdef DBG
			if( DBUG > DBGRUN ) printf( "Galilean Transformation to Rest Frame\n" );
		#endif
		galileantrans(SRDparticles, WALL, simMD, SP, in.KBT, AVV, GPOP, NBC, MD_mode, DIM );
		zeroExtraDims(SRDparticles, WALL, simMD, GPOP, NBC, MD_mode, DIM );
	}
	/* ****************************************** */
	/* *********** AVERAGES and OUTPUT ********** */
	/* ****************************************** */
	//Calculate the average scalar order parameter
	if( outFlag.AVSOUT>=OUT && runtime%outFlag.AVSOUT==0 ) {
		*AVS = avOrderParam( SRDparticles,in.LC,avDIR );
		*S4 = avS4( SRDparticles,in.LC,avDIR );
		avsout( outFiles.favs,time_now,*AVS,*S4,avDIR );
	}
	//Calculate density variation
	if( outFlag.DENSOUT>=OUT && runtime%outFlag.DENSOUT==0 ) {
		*stdN = stdNum( CL,GPOP,XYZ,XYZ_P1 );
		densSTDout( outFiles.fdensSTD,time_now,*stdN );
	}
	//Calculate centers of masses of multiphase fluids
	if( outFlag.AVCOMOUT>=OUT && runtime%outFlag.AVCOMOUT==0 ) {
		for (i=0;i<NSPECI;i++) for (j=0;j<_3D;j++) {
		com[i][j]=0.0;
		for (k=0;k<_3D;k++) sigcom[i][j][k]=0.0;
		}
		for (i=0;i<GPOP;i++){
			a=SRDparticles[i].SPID;
			for (j=0;j<DIM;j++) com[a][j]+=SRDparticles[i].Q[j];
		}
		for (i=0;i<NSPECI;i++) if(SP[i].POP>0) for (j=0;j<_3D;j++) com[i][j] /= (float) SP[i].POP;
		for (i=0;i<GPOP;i++){
			a=SRDparticles[i].SPID;
			for (j=0;j<DIM;j++) for (k=0;k<DIM;k++) sigcom[a][j][k]+=(com[a][j]-SRDparticles[i].Q[j])*(com[a][k]-SRDparticles[i].Q[k]);
		}
		for (i=0;i<NSPECI;i++) if(SP[i].POP>0) for (j=0;j<_3D;j++) for (k=0;k<_3D;k++) sigcom[i][j][k] /= (float) SP[i].POP;
		avcomout( outFiles.fmpcom,time_now,com,sigcom);
	}
	//Calculate binder cumulants
	if( outFlag.BINDER>=OUT && runtime%outFlag.BINDER==0 ) {
		UL=binderCumulant( CL,outFlag.BINDERBIN,in.LC );
		binderout( outFiles.fbinder,time_now,UL );
	}

	//Calculate average velocity and enstrophy
	if( (outFlag.AVVELOUT>=OUT && runtime%outFlag.AVVELOUT==0) || (outFlag.ENSTROPHYOUT>=OUT && runtime%outFlag.ENSTROPHYOUT==0) ) {
		if( outFlag.AVVELOUT>=OUT && runtime%outFlag.AVVELOUT==0 ) {
			for( i=0; i<_3D; i++ ) for( j=0; j<_3D; j++ ) avGradVel[i][j]=0.;
			for( a=0; a<XYZ[0]; a++ ) for( b=0; b<XYZ[1]; b++ ) for( c=0; c<XYZ[2]; c++ ) {
				for( i=0; i<_3D; i++ ) for( j=0; j<_3D; j++ ) avGradVel[i][j] += CL[a][b][c].E[i][j];
			}
			for( i=0; i<_3D; i++ ) for( j=0; j<_3D; j++ ) avGradVel[i][j] /= VOL;
			avveloutWithGradVel( outFiles.favvel,time_now,AVNOW,KBTNOW,avGradVel );
		}
		//Enstrophy
		if( outFlag.ENSTROPHYOUT>=OUT && runtime%outFlag.ENSTROPHYOUT==0 ) {
			wmf = avEnstrophy( CL );
			avenstrophyout( outFiles.fenstrophy,time_now,wmf );
		}
	}
	/* ****************************************** */
	/* *********** AVERAGE ORIENTATION ********** */
	/* ****************************************** */
	if( outFlag.AVORIOUT>=OUT && runtime%outFlag.AVORIOUT==0 ) {
		avoriout(outFiles.favori, time_now, AVORI);
	}
	/* ****************************************** */
	/* *************** TOTAL ENERGY ************* */
	/* ****************************************** */
	if( outFlag.ENOUT>=OUT && runtime%outFlag.ENOUT==0 ) {
		wmf = calcE_LC( CL,in.LC,SP );
		enout( outFiles.fenergy,SRDparticles,SP,WALL,time_now,KBTNOW,wmf );
	}
	if( outFlag.ENFIELDOUT>=OUT && runtime%outFlag.ENFIELDOUT==0 ) enfieldout( outFiles.fenergyfield,CL,SP,in.LC );
	if( outFlag.ENNEIGHBOURS>=OUT && runtime%outFlag.ENNEIGHBOURS==0 ) enneighboursout( outFiles.fenneighbours,time_now,CL,SP,in.LC );
	/* ****************************************** */
	/* ***** SWIMMERS' POSITONS/ORIENTATIONS **** */
	/* ****************************************** */
	if( outFlag.SWOUT>=OUT && runtime%outFlag.SWOUT==0 ) swimout( outFiles.fswimmers,swimmers,time_now );
	if( outFlag.SWORIOUT>=OUT && runtime%outFlag.SWORIOUT==0 ) swimoriout( outFiles.fswimmersOri,swimmers,time_now );
	/* ****************************************** */
	/* ********** SPATIAL CORRELATIONS ********** */
	/* ****************************************** */
	if( (outFlag.CVVOUT>=OUT && runtime%outFlag.CVVOUT==0) || (outFlag.CNNOUT>=OUT && runtime%outFlag.CNNOUT==0) || (outFlag.CWWOUT>=OUT && runtime%outFlag.CWWOUT==0) ) {
		#ifdef DBG
				if( DBUG >= DBGTITLE ) printf( "Calcualte spatial correlation functions.\n" );
		#endif
	}
	if( outFlag.CNNOUT>=OUT && runtime%outFlag.CNNOUT==0 ) {
		dirdirCorr( CL,maxXYZ,XYZ,corr,DIM );
		corrout( outFiles.fcorrNN,corr,time_now );
	}
	if( outFlag.CDDOUT>=OUT && runtime%outFlag.CDDOUT==0 ) {
		densdensCorr( CL,maxXYZ,XYZ,corr,DIM );
		corrout( outFiles.fcorrDD,corr,time_now );
	}
	if( outFlag.CSSOUT>=OUT && runtime%outFlag.CSSOUT==0 ) {
		if( in.LC ) {
			orderorderCorr( CL,maxXYZ,XYZ,corr,DIM );
			corrout( outFiles.fcorrSS,corr,time_now );
		}
	}
	if( outFlag.CPPOUT>=OUT && runtime%outFlag.CPPOUT==0 ) {
		if( in.MULTIPHASE!=MPHOFF ) {
			// phiphiCorr( CL,maxXYZ,XYZ,corr,DIM );
			printf("Warning: phiphiCorr() needs to be checked.\n");
			corrout( outFiles.fcorrPP,corr,time_now );
		}
	}
	//Velocity correlations and energy spectrum
	if( (outFlag.CVVOUT>=OUT && runtime%outFlag.CVVOUT==0) || (outFlag.ENERGYSPECTOUT>=OUT && runtime%outFlag.ENERGYSPECTOUT==0) ) {
		//Calculate the un-normalized correlation function
		velvelCorr( CL,maxXYZ,XYZ,corr,DIM );
		if( outFlag.ENERGYSPECTOUT>=OUT && runtime%outFlag.ENERGYSPECTOUT==0 ) {
			//FT into energy spectrum
			FTspectrum( corr,spect,maxXYZ,DIM );
			//Output spectrum
			spectout( outFiles.fenergyspect,spect,time_now );
		}
		if( outFlag.CVVOUT>=OUT && runtime%outFlag.CVVOUT==0 ) {
			//Normalize correlation function
			normCorr( corr,maxXYZ );
			//Output correlation function
			corrout( outFiles.fcorrVV,corr,time_now );
		}
	}
	//Vorticity correlations and energy spectrum
	if( (outFlag.CWWOUT>=OUT && runtime%outFlag.CWWOUT==0) || (outFlag.ENSTROPHYSPECTOUT>=OUT && runtime%outFlag.ENSTROPHYSPECTOUT==0) ) {
		//Calculate the un-normalized correlation function
		vortvortCorr( CL,maxXYZ,XYZ,corr,_3D );
		if( outFlag.ENSTROPHYSPECTOUT>=OUT && runtime%outFlag.ENSTROPHYSPECTOUT==0 ) {
			//FT into energy spectrum
			FTspectrum( corr,spect,maxXYZ,DIM );
			//Output spectrum
			spectout( outFiles.fenstrophyspect,spect,time_now );
		}
		if( outFlag.CWWOUT>=OUT && runtime%outFlag.CWWOUT==0 ) {
			//Normalize correlation function
			normCorr( corr,maxXYZ );
			//Output correlation function
			corrout( outFiles.fcorrWW,corr,time_now );
		}
	}
	/* ****************************************** */
	/* ************ WRITE COORDINATES *********** */
	/* ****************************************** */
	#ifdef DBG
		if( DBUG >= DBGTITLE ) printf( "Write Data Out.\n" );
	#endif
	if(outFlag.printSP>0) if( outFlag.TRAJOUT>=OUT  && runtime%outFlag.TRAJOUT==0 ) coordout( outFiles.fdetail,outFlag.printSP,time_now,SRDparticles,SP );
	if( outFlag.FLOWOUT>=OUT && runtime%outFlag.FLOWOUT==0 ) flowout( outFiles.fflow,CL,outFlag.FLOWOUT, time_now);
	if( outFlag.VELOUT>=OUT && runtime%outFlag.VELOUT==0 ) velout( outFiles.fvel, CL, time_now);
	if( outFlag.SWFLOWOUT>=OUT && runtime%outFlag.SWFLOWOUT==0 && runtime!=0) swflowout( outFiles.fswflow,CL,outFlag.SWFLOWOUT, time_now);
	if( outFlag.COAROUT>=OUT && runtime%outFlag.COAROUT==0 ) coarseout( outFiles.fcoarse,time_now,CL );
	if(in.LC!=ISOF) if( outFlag.ORDEROUT>=OUT && runtime%outFlag.ORDEROUT==0 ) orderout( outFiles.forder,time_now,CL,in.LC );
	if(in.LC!=ISOF) if( outFlag.QTENSOUT>=OUT && runtime%outFlag.QTENSOUT==0 ) orderQout( outFiles.forderQ,time_now,CL,in.LC );
	if(in.LC!=ISOF) if( outFlag.DISCLINOUT>=OUT && runtime%outFlag.DISCLINOUT==0 ) disclinationTensorOut( outFiles.fdisclination,time_now,CL,in.LC );
	if( outFlag.SPOUT>=OUT && runtime%outFlag.SPOUT==0 ) multiphaseout( outFiles.fmultiphase,time_now,CL );
	if( outFlag.PRESOUT>=OUT && runtime%outFlag.PRESOUT==0 ) pressureout( outFiles.fpressure,time_now,CL );
	if( outFlag.spPRESOUT>=OUT && runtime%outFlag.spPRESOUT==0 ) sppressureout( outFiles.fsppressure,time_now,CL );
	if(in.LC!=ISOF) if( outFlag.QKOUT>=OUT && runtime%outFlag.QKOUT==0 ) {
		#ifdef DBG
			if( DBUG >= DBGTITLE ) printf( "Calculate Q-tensor in reciprocal space.\n" );
		#endif
		orderQKout( outFiles.forderQK,time_now,SRDparticles,CL,in.LC );
	}
	/* ****************************************** */
	/* ************** TRACK DEFECTS ************* */
	/* ****************************************** */
	if(in.LC!=ISOF && DIM==_2D) if((outFlag.TOPOOUT>=OUT && runtime%outFlag.TOPOOUT==0)||(outFlag.DEFECTOUT>=OUT && runtime%outFlag.DEFECTOUT==0)) topoChargeAndDefectsOut( outFiles.ftopo, outFlag.TOPOOUT, outFiles.fdefects, outFlag.DEFECTOUT, time_now, CL, in.tolD);
}


///
/// @brief This function outputs histogram results.
///
/// This function calls relevent functions to calculate values and output the histogram results into data files.
///
/// @param CL This is a pointer to the co-ordinates and cell of each particle.
/// @param runtime This is the length of time the simulation runs.
/// @param in This is the list of inputs from input.json.
/// @param outFlag This is a flag for .dat files to be output.
/// @param outFiles This is the list of output files.
/// @see histVelout()
/// @see histSpeedout()
/// @see histEnstout()
/// @see histVortout()
/// @see histNout()
///
void outputHist( cell ***CL,int runtime, inputList in,outputFlagsList outFlag,outputFilesList outFiles ) {
	int a,b,c,i,j;
	double time_now = runtime*in.dt;
	double myVec[_3D];													//Velocity (etc) actual values for every MPCD cell
	double maxRange;														//Maximum for range for histograms
	int nc=VOL;
	int hist[_3D][BINS];												//Velocity (etc) histogram for each of the D3 components
	double myValues[_3D][XYZ[0]*XYZ[1]*XYZ[2]];	//Velocity (etc) actual values for every MPCD cell

	/* ****************************************** */
	/* ************ HISTOGRAM BINNNING ********** */
	/* ****************************************** */
	#ifdef DBG
		if( DBUG >= DBGTITLE ) printf( "Collect Distributions.\n" );
	#endif
	//Velocity
	if( outFlag.HISTVELOUT>=OUT && runtime%outFlag.HISTVELOUT==0 ) {
		#ifdef DBG
			if( DBUG >= DBGHIST ) printf("\tBin velocity\n");
		#endif
		//Zero the counter array
		zeroHISTVEC( hist );
		maxRange=0.0;
		j=0;
		//Sort the value array
		for( a=0; a<XYZ[0]; a++ ) for( b=0; b<XYZ[1]; b++ ) for( c=0; c<XYZ[2]; c++ ) {
			for( i=0; i<_3D; i++ ) {
				myValues[i][j]=CL[a][b][c].VCM[i];
				if( fabs(myValues[i][j])>maxRange ) maxRange=fabs(myValues[i][j]);
			}
			j++;
		}
		//Bin
		for( i=0; i<_3D; i++ ) histbin( myValues[i],hist[i],-1.*maxRange,maxRange,nc );
		histVelout( outFiles.fhistVel,hist,-1.*maxRange,maxRange,time_now );
	}
	//Speed
	if( outFlag.HISTSPEEDOUT>=OUT && runtime%outFlag.HISTSPEEDOUT==0 ) {
		#ifdef DBG
			if( DBUG >= DBGHIST ) printf("\tBin speed\n");
		#endif
		//Zero the counter array
		zeroHISTSCALAR( hist[0] );
		maxRange=0.0;
		j=0;
		//Sort the value array
		for( a=0; a<XYZ[0]; a++ ) for( b=0; b<XYZ[1]; b++ ) for( c=0; c<XYZ[2]; c++ ) {
			myValues[0][j]=length( CL[a][b][c].VCM,DIM );
			if( myValues[0][j]>maxRange ) maxRange=myValues[0][j];
			j++;
		}
		//Bin
		histbin( myValues[0],hist[0],0.0,maxRange,nc );
		histSpeedout( outFiles.fhistSpeed,hist[0],0.0,maxRange,time_now );
	}
	//Vorticity
	if( outFlag.HISTVORTOUT>=OUT && runtime%outFlag.HISTVORTOUT==0 ) {
		#ifdef DBG
			if( DBUG >= DBGHIST ) printf("\tBin vorticity\n");
		#endif
		//Zero the counter array
		zeroHISTVEC( hist );
		maxRange=0.0;
		j=0;
		//Sort the value array
		for( a=0; a<XYZ[0]; a++ ) for( b=0; b<XYZ[1]; b++ ) for( c=0; c<XYZ[2]; c++ ) {
			myValues[0][j]=(CL[a][b][c].E[2][1] - CL[a][b][c].E[1][2]);
			myValues[1][j]=(CL[a][b][c].E[0][2] - CL[a][b][c].E[2][0]);
			myValues[2][j]=(CL[a][b][c].E[1][0] - CL[a][b][c].E[0][1]);
			for( i=0; i<_3D; i++ ) if( fabs(myValues[i][j])>maxRange ) maxRange=fabs(myValues[i][j]);
			j++;
		}
		//Bin
		for( i=0; i<_3D; i++ ) histbin( myValues[i],hist[i],-1.*maxRange,maxRange,nc );
		histVortout( outFiles.fhistVort,hist,-1.*maxRange,maxRange,time_now );
	}
	//Enstrophy
	if( outFlag.HISTENSTROUT>=OUT && runtime%outFlag.HISTENSTROUT==0 ) {
		#ifdef DBG
			if( DBUG >= DBGHIST ) printf("\tBin enstrophy\n");
		#endif
		//Zero the counter array
		zeroHISTSCALAR( hist[0] );
		maxRange=0.0;
		j=0;
		//Sort the value array
		for( a=0; a<XYZ[0]; a++ ) for( b=0; b<XYZ[1]; b++ ) for( c=0; c<XYZ[2]; c++ ) {
			myVec[0]=(CL[a][b][c].E[2][1] - CL[a][b][c].E[1][2]);
			myVec[1]=(CL[a][b][c].E[0][2] - CL[a][b][c].E[2][0]);
			myVec[2]=(CL[a][b][c].E[1][0] - CL[a][b][c].E[0][1]);
			myValues[0][j]=0.5*dotprod( myVec,myVec,_3D );
			if( myValues[0][j]>maxRange ) maxRange=myValues[0][j];
			j++;
		}
		//Bin
		histbin( myValues[0],hist[0],0.0,maxRange,nc );
		histEnstrout( outFiles.fhistEnstr,hist[0],0.0,maxRange,time_now );
	}
	//Director
	if( outFlag.HISTDIROUT>=OUT && runtime%outFlag.HISTDIROUT==0 ) {
		#ifdef DBG
			if( DBUG >= DBGHIST ) printf("\tBin director\n");
		#endif
		//Zero the counter array
		zeroHISTVEC( hist );
		j=0;
		//Sort the value array
		for( a=0; a<XYZ[0]; a++ ) for( b=0; b<XYZ[1]; b++ ) for( c=0; c<XYZ[2]; c++ ) {
			for( i=0; i<_3D; i++ ) myValues[i][j]=fabs(CL[a][b][c].DIR[i]);
			j++;
		}
		//Bin
		for( i=0; i<_3D; i++ ) histbin( myValues[i],hist[i],0.0,1.0,nc );
		histVortout( outFiles.fhistDir,hist,0.0,1.0,time_now );
	}
	//Scalar order parameter
	if( outFlag.HISTSOUT>=OUT && runtime%outFlag.HISTSOUT==0 ) {
		#ifdef DBG
			if( DBUG >= DBGHIST ) printf("\tBin scalar order\n");
		#endif
		//Zero the counter array
		zeroHISTSCALAR( hist[0] );
		j=0;
		//Sort the value array
		for( a=0; a<XYZ[0]; a++ ) for( b=0; b<XYZ[1]; b++ ) for( c=0; c<XYZ[2]; c++ ) {
			myValues[0][j]=CL[a][b][c].S;
			j++;
		}
		//Bin
		histbin( myValues[0],hist[0],0.0,1.0,nc );
		histSout( outFiles.fhistS,hist[0],0.0,1.0,time_now );
	}
	//Number per cell
	if( outFlag.HISTNOUT>=OUT && runtime%outFlag.HISTNOUT==0 ) {
		#ifdef DBG
			if( DBUG >= DBGHIST ) printf("\tBin density\n");
		#endif
		//Zero the counter array
		zeroHISTSCALAR( hist[0] );
		maxRange=0.0;
		j=0;
		//Sort the value array
		for( a=0; a<XYZ[0]; a++ ) for( b=0; b<XYZ[1]; b++ ) for( c=0; c<XYZ[2]; c++ ) {
			myValues[0][j]=CL[a][b][c].POP;
			if( myValues[0][j]>maxRange ) maxRange=myValues[0][j];
			j++;
		}
		//Bin
		histbin( myValues[0],hist[0],0.0,maxRange,nc );
		histNout( outFiles.fhistDens,hist[0],0.0,maxRange,time_now );
	}
}

///
/// @brief This function closes output files after writing.
///
/// This function closes output files after writing.
///
/// @param SP This is the species-wide information about MPCD particles.
/// @param WALL This is a pointer to boundary position information.
/// @param outFlag This is a flag for .dat files to be output.
/// @param outFiles This is the list of output files.
///
void closeOutputFiles( spec *SP,bc WALL[],outputFlagsList outFlag,outputFilesList outFiles ) {
	int i;

	if( outFlag.TRAJOUT>=OUT ) for( i=0; i<NSPECI; i++ ) if( SP[i].POP>=1 ) fclose( outFiles.fdetail[i] );
	if( outFlag.COAROUT>=OUT ) fclose( outFiles.fcoarse );
	if( outFlag.AVVELOUT>=OUT ) fclose( outFiles.favvel );
	if( outFlag.AVORIOUT>=OUT ) fclose( outFiles.favori );
	if( outFlag.ORDEROUT>=OUT ) fclose( outFiles.forder );
	if( outFlag.QTENSOUT>=OUT ) fclose( outFiles.forderQ );
	if( outFlag.QKOUT>=OUT ) fclose( outFiles.forderQK );
	if( outFlag.AVSOUT>=OUT ) fclose( outFiles.favs );
	if( outFlag.DENSOUT>=OUT ) fclose( outFiles.fdensSTD );
	if( outFlag.ENSTROPHYOUT>=OUT ) fclose( outFiles.fenstrophy );
	if( outFlag.FLOWOUT>=OUT ) fclose( outFiles.fflow );
	if( outFlag.VELOUT>=OUT ) fclose( outFiles.fvel );
	if( outFlag.SWFLOWOUT>=OUT ) fclose( outFiles.fswflow);
	if( outFlag.ENOUT>=OUT ) fclose( outFiles.fenergy );
	if( outFlag.ENFIELDOUT>=OUT ) fclose( outFiles.fenergyfield );
	if( outFlag.ENNEIGHBOURS>=OUT ) fclose( outFiles.fenneighbours );
	if( outFlag.BINDER>=OUT ) fclose( outFiles.fbinder );
	if( outFlag.CVVOUT>=OUT ) fclose( outFiles.fcorrVV );
	if( outFlag.CNNOUT>=OUT ) fclose( outFiles.fcorrNN );
	if( outFlag.CWWOUT>=OUT ) fclose( outFiles.fcorrWW );
	if( outFlag.CDDOUT>=OUT ) fclose( outFiles.fcorrDD );
	if( outFlag.CSSOUT>=OUT ) fclose( outFiles.fcorrSS );
	if( outFlag.CPPOUT>=OUT ) fclose( outFiles.fcorrPP );
	if( outFlag.ENERGYSPECTOUT>=OUT ) fclose( outFiles.fenergyspect );
	if( outFlag.ENSTROPHYSPECTOUT>=OUT ) fclose( outFiles.fenstrophyspect );
	if( outFlag.HISTVELOUT>=OUT ) fclose( outFiles.fhistVel );
	if( outFlag.HISTSPEEDOUT>=OUT ) fclose( outFiles.fhistSpeed );
	if( outFlag.HISTVORTOUT>=OUT ) fclose( outFiles.fhistVort );
	if( outFlag.HISTENSTROUT>=OUT ) fclose( outFiles.fhistEnstr );
	if( outFlag.HISTDIROUT>=OUT ) fclose( outFiles.fhistDir );
	if( outFlag.HISTSOUT>=OUT ) fclose( outFiles.fhistS );
	if( outFlag.HISTNOUT>=OUT ) fclose( outFiles.fhistDens );
	if( outFlag.TOPOOUT>=OUT ) fclose( outFiles.ftopo );
	if( outFlag.DEFECTOUT>=OUT ) fclose( outFiles.fdefects );
	if( outFlag.DISCLINOUT>=OUT ) fclose( outFiles.fdisclination );
	if( outFlag.SPOUT>=OUT ) fclose( outFiles.fmultiphase );
	if( outFlag.PRESOUT>=OUT ) fclose( outFiles.fpressure );
	if( outFlag.spPRESOUT>=OUT ) fclose( outFiles.fsppressure );
	if( outFlag.SWOUT>=OUT ) fclose( outFiles.fswimmers );
	if( outFlag.SWORIOUT>=OUT ) fclose( outFiles.fswimmersOri );
	if( outFlag.RTOUT>=OUT ) fclose( outFiles.fruntumble );
	if( outFlag.AVCOMOUT>=OUT ) fclose( outFiles.fmpcom );
	if( outFlag.SOLOUT>=OUT ) for( i=0; i<NBC; i++ ) if( WALL[i].DSPLC ) fclose( outFiles.fsolids[i] );
}

///
/// @brief This function writes output files.
///
/// This function writes output files.
///
/// @param t This is the time step.
/// @param f This is a flag for .dat files to be output.
/// @param RFRAME This is a pointer to rest frame.
/// @param zeroNetMom This is momentum correction term to reset to the rest frame.
///
int writeOutput( int t,outputFlagsList f,int RFRAME,int zeroNetMom ) {
	if( ( RFRAME && t%zeroNetMom==0 ) || ( f.ENOUT>=OUT && t%f.ENOUT==0 ) || ( f.TRAJOUT>=OUT  && t%f.TRAJOUT==0 ) || ( f.AVVELOUT>=OUT && t%f.AVVELOUT==0 ) || ( f.AVORIOUT>=OUT && t%f.AVORIOUT==0 ) || ( f.QKOUT && t%f.QKOUT==0 ) || ( f.AVSOUT>=OUT && t%f.AVSOUT==0 ) || ( f.ENNEIGHBOURS>=OUT && t%f.ENNEIGHBOURS==0 ) || ( f.SOLOUT>=OUT && t%f.SOLOUT==0 ) || ( f.BINDER && t%f.BINDER==0 ) || ( f.SWOUT && t%f.SWOUT==0 ) || ( f.SWORIOUT && t%f.SWORIOUT==0 ) || ( f.spPRESOUT && t%f.spPRESOUT==0 ) || ( f.AVCOMOUT && t%f.AVCOMOUT==0 ) ) {
		return 1;
	}
	//Fields
	else if( ( f.FLOWOUT>=OUT && t%f.FLOWOUT==0 ) || ( f.VELOUT>=OUT && t%f.VELOUT==0 ) || ( f.SWFLOWOUT>=OUT && t%f.SWFLOWOUT==0 ) || ( f.COAROUT>=OUT && t%f.COAROUT==0 ) || ( f.ENFIELDOUT>=OUT && t%f.ENFIELDOUT==0 ) || ( f.ORDEROUT && t%f.ORDEROUT==0 ) || ( f.QTENSOUT && t%f.QTENSOUT==0 ) || ( f.TOPOOUT && t%f.TOPOOUT==0 ) || ( f.DEFECTOUT && t%f.DEFECTOUT==0 ) || ( f.DISCLINOUT && t%f.DISCLINOUT==0 ) || ( f.SPOUT && t%f.SPOUT==0 ) || ( f.PRESOUT && t%f.PRESOUT==0 ) ) {
		return 1;
	}
	//Correlation functions
	else if( ( f.CVVOUT && t%f.CVVOUT==0 ) || ( f.CNNOUT && t%f.CNNOUT==0 ) || ( f.CWWOUT && t%f.CWWOUT==0 ) || ( f.CDDOUT && t%f.CDDOUT==0 ) || ( f.CSSOUT && t%f.CSSOUT==0 ) || ( f.CPPOUT && t%f.CPPOUT==0 ) ) {
		return 1;
	}
	else return 0;
}

///
/// @brief This function writes histogram output files.
///
/// This function writes histogram output files
///
/// @param t This is the time step.
/// @param f This is momentum correction term to reset to the rest frame.
///
int writeHistograms( int t,outputFlagsList f ) {
	if( ( f.HISTVELOUT>=OUT && t%f.HISTVELOUT==0 ) || ( f.HISTSPEEDOUT>=OUT && t%f.HISTSPEEDOUT==0 ) || ( f.HISTVORTOUT>=OUT && t%f.HISTVORTOUT==0 ) || ( f.HISTENSTROUT>=OUT && t%f.HISTENSTROUT==0 ) || ( f.HISTDIROUT>=OUT && t%f.HISTDIROUT==0 ) || ( f.HISTSOUT>=OUT && t%f.HISTSOUT==0 ) || ( f.HISTNOUT>=OUT && t%f.HISTNOUT==0 ) ) {
		return 1;
	}
	else return 0;
}
