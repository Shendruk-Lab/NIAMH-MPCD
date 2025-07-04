///
/// @file
/// 
/// @brief This file contains the specific function for a multiphase twist on the NIAMH-MPCD algorithm.
///
/// The file contains two functions that interact with the regular collision operators to make them compatible with multiphase fluids.
///
/// @see andersenMPC()
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
# include "../headers/multiphase.h"

# include "../../md/mdsrd.h"


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
		// MPCD particles
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
	// MPCD particles
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

	if( multiphaseMode==MPHSURF ) {
	
		int i,j,k,d,id;
		double timestep = 0.1; //it would be better to get this value from the input file, as it must correspond to the dt set there. 
		smono *tsm;
		int ncoeff; 
		int orderApprox = (SP+1)->M[1]; 
		//for orderApprox, it may be better to set up new input variable. Currently the second element of "interMatr" for the second species in the input file 
		//can be set to 2 or 3 for parabolic and cubic approximations, respectively, of the order parameter field. 
		int incompressibility = 0; //0 for no incompressibility, 1 for ideal gas, 2 for non-linear 
		// the incompressibility method here needs revising, so the above variable should be set to zero
		//double incompStrength = 0.0001;
		//it is potentially better to have a dedicated input variable for setting the three free energy density pareameters below.
		double tau = (SP)->M[0]; // this is the reduced temperature of the system ( tau = (T-Tc)/Tc ) used in the free energy density "f", and should be negative to achieve phase separation. 
		double b = (SP)->M[1]; // the phi^4 coefficient in f. b>0 for stability.
		double kappa = (SP+1)->M[0]; // the coefficient of (grad phi)^2 in f, associated with surface tension. 
		double tempMag = 0.0;
		double FbulkMultiplier = 100.0; // scales the effect of bulk forces. 
		double FdenMultiplier = 0.0; // set to zero here as the Fden method is not yet complete.
		double FintMultiplier = 1.0; // scales the effect of interfacial forces. 
		double FswimMultiplier = 100000.0; //a very high value means that the chemotactic effect dominates whenever there is a finite order parameter gradient.
		double K = 0.1; // parameterises used in the chemotactic force on swimmers. 
		double Q[DIM]; 
		double Vtemp[DIM], velMag=0.0;
		double Fbulk[DIM], Fchem[DIM], Fint[DIM], Fden[DIM];
		double N=0.0,Na=0.0,Nb=0.0; //counts number of particles in cell, that of type A, and that of type B respectively, 
		particleMPC *tmpc;
		tmpc = CL->pp;
		id = tmpc->SPID;
		if (DIM == 2) // note that parabolic approx in 2D requires a specification of 6 coefficients, while cubic requires 10 coefficients. 
		{
			if (orderApprox==2){
				ncoeff = 6;
			}
			else if (orderApprox==3){
				ncoeff = 10;
			}
		}
		else if (DIM == 3)// note that parabolic approx in 3D requires a specification of 10 coefficients. Cubic for 3D has not yet been implemented. 
		{
			ncoeff = 10; 
			if (orderApprox == 3){
				printf("\nCubic approximation of order parameter field in 3D systems not yet implemented - switching to parabolic approximation.\n");
			}
		}
		double a[CL->POP][ncoeff]; //elements of matrix A - which contains products of position data for all particles according to approximation method used. 
		double phi[CL->POP]; //compositional order parameter values for each particle in MPCD cell.
		double C[ncoeff][ncoeff]; //elements of cofactor matrix of A^TA.
		double coeff[ncoeff]; //coefficients of the polynomial order parameter approximnation.
		// the following variables are the counterparts to a, phi, C, and coeff when finding density field approximation - but this method needs fixing so the associated code has been commented out.
		double a_rho[CL->POP+1][ncoeff];
		double rho[CL->POP+1];  
		double C_rho[ncoeff][ncoeff]; 
		double coeff_rho[ncoeff];

		double det = 0.0; //determinant of matrix A^TA.
		double chemotaxisFactor[2][DIM]; 
		//zero all necessary values
		for (i=0;i<DIM;i++) {
			Q[i]=0.0;
			Fbulk[i]=0.0;
			Fint[i]=0.0;
			Fden[i]=0.0;
			Fchem[i]=0.0;
			chemotaxisFactor[0][i] = 0.0;
			chemotaxisFactor[1][i] = 0.0;
		}
		for (i=0;i<CL->POP;i++) {
			phi[i]=0; 
			rho[i]=1; 
			for (j=0;j<ncoeff;j++) {
				a[i][j]=0.0; 
				a_rho[i][j]=0.0;
			}
		}
		rho[CL->POP] = -1.0*(CL->POP); //one heavy phantom particle placed at centre of cell with "negavtive" density 
		for (i=0;i<ncoeff;i++) {
			coeff[i]=0.0;
			coeff_rho[i]=0.0;
			for (j=0;j<ncoeff;j++) {
				C[i][j]=0.0; //elements of matrix A
				C_rho[i][j]=0.0;
			}
		}
		i=0;
		//MPC particles
		tmpc = CL->pp;
		N = (double)(CL->POP);
		while( tmpc!=NULL ) {
			//fetch position data for this particle, relative to the geometric centre of the cell
			for( d=0; d<DIM; d++ ) Q[d] = tmpc->Q[d] - 0.5 - (double)((int)tmpc->Q[d]);
			if (tmpc->SPID==0) 
			{
				//this is type A
				phi[i] = 1.0/(double)CL->POP; 
				Na+=1; 
			}
			else if (tmpc->SPID==1)
			{
				//this is type B
				phi[i] = -1.0/(double)CL->POP; 
				Nb+=1;
			}
			//compute A matrix elements
			if (DIM == 2)
			{
				switch(orderApprox)
				{
					case(2): // 2D parabolic
						a[i][0] = Q[0]*Q[0]; //x^2
						a[i][1] = Q[1]*Q[1]; //y^2
						a[i][2] = Q[0]*Q[1]; //xy
						a[i][3] = Q[0]; //x
						a[i][4] = Q[1]; //y
						a[i][5] = 1.0; //1
						//the code below is intended to set up a matrix a_rho for finding the density field approximation, however this method needs revising. 
						/* 
						a_rho[i][0] = Q[0]*Q[0]; //x^2
						a_rho[i][1] = Q[1]*Q[1]; //y^2
						a_rho[i][2] = Q[0]*Q[1]; //xy
						a_rho[i][3] = Q[0]; //x
						a_rho[i][4] = Q[1]; //y
						a_rho[i][5] = 1.0; //1
						if (i==(CL->POP-1)){//this is for the "phantom" particle at cell centre.
							double q0 = 0.0; 
							double q1 = 0.0;
							a_rho[CL->POP][0] = q0*q0; //x^2
							a_rho[CL->POP][1] = q1*q1; //y^2
							a_rho[CL->POP][2] = q0*q1; //xy
							a_rho[CL->POP][3] = q0; //x
							a_rho[CL->POP][4] = q1; //y
							a_rho[CL->POP][5] = 1.0; //1
						}
						*/
						break;
					case(3): // 2D cubic  
						a[i][0] = Q[0]*Q[0]*Q[0]; //x^3
						a[i][1] = Q[1]*Q[1]*Q[1]; //y^3
						a[i][2] = Q[0]*Q[0]*Q[1]; //x^2y
						a[i][3] = Q[1]*Q[1]*Q[0]; //y^2x
						a[i][4] = Q[0]*Q[0]; //x^2
						a[i][5] = Q[1]*Q[1]; //y^2
						a[i][6] = Q[0]*Q[1]; //xy
						a[i][7] = Q[0]; //x
						a[i][8] = Q[1]; //y
						a[i][9] = 1.0; //1
						// also set a_rho values here (when this method is fixed)
						break;
				}
			}
			else if (DIM == 3)
			{
				// 3D parabolic  
				a[i][0] = Q[0]*Q[0]; //x^2
				a[i][1] = Q[1]*Q[1]; //y^2
				a[i][2] = Q[2]*Q[2]; //z^2
				a[i][3] = Q[0]*Q[1]; //xy
				a[i][4] = Q[1]*Q[2]; //yz
				a[i][5] = Q[2]*Q[0]; //zx
				a[i][6] = Q[0]; //x
				a[i][7] = Q[1]; //y
				a[i][8] = Q[2]; //z
				a[i][9] = 1.0; //1
				// also set a_rho values here (when this method is fixed)
			}
			tmpc = tmpc->next;
			i++;
		}
		// note that if one particle population is absent, then it is not necessary (or possible) to phase separate in this cell, and swimmers should undergo usual run-tumble dynamics. 
		if (Na != 0 && Nb != 0) 
		{
			// regardless whether using 2D or 3D, find the C matrix and the determinant according to the number of coefficients.
			if (ncoeff==6){
				cofactors6x6(a,CL->POP,C,det);
			}
			else if (ncoeff==10){
				cofactors10x10(a,CL->POP,C,det);
			}
			//find the coefficients 
			for (i=0;i<ncoeff;i++)
			{
				for (j=0;j<ncoeff;j++)
				{
					for (k=0;k<CL->POP;k++)
					{
						coeff[i]+=C[i][j]*a[k][j]*phi[k]; 
						double cvec2[1];
						cvec2[0] = coeff[i]; 
						// note that for cells with five particles or fewer (likely to occur near interfaces), matrix A will be singular and the coefficients will be NAN.
						// There is a non-zero but quickly negligible likelihood for the same to occur for higher particle numbers. 
						// Where NAN values are returned, do not apply any phase separation force. 
						if (checkNAN_vec(cvec2,1)!=0)
							{
								coeff[i]=0.0;
								break;
							}
					}
				}
			}
			//The code below determines the continuous approximation of the density field given an appropriate a_rho matrix, and computes the appropriate 
			//force resulting from the equation of state (using one of two different methods). This code has been commented out here as the method for 
			//defining a_rho above needs revising in order to determine an appropriate incompressibility force. 

			/* 
			if (incompressibility!=0) 
			{
				switch(ncoeff)
				{
					case(6):
						cofactors6x6(a_rho,CL->POP+1,C_rho,det);
						break;
					case(10):
						cofactors10x10(a_rho,CL->POP+1,C_rho,det);
						break;
				}
				//find the coefficients 
				for (i=0;i<ncoeff;i++)
				{
					for (j=0;j<ncoeff;j++)
					{
						for (k=0;k<CL->POP+1;k++)
						{
							coeff_rho[i]+=C_rho[i][j]*a_rho[k][j]*rho[k]; 
							double cvec3[1];
							cvec3[0] = coeff_rho[i];
							if (checkNAN_vec(cvec3,1)!=0)
								{
									coeff_rho[i]=0.0;
									break;
								}
						}
					}
				}
				if (incompressibility==1)
				{
					//this is for ideal-gas-like pressure tensor so that f = -chi*grad(rho)
					Fden[0]= -1.0*incompStrength*coeff_rho[3];
					Fden[1]= -1.0*incompStrength*coeff_rho[4];
					//also do 3d and cubic cases 
				}
				if (incompressibility==2)
				{
					//this is for non-linear state function so that f = -chi*rhp*grad(rho)
					Fden[0]= -1.0*coeff_rho[5]*incompStrength*coeff_rho[3];
					Fden[1]= -1.0*coeff_rho[5]*incompStrength*coeff_rho[4];
					//also do 3d and cubic cases 
				}
			}
			*/
			//use coefficients to determine continuous approximation of phi and its derivatives 
			if (DIM == 2)
			{
				switch(orderApprox)
				{
					case(2): // 2D parabolic 
						Fbulk[0] = -fabs(coeff[5])*(tau+3*b*coeff[5]*coeff[5])*coeff[3];
						Fbulk[1] = -fabs(coeff[5])*(tau+3*b*coeff[5]*coeff[5])*coeff[4];
						Fint[0] = 0;
						Fint[1] = 0;
						Fchem[0]= coeff[3]; 
						Fchem[1]= coeff[4];
						//it is possible to insert other force components here as desired - eg Fden, when this method is corrected. 
						break;
					case(3): // 2D cubic
						Fbulk[0] = -fabs(coeff[9])*(tau+3*b*coeff[9]*coeff[9])*coeff[7];
						Fbulk[1] = -fabs(coeff[9])*(tau+3*b*coeff[9]*coeff[9])*coeff[8];
						Fint[0] = kappa*fabs(coeff[9])*(6*coeff[0]+2*coeff[3]);
						Fint[1] = kappa*fabs(coeff[9])*(6*coeff[1]+2*coeff[2]);
						Fchem[0]= coeff[7];
						Fchem[1]= coeff[8];
						//it is possible to insert other force components here as desired - eg Fden, when this method is corrected. 
						break;
				}
			}
			else if (DIM == 3)
			{
				//recall that in the code above we have ensured 3D systems only execute the parabolic approximation. 
				Fbulk[0] = -fabs(coeff[9])*(tau+3*b*coeff[9]*coeff[9])*coeff[6];
				Fbulk[1] = -fabs(coeff[9])*(tau+3*b*coeff[9]*coeff[9])*coeff[7];
				Fbulk[2] = -fabs(coeff[9])*(tau+3*b*coeff[9]*coeff[9])*coeff[8];
				Fint[0] = 0;
				Fint[1] = 0;
				Fint[2] = 0;
				Fchem[0]= coeff[6];
				Fchem[1]= coeff[7];
				Fchem[2]= coeff[8];
			}
			i=0;
			// MPC particles
			tmpc = CL->pp;
			while( tmpc!=NULL ) {
				id = tmpc->SPID;
				velMag=0.0;
				for(i=0;i<DIM;i++)
				{
					velMag+=tmpc->V[i]*tmpc->V[i];
					Vtemp[i] = tmpc->V[i]; 
				}
				velMag=sqrt(velMag);
				//note that we add velocity components to a temporary velocity vector, and later set the velocity to this temporary variable. 
				if (id==0) 
				{
					//this is type A
					for( j=0; j<DIM; j++ ) Vtemp[j] += (FbulkMultiplier*Fbulk[j]*(0.5*N/Na)+FintMultiplier*Fint[j]*(0.5*N/Na)+FdenMultiplier*Fden[j])*timestep;
				}
				else if (id==1)
				{
					//this is type B
					for( j=0; j<DIM; j++ ) Vtemp[j] += (-1.0*FbulkMultiplier*Fbulk[j]*(0.5*N/Nb)-FintMultiplier*Fint[j]*(0.5*N/Nb)+FdenMultiplier*Fden[j])*timestep;
				}
				// now normalise the velocity and give it the original velocity magnitude, to maintain the temperature of the system. 
				tempMag = dotprod(Vtemp,Vtemp,DIM);
				if (tempMag!=0)
				{
					norm(Vtemp,DIM);
				}
				else
				{
					for(i=0;i<DIM;i++) Vtemp[i]=0.0;
				}
				for(i=0;i<DIM;i++) tmpc->V[i]=velMag*Vtemp[i];
				//Increment link in list
				tmpc = tmpc->next;
				i++;
			}
			//Swimmer monomers
			switch(orderApprox)
			{
				case(2):
					chemotaxisFactor[0][0]=K/((K+0.5*(1.0+ N*coeff[5]))*(K+0.5*(1.0+ N*coeff[5])))*0.5*(N*coeff[3]);
					chemotaxisFactor[0][1]=K/((K+0.5*(1.0+ N*coeff[5]))*(K+0.5*(1.0+ N*coeff[5])))*0.5*( N*coeff[4]);
					chemotaxisFactor[1][0]=K/((K+0.5*(1.0- N*coeff[5]))*(K+0.5*(1.0- N*coeff[5])))*0.5*(- N*coeff[3]);
					chemotaxisFactor[1][1]=K/((K+0.5*(1.0- N*coeff[5]))*(K+0.5*(1.0- N*coeff[5])))*0.5*(- N*coeff[4]);
					break;
				case(3):
					chemotaxisFactor[0][0]=K/((K+0.5*(1.0+ N*coeff[9]))*(K+0.5*(1.0+ N*coeff[9])))*0.5*( N*coeff[7]);
					chemotaxisFactor[0][1]=K/((K+0.5*(1.0+ N*coeff[9]))*(K+0.5*(1.0+ N*coeff[9])))*0.5*( N*coeff[8]);
					chemotaxisFactor[1][0]=K/((K+0.5*(1.0- N*coeff[5]))*(K+0.5*(1.0- N*coeff[5])))*0.5*(- N*coeff[7]);
					chemotaxisFactor[1][1]=K/((K+0.5*(1.0- N*coeff[5]))*(K+0.5*(1.0- N*coeff[5])))*0.5*(- N*coeff[8]);
					break;
			}
			tsm = CL->sp;
			while( tsm!=NULL ) {
				if( tsm->HorM ) id = SS.MSPid;
				else id = SS.HSPid;
				velMag=0.0;
				for(i=0;i<DIM;i++)
				{
					velMag+=tsm->V[i]*tsm->V[i];
					Vtemp[i] = tsm->V[i]; 
				}
				velMag=sqrt(velMag);
				if (id==0) 
				{
					for( j=0; j<DIM; j++ ) Vtemp[j] += FswimMultiplier*chemotaxisFactor[0][j];
				}
				else if (id==1)
				{
					for( j=0; j<DIM; j++ ) Vtemp[j] += FswimMultiplier*chemotaxisFactor[1][j];
				}
				tempMag = dotprod(Vtemp,Vtemp,DIM);
				if (tempMag!=0)
				{
					norm(Vtemp,DIM);
				}
				else 
				{
					for(i=0;i<DIM;i++) Vtemp[i]=0.0;
				}
				for(i=0;i<DIM;i++) Vtemp[i]=velMag*Vtemp[i];
				for(i=0;i<DIM;i++) tsm->V[i]=Vtemp[i];
				//Increment link in list
				i++;
				//Increment link in list
				tsm = tsm->next;
				}
		}
    }
	else if( multiphaseMode==MPHPOINT ) multiphaseCollPoint(CL, SP, SS, KBT, MD_mode, CLQ, outP );
	else {
		printf( "Error: Multiphase interaction  technique unacceptable.\n" );
		exit( 1 );
	}
}