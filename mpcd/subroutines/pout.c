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
/* ************** VERSION NOTES ************* */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
void printVersionSummary( ) {
	printf( "Version 153\n\tComputeNemForces() was causing the simulation to crash --- rewrite\n" );
	printf( "\t\tAlso added another check in ComputeBendForces() and use sinc() in bendHarmonic() as it should be more stable\n" );
	printf( "\t\tPut a check in ComputeNemForces()\n" );
	printf( "\t\tPut a check in ComputeNemForces()\n" );

				printf( "\t\t(2) BC colloid that overlaps with a PBC feels a 'buoyancy force (swimmers can get a bit inside the colloid on the PBC side)'\n" );


	printf( "Version 152\n\tMany swimmers in PBC will start to drift\n" );
	printf( "\t\tIn swimmerForceDipole(), would sometimes get net force since only accelerated mpcd particles; now acc all\n" );
	printf( "\t\tAdded some checks in bendHarmonic() and bendNematic() in mssrd.c\n" );
	printf( "\t\tIn ComputeDispersionForcesSRDCell() and other MD subroutines associated with SRD put while loops on the weird PBCs since often crashes because this run outside of array size somehow\n" );
	printf( "\t\tIn ComputeForcesSRD() use ComputeDispersionForcesSRDCell() if MPCinMD otherwise use ComputeDispersionForces()\n" );
	printf( "Version 151\n\tThe polymer in a nematic is drifting and always aligns with x (never y nor z)\n" );
	printf( "\t\tRe-wrote bendNematic() to only apply two forces and to try to match the torque on the nematic to the torque on the polymer\n" );
	printf( "\t\tOnly enter bendNematic() if kNem*S>0\n" );
	printf( "\t\tModify to use the cosine angular potential instead of the harmonic angle\n" );
	printf( "\t\tAlso needed to include a check on the direction of the force (since nematic) and corresponding angle change\n" );
	printf( "\t\tTested that works on given set directors; HOWEVER, found polymer 'melts' nematic then found that INT fails\n" );
	printf( "\t\tTested that works on given set directors; HOWEVER, found polymer 'melts' nematic then found that INT fails (even without polymer)\n" );
	printf( "\t\tBacktracked from v144 and found added if-statement to LCcollision() in v150\n" );
	printf( "\t\tRemoved if( dotprod(DIR,DIR,_3D)<=TOL )\n" );
	printf( "Version 150\n\tRe-do bend polymer (last attempt in v110)\n");
	printf( "\t\tAdd ComputeBendForces(), ComputeNemForces(), bendHarmonic() and bendNematic() in mdsrd.c and mdsrd.h. BUT comment out content\n" );
	printf( "\t\tComputeNemForces() uses list2STD fene as a convenient pre-established list of nearest neighbors.\n" );
	printf( "\t\tAdd kBend, kNemMPC, bendE, nemE to simulation in mdtypes.h.\n" );
	printf( "\t\tIn mdsetup.c, create and call 	SetupBendList().\n" );
	printf( "\t\tIn mdmemory.c, FreeMemory() add bend. In mdmemory.c create GrowList3STD(), AddItem3STD() and ResetList3STD().\n" );
	printf( "\t\tCreate and add list3STD (a list of three beads for calculating bends) in mdtypes.h.\n" );
	printf( "\t\tIn mdmeasure.c, routine EvalProperties() add bendE and nemE.\n" );
	printf( "\t\tTo apply the force to the field, we commendeer magTorque_CL() and, in bendNematic(), give the force as the magnetic field (still use magnetic susceptibility)\n" );
	printf( "\t\tAdded if() to LCcollision() to only do LC collision if DIR isn't zeros\n" );
	printf( "Version 149\n\tCleaned up rotlet dipole and created run/tumble output\n");
	printf( "\t\tMade and fixed new input files\n");
	printf( "\t\tFixed segfault that happen when swimmers are present\n" );
	printf( "\t\tIn swimmerForceDipole(), swimmerDipole() and swimmerRotletDipole() the z-position of the temporary swimmer needed to be initialized to zero\n");
	printf( "\t\tNOTICE: Run/Tumble broken --- works but crazy slow!\n");
	printf( "Version 148\n\tModify the rotlet dipole\n");
	printf( "\t\tPut force and rotlet dipole in combined subroutine swimmerDipole()\n");
	printf( "\t\tIn swimmerDipole(), shift (and bin) all MPCD particls such that swimmer at centre of cell then shift back\n");
	printf( "\t\tMove swimmerDipole() to AFTER the MPCD collision and before streaming (but after grid shifting back)\n");
	printf( "\t\tAdd swimmer orientation output\n");
	printf( "\t\tAdd swimmer run/tumble times and angle change output\n");
	printf( "Version 148\n\tChange input/output so that md.inp MUST be in same directory as others\n");
	printf( "\t\tBy calling launchMD() in initializeSIM(), the simMD was not initialized in main(). Just moved launchMD() into main().\n");
	printf( "Version 147\n\tFix the rotlet dipole:\n");
	printf( "\t\tPrevious version did not conserve translational momentum\n");
	printf( "\t\tDon't apply torque about head but rather about CM of MPCD cell\n");
	printf( "\t\tAlso don't apply constant torque but rather do solid body rotation i.e. give all same change to angular speed\n");
	printf( "Version 146\n\tModify the swimmer:\n");
	printf( "\t\tSplit the propulsion force between the head and body monomers \n");
	printf( "\t\tApply a torque on the head and equal/opposite on the tail using swimmerRotletDipole()\n");
	printf( "\t\tAdd the ability to fix the swimmers in a plane\n");
	printf( "Version 145\n\tMake DUMBBELL_FIXED swimmer-type so that it is the dumbell swimmer but its movement isn't integrated \n");
	printf( "Version 144\n\tSolved North-East swimmer bias by re-binning and recalculating localPROP after GridShift Back and Boundary Translation steps\n");
	printf( "\tIncluded XYZPBC checking in binSwimmers\n");
	printf( "\tAdded initializing swimmer orientations within planes\n");
	printf( "Version 143\n\tNoticed that solving the eigen system uses some equality-expressions. Replace with feq() and fneq() from mtool\n" );
	printf( "\t\tcalcW(), calcW_BC(), calcW_swimmer, calcW_MD() and surf_func(), WALL.ROTSYMM[]==4.0\n" );
	printf( "\t\tcrosstime() and crosstimeReverse(), WALL.P[]==1.0, etc\n" );
	printf( "\t\tBC_BCcollision(), stillWall->MVN==1.0 && stillWall->MVT==1.0, etc\n" );
	printf( "\t\tnormal(), stillWall.ROTSYMM[]==4.0 and WALL.P[]==1.0, etc\n" );
	printf( "\t\tcheckSim(), SP[i].RFC==0. and SS.DS==0.0\n" );
	printf( "\t\tinitializeSIM(), (WALL+i)->O[0]==0.0, etc\n" );
	printf( "\t\tgenrand_maierSaupeMetropolis_2D(), ang0==-0.5*pi\n" );
	printf( "\t\toriBC(), WALL->MUN==0.0 and WALL->MUT==0.0\n" );
	printf( "\t\tIn vicsekAndersenMPC() and vicsekLangevinMPC(), ...\n" );
	printf( "\t\tIn atan2(), mominert() and dim_vol(), ...\n" );
	printf( "\t\tIn eigenvectors2x2() and eigenvalues3x3(), m[][]==m[][] or m[][]==eigenval[]\n" );
	printf( "\t\tIn bcin() and setPBC(),...\n" );
	printf( "\t\tIn shiftBC(), WALL->P[]!=1.0\n" );
	printf( "\t\tIn checkSim() and initializeSIM(), WALL[i].V[]!=0.0, etc\n" );
	printf( "\t\tIn dipoleAndersenROT_LC(), dipoleAndersenROT_LC(), activeSRD(), chateAndersenMPC(), chateLangevinMPC() and dipoleAndersenMPC(), ACT!= 0.\n" );
	printf( "\t\tIn norm() and normCopy(), l!= 0.\n" );
	printf( "\t\tIn eigenvalues3x3(), ...\n" );
	printf( "\t\tIn coordout(), MASS!=0.0\n" );
	printf( "\t\tIn bcin(), WALL->AINV[]!=0.0\n" );
	printf( "\t\tIn setswimmers() and runTumbleDynamics(), SS->sizeShrink!=1.0, etc\n" );
	printf( "\t\tIn scaleT(), VEL[j]!=0.0\n" );
	printf( "\tZero all mallocs to be super safe\n" );
	printf( "\tMy atan2 went atan2(x,y), whereas normally atan(y,x). Maybe icc uses math.h rather than mine. Switch my order.\n" );
	printf( "\tcrossprod() took a dimension but in 2D didn't give a result in z-dir. Force 3D input.\n" );
	printf( "\tsignedAngle() became wrong when atan2() got flipped.\n" );
	printf( "\tFinally solved by making DIR and dU 3D in LCcollision().\n" );
	printf( "Version 142\n\tTry to resolve all compiler warnings on SciGrid Macs and my desktop Mac\n" );
	printf( "\t\tIn setcoord(), gave filesuffix a size\n" );
	printf( "Version 141\n\tBabi Rae and I found active nematic wasn't working because PLANE was not setting ROTSYMM=4\n" );
	printf( "\t\tChanged in dipoleAndersenROT_LC() and dipoleAndersenMPC()\n" );
	printf( "\tAlso resolved warnings:\n" );
	printf( "\t\tIn oriBC(), UT0 was set but not used - delete\n" );
	printf( "\t\tIn rand.c, time(NULL) was replaced throughout by the timeval structure and microseconds used\n" );
	printf( "\t\tIn rand.c, <unistd.h> was included so getpid() works\n" );
	printf( "\t\tIn mdfiles.c, ReopenDataFiles() included a read variable as a check on ftruncate()\n" );
	printf( "\t\tIn swimmerDipole(), permanently removed the check to see if swimmer moving in direction of orientation\n" );
	printf( "Version 140\n\tSame problem as v138, v138 and v137\n" );
	printf( "\t\tSeem to have higher density, lower pressure and inward velocity of MPCD particles at walls\n" );
	printf( "\t\tHack the walls\n" );
	printf( "\t\tMade minor changes to get rid of warning when Wan compiles on her macbook\n" );
	printf( "\t\tI removed all 'inline' in mdsrd and in mdutils, which probably slows things down but hopefully compiles for her\n" );
	printf( "Version 139\n\tSame problem as v138 and v137\n" );
	printf( "\t\tMake a debug to print out all the swimmer details\n" );
	printf( "Version 138\n\tSame problem as v137\n" );
	printf( "\t\tTake swimmers and monomers our of localFLOW()\n" );
	printf( "\t\tMake two new kinds of swimmers (1) dumbbells with excluded volume between swimmer\n" );
	printf( "\t\t(2) Dumbbells without the dipolar counterforce\n" );
	printf( "\t\t(3) Remove the undulator-type swimmer (which wasn't implemented anyway)\n" );
	printf( "\t\tTo do this, I changed my velocity Verlet integrator\n" );
	printf( "\t\tIt appears that the swimming issue goes away if the random shift is turned off\n" );
	printf( "Version 137\n\tIn strong magnetic fields, swimmers tend to move to lower left corner\n" );
	printf( "\t\tThe torques can be very large so break into stepsMD steps\n" );
	printf( "\t\tHad to put swimmerMagTorque() into integrateSwimmers()\n" );
	printf( "\t\tTired of swimmers getting stuck between walls. Don't allow PBC initialization\n" );
	printf( "Version 136\n\tRarely a swimmer escapes from the system\n" );
	printf( "\t\tBug changes if change how often print position and/or energy to file\n" );
	printf( "\t\tintegrateSwimmers() changes the positions ever so slightly and so can break BCs" );
	printf( "\t\tI was doing the BCs *before* the verlet integration --- rather than after. Switching order fixes it" );
	printf( "Version 135\n\tThe last version worked BUT PBCs fail.\n" );
	printf( "\t\tStrategy: First improve algorithm THEN later fix PBCs\n" );
	printf( "\t\t(a) Make an interaction matrix for each species (which is as long as the number of species)\n" );
	printf( "\t\t    Leaving swimmers for now\n" );
	printf( "\t\t(b) Set loops in andersenMULTIPHASE()\n" );
	printf( "\t\t(c) Remove phi throughout\n" );
	printf( "\t\t(d) Put the swimmers and polymer back in\n" );
	printf( "\t\t(e) I put back in RVsum but that seems to generate momentum - take out BUT must check momentum cons\n" );
	printf( "Version 134\n\tWork on binary fluid with Romain in Oxford\n" );
	printf( "\t\tWe totally changed the algorithm. Now each species has self- and inter-species interaction tensor\n" );
	printf( "\t\tlike soluability parameters\n" );
	printf( "\t\tCurrently hacked together for a binary case in 2D but will be extended to multiphase fluids in 3D\n" );
	printf( "Version 132\n\tTry to fix the LC torque on a colloid again (having learned from the swimmer)\n" );
	printf( "\t\tIn oriBC() call a new routine torqueLCBC()\n" );
	printf( "Version 131\n\tAdd a multiplicative 'shrinking' or 'expanding' of swimmers' spring const during tumbling\n" );
	printf( "Version 130\n\tGetting rare escaped LC-MPC particles\n" );
	printf( "\t\tSeems to be only for homogeneous (planar) BC\n" );
	printf( "\t\tSeems to not occur if homeotropic (normal) or no anchoring BC\n" );
	printf( "\t\tSeems to not occur if alignment fixed in x, y or z\n" );
	printf( "\t\tFixed in 2D: Pairs of particles would get zero velocity then divide by zero\n" );
	printf( "\t\tStill get rare escape in 3D. ");
	printf( "\t\tStill only for homogeneous (not homeotropic, xyz, or no anchoring)\n" );
	printf( "\t\tOnly for homogeneous for BC with z normal\n" );
	printf( "\t\tEscape occurs when after a cell's director is zero.\n");
	printf( "\t\tThe cell's director is zero when the pop<=3 and z-comp of all is zero.\n");
	printf( "\t\tThe tensor order parameter LOOKS 2D and isn't being solved!\n" );
	printf( "\t\tGot eigenvalues correct BUT not eigenvectors\n" );
	printf( "\t\tIn eigenvectors3x3() the check for zero-eigenvectors was wrong.\n" );
	printf( "\t\tMake the swimmers magnetotaxic by giving them a magentic moment strength.\n" );
	printf( "\t\tAdd swimmerMagTorque(), which is different than for nematic particles because\n" );
	printf( "\t\tthe swimmers have a vectorial moment and align their head (rather than a susceptibility)\n" );
	printf( "Version 129\n\tGetting rare escaped MPC particles --- NOT SOLVED\n" );
	printf( "\t\tONLY happens when LC turned on! ---- NOT TRUE happens for iso too\n" );
	printf( "\t\tAdd error message (to terminal) and exit() if reading input fails.\n" );
	printf( "\t\tAdd error message (to terminal) and exit() if input files are too short.\n" );
	printf( "\t\tEnsure that can only run with LC on IF MFPOT>0\n" );
	printf( "Version 128\n\tThere is a density depletion at walls AND velocity towards walls (due to depletion?) I think it is due to the MPC and shift.\n" );
	printf( "\t\tPut the system shift into a routine gridShift_all()\n" );
	printf( "\t\tFix the ghost particles. Variance in Gaussian was wrong and took out BC velocity frame" );
	printf( "\t\tFix the colloid-LC issue - Would crash if torque happened to be zero. Fix with IF-statement" );
	printf( "\t\tColloid did not respect hard BCs. This is just because it was using the radius to power PR and not just the radius" );
	printf( "\t\tTo get colloid-LC fix, change oriBC() to model an effective traction force --- required adding an effective rod size to input.inp" );
	printf( "Version 127\n\tFlow still not fixed.\n" );
	printf( "\t\tMake routines for rudimentary PBCs, box and channels " );
	printf( "\t\tIn mpc.c (line 3410) replace MPC_BCcollision() with rudimentaryPBC_box() --- these still show a 'slip' at PBCs." );
	printf( "\t\tAdd RSHIFT to swimmers." );
	printf( "\t\tCheck that initial zeroing of momentum works (yes). Change flag to RFRAME & random shift flag name to GALINV." );
	printf( "\t\tWorking on flow: Move acceleration (acc_all()) BEFORE calculate VCM for collision. This FIXED the flow.\n" );
	printf( "Version 126\n\tLet MPC and BC masses be a double rather than int\n" );
	printf( "Version 125\n\tWhen channel flow occurs off axis, there is a depletion of MPC particles at corners\n" );
	printf( "\t\tAdd an input flag for whether or not to do the random shift(make RSHIFT=0 in ranshift())\n" );
	printf( "\t\tAdd an flag for each axis on whether or not a PBC\n" );
	printf( "\t\tIf flagged as PBC then wrap shifted MPC particles in binning\n" );
	printf( "\t\tThis is JUST binning so it might effect things like calculating CM of the cell\n" );
	printf( "Version 124\n\tGive an orientation to the BC surfaces.\n" );
	printf( "\t\tAllows us to set objects at arbitary orientation AND should allow rotational diffusion of non-symm objects.\n" );
	printf( "\t\tNOTICE: The current implementation is very wasteful:\n" );
	printf( "\t\tEvery particle is rotated about the centre of each BC - though simplest, there are very many particles\n" );
	printf( "\t\tThis is done through rotateBC() and rotatebackBC() (after and before, respectively, every shiftBC() shiftbackBC())\n" );
	printf( "\t\tASSUMES rotational matrices commute (ie order doesn't matter; only sign matters) but this may not be true (?)\n" );
	printf( "\t\tASSUMES moment of inertia tensor can be approximated as nearest ellipsoid but this MUST be improved in the future\n" );
	printf( "\t\tCurrently works for diffusing spheres or fixed anything else but general shapes can't diffuse so the problem is likely in BC_MPCcollision()\n" );
	printf( "\t\tColloids appear to have trouble diffusing across PBCs\n" );
	printf( "Version 123\n\tReturning to the fucking binary fluid.\n" );
	printf( "\t\tRevamp the BCs by adding more parameters can (obviously) get more shapes.\n" );
	printf( "\t\t(1) Make A[3] be for division (as is costumary in ellipse formuli).\n" );
	printf( "\t\t(2) Give each term its own power (not just one) and have the 'radius' be subject to a power too.\n" );
	printf( "\t\t(3) Have an absolute operator flag. If 0 then no abs() but if 1 then the absolute value of each term is taken.\n" );
	printf( "\t\t    Any power that is not unity or even should probably use the absolute operator.\n" );
	printf( "\t\t(4) Following Gielis (2003), introduce parameters m and l to allow m- and l-fold rotational symmetry.\n" );
	printf( "\t\tMake a debug mode that keeps track of any potential escaped particles.\n" );
	printf( "Version 122\n\tAdd an extension/shrinking phase to run/tumble dynamics.\n" );
	printf( "\t\tWas getting so much energy because the FENE spring is too sharp. Switched to Hookean always\n" );
	printf( "\t\tPut binary fluid back to about the CM to match Romain's advice\n" );
	printf( "Version 121\n\tMoving on while Romain works on binary algorithm.\n" );
	printf( "\t\tCleanup the debug options\n" );
	printf( "\t\tWith many BCs the algorithm gets stuck - in MPC_BCcollision() add a slow extra rewind if stuck then re-thermalize if really bad\n" );
	printf( "\t\tShrinking swimmers break the FENE spring so try to put hard wall on extension\n" );
	printf( "\t\t\tReveals that swimmer 'hops' when shrinks or extends --- because I integrated over more steps but kept time step the same\n" );
	printf( "Version 120\n\tRejected. Failed version. Went back to v119\n" );
	printf( "Version 119\n\tThis is a back-tracked version to try to reproduce the surface tension that in binary fluids, which we had and then lost\n" );
	printf( "\t\tv119_a - Just implement the new gradient - worked\n" );
	printf( "\t\tv119_b - Implement the FIRST new grad laplacian (incorrect; doesn't even include PHI!!!) - worked (with surface tension)\n" );
	printf( "\t\tv119_c - Add phi to the grad laplacian - FAILS\n" );
	printf( "\t\tv119_d - Try dividing by N - Doesn't help\n" );
	printf( "\t\tv119_e - Get rid of phi and N (ie back to v119_b) then use the geometic centre of the cell rather than the center of mass - FAILED\n" );
	printf( "\t\tv119_f - Add back in phi (with geometric center but no N yet) - FAILED\n" );
	printf( "\t\tv119_g - divide by N - FAILED\n" );
	printf( "\t\tv119_h - Use Romain's NEW gradPhi and grad lap derivatives - FAILED\n" );
	printf( "\t\tv119_i - Use Romain's NEW gradMu and VMU - Not crap but FAILED to get surface tension\n" );
	printf( "\t\tv119_j - Same as i BUT remove PHI from gradLaplacePhi (keep N) - \n" );
	printf( "\t\tv119_k - Return the code to i and send to Romain - \n" );
	printf( "Version 117\n\tMake the run and tumble times exponentially distributed\n" );
	printf( "Version 116\n\tThe Poisson distribution of swimmers shouldn't just depend on runTime/dt but should have a coefficient in front to control shape\n" );
	printf( "\t\tUse better functions for the theoretical viscosity and self-diffusivity in theory_trans()\n" );
	printf( "\t\t--- realized that I didn't have the angular momentum conserving Langevin version so add langevinROT()\n" );
	printf( "\t\tWhen I used Slurm to launch my jobs at Rockefeller, they all produced the same data --- use (process ID)*(clock) instead of just clock to seed random number generator\n" );
	printf( "Version 115\n\tWhen the swimmer velocity is VERY small then the swimmer moves backwards!\n" );
	printf( "\t\tHypothesis: Setting the velocity effectively to zero still applies a force to the fluid on average, which entrains the swimmer\n" );
	printf( "\t\tSo try to apply a constant force rather than constant velocity\n" );
	printf( "\tWhen the swimmer passes PBCs, often a strong flow occurs at random\n" );
	printf( "\t\t\n" );
	printf( "Version 114\n\tShrink the swimmers during the tumble phase\n" );
	printf( "\t\tHad to give LJ sigma and FENE ro to each swimmer instead of just swimmer type cuz it is different for each\n" );
	printf( "\t\tUse Normal distr approx for large lambda\n" );
	printf( "\t\tRandomly initiate as running or tumbling based on percentage of time in each\n" );
	printf( "\t\tUniform random dist of times through the initial run/tumble time\n" );
	printf( "\t\tOutput run/tumble phase flad in swimmer.dat\n" );
	printf( "Version 113\n\tAdd swimmer parameters for puller-type or pusher-type (by a dipole strength DS)\n" );
	printf( "\t\tAdd separate run and tumble average times\n" );
	printf( "\t\tWrite a random number generator for producing random integers belonging to a Poisson dist genrand_Poisson()\n" );
	printf( "\t\tThe phantom force position is now about the center of mass position of the swimmer (not mirrored about the mid-point)\n" );
	printf( "\t\tOnly do the force dipole and velocity if in the running phase\n" );
	printf( "Version 112\n\tFix the velocity averaging problem by moving the running sum of the velocity field from outputResults() to timestep()\n" );
	printf( "\t\tavVel.dat now agrees with the average of flowfield.dat by recalculating average vel just before output\n" );
	printf( "\t\tUnfortunately, it means that I need to rebin things --- I should have a flag to see if things need to be re-binned or not\n" );
	printf( "\t\tFix the swimmer/bc problem by only applying the phantom force to the fluid if head velocity and orientation point in same direction\n" );
	printf( "\t\tWould sometimes get stuck in first integrateSwimmers() --- didn't check shifted BCs. Fixed\n" );
	printf( "\t\tIt still sometimes gets stuck at later integrateSwimmers() instances --- found MPCD particles with NAN velocity coming form integrateSwimmers(), which doesn't use MPCD paticles\n" );
	printf( "\t\t\t---swimmerDipole() would divide by cell mass even if no particles present. Moved swimmer integration. Fixed\n" );
	printf( "Version 111\n\tMake the dumb-swimmers interact with hard obstacles\n" );
	printf( "\t\tMake sure that initial placement respects BCs\n" );
	printf( "\t\tDebug swimmer collision with BCs (in mdbc.c) --- BCs need to be checked in MD integrator\n" );
	printf( "Version 110\n\tAdd a bending potential to the MD polymer, a potential to align with the MPCD nematic and a weighting of the MPCD nematic to the MD polymer\n" );
	printf( "\t\tHold off on input; Follow fene\n\t\tAdd kBend, kNemMPC, weightNem, bendE, nemE to simulation in mdtypes.h.\n" );
	printf( "\t\tCreate list3STD (a list of three beads for calculating bends)\n" );
	printf( "\t\tAdd bendE and nemE in energy.dat in mdfilesdef.h\n" );
	printf( "\t\tIn mdmeasure.c and routine EvalProperties() add bendE and nemE. In mdmemory.c and FreeMemory() add bend. In mdmemory.c create GrowList3STD() and AddItem3STD()\n" );
	printf( "\t\tAdd kBend,kNemMPC and weightNem to SetupParameters() in mdsetup.c\n" );
	printf( "\t\tAdd ComputeBendForces() and ComputeNemForces() in mdsetup.c\n" );
	printf( "\t\tBending potential incomplete in this version.\n" );
	printf( "\t\tFixed dumb-swimmers (weren't dipolar) by making a head and a middle that are attached by springs, interact with the fluid and a phantom tail that just acts on the fluid at the mirrored position of the head about the middle.\n" );
	printf( "Version 109\n\tDebug the time problem\n\t\tmake all input times 'iterations' rather than simulation time in sim units\n" );
	printf( "\t\tFixed some typos in initiating velocity distribution (didn't see before when mass=1 always)\n" );
	printf( "Version 108\n\tWe've settled on a binary collision operator but now we implement better derivatives for grad(phi) and gard(laplacian(phi))\n" );
	printf( "Version 107\n\tChanged the binary collision operation again\n" );
	printf( "\t\tLet the swimmers' heads and tails have different phi and mass values\n" );
	printf( "\t\tSwimmers' needed to be debugged and are now far more robust\n" );
	printf( "\t\tModified correlation functions - will take longer to calculate but use more points and so should be more accurate\n" );
	printf( "\t\tGot segfault if many swimmers used: Initialized tails could be outside of control volume. Also fixed list problem\n" );
	printf( "\t\tSwimmers in andersenMULTIPHASE() go nuts (even if activity turned off): Wasn't summing the random velocities of each type\n" );
	printf( "\t\tChanged order of BCs, dipole kicks and MD integration of swimmers in order to produce more efficient swimmer\n" );
	printf( "Version 106\n\tSame as 105 except changed the chemical potential part to NOT include (1-phi^2)\n" );
	printf( "\t\tOutput instantaneous pressure (actually stress tensor) field - need for laplace pressure\n" );
	printf( "Version 105\n\tSame as 104 except now the collision operation has no 'interspecies' noise and all the particles feel the same chemical potential\n" );
	printf( "Version 104\n\tNew version of the binary fluid based on flux from chemical potential\n" );
	printf( "Version 103\n\tChanged the 'temperature' part of the binary SEP\n" );
	printf( "\t\tAdded energy difference.\n" );
	printf( "Version 102\n\tFirst implementation of a binary fluid\n" );
	printf( "\t\tIn localPROP() division by mass was inside particle loop ---moved out\n" );
	printf( "\t\tCURRENTLY have 3 surface tension version HARDCODED in ie v102a,v102b,v102c,v102d\n" );
	printf( "Version 101\n\tv100 fixed the sticking to the walls in 2D but crashes in 3D. Fix: \n" );
	printf( "\t\tProblem was that histograms took up too much memory in such huge 3D systems. Just separated from other output\n" );
	printf( "\t\tAdd an active version that uses the mean activity (times rho) rather than the sum\n" );
	printf( "Version 100\n\tSticking to the wall happens in isotropic fluids: FIX THIS\n" );
	printf( "\tRewrite the collisions with the walls entirely\n" );
	printf( "\t\tTranslate BCs at end of timestep rather than beginning\n" );
	printf( "\t\tAccelerate BCs after translating, rather than before\n" );
	printf( "\t\tCollision with BCs are added up and then the impulse is applied AFTER translating (keep track of change in velBC() but add in )\n" );
	printf( "\t\tDelete BConBC()\n" );
	printf( "\t\tDelete checkEscape_all()\n" );
	printf( "\t\tFor Berendsen thermostat changed the cutoff on TAU to be related to time step\n" );
	printf( "\t\tAll BC-BC collisions (except PBCs) are now bounce-back\n" );
	printf( "\t\tIn BC_MPCcollision() do collision in non-rotating frame ie subtract off angular velocity during collision\n" );
	printf( "Version 99\n\tRe-orientation of LC at mobile surface requires a torque at the surface that induces an impulse on the body's CM\n" );
	printf( "\tDoesn't work --- sticks to walls\n" );
	printf( "\tTo run isotropic fluid without crashing oriBC() needs an if-LC comment in front. Fixed.\n" );
	printf( "Version 98\n\tMake input file for swimmers (dumbbell only but undulator later)\n" );
	printf( "\t\tDIDN'T CHECK POINT THE SWIMMERS!!!\n" );
	printf( "\t\tAlso removed the mainmpc.c file since it was confusing and redundant to the code\n" );
	printf( "\t\t--- moved routines into more natural *.c files\n" );
	printf( "\t\tMade the active-MPCD algorithm conserve angular momentum in dipoleAndersenROT_LC()\n" );
	printf( "\t\tIf scalar order parameter is numerically found to be a bit >1 then set to unity\n" );
	printf( "Version 97\n\tDiscovered activity flips from extensive to contractile if alpha is small\n" );
	printf( "\t\tFix by making activity just alpha\n" );
	printf( "Version 96\n\tAdd additional output (enstrophy, etc)\n" );
	printf( "\t\tAdd additional histograms and modify how they're calculated\n" );
	printf( "\t\tOutput energy and enstrophy spectra, and output defect tracking\n" );
	printf( "Version 95\n\tAdd a damping/friction term to collision operator to go from wet to dry\n" );
	printf( "\t\tProblem with velocity gradient tensor\n" );
	printf( "Version 94\n\tDon't do any normalization of the correlation functions\n" );
	printf( "\t\tInclude (but don't use) cellVelSet() and cellVelForce() in mpc.c\n" );
	printf( "Version 93\n\tDebug nematic symmetry breaking that occurs when both HI coupling is on and tumbling parameter is large\n" );
	printf( "\t\tMust use gradients from lattice Boltzmann literature to 'ensure' isotropy --- i.e. D2Q9 and D3Q15\n" );
	printf( "Version 92\n\tModify the angular momentum conservation term in dipoleAndersenROT_LC() to account for force dipole\n\tOutput correlation functions (velocity, director and vorticity)\n\tOutput Binder cumulant for nematic-MPCD\n" );
	printf( "Version 91\n\tModifies Chate-like active MPCD to subtract off any residual force on the cell in chateAndersenMPC() and chateLangevinMPC()\n\tModifies dipole active MPCD to relax to a set velocity in either direction rather growing unconstrained\n" );
	printf( "\tNow subtract off mass weighted average of noise rather than just average in andersenMPC() and andersenROT().\n" );
	printf( "Version 90\n\tAdd a warmup time\n\tPut some huge routines into the main routine to capture the main beats of the program\n" );
	printf( "Version 89\n\tDouble checks torque on the MPCD side from the nematic side (no changes made)\n" );
	printf( "Version 88\n\tThis adds an output of the entire order tensor field AND the order tensor field in reciprocal space\n\tModify the output of the reciprocal space so that it is in the so-called 13-plane of Allen etal 1996 or O'Brien 2011\n" );
	printf( "\tAlso add output for orientational energy per cell\n\tAlso added a MUZ etc to the orientation BC so that can force planar orientation BCs to not have z-component\n" );
	printf( "Version 87\n\tIn strong magnetic fields the LC seems to melt\n\tFinally found the problem with the Freedericksz transition (or at least the magnetic relaxation time\n\t\t--- testing Freedericksz transition now).\n" );
	printf( "\toriBC() had a minor error when setting the orientation along a cartesian axis\n\tAlso minor changes to the action when the BCs can't be worked out.\n" );
	printf( "\tInstead of using replace(), now replacePos_WithCheck() is used\n\tand if the binning fails then the particles are replaced. \n" );
	printf( "Version 86\n\t\n\tTry to figure out the problem with the Freedericksz transition\n\tTry using local director instead of particle orientation\n\tFound another problem:\n" );
	printf( "\t\tWhen LC initialized with director all in y-direction, the LCcollision operation snaps the orientation down to x-direction\n\t\t- just do a slight offset during the initialization\n" );
	printf( "Version 85\n\tThe Chate algorithms are working right. Perhaps velocity needs to be normalized in calculating S?\n\tFixed. I think just the global average S was wrong (Chate needs large activities).\n" );
	printf( "\tAdd a dipole force in the direction of director\n\t\t(use the same routine as dipole force in direction of centre of mass vel)\n\tThe dipole force for the centre of mass shouldn't use VCM but the director calculated from the velocities\n" );
	printf( "Version 84\n\tAdd a checkpoint (both input and output)\n\tTo read the checkpoint set random seed to -1\n" );
	printf( "Version 83\n\tsaintillanRotRate() doesn't allow control of Leslie angle therefore, replace with larsonRotRate().\n\tThis requires having a 'tumbling parameter' for the nematogens and so need to change input\n" );
	printf( "\tUse tumbling parameter rather than aspect ratio so that can have tumbling parameter values greater than unity\n\t\t- If use aspect ratio, I can only get tumbling\n" );
	printf( "\t\t\tBUT if use too large a tumbling paameter then code crashes\n\tSmall (?) error in magTorque() caught.\n\t\tAngular velocity (w) needed to be normalized into unit vector before going into rodriguesRotation().\n" );
	printf( "\t\tDitto for jefferysTorque()\n\trotational diffusion coefficient times tumbling parameter must not be too large (put in check)\n\tOutput S4, in avS.dat\n\tAdd total nematic inteaction energy to energy.dat\n" );
	printf( "Version 82\n\tAdd a switch to use global S ---> LC=0 for isotropic, LC=1 for local S and LC=2 for global S\n\tWhen it is globalS, it should begin aligned rather than gain that alignment so a warning is added if\n\tNow MFPOT*2/DIM is passed to genrand_maierSaupe()\n" );
	printf( "\t\t---should make the isotropic-nematic transition occur at same place in all dimensions\n\tPut checks into readin() and init() when reading input files (just to get rid of the compile errors)\n" );
	printf( "Version 81\n\tWhen the particles are all given an initial orientation (0,1)  the director was still calculated to be (1,0), at least in 2D\n\tThe else-statement in eigenvectors2x2() was not always true - fixed\n" );
	printf( "\t\t---> I think this solved the defects from the walls!\n" );
	printf( "\tAlso eigenvectors3x3() failed for diagonal matrices - fixed\n\tMake the outputted S based on self, neighbours and next-nearest neighbours\n\tSomehow, these fixes lead to NANs in the 3D Q-order tensor\n" );
	printf( "\t\t--- The diagonal elements needed to be sorted\n\tGot rid of orientation in ghostPart()\n" );
	printf( "Version 80\n\tAppears to be bias in director direction to pi/2 and/or 5pi/2:\n\tInitialization of orientations only in +x quadrants\n\tIndex error in eigenvector choice (both in MPCD cell and global calculations) ---> FIXED\n" );
	printf( "\t\t---There was an error in version 79. Went back to 78 and it wasn't there so re-did\n\tInclude ghost particles in orientation collision through ghostPart()\n" );
	printf( "Version 78\n\tThe isotropic-to-nematic transition from simulations is not 1st order but the Maier-Saupe self-consistent theory is.\n\tTest the number MC annealing steps \n\t\t--- Not this --- leave ability to set MCINT and MCSLOPE\n" );
	printf( "\tUsing global S in orientation dist rather than local value makes it first order\n\tMade correction to oriBC()\n" );
	printf( "Version 77\n\tThe eigenvalues and vectors for calculating the scalar order parameter and director ALWAYS come in order\n\t\tso don't need to search.\n\tAlso S is better approximated as -2 times the smaller (negative) values than the large value\n" );
	printf( "\t\t(definitely true at low order since will then fluctuate about zero)\n\tAlso, I now feel confident enough to use the fact that the tensor order parameter is always symmetric\n\t\tin tensOrderParam(), tensOrderParamVel() and avOrderParam()\n" );
	printf( "\tUse low and high limits of US/KBT to approximate the orientation distribution (save computational time)\n" );
	printf( "Version 76\n\tThe order parameter and director is failing in 2D (works perfectly in 3D)\n\tIt WAS NOT THE ORDER PARAMETER CALCULATION. It was the initialization on a hypersphere.\n\tFixed by generating random ANGLES in 2D rather than magnitude along DIR.\n" );
	printf( "Version 75\n\tNow that the rotation is 'fixed' I don't get any nematic ordering\n\tSignificant debugging of genrand_maierSaupeMC() was required.\n\tProblems with random number generation AND WHEN the rotation occured\n" );
	printf( "Version 74\n\tWhen the alignment order parameter is large particles' orientation can diverge, which causes velocity to diverge\n\tProblem was due to orientation update. Now use angular velocity rather than rate of change of orientation\n" );
	printf( "\tAdded many different forms of Jeffery's equation for the alignment --- not sure which is best --- using Scaintillan's\n\tTensor dot product routine was only good for one order. Now have 2\n\tOutput director in avS.dat\n" );
	printf( "\tProblem with temperature in 2D for Andersen-type algorithms \n\t\t--- genrand_gaussScaled() was wrong but happened to be correct in 3D\n\tAlso orientation was always pointing in x-direction\n\t\t--- found Rodrigues' rotation formula requires UNIT vector for axis. Now done in routine\n" );
	printf( "\tAlso noticed could calculate rotation angle and axis for each cell rather than each particle\n" );
	printf( "Version 73\n\tTo couple the orientation and flow, I now calculate the torque from:\n\t\tthe collision operation, the shear and the magnetic field\n\tLCcollision must go first to set the torque\n\tWas failing in 2D because torque was in 3rd dimension but was only doing cross products in 2D\n\t\t--- made 2D and 3D variables\n" );
	printf( "\tAlso in order to get the twist measurement for the Frank coefficient, need in-plane orienational BCs to be different;\n\t\ttherefore, add to BC\n" );
	printf( "Version 72\n\tBack to the LC\n\tUse an MC model to generate the Maier-Saupe distribution \n\t\t--- I don't love this but seems like the best choice\n\tIt may currently use too many 'annealing steps'\n\tCreated andersenROT_LC() for the collision operation.\n\tIt adjusts for the angular momentum due to the rotation of rods and so couples orienation --> velocity\n" );
	printf( "\tIn localPROP() calculate the velocity gradient to use for coupling velocity --> orientation \n\tThe actual coupling is Jeffery's equation jefferysRotation()" );
	printf( "\n\tjefferysRotation() must FOLLOW LCcollision() just cuz LCcollision() sets dU and jefferysRotation() adds to that\n\tCode not working so debug:\n\t\t- orient() in setcoord() was not good in 2D because input was 2pi not pi \n" );
	printf( "\t\t--- was 'a' bug not 'the' bug\n\t\t- there was a problem with the input file printcom.inp and readpc()\n\t\t- found a minor error in calculating the tensor order parameter for LCs\n" );
	printf( "\t\t- Found the error in calculating moment of inertia. Used oppertunity to move calculation of CM and I to subroutines\n" );
	printf( "Version 71\n\tI still haven't got the LC collision operator correct BUT\n\tMoving on to dipolar active matter algorithm\n" );
	printf( "Version 70\n\tCreate a passive nematic liquid crystal MPCD algorithm\n\tThis first version is not coupled to the fluid flow BUT orientation is just carried on each SRD particle\n\tEach particle has an ortientation\n" );
	printf( "\tThe LC collision operator IS NOT CORRECT!!! It doesn't produce normalized vectors. It needs to be improved\n\tTo calcualte scalar order parameter (and director) must calculate tensor order parameter and solve eigensystem\n\tAdd a magnetic susceptibility to the particles and an external magnetic field to the system\n" );
	printf( "\tFix the MPCD-Chate models by outputting S as calculated by solving the eigensystem\n" );
	printf( "Version 69\n\tAdd a Chate model\n\tThis requires knowing the orientation order tensor for each cell\n\tand finding the first eigenvector of that matrix\n\tAlso threw in a recursive determinant algorithm rather than my dumb versions\n\tOutput tensor and scalar order parameters\n\tTHIS DID NOT SEEM TO PRODUCE AN ACTIVE FLUID...\n" );
	printf( "Version 68\n\tChange the directory names and where makefile is etc.\n\tAdd a literal Vicsek model\n" );
	printf( "Version 67\n\tChanged the name of the active Andersen and Langevin routines to reflect that they are the Vicsek versions\n\t\t--- flags={VICSEK_MPCAT,VICSEK_LANG] and routines={vicsekAndersenMPC(),vicsekLangevinMPC()}\n\tAdded the Chate variations on Vicsek model for 'nematic' or bipolar active fluids (in the Vicsek routines)\n" );
	printf( "\tI would like to output coarse-grained data separately for each species.\n\tAlso include the Langevin active MPCD --- active particles relax toward Vicsek speed rather than VCM\n\tThis meant that FRICCO could no longer double as the relaxation parameter so...\n" );
	printf( "Version 66\n\tMy second attemp at an MPCD-active fluid\n\tNow the 'activity' belongs to each particle rather than each cell\n\tI hope that this will encourage large density fluctuations and phase separation between active and passive particle populations\n" );
	printf( "\t\t- It doesn't appear to though\n\tAdded MPCD version of active fluid that relaxes towards a set velocity\n" );
	printf( "Version 65\n\tMy first attempt at an MPCD-active fluid algorithm.\n\tMy first attempt is a modification to SRD that DOES NOT conserve momentum but DOES conserve energy\n\tAlso modified the initialization so that temperature ALWAYS has EXACTLY the initialization value\n" );
	printf( "Version 64\n\tBC_BCcollision_putBack() didn't work.\n\tSo I put a 'putBack' if statement into BC_Bccollision()\n" );
	printf( "Version 63\n\tTHIS NEEDS TO BE TESTED AS SOON AS I GET HOME!!! TEST AGAINST THE FFF. TEST AGAINST THE ONES THAT FAIL.\n\tAdd a memory of last position of BC so can put back where it was.\n\tI'm not sure if this will work well. Kill it and go back to v062 if it doesn't.\n\tAlso added Lennard-Jones Tube to the MD-side of things\n" );
	printf( "Version 62\n\tAttempt to fix the first monomer getting cut off when initialized in a straight line (change made in v060)\n" );
	printf( "Version 61\n\tMartin had a situation which could use the probabilistic reflections/thermal walls\n\tBUT it required that every wall have it's own temperature\n" );
	printf( "Version 60\n\tThe MD portion is modified to put the polymer in a straight line\n" );
	printf( "Version 59\n\tAHH!!! BC problems again!!! What the hell?\n\tThe hack to fix this is now that if a hard sphere still overlaps with another BC after attempting to resolve the collision\n\tan ever stronger force is applied until the sphere moves out of the vicinity of the wall.\n" );
	printf( "Version 58\n\tVery minor changes made to DBG options.\n\tAlso included a compile option to force flush the buffer every timestep\n" );
	printf( "Version 57\n\tWhen there is a flow then the thermal energy in the other directions is surpressed.\n\t\tThis is trouble for the FFF simulations\n\tProblem doesn't exist if Andersen collison is used by itself\n\t\t- add warning if thermostat is on, ontop of Andersen collison\n" );
	printf( "\tAdded a new 'thermostat'\n\t\t--> MAXV which is like the Berendsen thermostat BUT relaxes the fluid toward the velocity given by GRAV[3]\n\tAlso added a replace() in the bin() rather then exitting\n" );
	printf( "Version 56\n\tAdd option to output global average MPCD velocity (will need to change all printcom.inp files)\n\tAlso at this time I modified the MD code so that the polymer could be placed as a rod (x-direction)\n" );
	printf( "\t\teither in the middle or near the top or near the bottom (y-direction): LAYOUT_MID, LAYOUT_TOP, LAYOUT_BOT\n\t\t	Add to mdtypes.h\n\t\t	Add to mdfiles.c\n" );
	printf( "Version 55\n\t2D-FFF simulations keep crashing!!!\n\t1) Tried just doing a floating bead with no MPCD particles - crashed.\n\t\tAdded a check to BC_MPCcollision() for the rare case when a BC doesn't collide with ANY MPCD particles\n\t2) Looked at BCs moving with no particles.\n" );
	printf( "\t\tThings worked until BC collided with another, nonplanar BC\n\t\tFaulty IF() statement fixed\n\t3) Added gain in angular momentum when collide with wall or other BC\n\t4) Put in MPCd particles with no flow but strong gravity downward - still had crashes." );
	printf( "\n\t\t The BC eventually 'broke through' the bottom wall. Once y<radius then the bead stopped moving in x.\n\t\tThe velocity kept fluctuating but the position was stuck in the wall\n\t\tBy checking whether the movingWall was approaching or receding from the stillWall was able to help\n\t\tBUT sometimes not enough\n" );
	printf( "\t5) If a moving BC is still violating a stillBC after the transformation then just let it stream out for as many timesteps as necessary.\n\t\tThis passes all the tests that I throw at it.\n" );
	printf( "Version 54\n\tAdded DVxyz[3] to BCs. In this way, a specific region of space can drive flow. This is needed for David's project.\n" );
	printf( "Version 53\n\tMoving BCs were getting stuck in walls. I finally got it! I was forgetting to to sum VN and VT in BC_BCcollision()\n\tPlus, I added a couple more output notes to synopsis.dat to let me know where in the program, the code gets stuck if it does\n" );
	printf( "Version 52\n\tStill restructuring the BC-BC collisions\n\tBCtrans_collision() now replaced with BC_BCcollision() followed by BC_MPCcollision()\n\t\t- therefore, delete BCtrans_collision() and BConBC()\n\tAlso - editted out resetting i=0 in if(tc < t_min) in chooseP()\n" );
	printf( "Version 51\n\tCorrect to deciding if BCs in contact\n\tRecreated BConBC. Colliding BC always moves AWAY from fixed BC. Moves away from wall too.\n" );
	printf( "Version 50\n\tNot actually a new version of MPC but the video output (MD positions) is altered from scenic-format to VMD-format\n\tAdd harmonic potential for MD - for Martin project\n\tClean SRD arg handling so that doesn't matter if input or output directories end with '/'\n" );
	printf( "Version 49\n\tIn order to get shear flow from a moving plane, we needed to add the wall's velocity to the ghost particles\n\tUSE PHANTOM PARTICLES!\n" );
	printf( "Version 48\n\tIt appears we weren't giving the polymers the proper charge because we weren't always passing integrateMD the SRD particles\n\tActually problem was bigger:\n" );
	printf( "Version 47\n\tshiftBC() is called often and so is expensive even though fast. Changed MAYBE faster now...\n\tGave each BC a PLANE flag. If the BC is a simple plane then it skips some of the power calcs\n\tCreated specialized calcW and crosstime (called calcW_PLANE and crosstime_PLANE) for planes i.e. when PLANE==1\n" );
	printf( "Version 46\n\tProfiled the code and found localPROP was taking lots of time.\n\tRemoved calc of KBT and added  IF-statement so only calc MOMINERT and CM\n" );
	printf( "Version 45\n\tWhen shifting BCs in chooseP() don't shift back every time: keep track of the total shift.\n" );
	printf( "Version 44\n\tI improved the way multiple bounces are handles (important for systems like nanopits (see NBOUNCE))\n\tPut an if-statement around output so that don't call bin() twice unless need to\n\tThere are if-statements everywhere. Let's make them depend on compiled version (in definitions.h)\n" );
	printf( "Version 43\n\tFuck me. A new error for when BCs cross PBCs. Fixed\n\tacc() was truncating time to an integer. Casting error now fixed.\n\tOwen's the BEST he found the symmetry problem with the MPCinMD coupling\n" );
	printf( "Version 42\n\tRemoved redundant bin()\n\tAdd forgotten binMD()\n\tCreated zeroExtraDims() so not in main routine\n\tghostPart() was not adding phantom particles correctly. Fixed it.\n\tIn localPROP() there was an && (instead of an ||) for checking if there were MPC or MD particles\n\t\t---screwed up things like flowout()\n" );
	printf( "\tFixed when GRAV_FLAG was turned on (what was I on?)\n\tMoved recording flow out of localRPOP() so that doesn't have shiftR bias.\n\tIn Andersen algorithm, there were two times (once after MPC particles then again after MD) that RS was divided by pop\n\tIn Andersen algorithm, RV[][] was keeping values from last time in routine. Needed to be zeroed at start.\n" );
	printf( "Version 41\n\tTook STDY_FLAG out of printcom.inp. It was dumb\n\tAdd double the tolerance to the periodic BCs to make sure they are inside the control box\n\tGive separate acceleration to each bc.\n" );
	printf( "Version 40\n\tTook out 'thermostating' of BC-MPC collisions since destroys conservation of momentum\n\tAdds spherical-BC on BC collisions\n\tCleaned up velBC() routine\n" );
	printf( "Version 39\n\tZeros the net momentum periodically if no acceleration\n\tPeriodically zeros pos and vel of unused dimensions\n\tFixes user print options to be more intuitive and useful\n\tAdds Heyes cell thermostating\n" );
	printf( "Version 38\n\tTries to get rid of my BC bug\n\tFINALLY DOES IT!!!\n" );
	printf( "Version 37\n\tThis version debugs the D3D mode\n\tAnd couples to the MD\n" );
	printf( "Version 36\n\tThis version double checks all the gravity - especially in crosstime()\n\tUses secant method to evaluate if particle is on surface\n\tThe above point means it can do spheres in a field and can do higher power shapes\n\tcast time step as a double so that don't have to every time\n\tDIDN't Solve!!!\n" );
	printf( "Version 35\n\tThis version makes some variables global\n\tIt removes the excessive amount of DBUG IF-statements\n\tand comments to myself\n" );
	printf( "Version 34\n\tThis version debugs a sparce error in collision time\n" );
	printf( "Version 33\n\tThis version uses dynamic memory allocation on cell CL\n" );
	printf( "Version 32\n\tCoupled MPC into MD algorithm\n\tNow define new structures (like particleMPC, bc and spec)\n" );
	printf( "Version 31\n\tWorks towards a MPC-MD hybrid\n\tAdds Brownian thermostat and Langevin thermostat\n\tFix angular momentum conservation\n" );
	printf( "Version 30\n\tImplements angular momentum conservation\n" );
	printf( "Version 029\n\tImplements the Andersen thermostat version of MPC\n" );
	printf( "Version 028\n\tImplements phantom particles at walls to reduce slip.\n" );
	printf( "Version 027\n\tDebug anisotropy in diffusion\n" );
	printf( "Version 026\n\tTakes arguments for path to input files and to saving output. All other input still through input files.\n\tThermal BCs up and running\n" );
	printf( "Version 025\n\tThis version uses dynamic memory allocation by the malloc() function\n" );
	printf( "Version 024\n\tDeals with the energy conservation issue when moving beads are present\n" );
	printf( "Version 023\n\tUse a linked list instead of an array for listing the particles in the cells.\n" );
	printf( "Version 022\n\tClean the code, divide into headers to use make file\n\tCreated first_col to find the first collision outside the collision IF\n" );
	printf( "Version 021\n\tThis version includes different (non-energy conserving i.e. need thermostat) boundary conditions\n" );
	printf( "Version 020\n\tEverything seems to be working EXCEPT bounce back with multiple particles and a movable BC object\n\tZeroed the system's initial net momentum by a galilean transformation to rest frame.\n\tThere's still an error when lots of particles are INSIDE a very light shell\n" );
	printf( "Version 019\n\tContinued debugging impulse proceedure. Problem lies in translation and collision of BC objects\n\tCROSSTIME MUST BE MODIFIED BEFORE CAN DO ANYTHING BETTER THAN A SPHERE!!!\n\tAND MOMINERT!!!\n" );
	printf( "Version 018\n\tHawaiian addition. This version cleans up the code.\n\tFirst made sure that there are no more problems when PLACING particles while obstacles are around CHECK\n\tSecond made the debugger topic specific rather than steps in detail\n\tThird made so that multiple BC collisions per timestep were handled properly\n\tFourth put BC routines into functions\n" );
	printf( "\tFifth made many of the elements in structures into vectors\n\tSixth: impulse method used to handle collisions\n\tDID NOT FINISH DEBUGGING - Moved to Version 019\n" );
	printf( "Version 017\n\tThis version reconsiders what was done in v016 and tries to make it cleaner and well work properly\n" );
	printf( "Version 013\n\tMore complex boundaries such as circles and squircles work\n\tObjects defined by BCs can move and transfer momentum\n" );
	printf( "Version 012\n\t2D Mode works\n\tNow checks if the placement of particles are acceptable.\n" );
	printf( "Version 011\n\tFixed the BC. They work for planar surfaces at least.\n" );
	printf( "Version 010\n\tAllows the user to state what data sets should be printed out (in some cases saves much time)\n" );
	printf( "Version 009\n\tCalculates the boundary conditions in a much more universal way\n\tby considering components normal and tangential to the surfaces.\n" );
	printf( "Version 008\n\tUses the Mersenne Twister random number generator and specifies\n\tinitial velocities by temperature.\n\tOutputs velocity distributions.\n\tIncludes velocity scaling and Berendsen thermostats.\n" );
	printf( "Version 007\n\tBy randomly shifting the arrays (actually everything else)\n\tthe code is now galillean invariant.\n\tRemoved LX,LY,LZ since not using them yet (or ever?)\n" );
	printf( "Version 006\n\tThis version uses an array of boundary conditions so that the number\n\tof them is easily variable and the code is prettier. It also adds\n\tterms to the bc structure that determine how a particle behaves if\n\tit crosses a boundary. This is nicer and involves less code than using\n" );
	printf( "\tglobal flags BUT it means that the user must be a little more on top\n\tof things. Still galillean variant.\n\tIt also outputs some macroscopic values.\n" );
	printf( "Version 005\n\tApplies what was learnt from v004 to the SRD code. It is\n\tvery different from v003. It is still gallelian variant!\n" );
	printf( "Version 004\n\tIS NOT A SRD CODE. This version is just practice for\n\tusing/passing pointers instead of large structures.\n\tIt passes pointers to letters to spell things and\n\treorganizes them\n" );
	printf( "Version 003\n\tDeveloping Code. Doesn't get to collisions just gets to binning\n" );
	printf( "Version 001 is the most basic SRD imaginable\n\tOnly supports rectangular BC\n\tNo angular momentum\n\tNo thermostating\n" );
	printf( "Version 000 had two separate lists of all the particles\n\tA list of species and an array of cells.\n\tIn this version we will try to free up some memory by only one\n" );
	printf( "\tI don't want to use species because I would have to search the\n\twhole list once for every cell whereas if I use cells I will\n\thave to pass particles from cell to cell.\n" );
}

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* **************** PRINTING **************** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
void outheader( FILE *fout,int SP ) {
/* Simple header output for output files */
	fprintf( fout," **********************************************\n" );
	fprintf( fout," ******** Stochastic Rotation Dynamics ********\n" );
	fprintf( fout," ************* By Tyler Shendruk **************\n\n" );
	fprintf( fout," **********************************************\n" );
	fprintf( fout," *** For the Polymer Physics Research Group ***\n" );
	fprintf( fout," ******** At the University of Ottawa *********\n\n" );
	fprintf( fout," **********************************************\n" );
	fprintf( fout," ************ Code begun Oct. 2008 ************\n" );
	fprintf( fout," **********************************************\n\n" );
}
void coordheader( FILE *fout ) {
/* Simple header for output columns */
	fprintf( fout,"t\t QX\t\t QY\t\t QZ\t\tVX\t\tVY\t\tVZ\t\t|V|\t\tUX\t\tUY\t\tUZ\n" );
}
void coarseheader( FILE *fout ) {
/* Simple header for output columns */
	int n;
	fprintf( fout,"t\t\t\t\tQX\t\tQY\t\tQZ\tVcmX\t\t\tVcmY\t\t\tVcmZ\t\t\t\tPOP" );
	for( n=0; n<NSPECI; n++ ) fprintf( fout,"\t\tSP%d",n );
	fprintf( fout,"\n" );
}
void orderheader( FILE *fout ) {
/* Simple header for output columns */
	fprintf( fout,"t\t\tQX\t\tQY\t\tQZ\t\tNX\t\tNY\t\tNZ\t\tS\n" );
}
void orderQheader( FILE *fout ) {
/* Simple header for output columns */
	fprintf( fout,"X\tY\tZ\tQXX\tQXY\tQXZ\tQYX\tQYY\tQYZ\tQZX\tQZY\tQZZ\n" );
}
void orderQKheader( FILE *fout ) {
/* Simple header for output columns */
	fprintf( fout,"t\tK123_X\ttK123_Y\ttK123_Z\t|QXX|2\t|QXY|2\t|QXZ|2\t|QYX|2\t|QYY|2\t|QYZ|2\t|QZX|2\t|QZY|2\t|QZZ|2\n" );
}
void avvelheader( FILE *fout ) {
/* Simple header for output columns */
	fprintf( fout,"t\t VcmX\t\tVcmY\t\tVcmZ\t\tKBT\n" );
}
void avvelWithGradVelheader( FILE *fout ) {
/* Simple header for output columns */
	fprintf( fout,"t\t VcmX\t\tVcmY\t\tVcmZ\t\tKBT\t\tdVXX\t\tdVXY\t\tdVXZ\t\tdVYX\t\tdVYY\t\tdVYZ\t\tdVZX\t\tdVZY\t\tdVZZ\t\n" );
}
void corrheader( FILE *fout ) {
/* Simple header for output columns */
	fprintf( fout,"t\tdr\t C\n" );
}
void energyspectheader( FILE *fout ) {
/* Simple header for output columns */
	fprintf( fout,"t\t k\t\t E\n" );
}
void enstrophyspectheader( FILE *fout ) {
/* Simple header for output columns */
	fprintf( fout,"t\t k\t\t Omega\n" );
}
void defectheader( FILE *fout ) {
/* Simple header for output columns */
	fprintf( fout,"t\t QX\t\t QY\t\t QZ\t\t charge\t\t angle\n" );
}
void multiphaseheader( FILE *fout ) {
/* Simple header for output columns */
	int i;
	fprintf( fout,"t\t\tQX\t\tQY\t\tQZ" );
	for( i=0; i<NSPECI; i++ ) fprintf( fout,"\t\tN_%d",i );
	fprintf( fout,"\n" );
}
void pressureheader( FILE *fout ) {
/* Simple header for output columns */
	fprintf( fout,"t\t\tQX\t\tQY\t\tQZ\tPxx\tPxy\tPxz\tPyx\tPyy\tPyz\tPzx\tPzy\tPzz\n" );
}
void binderheader( FILE *fout,int binSize ) {
/* Simple header for output columns */
	fprintf( fout,"Bin Size:\t%d\n",binSize );
	fprintf( fout,"t\tBinderCumulant\n" );
}
void avsheader( FILE *fout ) {
/* Simple header for output columns */
	fprintf( fout,"t\t\t S\t\t S4\t\t nX\t\t nY\t\t nZ\n" );
}
void densheader( FILE *fout ) {
/* Simple header for output columns */
	fprintf( fout,"t\t\t densSTD\n" );
}
void avenstrophyheader( FILE *fout ) {
/* Simple header for output columns */
	fprintf( fout,"t\t\t enstrophy\n" );
}
void flowheader( FILE *fout ) {
/* Simple header for output columns */
	fprintf( fout,"   QX\t   QY\t   QZ\tVcmX\t\tVcmY\t\tVcmZ\n" );
}
void solidsheader( FILE *fout ) {
/* Simple header for output columns */
	fprintf( fout,"t \tQX\t\tQY\t\tQZ\t\tVX\t\tVY\t\tVZ\t\tOX\t\tOY\t\tOZ\t\tLX\t\tLY\t\tLZ\n" );
}
void histVelheader( FILE *fout ) {
/* Simple header for output columns */
	fprintf( fout,"t\t V\t\t\tPX\t\tPY\t\tPZ\n" );
}
void histVortheader( FILE *fout ) {
/* Simple header for output columns */
	fprintf( fout,"t\t w\t\t\tPX\t\tPY\t\tPZ\n" );
}
void histDirheader( FILE *fout ) {
/* Simple header for output columns */
	fprintf( fout,"t\t n\t\t\tPX\t\tPY\t\tPZ\n" );
}
void histSpeedheader( FILE *fout ) {
/* Simple header for output columns */
	fprintf( fout,"t\t |V|\t\tP\n" );
}
void histEnstrheader( FILE *fout ) {
/* Simple header for output columns */
	fprintf( fout,"t\t |w|\t\tP\n" );
}
void histNheader( FILE *fout ) {
/* Simple header for output columns */
	fprintf( fout,"t\t stdN\t\tP\n" );
}
void histSheader( FILE *fout ) {
/* Simple header for output columns */
	fprintf( fout,"t\t S\t\tP\n" );
}
void energyheader( FILE *fout ) {
	/* Simple header for output columns */
	fprintf( fout,"t\t\tMPC_kin\t\tMPC_nem\t\tBC_kin\t\tBC_rot\t\tTotal\t\tKBT\n" );
}
void energyfieldheader( FILE *fout ) {
	/* Simple header for output columns */
	fprintf( fout,"QX\tQY\tQZ\tMPC_kin\ttMPC_nem\n" );
}
void swimmerheader( FILE *fout ) {
	/* Simple header for output columns */
	fprintf( fout,"t\t\tHX\tHY\tHZ\tHVX\tHVY\tHVZ\tMX\tMY\tMZ\tMVX\tMVY\tMVZ\tRTphase\n" );
}
void swimmeroriheader( FILE *fout ) {
	/* Simple header for output columns */
	fprintf( fout,"t\t\tnX\tnY\tnZ\tRTphase\n" );
}
void runtumbleheader( FILE *fout ) {
	/* Simple header for output columns */
	fprintf( fout,"RTphase\t\tdt_cnt\tdAng\n" );
}
void energyneighboursheader( FILE *fout ) {
	/* Simple header for output columns */
	fprintf( fout,"t\ttMPC_nem\n" );
}
void coordout( FILE *fout[MAXSPECI],int pr,double T,particleMPC p[],spec SP[] ) {
/* Prints particleMPC coordinate data to output file */
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
void coarseout( FILE *fout,double t,cell ***CL ) {
	/*
	 Print coarse data to file
	 */
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
void cellout( cell ***CL ) {
/*
   This function states the coordinates of each
   particleMPC in the array of lists. It states the
   cell (and population of the cell) that the
   particleMPC resides in.
*/
	int i,j,k,l;
	particleMPC *pMPC;	//Temporary pointer to MPC particles
	particleMD *pMD;		//Temporary pointer to MD particles
	smono *pSW;					//Temporary pointer to swimmer monomers

	for( i=0; i<XYZ_P1[0]; i++ ) for( j=0; j<XYZ_P1[1]; j++ ) for( k=0; k<XYZ_P1[2]; k++ ) {
		l=0;
		// MPC
		pMPC = CL[i][j][k].pp;
		pMD = CL[i][j][k].MDpp;
		pSW = CL[i][j][k].sp;

		if( pMPC != NULL || pMD != NULL || pSW != NULL ) printf( "In Cell [%d,%d,%d]:\n",i,j,k );
		while( pMPC != NULL ) {
			l++;
			printf( "\tMPC Particle %d:",l );
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
void listout( cell ***CL,int XYZ_P1[_3D] ) {
/*
    Output entire list
*/
	int a,b,c,d;
	particleMPC *pMPC;
	particleMD *pMD;
	smono *pSW;

	printf( "Local properties:\n" );
	for( a=0; a<XYZ_P1[0]; a++ ) for( b=0; b<XYZ_P1[1]; b++ ) for( c=0; c<XYZ_P1[2]; c++ ) {
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
void pcoord( particleMPC p ) {
/*
    Prints the particleMPC's coordinates (position
    and velocity) to the terminal
*/
	printf( "\tQ=(%6.12e,%6.12e,%6.12e)\n",p.Q[0],p.Q[1],p.Q[2] );
	printf( "\tV=(%6.12e,%6.12e,%6.12e)\n",p.V[0],p.V[1],p.V[2] );
	printf( "\tU=(%6.12e,%6.12e,%6.12e)\n",p.U[0],p.U[1],p.U[2] );
}
void bccoord( bc WALL ) {
/*
    Prints the BC's coordinates (position, velocity
    and angular velocity) to the terminal
*/
	printf( "\tQ=(%6.12e,%6.12e,%6.12e)\n",WALL.Q[0],WALL.Q[1],WALL.Q[2] );
	printf( "\tV=(%6.12e,%6.12e,%6.12e)\n",WALL.V[0],WALL.V[1],WALL.V[2] );
	printf( "\tO=(%6.12e,%6.12e,%6.12e)\n",WALL.O[0],WALL.O[1],WALL.O[2] );
	printf( "\tL=(%e,%e,%e)\n",WALL.L[0],WALL.L[1],WALL.L[2] );
}
void mdcoord( particleMD p ) {
/*
    Prints the particleMPC's coordinates (position
    and velocity) to the terminal
*/
	printf( "\tQ=(%lf,%lf,%lf)\n",p.rx,p.ry,p.rz );
	printf( "\tV=(%lf,%lf,%lf)\n",p.vx,p.vy,p.vz );
}
void swcoord( swimmer sw ) {
/*
    Prints the swimmers's coordinates (position
    and velocity of head and middle) to the terminal
*/
	printf( "\tH Q=(%lf,%lf,%lf) ",sw.H.Q[0],sw.H.Q[1],sw.H.Q[2] );
	printf( "\tV=(%lf,%lf,%lf)\n",sw.H.V[0],sw.H.V[1],sw.H.V[2] );
	printf( "\tM Q=(%lf,%lf,%lf) ",sw.M.Q[0],sw.M.Q[1],sw.M.Q[2] );
	printf( "\tV=(%lf,%lf,%lf)\n",sw.M.V[0],sw.M.V[1],sw.M.V[2] );
}
void pvec( double VEC[],int dimension ) {
/*
    Prints any vector
*/
	int i;
	printf( " (" );
	for( i=0; i<(dimension-1); i++ ) printf( "%lf,",VEC[i] );
	printf( "%lf)\n",VEC[dimension-1] );
}
void ptens3D( double TENS[][_3D] ) {
/*
    Prints any vector
*/
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
void ptens2D( double TENS[][_2D] ) {
/*
    Prints any vector
*/
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
void pvcoord( double POS[_3D],double VEL[_3D],double ANG[_3D],int dimension ) {
/*
    Prints the particleMPC's coordinates (position
    and velocity) to the terminal
*/
	printf( "\tQ=" );
	pvec( POS,dimension );
	printf( "\tV=" );
	pvec( VEL,dimension );
	printf( "\tW=" );
	pvec( ANG,dimension );

}
void pall( particleMPC p[] ) {
	/* Prints all particleMPC coordinate data to terminal */
	int i;
	for( i=0; i<GPOP; i++ ) {
		printf( "Particle %i\t",i );
		pcoord( p[i] );
	}
}
void listinput( inputList in,double AVVEL,spec SP[],kinTheory theory ) {
/*
    Lists the input data
*/
	int i,n;
	#ifdef DBG
		if( DBUG >= DBGINIT ){
			printf( "\tDebug Mode: %i\n",DBUG );
			printf( "\tDimensions: %i\n",DIM );
			printf( "\tRotation Technique: %i\n",in.RTECH );
			printf( "\tLiquid Crystal: %i\n",in.LC );
			printf( "\tLiquid Crystal Mean-Field Potential: %lf\n",in.MFPOT );
			printf( "\tRotation angle: %lf\n",in.RA );
			printf( "\tSystem Size: [%i,%i,%i]\n",XYZ[0],XYZ[1],XYZ[2] );
			printf( "\tConstant external acceleration:" );
			pvec( in.GRAV,DIM );
			printf( "\tConstant external magnetic field:" );
			pvec( in.MAG,DIM );
			printf( "\tAverage speed in any direction: %lf\n",AVVEL );
			printf( "\tNumber of Species: %d\n",NSPECI );
			printf( "\tSpecies population: " );
			printf( "%i",SP[0].POP );
			for( i=1; i<NSPECI; i++ ) printf( ",%i",SP[i].POP );
			printf( "\n\tParticle Number Density: %lf\n",nDNST );
			printf( "\n\tMass Density: %lf\n",mDNST );
			printf( "\tWarmup Iterations: %d\n",in.warmupSteps );
			printf( "\tSimulation Iterations: %d\n",in.simSteps );
			printf( "\tTime Step: %lf\n",in.dt );
			printf( "\tSystem population: %i\n",GPOP );
			printf( "\tTotal mass: %lf\n",theory.sumM );
			printf( "\tThe number of boundaries: %i\n",NBC );
			printf( "\tInputted System Temperature (units of KB): %lf\n",in.KBT );
			printf( "\tRemove system's net momentum (1=yes, 0=no): %i\n",in.RFRAME );
			printf( "\tThermostat Method: %i\n",in.TSTECH );
			printf( "\tThermal Relaxation Scale: %lf\n",in.TAU );
			printf( "\tMean Free Path: %lf\n",theory.MFP );
			printf( "\tKinematic Viscosity: %lf\n",theory.VISC );
			printf( "\tSelf Diffusion Coefficient: %lf\n",theory.SDIFF );
			printf( "\tSchmidt number: %lf\n",theory.VISC/theory.SDIFF );
			printf( "\tSpeed of sound: %lf\n",theory.SPEEDOFSOUND );
			printf( "\tThermal Diffusion Coefficient: %lf\n",theory.THERMD );
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
void stateinput( inputList in,spec SP[],bc WALL[],specSwimmer SS,outputFlagsList out,kinTheory theory,FILE *fsynopsis ) {
/*
    Lists input data
*/
	int i;

	if( out.SYNOUT == OUT ) {
		fprintf( fsynopsis,"\nUser defined variables:\n" );
		fprintf( fsynopsis,"Dimensionality: %i\n",DIM );
		fprintf( fsynopsis,"System dimensions: (%i,%i,%i)\n",XYZ[0],XYZ[1],XYZ[2] );
		fprintf( fsynopsis,"Rotation technique: %i\n",in.RTECH );
		fprintf( fsynopsis,"Nematic Liquid Crystal: ");
		if(in.LC) fprintf( fsynopsis,"YES\n" );
		else fprintf( fsynopsis,"NO\n" );
		fprintf( fsynopsis,"Thermal energy: %lf\n",in.KBT );
		fprintf( fsynopsis,"Remove system's net momentum (1=yes, 0=no): %i\n",in.RFRAME );
		fprintf( fsynopsis,"Randomly shift the MPC cells for Galilean invariance (1=yes, 0=no): %i\n",in.GALINV );
		fprintf( fsynopsis,"Thermostating method: %i\n",in.TSTECH );
		fprintf( fsynopsis,"Thermal relaxation time scale: %lf\n",in.TAU );
		fprintf( fsynopsis,"Rotation angle: %lf\n",in.RA );
		fprintf( fsynopsis,"Langevin friction coefficient: %lf\n",in.FRICCO );
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
			fprintf( fsynopsis,"\tRotational Friction Coefficient: %lf\n",SP[i].RFC);
			fprintf( fsynopsis,"\tEffective rod-length to couple MPC torque to BC force: %lf\n",SP[i].LEN);
			fprintf( fsynopsis,"\tTumbling parameter: %lf\n",SP[i].TUMBLE);
			fprintf( fsynopsis,"\tMagnetic Sysceptibility: %lf\n\tShear Sysceptibility: %lf\n",SP[i].CHIA,SP[i].CHIHI );
			fprintf( fsynopsis,"\tActivity: %lf\n",SP[i].ACT );
			fprintf( fsynopsis,"\tDamping friction: %lf\n",SP[i].DAMP );
			fprintf( fsynopsis,"\tPos. dist.: %i\n\tVel. dist.: %i\n\tOri. dist.: %i\n",SP[i].QDIST,SP[i].VDIST,SP[i].ODIST );
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
		}
		fprintf( fsynopsis,"\nSwimmer variables:\n" );
		fprintf( fsynopsis,"\tTyper: %d\n",SS.TYPE );
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
		fprintf( fsynopsis,"Print averaged flow data: %i\n",out.AVVELOUT );
		fprintf( fsynopsis,"Print energy data: %i\n",out.ENOUT );
		fprintf( fsynopsis,"Print director and scalar order parameter fields: %i\n",out.ORDEROUT );
		fprintf( fsynopsis,"Print tensor order parameter data: %i\n",out.QTENSOUT );
		fprintf( fsynopsis,"Print reciprocal space order parameter data: %i\n",out.QKOUT );
		fprintf( fsynopsis,"Print averaged order parameter data: %i\n",out.AVSOUT );
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
		fprintf( fsynopsis,"\tVel: %d\n\tSpeed: %d\n\tVorticity: %d\\n\tEnstrophy: %d\n\tDirector: %d\n\tScalar order parameter: %d\n\tDensity: %d\n",out.HISTVELOUT,out.HISTSPEEDOUT,out.HISTVORTOUT,out.HISTENSTROUT,out.HISTDIROUT,out.HISTSOUT,out.HISTNOUT );
		fprintf( fsynopsis,"Synopsis of Simulation: %i\n",out.SYNOUT );
	}
}
void histVelout( FILE *fout,int vel[_3D][BINS],double minRange,double maxRange,double t ) {
/*
    Print velocity distribution data to file
*/
	int i;
	double dv = (maxRange-minRange)/((float)BINS - 1.0);

	fprintf( fout,"%12.5e\n",t );
	for( i=0; i<BINS; i++ ) fprintf( fout,"\t%12.5e\t%7i\t%7i\t%7i\n",minRange+i*dv,vel[0][i],vel[1][i],vel[2][i] );

	fprintf( fout,"\n" );
	#ifdef FFLSH
		fflush(fout);
	#endif
}
void histSpeedout( FILE *fout,int speed[BINS],double minRange,double maxRange,double t ) {
/*
    Print speed distribution data to file
*/
	int i;
	double dv = (maxRange-minRange)/((float)BINS - 1.0);

	fprintf( fout,"%12.5e\n",t );
	for( i=0; i<BINS; i++ ) fprintf( fout,"\t%12.5e\t%7i\n",minRange+i*dv,speed[i] );

	fprintf( fout,"\n" );
	#ifdef FFLSH
		fflush(fout);
	#endif
}
void histVortout( FILE *fout,int vort[_3D][BINS],double minRange,double maxRange,double t ) {
/*
    Print vorticity distribution data to file
*/
	int i;
	double dw = (maxRange-minRange)/((float)BINS - 1.0);

	fprintf( fout,"%12.5e\n",t );
	for( i=0; i<BINS; i++ ) fprintf( fout,"\t%12.5e\t%7i\t%7i\t%7i\n",minRange+i*dw,vort[0][i],vort[1][i],vort[2][i] );

	fprintf( fout,"\n" );
	#ifdef FFLSH
		fflush(fout);
	#endif
}
void histEnstrout( FILE *fout,int enstrophy[BINS],double minRange,double maxRange,double t ) {
/*
    Print enstrophy distribution data to file
*/
	int i;
	double dw2 = (maxRange-minRange)/((float)BINS - 1.0);

	fprintf( fout,"%12.5e\n",t );
	for( i=0; i<BINS; i++ ) fprintf( fout,"\t%12.5e\t%7i\n",minRange+i*dw2,enstrophy[i] );

	fprintf( fout,"\n" );
	#ifdef FFLSH
		fflush(fout);
	#endif
}
void histDirout( FILE *fout,int dir[_3D][BINS],double minRange,double maxRange,double t ) {
/*
    Print vorticity distribution data to file
*/
	int i;
	double dn = (maxRange-minRange)/((float)BINS - 1.0);

	fprintf( fout,"%12.5e\n",t );
	for( i=0; i<BINS; i++ ) fprintf( fout,"\t%12.5e\t%7i\t%7i\t%7i\n",minRange+i*dn,dir[0][i],dir[1][i],dir[2][i] );

	fprintf( fout,"\n" );
	#ifdef FFLSH
		fflush(fout);
	#endif
}
void histSout( FILE *fout,int S[BINS],double minRange,double maxRange,double t ) {
/*
    Print enstrophy distribution data to file
*/
	int i;
	double dS = (maxRange-minRange)/((float)BINS - 1.0);

	fprintf( fout,"%12.5e\n",t );
	for( i=0; i<BINS; i++ ) fprintf( fout,"\t%12.5e\t%7i\n",minRange+i*dS,S[i] );

	fprintf( fout,"\n" );
	#ifdef FFLSH
		fflush(fout);
	#endif
}
void histNout( FILE *fout,int dens[BINS],double minRange,double maxRange,double t ) {
/*
    Print enstrophy distribution data to file
*/
	int i;
	double dp = (maxRange-minRange)/((float)BINS - 1.0);

	fprintf( fout,"%12.5e\n",t );
	for( i=0; i<BINS; i++ ) fprintf( fout,"\t%12.5e\t%7i\n",minRange+i*dp,dens[i] );

	fprintf( fout,"\n" );
	#ifdef FFLSH
		fflush(fout);
	#endif
}
void enout( FILE *fout,particleMPC *pp,spec *pSP,bc WALL[],double t,double KBT,double wmf ) {
/*
    Print distribution data to file
*/
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
void enfieldout( FILE *fout,cell ***CL,spec *SP,double MFPOT,int LC ) {
/*
    Calculate and print the average orientation interaction energy with neighbouring cells
*/
	int a,b,c,d,id;
	double enK,wmf,S,un,DIR[_3D],u[_3D],m;
	double invdim=1./((double)DIM);
	particleMPC *tmpc;	//Temporary pointer to MPC particles

	for( a=0; a<XYZ[0]; a++ ) for( b=0; b<XYZ[1]; b++ ) for( c=0; c<XYZ[2]; c++ ) {
		wmf=0.;
		enK=0.;
		if( CL[a][b][c].POP > 1 ) {
			S = CL[a][b][c].S;
			for( d=0; d<DIM; d++ ) DIR[d] = CL[a][b][c].DIR[d];
			tmpc = CL[a][b][c].pp;
			while( tmpc != NULL ) {
				id = tmpc->SPID;
				//Nematic energy
				if( LC ) {
					for( d=0; d<DIM; d++ ) u[d] = tmpc->U[d];
					un = dotprod( u,DIR,DIM );
					wmf += S*un*un  + (1.-S)*invdim;
				}
				//Kinetic energy
				m = (SP+id)->MASS;
				for( d=0; d<DIM; d++ ) u[d] = tmpc->V[d];
				un = dotprod( u,u,DIM );
				enK += 0.5*m*un;

				tmpc = tmpc->next;
			}
		}
		wmf*=MFPOT;
		fprintf( fout, "%5i\t%5i\t%5i\t%e\t%e\n",a,b,c,enK,wmf );
	}
}
void enneighboursout( FILE *fout,double t,cell ***CL,double MFPOT,int LC ) {
/*
    Calculate and print the average orientation interaction energy between neighbouring cells.
    To do this, calculate the total director of all the cells under consideration and
    calculate the energy based on the local director and this value.
*/
	int a,b,c,d;
	double wmf,un,sumWMF;
	double local_DIR[DIM],nnn_DIR[DIM],local_S,nnn_S;
	double **Q,eigval[DIM];
	//double invDIM=1.0/((double)DIM);

	sumWMF=0.;
	//Allocate memory for tensor order parameter Q
	Q = malloc ( DIM * sizeof( *Q ) );
	for( a=0; a<DIM; a++ ) Q[a] = malloc ( DIM * sizeof( *Q[a] ) );
	for( a=0; a<DIM; a++ ) for( b=0; b<DIM; b++ ) Q[a][b] = 0.0;

	if( LC ) for( a=0; a<XYZ[0]; a++ ) for( b=0; b<XYZ[1]; b++ ) for( c=0; c<XYZ[2]; c++ ) if( CL[a][b][c].POP>1 ) {
		//Local values
		for( d=0; d<DIM; d++ ) local_DIR[d]=CL[a][b][c].DIR[d];
		local_S=CL[a][b][c].S;
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
		wmf*=MFPOT;
		sumWMF+=wmf;
	}
	fprintf( fout, "%12.5e\t%12.5e\n",t,sumWMF );

	for( d=0; d<DIM; d++ ) free( Q[d] );
	free( Q );
}
void avvelout( FILE *fout,double t,double vel[_3D],double KBT ) {
/*
    Print coarse data to file
*/
	fprintf( fout, "%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\n",t,vel[0],vel[1],vel[2],KBT );
	#ifdef FFLSH
		fflush(fout);
	#endif
}
void avveloutWithGradVel( FILE *fout,double t,double vel[_3D],double KBT,double gradVel[_3D][_3D] ) {
/*
    Print coarse data to file
*/
	fprintf( fout, "%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t",t,vel[0],vel[1],vel[2],KBT );
	fprintf( fout, "%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\n",gradVel[0][0],gradVel[0][1],gradVel[0][2],gradVel[1][0],gradVel[1][1],gradVel[1][2],gradVel[2][0],gradVel[2][1],gradVel[2][2] );
	#ifdef FFLSH
		fflush(fout);
	#endif
}
void avsout( FILE *fout,double t,double S,double S4,double DIR[] ) {
/*
    Print coarse data to file
*/
	fprintf( fout, "%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\n",t,S,S4,DIR[0],DIR[1],DIR[2] );
	#ifdef FFLSH
		fflush(fout);
	#endif
}
void densSTDout( FILE *fout,double t,double stdN ) {
/*
    Print coarse data to file
*/
	fprintf( fout, "%12.5e\t%12.5e\n",t,stdN );
	#ifdef FFLSH
		fflush(fout);
	#endif
}
void avenstrophyout( FILE *fout,double t,double E ) {
/*
    Print coarse data to file
*/
	fprintf( fout, "%12.5e\t%12.5e\n",t,E );
	#ifdef FFLSH
		fflush(fout);
	#endif
}
void binderout( FILE *fout,double t,double UL ) {
/*
    Print coarse data to file
*/
	fprintf( fout, "%12.5e\t%12.5e\n",t,UL );
	#ifdef FFLSH
		fflush(fout);
	#endif
}
void flowout( FILE *fout,cell ***CL,int interval ) {
/*
    Turns sum of cells' vcm into average, prints average to
    file, zeros sums to start anew
*/
	int h,i,j,k;
	double av[_3D];
	// for( i=0; i<_3D; i++ ) av[i] = 0.0;
	double dint = (double)interval;

	for( i=0; i<XYZ[0]; i++ ) for( j=0; j<XYZ[1]; j++ ) for( k=0; k<XYZ[2]; k++ ) {
	// for( i=0; i<XYZ_P1[0]; i++ ) for( j=0; j<XYZ_P1[1]; j++ ) for( k=0; k<XYZ_P1[2]; k++ ) {
		for( h=0; h<DIM; h++ ) av[h] = CL[i][j][k].FLOW[h]/dint;		//Normalize the sum to get the average
		fprintf( fout, "%5d\t%5d\t%5d\t",i,j,k );
		fprintf( fout, "%12.5e\t%12.5e\t%12.5e\n",av[0],av[1],av[2] );
		for( h=0; h<DIM; h++ ) CL[i][j][k].FLOW[h] = 0.0;		//Reset sum
	}
	#ifdef FFLSH
		fflush(fout);
	#endif
}
void velout( FILE *fout,cell ***CL, double t) {
/*
    Turns sum of cells' vcm into average, prints average to
    file, zeros sums to start anew
*/
    int h,i,j,k;
    double vel[_3D];
    // for( i=0; i<_3D; i++ ) av[i] = 0.0;

    for( i=0; i<XYZ[0]; i++ ) for( j=0; j<XYZ[1]; j++ ) for( k=0; k<XYZ[2]; k++ ) {
                // for( i=0; i<XYZ_P1[0]; i++ ) for( j=0; j<XYZ_P1[1]; j++ ) for( k=0; k<XYZ_P1[2]; k++ ) {
                for( h=0; h<DIM; h++ ) vel[h] = CL[i][j][k].VCM[h];
                fprintf( fout,"%12.5e\t", t); // print time
                fprintf( fout, "%5d\t%5d\t%5d\t",i,j,k );
                fprintf( fout, "%12.5e\t%12.5e\t%12.5e\n",vel[0],vel[1],vel[2] );
                for( h=0; h<DIM; h++ ) CL[i][j][k].FLOW[h] = 0.0;		//Reset sum
            }
#ifdef FFLSH
    fflush(fout);
#endif
}
void solidout( FILE *fout,bc WALL,double t ) {
/*
    Print solids data to file
*/
	fprintf( fout,"%12.5e\t",t );
	fprintf( fout,"%12.5e\t%12.5e\t%12.5e\t",WALL.Q[0],WALL.Q[1],WALL.Q[2] );
	fprintf( fout,"%12.5e\t%12.5e\t%12.5e\t",WALL.V[0],WALL.V[1],WALL.V[2] );
	fprintf( fout,"%12.5e\t%12.5e\t%12.5e\t",WALL.O[0],WALL.O[1],WALL.O[2] );
	fprintf( fout,"%12.5e\t%12.5e\t%12.5e\n",WALL.L[0],WALL.L[1],WALL.L[2] );
	#ifdef FFLSH
		fflush(fout);
	#endif
}
void topochargeout( FILE *fout,double t,cell ***CL ) {
	/*
	 Print topological charge data to file
	 */
	//FIXME: only designed to work for 2D!
	int i,j,k;

	double topoC[XYZ[0]][XYZ[1]]; //init topo charge array
	for( i=0; i<XYZ[0]; i++ ) for( j=0; j<XYZ[1]; j++ ) topoC[i][j] = .0;
	double topoAngle[XYZ[0]][XYZ[1]]; //init topo angle array
	for( i=0; i<XYZ[0]; i++ ) for( j=0; j<XYZ[1]; j++ ) topoAngle[i][j] = .0;

	//loop through non-CB boundary cells and calculate topo charge
	for( i=1; i<XYZ[0]-1; i++ ) for( j=1; j<XYZ[1]-1; j++ ) topoC[i][j] = topoChargeLocal(CL, i, j, 0); 
	//loop through non-CB boundary cells and calculate topo angle
	for( i=2; i<XYZ[0]-2; i++ ) for( j=2; j<XYZ[1]-2; j++ ){
		///FIXME: Too lazy to handle derivatives properly at the boundaries, so we just ignoring another layer there instead. Oopsies. Same goes for the loop above. 
		if (fabs(topoC[i][j]) > TOL) topoAngle[i][j] = topoAngleLocal(CL, i, j, 0, topoC[i][j]); 
	} 
	
	for( i=0; i<XYZ[0]; i++ ) for( j=0; j<XYZ[1]; j++ ) for( k=0; k<XYZ[2]; k++ ) {
		//Output
		fprintf( fout,"%.2f\t",t );
		fprintf( fout,"%5d\t%5d\t%5d\t",i,j,k );
		if( CL[i][j][k].POP == 0 ) fprintf( fout, "%06.3f\t%12.5e\n", 0.0, 0.0);
		else fprintf( fout, "%06.3f\t%12.5e\n",topoC[i][j], topoAngle[i][j]);
	}
	#ifdef FFLSH
		fflush(fout);
	#endif
}
void orderout( FILE *fout,double t,cell ***CL,int LC ) {
	/*
	 Print scalar order parameter and director field data to file
	 */
	int i,j,k;
	for( i=0; i<XYZ[0]; i++ ) for( j=0; j<XYZ[1]; j++ ) for( k=0; k<XYZ[2]; k++ ) {
		//Output
		fprintf( fout,"%.2f\t",t );
		fprintf( fout,"%5d\t%5d\t%5d\t",i,j,k );
		if( CL[i][j][k].POP == 0 ) fprintf( fout, "%12.5e\t%12.5e\t%12.5e\t%12.5e\n",0.0,0.0,0.0,0.0 );
		else fprintf( fout, "%12.5e\t%12.5e\t%12.5e\t%12.5e\n",CL[i][j][k].DIR[0],CL[i][j][k].DIR[1],CL[i][j][k].DIR[2],CL[i][j][k].S );
	}
	#ifdef FFLSH
		fflush(fout);
	#endif
}
void multiphaseout( FILE *fout,double t,cell ***CL ) {
	/*
	 Print phi/colour/species-type field data to file
	 */
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
void pressureout( FILE *fout,double t,cell ***CL ) {
	/*
	 Calculate and print pressure field data to file
	 */
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
void orderQout( FILE *fout,double t,cell ***CL,int LC ) {
	/*
	 Print scalar order parameter data to file
	 */
	int i,j,k;
	double **Q;

	//Allocate memory for tensor order parameter Q
	Q = malloc ( _3D * sizeof( *Q ) );
	for( i=0; i<_3D; i++ ) Q[i] = malloc ( _3D * sizeof( *Q[i] ) );
	for( i=0; i<_3D; i++ ) for( j=0; j<_3D; j++ ) Q[i][j] = 0.0;

	for( i=0; i<XYZ[0]; i++ ) for( j=0; j<XYZ[1]; j++ ) for( k=0; k<XYZ[2]; k++ ) {
		//// Find the tensor order parameter based on self, neighbours and next-nearest neighbours (NNN)
		//tensOrderParamNNN( CL,Q,LC,i,j,k );
		// Find the local tensor order parameter
		tensOrderParam( &CL[i][j][k],Q,LC );
		//Output
		fprintf( fout,"%5d\t%5d\t%5d\t",i,j,k );
		if( CL[i][j][k].POP == 0 ) fprintf( fout, "%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\n" ,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 );
		else fprintf( fout, "%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\n", Q[0][0],Q[0][1],Q[0][2],Q[1][0],Q[1][1],Q[1][2],Q[2][0],Q[2][1],Q[2][2] );
	}

	for( i=0; i<_3D; i++ ) free( Q[i] );
	free( Q );
	#ifdef FFLSH
		fflush(fout);
	#endif
}
void orderQKout( FILE *fout,double t,particleMPC pMPC[],cell ***CL,int LC ) {
	/*
	 Print scalar order parameter data to file
	 */
	int i,j,k,a,b,n;
	double S,ReQ[_3D][_3D],ImQ[_3D][_3D],temp_ReQ[_3D][_3D],temp_ImQ[_3D][_3D],modQ2[_3D][_3D];
	double rotMat[_3D][_3D],rotMatTranspose[_3D][_3D];
	double invL[_3D],DIR[_3D],xhat[_3D],zhat[_3D];
	double waveNum[_3D],K123[_3D],Kprime[_3D],k3;
	double U[_3D],pos[_3D],kr,ckr,skr;
	double c1=1./((double)DIM-1.);
	double c2=((double)XYZ[0]*XYZ[1]*XYZ[2])/((double)GPOP);
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
void corrout( FILE *fout,double corr[],double t ) {
	/*
	 Print scalar order parameter data to file
	 */
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
void spectout( FILE *fout,double spect[],double t ) {
	/*
	 Print scalar order parameter data to file
	 */
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
void checkpoint( FILE *fout,inputList in,spec *SP,particleMPC *pSRD,int MDmode,bc *WALL,outputFlagsList outFlag,int runtime,int warmtime,double AVVEL,double AVS,double avDIR[_3D],double S4,double stdN,double KBTNOW,double AVV[_3D],double AVNOW[_3D],kinTheory theory,specSwimmer specS,swimmer *sw ) {
	/*
	 Checkpoint the entire simulation
	 */
	int i,j;

	fprintf( fout,"%d\n",in.simSteps );		//total time (or number of iterations)
	fprintf( fout,"%d %lf\n",in.warmupSteps,in.dt );		//Warmup time iterations, and time step

	fprintf( fout,"%ld\n",in.seed );				//Random seed (0 if read from time)
	fprintf( fout,"%d %d %d %d %lf %lf\n",DIM,XYZ[0],XYZ[1],XYZ[2],in.KBT,KBTNOW );
	fprintf( fout,"%d %d %d %d %d %d\n",in.RFRAME,in.zeroNetMom,in.GALINV,in.TSTECH,in.RTECH,in.LC );
	fprintf( fout,"%lf %lf %lf %lf\n",in.TAU,in.RA,in.FRICCO,in.MFPOT );
	fprintf( fout,"%lf %lf %lf\n",in.GRAV[0],in.GRAV[1],in.GRAV[2] );		//Acceleration (external force)
	fprintf( fout,"%lf %lf %lf\n",in.MAG[0],in.MAG[1],in.MAG[2] );			//External magnetic field
	fprintf( fout,"%d %d\n",MDmode,in.stepsMD );		//MD coupling mode and number of MD steps per SRD step
	fprintf( fout,"%d %d\n",GPOP,NSPECI);			//Total number of particles and number of species

	fprintf( fout,"%d %d %lf %lf %d %d\n",runtime,warmtime,in.C,in.S,in.GRAV_FLAG,in.MAG_FLAG );
	fprintf( fout,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", theory.MFP, theory.VISC, theory.THERMD, theory.SDIFF, theory.SPEEDOFSOUND, theory.sumM, AVVEL, AVS, avDIR[0], avDIR[1], avDIR[2], S4, stdN, nDNST, mDNST );
	fprintf( fout,"%lf %lf %lf %lf %lf %lf\n",AVV[0], AVV[1], AVV[2], AVNOW[0], AVNOW[1], AVNOW[2] );

	//Output variables
	fprintf( fout,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",DBUG,outFlag.TRAJOUT,outFlag.printSP,outFlag.COAROUT,outFlag.FLOWOUT,outFlag.VELOUT,outFlag.AVVELOUT,outFlag.ORDEROUT,outFlag.QTENSOUT,outFlag.QKOUT,outFlag.AVSOUT,outFlag.SOLOUT,outFlag.ENOUT,outFlag.ENFIELDOUT,outFlag.ENNEIGHBOURS,outFlag.ENSTROPHYOUT,outFlag.DENSOUT,outFlag.CVVOUT,outFlag.CNNOUT,outFlag.CWWOUT,outFlag.CDDOUT,outFlag.CSSOUT,outFlag.CPPOUT,outFlag.BINDER,outFlag.BINDERBIN,outFlag.SYNOUT,outFlag.CHCKPNT,outFlag.CHCKPNTrcvr );
	fprintf( fout,"%d %d\n",outFlag.SPOUT,outFlag.PRESOUT );
	fprintf( fout,"%d %d %d %d %d %d %d\n",outFlag.HISTVELOUT,outFlag.HISTSPEEDOUT,outFlag.HISTVORTOUT,outFlag.HISTENSTROUT,outFlag.HISTDIROUT,outFlag.HISTSOUT,outFlag.HISTNOUT );
	fprintf( fout,"%d %d %d\n",outFlag.ENERGYSPECTOUT,outFlag.ENSTROPHYSPECTOUT,outFlag.DEFECTOUT );
	fprintf( fout,"%d %d %d\n",outFlag.SWOUT,outFlag.SWORIOUT,outFlag.RTOUT );

	//Species of MPCD particles
	for( i=0; i<NSPECI; i++ ) {
		fprintf( fout,"%lf %i %i %i %i ",(SP+i)->MASS,(SP+i)->POP,(SP+i)->QDIST,(SP+i)->VDIST,(SP+i)->ODIST );
		fprintf( fout,"%lf %lf %lf %lf %lf %lf %lf\n",(SP+i)->RFC, (SP+i)->LEN, (SP+i)->TUMBLE, (SP+i)->CHIHI, (SP+i)->CHIA, (SP+i)->ACT, (SP+i)->DAMP );
		for( j=0; j<NSPECI; j++ ) fprintf( fout,"%lf ",(SP+i)->M[j] );			//Binary fluid control parameters
		fprintf( fout,"\n" );
	}
	//BCs
	fprintf( fout,"%d\n",NBC );
	for( i=0; i<NBC; i++ ) {
		fprintf( fout,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",WALL[i].COLL_TYPE, WALL[i].PHANTOM, WALL[i].E, WALL[i].Q[0], WALL[i].Q[1], WALL[i].Q[2], WALL[i].V[0], WALL[i].V[1], WALL[i].V[2], WALL[i].O[0], WALL[i].O[1], WALL[i].O[2] );
		fprintf( fout,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", WALL[i].L[0], WALL[i].L[1], WALL[i].L[2], WALL[i].G[0], WALL[i].G[1], WALL[i].G[2], WALL[i].A[0], WALL[i].A[1], WALL[i].A[2], WALL[i].AINV[0], WALL[i].AINV[1], WALL[i].AINV[2],WALL[i].P[0],WALL[i].P[1],WALL[i].P[2],WALL[i].P[3], WALL[i].R );
		fprintf( fout,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", WALL[i].DN, WALL[i].DT, WALL[i].DVN, WALL[i].DVT, WALL[i].DVxyz[0], WALL[i].DVxyz[1], WALL[i].DVxyz[2], WALL[i].MVN, WALL[i].MVT, WALL[i].MUN, WALL[i].MUT, WALL[i].MUxyz[0], WALL[i].MUxyz[1], WALL[i].MUxyz[2] );
		fprintf( fout,"%lf %lf %lf %lf %d %d %lf\n", WALL[i].DUxyz[0], WALL[i].DUxyz[1], WALL[i].DUxyz[2], WALL[i].KBT, WALL[i].DSPLC, WALL[i].INV, WALL[i].MASS );
		fprintf( fout,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", WALL[i].W, WALL[i].VOL, WALL[i].Q_old[0], WALL[i].Q_old[1], WALL[i].Q_old[2], WALL[i].O_old[0], WALL[i].O_old[1], WALL[i].O_old[2], WALL[i].I[0][0], WALL[i].I[0][1], WALL[i].I[0][2], WALL[i].I[1][0], WALL[i].I[1][1], WALL[i].I[1][2], WALL[i].I[2][0], WALL[i].I[2][1], WALL[i].I[2][2] );
		fprintf( fout,"%d %d %d %lf %lf\n", WALL[i].PLANAR, WALL[i].REORIENT, WALL[i].ABS, WALL[i].ROTSYMM[0], WALL[i].ROTSYMM[1] );
		fprintf( fout,"%lf %lf %lf %lf %lf %lf\n", WALL[i].dV[0], WALL[i].dV[1], WALL[i].dV[2], WALL[i].dL[0], WALL[i].dL[1], WALL[i].dL[2] );
	}

	//MPCD particles
	for( i=0; i<GPOP; i++ ) {
		fprintf( fout,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",pSRD[i].S_flag, pSRD[i].SPID, pSRD[i].q, pSRD[i].Q[0], pSRD[i].Q[1], pSRD[i].Q[2], pSRD[i].V[0], pSRD[i].V[1], pSRD[i].V[2], pSRD[i].U[0], pSRD[i].U[1], pSRD[i].U[2], pSRD[i].T[0], pSRD[i].T[1], pSRD[i].T[2] );
	}

	//Swimmers 
	fprintf( fout,"%d %d %d %d %d %d %lf %lf %d %d\n", NS,specS.TYPE, specS.QDIST, specS.ODIST, specS.headM, specS.middM, specS.iheadM, specS.imiddM, specS.HSPid, specS.MSPid );
	fprintf( fout,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf\n", specS.FS, specS.TS, specS.DS, specS.sizeShrink, specS.springShrink, specS.fixDist, specS.k, specS.ro, specS.iro, specS.sig, specS.isig, specS.eps, specS.runTime, specS.tumbleTime, specS.shrinkTime, specS.MAGMOM );

	for( i=0; i<NS; i++ ) {
		fprintf( fout,"%d %lf %lf %lf %d %d %lf %lf %lf %lf %lf\n",(sw+i)->RT,(sw+i)->n0[0],(sw+i)->n0[1],(sw+i)->n0[2],(sw+i)->timeCNT,(sw+i)->timeRND,(sw+i)->ro,(sw+i)->iro,(sw+i)->sig,(sw+i)->isig,(sw+i)->k );
		fprintf( fout,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",(sw+i)->H.HorM,(sw+i)->H.Q[0],(sw+i)->H.Q[1],(sw+i)->H.Q[2],(sw+i)->H.V[0],(sw+i)->H.V[1],(sw+i)->H.V[2],(sw+i)->H.A[0],(sw+i)->H.A[1],(sw+i)->H.A[2] );
		fprintf( fout,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",(sw+i)->M.HorM,(sw+i)->M.Q[0],(sw+i)->M.Q[1],(sw+i)->M.Q[2],(sw+i)->M.V[0],(sw+i)->M.V[1],(sw+i)->M.V[2],(sw+i)->M.A[0],(sw+i)->M.A[1],(sw+i)->M.A[2] );
	}

	#ifdef FFLSH
		fflush(fout);
	#endif
}

void outputResults( cell ***CL,particleMPC *SRDparticles,spec SP[],bc WALL[],simptr simMD,specSwimmer SS, swimmer swimmers[],double AVNOW[_3D],double AVV[_3D],double avDIR[_3D], int runtime, inputList in, double AVVEL, double KBTNOW,double *AVS,double *S4,double *stdN,int MDmode,outputFlagsList outFlag,outputFilesList outFiles ) {
	/*
	 Output the results (except histograms of distributions, which are separate since they take lots of memory)
	 */
	int a,b,c,i,j;
	double time_now = runtime*in.dt;					//Simulation time
	double wmf;
	double corr[maxXYZ],spect[maxXYZ];				//Correlation functions and energy spectra
	double UL;																//Binder cumulant
	double avGradVel[_3D][_3D];								//Velocity gradient

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
	if( MDmode ) binMD( CL );
	//Calculate the local properties of each cell (VCM,in.KBT,POPulation,Mass)
	localPROP( CL,SP,SS,in.RTECH,in.LC );
	avVel( CL,AVNOW );
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
		galileantrans( SRDparticles,WALL,simMD,SP,in.KBT,AVV,GPOP,NBC,MDmode,DIM );
		zeroExtraDims( SRDparticles,WALL,simMD,GPOP,NBC,MDmode,DIM );
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
			for( i=0; i<_3D; i++ ) for( j=0; j<_3D; j++ ) avGradVel[i][j] /= (double)(XYZ[0]*XYZ[1]*XYZ[2]);
			avveloutWithGradVel( outFiles.favvel,time_now,AVNOW,KBTNOW,avGradVel );
		}
		//Enstrophy
		if( outFlag.ENSTROPHYOUT>=OUT && runtime%outFlag.ENSTROPHYOUT==0 ) {
			wmf = avEnstrophy( CL );
			avenstrophyout( outFiles.fenstrophy,time_now,wmf );
		}
	}
	/* ****************************************** */
	/* *************** TOTAL ENERGY ************* */
	/* ****************************************** */
	if( outFlag.ENOUT>=OUT && runtime%outFlag.ENOUT==0 ) {
		wmf = calcE_LC( CL,in.LC,in.MFPOT );
		enout( outFiles.fenergy,SRDparticles,SP,WALL,time_now,KBTNOW,wmf );
	}
	if( outFlag.ENFIELDOUT>=OUT && runtime%outFlag.ENFIELDOUT==0 ) enfieldout( outFiles.fenergyfield,CL,SP,in.MFPOT,in.LC );
	if( outFlag.ENNEIGHBOURS>=OUT && runtime%outFlag.ENNEIGHBOURS==0 ) enneighboursout( outFiles.fenneighbours,time_now,CL,in.MFPOT,in.LC );
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
		if( in.RTECH==MULTIPHASE ) {
			// phiphiCorr( CL,maxXYZ,XYZ,corr,DIM );
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
	if( outFlag.FLOWOUT>=OUT && runtime%outFlag.FLOWOUT==0 && runtime!=0 ) flowout( outFiles.fflow,CL,outFlag.FLOWOUT );
	if( outFlag.VELOUT>=OUT && runtime%outFlag.VELOUT==0 ) velout( outFiles.fvel,CL, time_now );
	if( outFlag.COAROUT>=OUT && runtime%outFlag.COAROUT==0 ) coarseout( outFiles.fcoarse,time_now,CL );
	if(in.LC!=ISOF) if( outFlag.ORDEROUT>=OUT && runtime%outFlag.ORDEROUT==0 ) orderout( outFiles.forder,time_now,CL,in.LC );
	if(in.LC!=ISOF) if( outFlag.QTENSOUT>=OUT && runtime%outFlag.QTENSOUT==0 ) orderQout( outFiles.forderQ,time_now,CL,in.LC );
	if( outFlag.SPOUT>=OUT && runtime%outFlag.SPOUT==0 ) multiphaseout( outFiles.fmultiphase,time_now,CL );
	if( outFlag.PRESOUT>=OUT && runtime%outFlag.PRESOUT==0 ) pressureout( outFiles.fpressure,time_now,CL );
	if(in.LC!=ISOF) if( outFlag.QKOUT>=OUT && runtime%outFlag.QKOUT==0 ) {
		#ifdef DBG
			if( DBUG >= DBGTITLE ) printf( "Calculate Q-tensor in reciprocal space.\n" );
		#endif
		orderQKout( outFiles.forderQK,time_now,SRDparticles,CL,in.LC );
	}
	/* ****************************************** */
	/* ************** TRACK DEFECTS ************* */
	/* ****************************************** */
	if (outFlag.DEFECTOUT>=OUT && runtime%outFlag.DEFECTOUT==0 && DIM==_2D) topochargeout( outFiles.fdefects, time_now, CL);
}

void outputHist( cell ***CL,int runtime, inputList in,outputFlagsList outFlag,outputFilesList outFiles ) {
	/*
	 Output histograms of distributions
	 */
	int a,b,c,i,j;
	double time_now = runtime*in.dt;
	double myVec[_3D];													//Velocity (etc) actual values for every MPCD cell
	double maxRange;														//Maximum for range for histograms
	int nc=XYZ[0]*XYZ[1]*XYZ[2];
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

void closeOutputFiles( spec *SP,bc WALL[],outputFlagsList outFlag,outputFilesList outFiles ) {
	int i;

	if( outFlag.TRAJOUT>=OUT ) for( i=0; i<NSPECI; i++ ) if( SP[i].POP>=1 ) fclose( outFiles.fdetail[i] );
	if( outFlag.COAROUT>=OUT ) fclose( outFiles.fcoarse );
	if( outFlag.AVVELOUT>=OUT ) fclose( outFiles.favvel );
	if( outFlag.ORDEROUT>=OUT ) fclose( outFiles.forder );
	if( outFlag.QTENSOUT>=OUT ) fclose( outFiles.forderQ );
	if( outFlag.QKOUT>=OUT ) fclose( outFiles.forderQK );
	if( outFlag.AVSOUT>=OUT ) fclose( outFiles.favs );
	if( outFlag.DENSOUT>=OUT ) fclose( outFiles.fdensSTD );
	if( outFlag.ENSTROPHYOUT>=OUT ) fclose( outFiles.fenstrophy );
	if( outFlag.FLOWOUT>=OUT ) fclose( outFiles.fflow );
	if( outFlag.VELOUT>=OUT ) fclose( outFiles.fvel );
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
	if( outFlag.DEFECTOUT>=OUT ) fclose( outFiles.fdefects );
	if( outFlag.SPOUT>=OUT ) fclose( outFiles.fmultiphase );
	if( outFlag.PRESOUT>=OUT ) fclose( outFiles.fpressure );
	if( outFlag.SWOUT>=OUT ) fclose( outFiles.fswimmers );
	if( outFlag.SWORIOUT>=OUT ) fclose( outFiles.fswimmersOri );
	if( outFlag.RTOUT>=OUT ) fclose( outFiles.fruntumble );
	if( outFlag.SOLOUT>=OUT ) for( i=0; i<NBC; i++ ) if( WALL[i].DSPLC ) fclose( outFiles.fsolids[i] );
}

int writeOutput( int t,outputFlagsList f,int RFRAME,int zeroNetMom ) {
	if( ( RFRAME && t%zeroNetMom==0 ) || ( f.ENOUT>=OUT && t%f.ENOUT==0 ) || ( f.TRAJOUT>=OUT  && t%f.TRAJOUT==0 ) || ( f.AVVELOUT>=OUT && t%f.AVVELOUT==0 ) || ( f.QKOUT && t%f.QKOUT==0 ) || ( f.AVSOUT>=OUT && t%f.AVSOUT==0 ) || ( f.ENNEIGHBOURS>=OUT && t%f.ENNEIGHBOURS==0 ) || ( f.SOLOUT>=OUT && t%f.SOLOUT==0 ) || ( f.BINDER && t%f.BINDER==0 ) || ( f.SWOUT && t%f.SWOUT==0 ) || ( f.SWORIOUT && t%f.SWORIOUT==0 ) ) {
		return 1;
	}
	//Fields
	else if( ( f.FLOWOUT>=OUT && t%f.FLOWOUT==0 ) || (f.VELOUT>=OUT && t%f.VELOUT==0 ) || ( f.COAROUT>=OUT && t%f.COAROUT==0 ) || ( f.ENFIELDOUT>=OUT && t%f.ENFIELDOUT==0 ) || ( f.ORDEROUT && t%f.ORDEROUT==0 ) || ( f.QTENSOUT && t%f.QTENSOUT==0 ) || ( f.SPOUT && t%f.SPOUT==0 ) || ( f.PRESOUT && t%f.PRESOUT==0 ) ) {
		return 1;
	}
	//Correlation functions
	else if( ( f.CVVOUT && t%f.CVVOUT==0 ) || ( f.CNNOUT && t%f.CNNOUT==0 ) || ( f.CWWOUT && t%f.CWWOUT==0 ) || ( f.CDDOUT && t%f.CDDOUT==0 ) || ( f.CSSOUT && t%f.CSSOUT==0 ) || ( f.CPPOUT && t%f.CPPOUT==0 ) ) {
		return 1;
	}
	else return 0;
}
int writeHistograms( int t,outputFlagsList f ) {
	if( ( f.HISTVELOUT>=OUT && t%f.HISTVELOUT==0 ) || ( f.HISTSPEEDOUT>=OUT && t%f.HISTSPEEDOUT==0 ) || ( f.HISTVORTOUT>=OUT && t%f.HISTVORTOUT==0 ) || ( f.HISTENSTROUT>=OUT && t%f.HISTENSTROUT==0 ) || ( f.HISTDIROUT>=OUT && t%f.HISTDIROUT==0 ) || ( f.HISTSOUT>=OUT && t%f.HISTSOUT==0 ) || ( f.HISTNOUT>=OUT && t%f.HISTNOUT==0 ) ) {
		return 1;
	}
	else return 0;
}
