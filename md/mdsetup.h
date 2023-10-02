
//================================================================================
//
// name:   mdsetup.h
// author: ftessier
// date:   2005-05-03 @ 11:04:36
//
//================================================================================


#ifndef MDSETUP_H
#define MDSETUP_H


//================================================================================
// Prototypes
//================================================================================

simptr			SetupSimulation				(int argc, char *argv[]);
void 			SetupNewWorld				(simptr sim);
void 			SetupCheckpointWorld		(simptr sim);
void 			SetupRandomGenerator 		(simptr sim);
void 			SetupParameters				(simptr sim);
void 			SetupParticles				(simptr sim);
void 			InitRandomCoord 			(simptr sim, int type, int layout);
void 			InitLatticeCoord			(simptr sim, int type, int layout);
void 			InitCapillaryRadius			(simptr sim);
void 			InitCapillaryPore			(simptr sim);
void 			InitPolymers 				(simptr sim);
void    		InitCharges					(simptr sim);
void 			InitVelocities				(simptr sim);
void 			SetupLayoutLists 			(simptr sim);
void 			SetupPolymerList 			(simptr sim);
void 			SetupChargeList 			(simptr sim);
void 			SetupNeighborList			(simptr sim);
void 			SetupAnchorList				(simptr sim);
void 			SetupFeneList				(simptr sim);
void 			SetupBendList				(simptr sim);
void 			SetupGroups 				(simptr sim);
void 			SetupMeasurements 			(simptr sim);
void 			SetupStepCounters			(simptr sim);
void			RelaxParticles				(simptr sim);
void			RelaxStep					(simptr sim);
void			JiggleParticles				(simptr sim);
void			FinishSimulation			(simptr sim);

particleMD*		GrowLinearChain 			(simptr sim, int type, int layout, int n, particleMD *p0, int *status);
particleMD*		GrowRodChain 			(simptr sim, int type, int layout, int n, particleMD *p0, int *status, int dir);
particleMD*		GrowUChain 			(simptr sim, int type, int layout, int n, particleMD *p0, int *status);
particleMD*       AtomInsert 					(simptr sim, int type, int layout, particleMD *p0, int chkolap, int chkptrs);
int       		AtomRemove 					(simptr sim, int type, int n, particleMD *p0);
void 			OffsetParticlePointers 		(simptr sim, int n, int offset);
void 			OffsetPointer 				(particleMD **p, particleMD *p1, particleMD *p2, int offset);
int 			LayoutRule 					(simptr sim, int layout, real x, real y, real z);
int 			AtomCheckOverlap 			(simptr sim, int type, real x, real y, real z);
void 			NewCounter					(listStep *list, int *counter, void	(*action)());


#endif
