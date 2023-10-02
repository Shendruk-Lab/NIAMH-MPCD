
//================================================================================
//
// name:   mdhistogramdef.h
// author: ftessier
// date:   2005-05-03 @ 11:04:25
//
// Definition of histograms for molecular dynamics program
//
//================================================================================

//----------------------------------------------------------------------------
// velocity distribution
//----------------------------------------------------------------------------
{
    // create histogram
    h = HistogramNew (&sim->histograms);
    snprintf (h->label,			STRMAX, "fv");
    snprintf (h->histfuncName,  STRMAX, "HistogramFunction");
    snprintf (h->funcName,		STRMAX, "HistogramFunction_fv");
    h->active 	= 1;
    h->file 	= 1;
    h->groupInc	= GROUP_ALL;
    h->normMode = NORM_INTEGRAL;
    h->histLoop = LOOP_ATOM;

    // space domain (only consider atoms in this region)
    h->sdomain	= DOMAIN_ALL;

    // histogram domain
    h->coord	= SPHERICAL;
    h->d	 	= 1;
    h->domain	= r_;
    {
		h->min1 = 0;
		h->max1 = 5*sqrt(0.5*DIM_MD*sim->kT[sim->phase]);
		h->bin1 = 0.01;
    }

    // histogram step counters
    CopyStepCounters (h->stepCollect, sim->stepHistCollect);
    CopyStepCounters (h->stepPrint,   sim->stepHistPrint);

    // allocate memory
    HistogramAllocateBins (h);
}


//----------------------------------------------------------------------------
// ALL: height particle density profile
//----------------------------------------------------------------------------
{
    // create histogram
    h = HistogramNew (&sim->histograms);
    snprintf (h->label,     	STRMAX, "rho-h");
    snprintf (h->histfuncName,  STRMAX, "HistogramFunction");
    h->active 	= 1;
    h->file 	= 1;
    h->groupInc	= GROUP_ALL;
    h->normMode = NORM_VOLUME;
    h->histLoop = LOOP_ATOM;

    // space domain (only consider atoms in this region)
    h->sdomain	= DOMAIN_ALL;

    // histogram domain
    h->coord	= CARTESIAN;
    h->d		= 3;
    h->domain	= xyz_;
    {
		h->min1 = 0;
		h->max1 = sim->box[x_];
		h->n1 = 1;
		h->bin1 = 0.05;

		h->min2 = 0;
		h->max2 = sim->box[y_];
		h->bin2 = 0.05;

		h->min3 = 0;
		h->max3 = sim->box[z_];
		h->n3   = 1;
    }

    // histogram step counters
    CopyStepCounters (h->stepCollect, sim->stepHistCollect);
    CopyStepCounters (h->stepPrint,   sim->stepHistPrint);

    // allocate memory
    HistogramAllocateBins (h);
}


//----------------------------------------------------------------------------
// ALL: radial particle density profile
//----------------------------------------------------------------------------
{
    // create histogram
    h = HistogramNew (&sim->histograms);
    snprintf (h->label,     	STRMAX, "rho-r");
    snprintf (h->histfuncName,  STRMAX, "HistogramFunction");
    h->active 	= 1;
    h->file 	= 1;
    h->groupInc	= GROUP_ALL;
    h->normMode = NORM_VOLUME;
    h->histLoop = LOOP_ATOM;

    // space domain (only consider atoms in this region)
    h->sdomain	= DOMAIN_ALL;

    // histogram domain
    h->coord	= CYLINDRICAL | x_;
    h->d		= 3;
    h->domain	= rqz_;
    {
		h->min1 = 0;
		h->max1 = sim->caprOut;
		h->bin1 = 0.05;

		h->min2 = 0;
		h->max2 = 2*pi;
		h->n2   = 1;

		h->min3 = 0;
		h->max3 = sim->box[x_];
		h->n3   = 1;
    }

    // histogram step counters
    CopyStepCounters (h->stepCollect, sim->stepHistCollect);
    CopyStepCounters (h->stepPrint,   sim->stepHistPrint);

    // allocate memory
    HistogramAllocateBins (h);
}


//----------------------------------------------------------------------------
// MONOMER: brush profile
//----------------------------------------------------------------------------
{
    // create histogram
    h = HistogramNew (&sim->histograms);
    snprintf (h->label,     	STRMAX, "MONOMER-rho-rz");
    snprintf (h->histfuncName,  STRMAX, "HistogramFunction");
    h->active 	= 1;
    h->file 	= 1;
    h->groupInc	= GROUP_MONOMER;
    h->normMode = NORM_MEASURE;
    h->histLoop = LOOP_ATOM;

    // space domain (only consider atoms in this region)
    h->sdomain	= DOMAIN_ALL;

    // histogram domain
    h->coord	= CYLINDRICAL | x_;
    h->d		= 2;
    h->domain	= rz_;
    {
		h->min1 = 0;
		h->max1 = sim->caprIn;
		h->bin1 = 1;

		h->min2 = -sim->boxHalf[x_];
		h->max2 = sim->boxHalf[x_];
		h->bin2 = 1;
    }

    // histogram step counters
    CopyStepCounters (h->stepCollect, sim->stepHistCollect);
    CopyStepCounters (h->stepPrint,   sim->stepHistPrint);

    // allocate memory
    HistogramAllocateBins (h);
}


//----------------------------------------------------------------------------
// MONOMER: azimutal density profile
//----------------------------------------------------------------------------
{
    // create histogram
    h = HistogramNew (&sim->histograms);
    snprintf (h->label,     	STRMAX, "MONOMER-rho-qz");
    snprintf (h->histfuncName,  STRMAX, "HistogramFunction");
    h->active 	= 1;
    h->file 	= 1;
    h->groupInc	= GROUP_MONOMER;
    h->normMode = NORM_MEASURE;
    h->histLoop = LOOP_ATOM;

    // space domain (only consider atoms in this region)
    h->sdomain	= DOMAIN_ALL;

    // histogram domain
    h->coord	= CYLINDRICAL | x_;
    h->d		= 2;
    h->domain	= qz_;
    {
		h->min1 = 0;
		h->max1 = 2*pi;
        if (sim->geometry == GEOM_CYLINDER) // fix for dumb FPE
            h->bin1 = 1/sim->caprIn;
        else
		    h->bin1 = 0;

		h->min2 = -sim->boxHalf[x_];
		h->max2 = sim->boxHalf[x_];
		h->bin2 = 1;
    }

    // histogram step counters
    CopyStepCounters (h->stepCollect, sim->stepHistCollect);
    CopyStepCounters (h->stepPrint,   sim->stepHistPrint);

    // allocate memory
    HistogramAllocateBins (h);
}


//----------------------------------------------------------------------------
// MONOMER: radial density profile
//----------------------------------------------------------------------------
{
    // create histogram
    h = HistogramNew (&sim->histograms);
    snprintf (h->label,     	STRMAX, "MONOMER-rho-r");
    snprintf (h->histfuncName,  STRMAX, "HistogramFunction");
    h->active 	= 1;
    h->file 	= 1;
    h->groupInc	= GROUP_MONOMER;
    h->normMode = NORM_VOLUME;
    h->histLoop = LOOP_ATOM;

    // space domain (only consider atoms in this region)
    h->sdomain	= DOMAIN_ALL;

    // histogram domain
    h->coord	= CYLINDRICAL | x_;
    h->d		= 3;
    h->domain	= rqz_;
    {
		h->min1 = 0;
		h->max1 = sim->caprOut;
		h->bin1 = 0.05;

		h->min2 = 0;
		h->max2 = 2*pi;
		h->n2   = 1;

		h->min3 = -sim->boxHalf[x_];
		h->max3 = sim->boxHalf[x_];
		h->n3   = 1;
    }

    // histogram step counters
    CopyStepCounters (h->stepCollect, sim->stepHistCollect);
    CopyStepCounters (h->stepPrint,   sim->stepHistPrint);

    // allocate memory
    HistogramAllocateBins (h);
}


//----------------------------------------------------------------------------
// MONOMER: cross-sectional density profile
//----------------------------------------------------------------------------
{
    // create histogram
    h = HistogramNew (&sim->histograms);
    snprintf (h->label,    		STRMAX, "MONOMER-rho-yz");
    snprintf (h->histfuncName, 	STRMAX, "HistogramFunction");
    h->active   = 1;
    h->file 	= 1;
    h->groupInc	= GROUP_MONOMER;
    h->normMode = NORM_VOLUME;
    h->histLoop = LOOP_ATOM;

    // space domain (only consider atoms in this region)
    h->sdomain  = DOMAIN_ALL;

    // histogram domain
    h->coord	= CARTESIAN;
    h->d		= 3;
    h->domain	= xyz_;
    {
		h->min1 = -sim->boxHalf[x_];
		h->max1 = sim->boxHalf[x_];
		h->n1   = 1;

		h->min2 = -sim->boxHalf[y_];
		h->max2 = sim->boxHalf[y_];
		h->n2   = 32;

		h->min3 = -sim->boxHalf[z_];
		h->max3 = sim->boxHalf[z_];
		h->n3   = 32;
    }

    // histogram step counters
    CopyStepCounters (h->stepCollect, sim->stepHistCollect);
    CopyStepCounters (h->stepPrint,   sim->stepHistPrint);

    // allocate memory
    HistogramAllocateBins (h);
}


//----------------------------------------------------------------------------
// MONOMER: bond length distribution
//----------------------------------------------------------------------------
{
    // create histogram
    h = HistogramNew (&sim->histograms);
    snprintf (h->label,			STRMAX, "MONOMER-bond");
    snprintf (h->histfuncName,  STRMAX, "HistogramFunction");
    snprintf (h->funcName,		STRMAX, "HistogramFunction_bond");
    h->active 	= 1;
    h->file 	= 1;
    h->groupInc	= GROUP_MONOMER;
    h->normMode = NORM_INTEGRAL;
    h->histLoop = LOOP_ATOM;

    // space domain (only consider atoms in this region)
    h->sdomain	= DOMAIN_ALL;

    // histogram domain
    h->coord	= SPHERICAL;
    h->d		= 1;
    h->domain	= r_;
    {
		h->min1 = 0;
		h->max1 = 1.2;
		h->bin1 = 0.05;
    }

    // histogram step counters
    CopyStepCounters (h->stepCollect, sim->stepHistCollect);
    CopyStepCounters (h->stepPrint,   sim->stepHistPrint);

    // allocate memory
    HistogramAllocateBins (h);
}


//----------------------------------------------------------------------------
// MONOMER: bond length distribution
//----------------------------------------------------------------------------
{
    // create histogram
    h = HistogramNew (&sim->histograms);
    snprintf (h->label,			STRMAX, "MONOMER-bond");
    snprintf (h->histfuncName,  STRMAX, "HistogramFunction");
    snprintf (h->funcName,		STRMAX, "HistogramFunction_bond");
    h->active 	= 1;
    h->file 	= 1;
    h->groupInc	= GROUP_MONOMER;
    h->normMode = NORM_INTEGRAL;
    h->histLoop = LOOP_POLYMER;

    // space domain (only consider atoms in this region)
    h->sdomain	= DOMAIN_ALL;

    // histogram domain
    h->coord	= SPHERICAL;
    h->d		= 1;
    h->domain	= r_;
    {
		h->min1 = 0.8;
		h->max1 = 1.2;
		h->bin1 = 0.001;
    }

    // histogram step counters
    CopyStepCounters (h->stepCollect, sim->stepHistCollect);
    CopyStepCounters (h->stepPrint,   sim->stepHistPrint);

    // allocate memory
    HistogramAllocateBins (h);
}


//----------------------------------------------------------------------------
// MONOMER: bond forces as a function of monomer number
//----------------------------------------------------------------------------
{
    // create histogram
    h = HistogramNew (&sim->histograms);
    snprintf (h->label,			STRMAX, "MONOMER-fenex-N");
    snprintf (h->histfuncName,  STRMAX, "HistogramFunction");
    snprintf (h->funcName,		STRMAX, "HistogramFunction_fenex_N");
    h->active 	= 1;
    h->file 	= 1;
    h->groupInc	= GROUP_ALL;
    h->normMode = NORM_AVERAGE;
    h->histLoop = LOOP_POLYMER;

    // space domain (only consider atoms in this region)
    h->sdomain	= DOMAIN_ALL;

    // histogram domain
    h->coord	= CARTESIAN;
    h->d		= 1;
    h->domain	= x_;
    {
		h->min1 = 0;
		h->max1 = 30;
		h->n1   = 30;
    }

    // histogram step counters
    CopyStepCounters (h->stepCollect, sim->stepHistCollect);
    CopyStepCounters (h->stepPrint,   sim->stepHistPrint);

    // allocate memory
    HistogramAllocateBins (h);
}


//----------------------------------------------------------------------------
// MONOMER: distance to wall as a function of monomer number
//----------------------------------------------------------------------------
{
    // create histogram
    h = HistogramNew (&sim->histograms);
    snprintf (h->label,			STRMAX, "MONOMER-hwall-N");
    snprintf (h->histfuncName,  STRMAX, "HistogramFunction");
    snprintf (h->funcName,		STRMAX, "HistogramFunction_hwall_N");
    h->active 	= 1;
    h->file 	= 1;
    h->groupInc	= GROUP_ALL;
    h->normMode = NORM_AVERAGE;
    h->histLoop = LOOP_POLYMER;

    // space domain (only consider atoms in this region)
    h->sdomain	= DOMAIN_ALL;

    // histogram domain
    h->coord	= CARTESIAN;
    h->d		= 1;
    h->domain	= x_;
    {
		h->min1 = 0;
		h->max1 = 30;
		h->n1   = 30;
    }

    // histogram step counters
    CopyStepCounters (h->stepCollect, sim->stepHistCollect);
    CopyStepCounters (h->stepPrint,   sim->stepHistPrint);

    // allocate memory
    HistogramAllocateBins (h);
}


//----------------------------------------------------------------------------
// MONOMER: Effective charge as a function of monomer number
//----------------------------------------------------------------------------
{
    // create histogram
    h = HistogramNew (&sim->histograms);
    snprintf (h->label,			STRMAX, "MONOMER-qeff-N");
    snprintf (h->histfuncName,  STRMAX, "HistogramFunction");
    snprintf (h->funcName,		STRMAX, "HistogramFunction_qeff_N");
    h->active 	= 1;
    h->file 	= 1;
    h->groupInc	= GROUP_ALL;
    h->normMode = NORM_AVERAGE;
    h->histLoop = LOOP_POLYMER;

    // space domain (only consider atoms in this region)
    h->sdomain	= DOMAIN_ALL;

    // histogram domain
    h->coord	= CARTESIAN;
    h->d		= 1;
    h->domain	= x_;
    {
		h->min1 = 0;
		h->max1 = 100;
		h->n1   = 100;
    }

    // histogram step counters
    CopyStepCounters (h->stepCollect, sim->stepHistCollect);
    CopyStepCounters (h->stepPrint,   sim->stepHistPrint);

    // allocate memory
    HistogramAllocateBins (h);
}



//----------------------------------------------------------------------------
// MONOMER: monomer ellipsoid
//----------------------------------------------------------------------------
{
    // create histogram
    h = HistogramNew (&sim->histograms);
    snprintf (h->label,			STRMAX, "MONOMER-ellipsoid");
    snprintf (h->histfuncName,  STRMAX, "HistogramFunction");
    snprintf (h->funcName,		STRMAX, "HistogramFunction_ellipsoid");
    h->active 	= 1;
    h->file 	= 1;
    h->groupInc	= GROUP_ALL;
    h->normMode = NORM_AVERAGE;
    h->histLoop = LOOP_POLYMER;

    // space domain (only consider atoms in this region)
    h->sdomain	= DOMAIN_ALL;

    // histogram domain
    h->coord	= CARTESIAN;
    h->d		= 3;
    h->domain	= xyz_;
    {
		h->min1 = -5;
		h->max1 = 25;
		h->n1   = 30;

		h->min2 = -5;
		h->max2 = 5;
		h->n2   = 10;

		h->min3 = 0;
		h->max3 = 10;
		h->n3   = 10;
    }

    // histogram step counters
    CopyStepCounters (h->stepCollect, sim->stepHistCollect);
    CopyStepCounters (h->stepPrint,   sim->stepHistPrint);

    // allocate memory
    HistogramAllocateBins (h);
}


// //----------------------------------------------------------------------------
// // QPOS: qz density profile
// //----------------------------------------------------------------------------
// {
//     // create histogram
//     h = HistogramNew (&sim->histograms);
//     snprintf (h->label,     	STRMAX, "QPOS-rho-qz");
//     snprintf (h->histfuncName,  STRMAX, "HistogramFunction");
//     h->active 	= 1;
//     h->file 	= 1;
//     h->groupInc	= GROUP_ION_POS;
//     h->normMode = NORM_MEASURE;
//     h->histLoop = LOOP_ATOM;
//
//     // space domain (only consider atoms in this region)
// 	h->scoord 	= CYLINDRICAL | x_;
// 	h->sd		= 1;
//     h->sdomain	= r_;
// 	{
// 		h->smin1 = 0;
// 		h->smax1 = sim->caprIn;
// 	}
//
//     // histogram domain
//     h->coord	= CYLINDRICAL | x_;
//     h->d		= 2;
//     h->domain	= qz_;
//     {
// 		h->min1 = 0;
// 		h->max1 = 2*pi;
// 		h->bin1 = 1/sim->caprIn;
//
// 		h->min2 = -sim->boxHalf[x_];
// 		h->max2 = sim->boxHalf[x_];
// 		h->bin2 = 1;
//     }
//
//     // histogram step counters
// 	FillStepCounters (h->stepCollect, 100);
// 	FillStepCounters (h->stepPrint,   100);
//
//     // allocate memory
//     HistogramAllocateBins (h);
// }


//----------------------------------------------------------------------------
// QPOS: radial charge density profile
//----------------------------------------------------------------------------
{
    // create histogram
    h = HistogramNew (&sim->histograms);
    snprintf (h->label,    		STRMAX, "QPOS-r");
    snprintf (h->histfuncName,	STRMAX, "HistogramFunction");
    snprintf (h->funcName,		STRMAX, "HistogramFunction_q");
    h->active   = 1;
    h->file 	= 1;
    h->groupInc	= GROUP_ION_POS;
    h->normMode = NORM_VOLUME;
    h->histLoop = LOOP_ATOM;

    // space domain (only consider atoms in this region)
    h->sdomain  = DOMAIN_ALL;

    // histogram domain
    h->coord	= CYLINDRICAL | x_;
    h->d		= 3;
    h->domain	= rqz_;
    {
		h->min1 = 0;
		h->max1 = sim->caprOut;
		h->bin1 = 0.05;

		h->min2 = 0;
		h->max2 = 2*pi;
		h->n2   = 1;

		h->min3 = -sim->boxHalf[x_];
		h->max3 = sim->boxHalf[x_];
		h->n3   = 1;
    }

    // histogram step counters
    CopyStepCounters (h->stepCollect, sim->stepHistCollect);
    CopyStepCounters (h->stepPrint,   sim->stepHistPrint);

    // allocate memory
    HistogramAllocateBins (h);
}


//----------------------------------------------------------------------------
// QPOS; radial profile of the axial velocity
//----------------------------------------------------------------------------
{
    // create histogram
    h = HistogramNew (&sim->histograms);
    snprintf (h->label,         STRMAX, "QPOS-vx-r");
    snprintf (h->histfuncName,  STRMAX, "HistogramFunction");
    snprintf (h->funcName,      STRMAX, "HistogramFunction_vx");
    h->active   = 1;
    h->file     = 1;
    h->groupInc = GROUP_ION_POS;
    h->normMode = NORM_AVERAGE;
    h->histLoop = LOOP_ATOM;

    // space domain (only consider atoms in this region)
    h->sdomain  = DOMAIN_ALL;

    // histogram domain
    h->coord    = CYLINDRICAL | x_;
    h->d        = 1;
    h->domain   = r_;
    {
    	h->min1 = 0;
    	h->max1 = sim->caprIn;
    	h->bin1 = 0.05;
    }

    // histogram step counters
    CopyStepCounters (h->stepCollect, sim->stepHistCollect);
    CopyStepCounters (h->stepPrint,   sim->stepHistPrint);

    // allocate memory
    HistogramAllocateBins (h);
}


//----------------------------------------------------------------------------
// QPOS: axial charge density profile
//----------------------------------------------------------------------------
{
    // create histogram
    h = HistogramNew (&sim->histograms);
    snprintf (h->label,    		STRMAX, "QPOS-x");
    snprintf (h->histfuncName,	STRMAX, "HistogramFunction");
    snprintf (h->funcName,		STRMAX, "HistogramFunction_q");
    h->active   = 1;
    h->file 	= 1;
    h->groupInc	= GROUP_ION_POS;
    h->normMode = NORM_VOLUME;
    h->histLoop = LOOP_ATOM;

    // space domain (only consider atoms in this region)
    h->sdomain  = DOMAIN_ALL;

    // histogram domain
    h->coord	= CYLINDRICAL | x_;
    h->d		= 3;
    h->domain	= rqz_;
    {
		h->min1 = 0;
		h->max1 = sim->caprOut;
		h->n1   = 1;

		h->min2 = 0;
		h->max2 = 2*pi;
		h->n2   = 1;

		h->min3 = -sim->boxHalf[x_];
		h->max3 = sim->boxHalf[x_];
		h->bin3 = 0.1;
    }

    // histogram step counters
    CopyStepCounters (h->stepCollect, sim->stepHistCollect);
    CopyStepCounters (h->stepPrint,   sim->stepHistPrint);

    // allocate memory
    HistogramAllocateBins (h);
}


//----------------------------------------------------------------------------
// QNEG radial charge density profile
//----------------------------------------------------------------------------
{
    // create histogram
    h = HistogramNew (&sim->histograms);
    snprintf (h->label,    		STRMAX, "QNEG-r");
    snprintf (h->histfuncName, 	STRMAX, "HistogramFunction");
    snprintf (h->funcName,		STRMAX, "HistogramFunction_q");
    h->active   = 1;
    h->file 	= 1;
    h->groupInc	= GROUP_ION_NEG;
	h->groupExc = GROUP_WALL;
    h->normMode = NORM_VOLUME;
    h->histLoop = LOOP_ATOM;

    // space domain (only consider atoms in this region)
    h->sdomain  = DOMAIN_ALL;

    // histogram domain
    h->coord	= CYLINDRICAL | x_;
    h->d		= 3;
    h->domain	= rqz_;
    {
		h->min1 = 0;
		h->max1 = sim->caprOut;
		h->bin1 = 0.05;

		h->min2 = 0;
		h->max2 = 2*pi;
		h->n2   = 1;

		h->min3 = -sim->boxHalf[x_];
		h->max3 = sim->boxHalf[x_];
		h->n3   = 1;
    }

    // histogram step counters
    CopyStepCounters (h->stepCollect, sim->stepHistCollect);
    CopyStepCounters (h->stepPrint,   sim->stepHistPrint);

    // Allocate memory
    HistogramAllocateBins (h);
}


//----------------------------------------------------------------------------
// QNEG radial profile of the axial velocity
//----------------------------------------------------------------------------
{
    // create histogram
    h = HistogramNew (&sim->histograms);
    snprintf (h->label,         STRMAX, "QNEG-vx-r");
    snprintf (h->histfuncName,  STRMAX, "HistogramFunction");
    snprintf (h->funcName,		STRMAX, "HistogramFunction_vx");
    h->active   = 1;
    h->file     = 1;
    h->groupInc = GROUP_ION_NEG;
    h->groupExc = GROUP_WALL;
    h->normMode = NORM_AVERAGE;
    h->histLoop = LOOP_ATOM;

    // space domain (only consider atoms in this region)
    h->sdomain  = DOMAIN_ALL;

    // histogram domain
    h->coord    = CYLINDRICAL | x_;
    h->d        = 1;
    h->domain   = r_;
    {
    	h->min1 = 0;
    	h->max1 = sim->caprIn;
     	h->bin1 = 0.05;
	}

    // histogram step counters
    CopyStepCounters (h->stepCollect, sim->stepHistCollect);
    CopyStepCounters (h->stepPrint,   sim->stepHistPrint);

    // allocate memory
    HistogramAllocateBins (h);
}


// //----------------------------------------------------------------------------
// // QNEG: qz density profile
// //----------------------------------------------------------------------------
// {
//     // create histogram
//     h = HistogramNew (&sim->histograms);
//     snprintf (h->label,     	STRMAX, "QNEG-rho-qz");
//     snprintf (h->histfuncName,  STRMAX, "HistogramFunction");
//     h->active 	= 0;
//     h->file 	= 1;
//     h->groupInc	= GROUP_ION_NEG;
//     h->normMode = NORM_MEASURE;
//     h->histLoop = LOOP_ATOM;
//
//     // space domain (only consider atoms in this region)
// 	h->scoord 	= CYLINDRICAL | x_;
// 	h->sd		= 1;
//     h->sdomain	= r_;
// 	{
// 		h->smin1 = sim->caprIn;
// 		h->smax1 = sim->caprOut;
// 	}
//
//     // histogram domain
//     h->coord	= CYLINDRICAL | x_;
//     h->d		= 2;
//     h->domain	= qz_;
//     {
// 		h->min1 = 0;
// 		h->max1 = 2*pi;
// 		h->bin1 = 1/sim->caprIn;
//
// 		h->min2 = -sim->boxHalf[x_];
// 		h->max2 = sim->boxHalf[x_];
// 		h->bin2 = 1;
//     }
//
//     // histogram step counters
// 	FillStepCounters (h->stepCollect, 100);
// 	FillStepCounters (h->stepPrint,   100);
//
//     // allocate memory
//     HistogramAllocateBins (h);
// }


