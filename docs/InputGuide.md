# Guide to MPCD JSON Input

## Contents
1. [Introduction](#introduction-in)
2. [Tag Tables](#tag-tables)
    - [Main Input Tag Table](#main-input-tag-table)
        - [Overrides](#overrides)
    - [Species Tag Table](#species-tag-table)
        - [Species Overrides](#species-overrides)
    - [BC Tag Table](#bc-tag-table)
        - [BC Overrides](#bc-overrides)
3. [Example](#example)
4. [Guide to adding new input parameters](#guide-to-adding-new-input-parameters)

## Introduction         {#introduction-in}

This is just meant to be a quick and dirty guide to the new MPCD JSON input system. It is not meant to be a complete reference. For basics on the JSON file format take a look at the tutorial series starting on [this page](https://www.w3schools.com/js/js_json_intro.asp).

**Note:** The MPCD JSON input system is unverified with MD simulations as of yet.

Everything listed below should be in a single JSON file. 
Every parameter in the simulation can be fed in with a corresponding name/value pair. 
The name, or "tag", uniquely represents each parameter, is case sensitive, but can be given in any order within the json file. 

There are some exceptions to this rule, namely overrides.
Overrides do not correspond to simulation input parameters, but instead will override some other behaviour in a helpful manner.
These are given in their own table, and each override explained in more detail.

If a tag is not specified, then that parameter will assume the default value in the table below --- This means that, with the exception of some parameters used when declaring a BC, you may only need to specify several parameters in your input file.
Some defaults will instead set an override.

The table below lists all tags, their types, their default values, and their description.
Tags will generally be named the same as they are named in the code (in the various input structs), unless there is a good reason to change it for legibility or clarity.
When in doubt about which tags correspond to which variable names in the code, please consult SRDclss.h, definitions.h and the subroutine readJson() in read.c.

Some tags will take an array of custom objects (for example, species and boundary conditions). 
These are each defined in their own table following the main one.

We have added support for comments in the JSON by using comment tags.
If your tag is one of `"c"`, `"comment"`, `"//"`, or `"#"` then it will be ignored by the parser.
Comment tags can be repeated (despite it being invalid JSON) and not cause issues with the simulation.

## Tag Tables           {#tag-tables}

### Main Input Tag Table            {#main-input-tag-table}

All tags suffixed with "out", unless otherwise specified, take a value representing how frequently (in timesteps) you want that quantity to be dumped to file. 

Tag             | Type          | Default Value | Description
---             | ---           | ---           | ---
`domain`        | array(int)    | [30, 30]      | The system size, given as [X, Y] or [X, Y, Z]. Array dimensions correspond to DIM. If attempting to run in 1D or D>3, sim won't run
`kbt`           | double        | 1             | The thermal energy of the system. Should **always** be 1 since sets energy scale 
`dt`            | double        | 0.1           | The size of the simulation timestep in units of MPCD time (cell size/sqrt(particle mass*thermal energy))
`simSteps`      | int           | 2000          | How many total timesteps to run the simulation for
`warmUp`        | int           | 0             | How much warmup time to run before starting the simulation. Note that during warmup no output is written to file
`rFrame`        | int           | 1             | Takes of net non-zero average momentum post initialisation to remove any drift. 1 = on, 0 = off
`zeroNetMom`    | int           | 0             | This substracts any excess momentum EVERY time step (generally have this off)
`galInv`        | int           | 1             | Whether to enable the random shift of particles to counter Gallilean variance of grid-based algorithm
`tsTech`        | int           | 0             | Which thermostat technique to use. 0 = off. Note some collsion operators have thermostats built in, see definitions.h for list
`rTech` OR `collOp`| int           | 2             | Which collision operator to use. See definitions.h for list. `collOp` is prioritised
`lc`            | int           | 0             | Liquid crystal mode. 2 = global S, 1 = local S, 0 = off
`tau`           | double        | 0.5           | Thermal relaxation time scale
`rotAng`        | double        | 1.570796      | This is the angle used in the original SRD collision operator
`fricCoef`      | double        | 1.0           | Friction coefficient for langevin thermostat
`mfpot`         | double        | 10            | Liquid crystal mean field potential in units of thermal energy
`noHI`          | int           | 0             | Enable no HI mode. 1 = on, 0 = off
`incomp`        | int           | 0             | Enable incompressibility correction. 1 = on, 0 = off
`multiphase`    | int           | 0             | Enable multiphase mode, applying multiphase interactions. 1 = on, 0 = off
`grav`          | array(double) | [0,0,0]       | Constant acceleration due to external force. **Must** be 3D
`mag`           | array(double) | [0,0,0]       | Constant external magnetic field. **Must** be 3D
`seed`          | int           | 0             | Seed for random number generator. 0 for pseudorandom seed. Set to -1 to load a checkpoint
`mdIn`          | string        | ""            | Path to the MD input file. This also acts as the switch for enabling MD --- If set to `""`, then MD is disabled, otherwise MD is enabled with the corresponding input file
`mdCoupleMode`  | int           | 1             | Coupling mode for MD. Only matters if `mdIn` is set. Set to 1 for MD particles to be treated as MPCD particles, 2 for MD particles to be treated as MD particles in the MPCD simulator
`stepsMD`       | int           | 20            | MD time steps per MPCD time step
`mfpLayerH`     | int           | 0             | Height above which MFP goes to 0. For simulating thin ordered films below disordered fluids. If 0, disable this functionality entirely.
`species`       | array(species)| 1 default spec  | An array of species objects.  See the [species table](#species-tag-table) for species tags
---             | ---           | ---           | ---
`debugOut`      | int           | 3             | Debug (verbosity) level. See definitions.h for list
`trajOut`       | int           | 0             | Detailed particle trajectories for every particle of species type given by `trajSpecOut`
`trajSpecOut`   | int           | 0             | Which number of species whose detailed trajectories to output
`coarseOut`     | int           | 0             | Coarse grain data (cell velocity, densities, density of each species) field
`flowOut`       | int           | 0             | Flow field averaged between output times
`velOut`        | int           | 0             | Instantaneous velocity field.
`avVelOut`      | int           | 0             | Total average MPCD velocity. System-averaged single value
`dirSOut`       | int           | 0             | Director and scalar order parameter fields
`qTensOut`      | int           | 0             | Q tensor field
`qkTensOut`     | int           | 0             | Reciprocal Q tensor field
`oriEnOut`      | int           | 0             | Orientational energy field
`colourOut`     | int           | 0             | Colour/ phi/ species-type field
`pressureOut`   | int           | 0             | Pressure field
`neighbourEnOut`| int           | 0             | Orientational energy from neighbours. System-averaged single value
`avSOut`        | int           | 0             | Total average scalar order parameter. System-averaged single value
`densSDOut`     | int           | 0             | SD of the number per cell. System-averaged single value
`enstrophyOut`  | int           | 0             | Enstrophy field
`histVelOut`    | int           | 0             | Velocity distribution
`histSpeedOut`  | int           | 0             | Speed distribution
`histVortOut`   | int           | 0             | Vorticity distribution
`histEnsOut`    | int           | 0             | Enstrophy distribution
`histDirOut`    | int           | 0             | Director distribution
`histSOut`      | int           | 0             | Scalar order parameter distribution
`histNOut`      | int           | 0             | Number per cell distribution
`solidTrajOut`  | int           | 0             | Solid BC trajectories
`topoFieldOut`  | int           | 0             | Topological charge field
`defectsOut`    | int           | 0             | Defect positions and orientations
`disclinOut`    | int           | 0             | Disclination tensor field
`energyOut`     | int           | 0             | System energy field
`velCorrOut`    | int           | 0             | Velocity autocorrelation (radial function)
`dirCorrOut`    | int           | 0             | Director autocorrelation (radial function)
`vortCorrOut`   | int           | 0             | Vorticity autocorrelation (radial function)
`densCorrOut`   | int           | 0             | Densty autocorrelation (radial function)
`orderCorrOut`  | int           | 0             | Scalar order parameter autocorrelation (radial function)
`phaseCorrOut`  | int           | 0             | Binary fluid correlation (radial function)
`energySpecOut` | int           | 0             | Energy spectra  (radial function)
`enstrophySpecOut`| int           | 0             | Enstrophy spectra (radial function)
`binderOut`     | int           | 0             | Binder cumulant
`binderBin`     | int           | 0             | Binder cumulant bin size
`swimQOut`      | int           | 0             | Swimmer positions
`swimOOut`      | int           | 0             | Swimmer orientations
`swimROut` OR `swimRTOut`| int           | 0             | Swimmer run/ tumble. `swimRTOut` is prioritised
`synopsisOut`   | int           | 1             | Synopsis output. Highly recommended to be on. 1 = on, 0 = off
`checkpointOut` | int           | 0             | Simulation checkpointing. Unlike all other `out` tags this **will** work during the warmup stage, and will overwrite the file every time it runs. Note that this just controls dump rate --- To actually load a saved checkpoint, set `seed` to -1
---             | ---           | ---           | ---
`BC`            | array(BC)     | PBCs around domain | The array of boundary objects. See the [BC table](#bc-tag-table) for BC tags
---             | ---           | ---                | ---
`typeSwim`      | int           | 2                  | Swimmer type. Default is a dumbell swimmer with excluded volume interactions
`nSwim`         | int           | 0                  | Swimmer population
`qDistSwim`     | int           | 0                  | Initial distribution of swimmer positions
`oDistSwim`     | int           | 0                  | Initial distribution of swimmer orientations
`headMSwim`     | int           | 20                 | Mass of head monomer
`midMSwim`      | int           | 20                 | Mass of middle monomer
`hspIdSwim`     | int           | 1                  | Multiphase fluid particle type of the head monomer in the swimmer
`mspIdSwim`     | int           | 1                  | Multiphase fluid particle type of the middle monomer in the swimmer
`fsSwim`        | double        | 20                 | Magnitude of swimmer propulsion force
`dsSwim`        | double        | 1                  | Dipole strength
`tsSwim`        | double        | 0                  | Magnitude of swimmer torque
`sizeShrinkSwim`| double        | 0.1                | How much do LJ sigma ro shrink when tumbling
`springShrinkSwim`| double        | 0.1                | How much the spring constant is shrunk when tumbling
`kSwim`         | double        | 30                 | Spring constant
`roSwim`        | double        | 4                  | Spring seperation
`sigSwim`       | double        | 4                  | Diameter approx sigma
`epsSwim`       | double        | 1                  | Interaction energy
`runTSwim`      | double        | 0                  | Average run time in units of MPCD timesteps dt (iterations, not MPCD time units)
`tumTSwim`      | double        | 0                  | Average tumble time in units of MPCD timesteps dt (iterations, not MPCD time units)
`shrTSwim`      | double        | 2                  | Set time to shrink/ extend in units of MPCD timesteps dt (iterations, not MPCD time units)
`magMomSwim`    | double        | 1                  | Magnetic moment/ strength
`fixDistSwim`   | double        | 0                  | The fixed distance from the wall for `DUMBELL_NEARWALL` mode

#### Overrides          {#overrides}
Override Tag    | Type  | Override param | Description
---             |-------|--------------| ---
`domainWalls`   | int   | `BC`         | This override will add extra BCs to the simulation, on top of the declared ones, on the domain walls. If set to 1, it places PBCs, and if set to 0 it places solid walls. 
`checkpointTimerOut`| float| `checkpointOut`| This override enables checkpointing, but puts it on a timer. It will checkpoint every X **hours**, where X is specified by this parameter. 

### Species Tag Table           {#species-tag-table}
Tag             | Type          | Default Value | Description
---             | ---           | ---           | ---
`mass`          | double        | 1             | Mass of this species of particles. Should **always** be 1 for at least one species since sets mass scale 
`pop`           | int           | 18000         | Number of particles of this species
`qDist`         | int           | 0             | Initial positional distribution function for the species. See definitions.h for a list
`vDist`         | int           | 0             | Initial velocity distribution function for the species. See definitions.h for a list
`oDist`         | int           | 2             | Initial orientation distribution function for the species. See definitions.h for a list
`interMatr`     | array(double) | [0,...]       | Interaction matrix for this species, against other species. **Must** be of the same length as the number of species. Default will autopopulate with 0s
`rfc`           | double        | 0.01          | Nematogen rotational friction coefficient
`len`           | double        | 0.007         | Effective rod length (specifically for solid boundary interactions)
`tumble`        | double        | 2             | Tumbling parameter
`shearSusc`     | double        | 0.5           | Shear susceptibility
`magnSusc`      | double        | 0.001         | Magnetic susceptibility
`act`           | double        | 0.05          | Species activity
`sigWidth`      | double        | 1.0           | Sigmoid width for activity fall-off, specifically for collision operator #20. **Cannot be 0**
`sigPos`        | double        | `sigWidth`    | Sigmoid position for activity fall-off, specifically for collision operator #20
`minActRatio`   | double        | 0             | Minimum ratio of particles (of mean density) of this species to allow activity to be calculated in this cell. If 0, then this is ignored
`damp`          | double        | 0             | Damping friction to kill hydrodynamics. Between 0 and 1

#### Species Overrides          {#species-overrides}
Override Tag    | Type          | Override param| Description
---             | ---           | ---           | ---
`dens`          | double        | `pop`         | This override will set a given species population to correspond to a cell density. This sets pop to the volume of the **whole system** (i.e. volume of the domain, not excluding any excluded regions due to BCs) times this given value

### BC Tag Table            {#bc-tag-table}

Unlike the previous tags, if you declare a new BC object then there are some tags that **must** be declared. 
These are noted as NECESSARY in the default value column.

Tag             | Type          | Default Value | Description
---             | ---           | ---           | ---
`colType`       | int           | 1             | Which BC type do you want to use. See definitions.h for a list
`phantom`       | int           | 0             | Use phantom particles if 1, 0 otherwise
`E`             | double        | -1            | Coefficient of restitution
`Q`             | array(double) | [0,0,0]       | Position of the boundary, **must** be 3D
`V`             | array(double) | [0,0,0]       | Velocity of the boundary, **must** be 3D
`O`             | array(double) | [0,0,0]       | Orientation of the boundary about (x,y,z), **must** be 3D
`L`             | array(double) | [0,0,0]       | Angular velocity of the boundary, **must** be 3D
`G`             | array(double) | [0,0,0]       | External acceleration (ie, due to gravity) of the boundary, **must** be 3D
`aInv`          | array(double) | NECESSARY     | Sets the geometry of the surface. Principal semi-axes of the ellipsoid (see SRDClass for explanation). **Must** be 3D. If you declare a BC then this **must** be given, or the simulation will not run
`rotSym`        | array(double) | [4,4]         | Sets rotational symmetry of shapes. **Must** be of form (a,b)
`abs`           | int           | 0             | Flags if each term should be absolute only. 1 = yes, 0 = no
`P`             | array(double) | NECESSARY     | Sets the geometry of the surface. **Must** be of form (a,b,c,d). See SRDClass for explanation. If you declare a BC then this **must** be given, or the simulation will not run
`R`             | double        | NECESSARY     | Sets the geometry of the surface. See SRDClass for explanation. If you declare a BC then this **must** be given, or the simulation will not run
`DN`            | double        | NECESSARY     | Displacement of the particle in the normal direction. If you declare a BC then this **must** be given, or the simulation will not run
`DT`            | double        | 0             | Displacement of the particle in the tangential direction
`DVN`           | double        | 0             | Add velocity in the normal direction on contact
`DVT`           | double        | 0             | Add velocity in the tangential direction on contact
`DVxyz`         | array(double) | [0,0,0]       | Add velocity in the cartesian directions on contact. **Must** be 3D
`MVN`           | double        | NECESSARY     | Multiplies the velocity in the normal direction on contact. If you declare a BC then this **must** be given, or the simulation will not run
`MVT`           | double        | NECESSARY     | Multiplies the velocity in the tangential direction on contact. If you declare a BC then this **must** be given, or the simulation will not run
`MUN`           | double        | 1             | Multiplies the orientation in the normal direction on contact
`MUT`           | double        | 1             | Multiplies the orientation in the tangential direction on contact
`MUxyz`         | array(double) | [1,1,1]       | Multiplies the orientation in the cartesian direction on contact. **Must** be 3D
`DUxyz`         | array(double) | [0,0,0]       | Add the orientation in the cartesian direction on contact. **Must** be 3D
`kbt`           | double        | 1             | Temperature of the wall
`dsplc`         | int           | 0             | Whether the wall can displace/ is mobile. 0 = no, 1 = yes
`inv`           | int           | 0             | Whether to invert the bc (ie, multiply the A's by -1). 0 = no, 1 = yes
`mass`          | double        | 1             | Mass of the wall in MPCD units. Should be the same density as the fluid if its displaceable
`wavy`          | array(double) | [0,0,0]       | Generalized amplitudes and frequencies for wavy-walls. **Must** be 3D

#### BC Overrides           {#bc-overrides}
Override Tag    | Type          | Override param| Description
---             | ---           | ---           | ---
`homeotropic`   | int           | `MUN`, `MUT`  | Setting this override with a value of 1 will give the wall homeotropic anchoring conditions
`planar`        | int           | `MUN`, `MUT`  | Setting this override with a value of 1 will give the wall planar anchoring conditions

## Example          {#example}
Below is an example JSON input file that uses all tags but sets them to the defaults.
As a reminder, if you wish to use the default value for a tag, you can leave it out of the JSON file, so this is primarily for illustrative purposes.

```json
{
    "domain":           [30, 30],
    "kbt":              1,
    "dt":               0.1,
    "simSteps":         2000,
    "warmUp":           0,
    "rFrame":           1,
    "zeroNetMom":       0,
    "galInv":           1,
    "tsTech":           0,
    "collOp":            2,
    "lc":               0,
    "tau":              0.5,
    "rotAng":           1.570796,
    "fricCoef":         1.0,
    "mfpot":            10,
    "grav":             [0, 0, 0],
    "mag":              [0, 0, 0],
    "seed":             0,
    "mdIn":             "",
    "stepsMD":          20,
    "mfpLayerH":        0,
    "species":
    [
        {
            "mass":         1,
            "pop":          18000,
            "qDist":        0,
            "vDist":        0,
            "oDist":        2,
            "interMatr":    [0],
            "rfc":          0.01,
            "len":          0.007,
            "tumble":       2,
            "shearSusc":    0.5,
            "magnSusc":     0.001,
            "act":          0.05,
            "damp":         0,
            "sigWidth":     1.0,
            "sigPos":       0.0,
            "minActRatio":  0.0
        }
    ],
    "debugOut":         3,
    "trajOut":          0,
    "trajSpecOut":      0,
    "coarseOut":        0,
    "flowOut":          0,
    "velOut":           0,
    "avVelOut":         0,
    "dirSOut":          0,
    "qTensOut":         0,
    "qkTensOut":        0,
    "oriEnOut":         0,
    "colourOut":        0,
    "pressureOut":      0,
    "neighbourEnOut":   0,
    "avSOut":           0,
    "densSDOut":        0,
    "enstrophyOut":     0,
    "histVelOut":       0,
    "histSpeedOut":     0,
    "histVortOut":      0,
    "histEnsOut":       0,
    "histDirOut":       0,
    "histSOut":         0,
    "histNOut":         0,
    "solidTrajOut":     0,
    "topoFieldOut":     0,
    "defectsOut":       0,
    "disclinOut":       0,
    "energyOut":        0,
    "velCorrOut":       0,
    "dirCorrOut":       0,
    "vortCorrOut":      0,
    "densCorrOut":      0,
    "orderCorrOut":     0,
    "phaseCorrOut":     0,
    "energySpecOut":    0,
    "enstrophySpecOut": 0,
    "binderOut":        0,
    "binderBin":        0,
    "swimQOut":         0,
    "swimOOut":         0,
    "swimRTOut":         0,
    "synopsisOut":      1,
    "checkpointOut":    0,
    "BC":
    [
        {
            "colType":      1,
            "phantom":      0,
            "E":            -1,
            "Q":            [0, 0, 0],
            "V":            [0, 0, 0],
            "O":            [0, 0, 0],
            "L":            [0, 0, 0],
            "G":            [0, 0, 0],
            "aInv":         [1, 0, 0],
            "rotSym":       [4, 4],
            "abs":          0,
            "P":            [1, 1, 1, 1],
            "R":            0,
            "DN":           30,
            "DT":           0,
            "DVN":          0,
            "DVT":          0,
            "DVxyz":        [0, 0, 0],
            "MVN":          1,
            "MVT":          1,
            "MUN":          1,
            "MUT":          1,
            "MUxyz":        [1, 1, 1],
            "DUxyz":        [0, 0, 0],
            "kbt":          1,
            "dsplc":        0,
            "inv":          0,
            "mass":         1
        },
        {
            "colType":      1,
            "phantom":      0,
            "E":            -1,
            "Q":            [0, 0, 0],
            "V":            [0, 0, 0],
            "O":            [0, 0, 0],
            "L":            [0, 0, 0],
            "G":            [0, 0, 0],
            "aInv":         [-1, 0, 0],
            "rotSym":       [4, 4],
            "abs":          0,
            "P":            [1, 1, 1, 1],
            "R":            -30,
            "DN":           30,
            "DT":           0,
            "DVN":          0,
            "DVT":          0,
            "DVxyz":        [0, 0, 0],
            "MVN":          1,
            "MVT":          1,
            "MUN":          1,
            "MUT":          1,
            "MUxyz":        [1, 1, 1],
            "DUxyz":        [0, 0, 0],
            "kbt":          1,
            "dsplc":        0,
            "inv":          0,
            "mass":         1
        },
        {
            "colType":      1,
            "phantom":      0,
            "E":            -1,
            "Q":            [0, 0, 0],
            "V":            [0, 0, 0],
            "O":            [0, 0, 0],
            "L":            [0, 0, 0],
            "G":            [0, 0, 0],
            "aInv":         [0, 1, 0],
            "rotSym":       [4, 4],
            "abs":          0,
            "P":            [1, 1, 1, 1],
            "R":            0,
            "DN":           30,
            "DT":           0,
            "DVN":          0,
            "DVT":          0,
            "DVxyz":        [0, 0, 0],
            "MVN":          1,
            "MVT":          1,
            "MUN":          1,
            "MUT":          1,
            "MUxyz":        [1, 1, 1],
            "DUxyz":        [0, 0, 0],
            "kbt":          1,
            "dsplc":        0,
            "inv":          0,
            "mass":         1
        },
        {
            "colType":      1,
            "phantom":      0,
            "E":            -1,
            "Q":            [0, 0, 0],
            "V":            [0, 0, 0],
            "O":            [0, 0, 0],
            "L":            [0, 0, 0],
            "G":            [0, 0, 0],
            "aInv":         [0, -1, 0],
            "rotSym":       [4, 4],
            "abs":          0,
            "P":            [1, 1, 1, 1],
            "R":            -30,
            "DN":           30,
            "DT":           0,
            "DVN":          0,
            "DVT":          0,
            "DVxyz":        [0, 0, 0],
            "MVN":          1,
            "MVT":          1,
            "MUN":          1,
            "MUT":          1,
            "MUxyz":        [1, 1, 1],
            "DUxyz":        [0, 0, 0],
            "kbt":          1,
            "dsplc":        0,
            "inv":          0,
            "mass":         1
        }
    ],
    "typeSwim":         2,
    "nSwim":            0,
    "qDistSwim":        0,
    "oDistSwim":        0,
    "headMSwim":        20,
    "midMSwim":         20,
    "hspIdSwim":        1,
    "mspIdSwim":        1,
    "fsSwim":           20,
    "dsSwim":           1,
    "tsSwim":           0,
    "sizeShrinkSwim":   0.1,
    "springShrinkSwim": 0.1,
    "kSwim":            30,
    "roSwim":           4,
    "sigSwim":          4,
    "epsSwim":          1,
    "runTSwim":         0,
    "tumTSwim":         0,
    "shrTSwim":         2,
    "magMomSwim":       1,
    "fixDistSwim":      0
}
```

## Guide to Adding New Input Parameters         {#guide-to-adding-new-input-parameters}
Adding new parameters needs to be done within /mpcd/subroutines/read.c . 

First, ensure that the parameter you wish to "write to" exists in the code somewhere: 
This can be in the inputList struct or a global, but it needs to be accessible and passed to/ from read.c . 

Parameters are read from the JSON in the method `readJson()`.
The new parameter needs to be read **after** `initLL()` is called, but **before** `verifyJson()` is called.
As the order of parameters is irrelevent with the JSON, you can just append the new parameter to the end of the list.

The API for reading "tags" from the JSON makes use of three methods, one for varying tag types:
- `int getJObjInt(cJSON *cJSONRoot, const char* jsonTag, int d, linkedList *head)`:
- `double getJObjDou(cJSON *cJSONRoot, const char* jsonTag, double d, linkedList *head)`
- `void getJObjStr(cJSON *cJSONRoot, const char* jsonTag, const char* d, char **toReturn, linkedList *head)`

Each method has nearly identical parameters:
- `cJSONRoot`: the root of the JSON object, should be `jObj`
- `jsonTag`: the tag to search for
- `d`: the default value to return if the tag is not found
- `head`: the linked list to append the new parameter to, should be `jsonTagList`
Strings work a bit differently and you need to pass a reference to where you want to store the string in arg `toReturn`.

Setting up arrays is more complicated, a decent example is that of system dimensionality & size (controlled by tag `domain`).
Arrays of custom objects, such as species or BCs, are very complicated and not recommended to be done unless you're sure about what you're doing.
Overrides are, by definition, hacks, so how you implement them is very fast and loose.
