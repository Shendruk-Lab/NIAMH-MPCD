# Guide to MPCD Outputs

## Contents
1. [Introduction](#introduction-out)
2. [Output Tables](#output-tables)
    - [System Information](#system-information)
    - [Scalar Outputs](#scalar-outputs)
    - [Trajectory Outputs](#trajectory-outputs)
    - [Field Outputs](#field-outputs)
    - [Histograms](#histograms)
    - [Correlation Functions](#correlation-functions)
    - [Swimmers](#swimmers)


## Introduction         {#introduction-out}

This is a comprehensive reference guide to MPCD output data files and what one might find in them. 

For MPCD to output data, the corresponding flag must be set to non-zero in the `input.json`. These flags are listed here in the input tag column.

With the exception of `synopsis.dat`, the values in the input are the frequency (in units of MPCD time steps) with which the relevant data will be written to the outputs for the corresponding tag. 
All outputs are written to `.dat` files, which are uncompressed text files.

## Output Tables            {#output-tables}

### System Information          {#system-information}

| Input tag       | Output file      | Description                                                                                                                                                                                                      |
|-----------------|------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `synopsisOut`   | `synopsis.dat`   | Synopsis output. Highly recommended to be on. 1 = on, 0 = off. Provides a detailed documentation of simulation information as the code runs.                                                                     |
| `checkpointOut` | `checkpoint.dat` | Simulation checkpointing. Provides data for re-populating the system and restarting another simulation from this point. Can be set to be on a physical timer, rather than based on time-step, using an override. |

### Scalar Outputs          {#scalar-outputs}

| Input tag        | Output file          | Description                                                             | Outputs                       | Column Headers   |
|------------------|----------------------|-------------------------------------------------------------------------|-------------------------------|------------------|
| `neighbourEnOut` | `enneighbours.dat`   | Orientational energy from neighbours. System-averaged single value      | Time                          | `t`              |
|                  |                      |                                                                         | Nematic energy                | `tMPC_nem`       |
| `avSOut`         | `avS.dat`            | Total average scalar order parameter. System-averaged single value      | Time                          | `t`              |
|                  |                      |                                                                         | Scalar order parameter        | `S`              |
|                  |                      |                                                                         | 4th moment of Scalar order    | `S4`             |
|                  |                      |                                                                         | Director orientation          | `nX`,`nY`,`nZ`   |
| `densSDOut`      | `densSTD.dat `       | Standard deviation of the number per cell. System-averaged single value | Time                          | `t`              |
|                  |                      |                                                                         | Standard Deviation of density | `densSTD`        |
| `binderOut`      | `binderCumulant.dat` | Binder cumulant, bin size must be set in `binderBin`                    | Time                          | `t`              |
|                  |                      |                                                                         | Binder cumulant               | `BinderCumulant` |

### Trajectory Outputs          {#trajectory-outputs}
| Input tag          | Output file              | Description                                                                              | Outputs                        | Column Headers                                                    |
|--------------------|--------------------------|------------------------------------------------------------------------------------------|--------------------------------|-------------------------------------------------------------------|
| `trajOut`          | `detailedSP0.dat`        | Detailed particle trajectories for every particle of species type given by `trajSpecOut` | Time                           | `t`                                                               |
|                    |                          |                                                                                          | X, Y, Z co-ordinates           | `QX`,`QY`,`QZ`                                                    |
|                    |                          |                                                                                          | Cell velocities                | `VX`,`VY`,`VZ`                                                    |
|                    |                          |                                                                                          | Speed                          | \|`V`\|                                                  |
|                    |                          |                                                                                          | Species velocities             | `UX`,`UY`,`UZ`                                                    |

### Field Outputs           {#field-outputs}

| Input tag          | Output file              | Description                                                                              | Outputs                        | Column Headers                                                    |
|--------------------|--------------------------|------------------------------------------------------------------------------------------|--------------------------------|-------------------------------------------------------------------|
| `coarseOut`        | `coarsegrain.dat`        | Coarse grain data (cell velocity, densities, density of each species) field              | Time                           | `t`                                                               |
|                    |                          |                                                                                          | X, Y, Z co-ordinates           | `QX`,`QY`,`QZ`                                                    |
|                    |                          |                                                                                          | Centre of mass Velocities      | `VcmX`,`VcmY`,`VcmZ`                                              |
|                    |                          |                                                                                          | Cell population                | `POP`                                                             |
|                    |                          |                                                                                          | Species                        | `SP0`                                                             |
| `flowOut`          | `flowfield.dat`          | Flow field averaged between output times                                                 | Time                           | `t`                                                               |
|                    |                          |                                                                                          | X, Y, Z co-ordinates           | `QX`,`QY`,`QZ`                                                    |
|                    |                          |                                                                                          | Centre of mass Velocities      | `VcmX`,`VcmY`,`VcmZ`                                              |
| `avVelOut`         | `avVel.dat`              | Total average MPCD velocity. System-averaged single value and velocity gradient tensor   | Time                           | `t`                                                               |
|                    |                          |                                                                                          | Centre of mass Velocities      | `VcmX`,`VcmY`,`VcmZ`                                              |
|                    |                          |                                                                                          | Thermal energy                 | `KBT`                                                             |
|                    |                          |                                                                                          | X derivatives of velocity      | `dVXX`,`dVYX`,`dVZX`                                              |
|                    |                          |                                                                                          | Y derivatives of velocity      | `dVXY`,`dVYY`,`dVZY`                                              |
|                    |                          |                                                                                          | Z derivatives of velocity      | `dVXZ`,`dVYZ`,`dVZZ`                                              |
| `dirSOut`          | `directorfield.dat`      | Director and scalar order parameter fields                                               | Time                           | `t`                                                               |
|                    |                          |                                                                                          | X, Y, Z co-ordinates           | `QX`,`QY`,`QZ`                                                    |
|                    |                          |                                                                                          | Director orientation           | `NX`,`NY`,`NZ`                                                    |
| `qTensOut`         | `ordertens.dat`          | Q tensor field                                                                           | X, Y, Z co-ordinates           | `X`,`Y`,`Z`                                                       |
|                    |                          |                                                                                          | Q tensor components            | `QXX`,`QXY`...`QZY`,`QZZ`                                         |
| `qkTensOut`        | `recipOrder.dat`         | Reciprocal Q tensor field                                                                | Time                           | `t`                                                               |
|                    |                          |                                                                                          | Wave Vectors                   | `K123_X`,`K123_Y`,`K123_Z`                                        |
|                    |                          |                                                                                          | Squared-modulus of Q tensor    | \|`QXX`\|`2`,\|`QXY`\|`2`...,\|`QZZ`\|`2` |
| `oriEnOut`         | `enfield.dat`            | Orientational energy field                                                               | X, Y, Z co-ordinates           | `QX`,`QY`,`QZ`                                                    |
|                    |                          |                                                                                          | Kinetic energy                 | `MPC_kin`                                                         |
|                    |                          |                                                                                          | Nematic energy                 | `tMPC_nem`                                                        |
| `colourOut`        | `multiphase.dat`         | Colour/ phi/ species-type field                                                          | Time                           | `t`                                                               |
|                    |                          |                                                                                          | X, Y, Z co-ordinates           | `QX`,`QY`,`QZ`                                                    |
|                    |                          |                                                                                          | multiphase                     | `N_0`,`N_1`, ...                                                             |
| `pressureOut`      | `pressure.dat`           | Pressure field                                                                           | Time                           | `t`                                                               |
|                    |                          |                                                                                          | X, Y, Z co-ordinates           | `QX`,`QY`,`QZ`                                                    |
|                    |                          |                                                                                          | Pressure tensor components     | `Pxx`,`Pxy`...`Pzy`,`Pzz`                                         |
| `enstrophyOut`     | `avEnstrophy.dat`        | Enstrophy field                                                                          | Time                           | `t`                                                               |
|                    |                          |                                                                                          | Average enstrophy              | `enstrophy`                                                       |
| `topoFieldOut`     | `topochargefield.dat`    | Topological charge field and angles                                                      | Time                           | `t`                                                               |
|                    |                          |                                                                                          | X, Y, Z co-ordinates           | `QX`,`QY`,`QZ`                                                    |
|                    |                          |                                                                                          | Topological charge             | `charge`                                                          |
|                    |                          |                                                                                          | Angle of defect                | `angle`                                                           |
| `energyOut`        | `energy.dat`             | System energy field                                                                      | Time                           | `t`                                                               |
|                    |                          |                                                                                          | Kinetic energy                 | `MPC_kin`                                                         |
|                    |                          |                                                                                          | Nematic energy                 | `MPC_nem`                                                         |
|                    |                          |                                                                                          | Boundary kinetic energy        | `BC_kin`                                                          |
|                    |                          |                                                                                          | Boundary rotational energy     | `BC_rot`                                                          |
|                    |                          |                                                                                          | Total energy                   | `Total`                                                           |
|                    |                          |                                                                                          | Thermal energy                 | `KBT`                                                             |
| `energySpecOut`    | `energySpectrum.dat`     | Energy spectra  (radial function)                                                        | Time                           | `t`                                                               |
|                    |                          |                                                                                          | Wave number                    | `k`                                                               |
|                    |                          |                                                                                          | Energy                         | `E`                                                               |
| `defectsOut`       | `defects.dat`            | Defect positions and orientations                                                        | Time                           | `t`                                                               |
|                    |                          |                                                                                          | Number of defects at `t`       | `numDefects`                                                      |
|                    |                          |                                                                                          | X, Y co-ordinates              | `QX`,`QY`,`QZ`                                                    |
|                    |                          |                                                                                          | Topological charge             | `charge`                                                          |
|                    |                          |                                                                                          | Angle of defect                | `angle`                                                           |
| `disclinOut`       | `disclinTensorfield.dat` | Disclination tensor field                                                                | X, Y, Z co-ordinates           | `QX`,`QY`,`QZ`                                                    |
|                    |                          |                                                                                          | Disclination tensor components | `DXX`,`DXY`...`DZY`,`DZZ`                                         |
| `enstrophySpecOut` | `enstrophySpectrum.dat`  | Enstrophy spectra  (radial function)                                                     | Time                           | `t`                                                               |
|                    |                          |                                                                                          | Wave number                    | `k`                                                               |
|                    |                          |                                                                                          | Enstrophy                      | `Omega`                                                           |

### Histograms          {#histograms}

| Input tag      | Output file         | Description                                                             | Outputs                    | Column Headers  |
|----------------|---------------------|-------------------------------------------------------------------------|----------------------------|-----------------|
| `histVelOut`   | `histVel.dat`       | Velocity probability distribution in x, y, and z directions             | Time                       | `t`             |
|                |                     |                                                                         | Bin velocity               | `V`             |
|                |                     |                                                                         | Bin probability            | `PX`,`PY`,`PZ`  |
| `histSpeedOut` | `distSpeed.dat`     | Speed probability distribution                                          | Time                       | `t`             |
|                |                     |                                                                         | Bin speeds                 | \|`V`\| |
|                |                     |                                                                         | Bin probability            | `P`             |
| `histVortOut`  | `distVort.dat`      | Vorticity probability distribution in x, y, and z directions            | Time                       | `t`             |
|                |                     |                                                                         | Bin vorticity              | `W`             |
|                |                     |                                                                         | Bin probability            | `PX`,`PY`,`PZ`  |
| `histEnsOut`   | `distEnstrophy.dat` | Enstrophy probability distribution                                      | Time                       | `t`             |
|                |                     |                                                                         | Bin enstrophy              | \|`w`\| |
|                |                     |                                                                         | Bin probability            | `P`             |
| `histDirOut`   | `distDir.dat`       | Director orientation probability distribution in x, y, and z directions | Time                       | `t`             |
|                |                     |                                                                         | Bin orientation            | `n`             |
|                |                     |                                                                         | Bin probability            | `PX`,`PY`,`PZ`  |
| `histSOut`     | `distS.dat`         | Scalar order parameter probability distribution                         | Time                       | `t`             |
|                |                     |                                                                         | Bin scalar order parameter | `S`             |
|                |                     |                                                                         | Bin probability            | `P`             |
| `histNOut`     | `distN.dat`         | Number per cell probability distribution                                | Time                       | `t`             |
|                |                     |                                                                         | Bin density                | `stdN`          |
|                |                     |                                                                         | Bin probability            | `P`             |

### Correlation Functions           {#correlation-functions}

| Input tag      | Output file          | Description                                              | Outputs           | Column Headers |
|----------------|----------------------|----------------------------------------------------------|-------------------|----------------|
| `velCorrOut`   | `corrVelVel.dat`     | Velocity autocorrelation (radial function)               | Time              | `t`            |
|                |                      |                                                          | Separation        | `dr`           |
|                |                      |                                                          | Correlation value | `C`            |
| `dirCorrOut`   | `corrDirDir.dat`     | Director autocorrelation (radial function)               | Time              | `t`            |
|                |                      |                                                          | Separation        | `dr`           |
|                |                      |                                                          | Correlation value | `C`            |
| `vortCorrOut`  | `corrVortVort.dat`   | Vorticity autocorrelation (radial function)              | Time              | `t`            |
|                |                      |                                                          | Separation        | `dr`           |
|                |                      |                                                          | Correlation value | `C`            |
| `orderCorrOut` | `corrOrderOrder.dat` | Scalar order parameter autocorrelation (radial function) | Time              | `t`            |
|                |                      |                                                          | Separation        | `dr`           |
|                |                      |                                                          | Correlation value | `C`            |
| `phaseCorrOut` | `corrPhiPhi.dat`     | Binary fluid correlation (radial function)               | Time              | `t`            |
|                |                      |                                                          | Separation        | `dr`           |
|                |                      |                                                          | Correlation value | `C`            |

### Swimmers            {#swimmers}

| Input tag                 | Output file       | Description                                                               | Outputs                    | Column Headers    |
|---------------------------|-------------------|---------------------------------------------------------------------------|----------------------------|-------------------|
| `swimQOut`                | `swimmers.dat`    | Swimmer head(H) and middle(M) positions, velocities, and run-tumble phase | Time                       | `t`               |
|                           |                   |                                                                           | Head particle positions    | `HX`,`HY`,`HZ`    |
|                           |                   |                                                                           | Head particle velocities   | `HVX`,`HVY`,`HVZ` |
|                           |                   |                                                                           | Middle particle positions  | `MX`,`MY`,`MZ`    |
|                           |                   |                                                                           | Middle particle velocities | `MVX`,`MVY`,`MVZ` |
|                           |                   |                                                                           | Run-tumble phase           | `RTphase`         |
| `swimOOut`                | `swimmersOri.dat` | Swimmer orientations and run-tumble phase                                 | Time                       | `t`               |
|                           |                   |                                                                           | Swimmer orientation        | `nX`,`nY`,`nZ`    |
|                           |                   |                                                                           | Run-tumble phase           | `RTphase`         |
| `swimROut` OR `swimRTOut` | `runtumble.dat`   | Swimmer run/ tumble. `swimRTOut` is prioritised                           | Run-tumble phase           | `RTphase`         |
|                           |                   |                                                                           | Time interval              | `dt_cnt`          |
|                           |                   |                                                                           | Change in angle            | `dAng`            |
