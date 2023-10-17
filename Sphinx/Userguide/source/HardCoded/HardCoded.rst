Hard-coded options
##################

While MPCD is almost entirely configureable from its input files, there are a several options that are hard coded. 
These hard-coded options are in the form of preprocessor defines (such as ``# define DBG``), which are read by the compiler at compile time and are within the ``mpcd/headers/definitions.h`` file. 
The majority of these options relate to the operation of the simulator, rather than the simulation of the physics, and are detailed on this page.

There are two variants of these options, which are configured in different ways:

* Toggle options, which are either on or off, are set or unset by commenting them out. This is done by placing a ``//`` at the start of the line, which prevents it from being read by the compiler. 
* Numbered options, which have either a decimal or integer number following their definition (ie, ``# define NN 624``). These are set by changing the number to the desired value.

.. warning:: 
    Changing these options will require the simulator to be re-compiled to enable them. 

Random number generation options
********************************

MPCD uses two types of random number generators: 

1. A Xoshiro128++ generator, derived from `this example <https://prng.di.unimi.it/xoshiro128plusplus.c>`_.
2. A Mersenne-Twister, taken from `this source <http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/CODES/mt19937ar.c>`_.

By default, the Xoroshiro generator is used, as it is significantly faster while still retaining a sufficiently large period for MPCD simulations.
To enable the Mersenne-Twister method, uncomment the ``// # define RNG_MERSENNE`` line in ``mpcd/headers/definitions.h``.

Mersenne-Twister options
------------------------

.. note:: 
    These options only take effect if the Mersenne-Twister random number generator is enabled.

Of the two random number generator options, only the `Mersenne-Twister <http://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/emt.html>`_ has configureable options.

.. warning:: 
    These options are not designed to be configurable by a typical user. 
    We strongly discourage altering these options unless you are familiar with the Mersenne-Twister algorithm, and the implications of changing these options.

    This is especially true due to MPCD's dependence on random number generation.

These options include:

* ``# define NN 624``: 
    The size of the state vector.
* ``# define MM 397``: 
    A middle element of the state vector. 
    Used in the process of shifting the state vector between states, such that generating a new integer performs ``NN-MM`` iterations.
* ``# define MATRIX_A 0x9908b0dfUL``: 
    A constant vector used in the process of shifting the state vector.
* ``# define UPPER_MASK 0x80000000UL``: 
    The most significant w-r bits of the state vector.
* ``# define LOWER_MASK 0x7fffffffUL``: 
    The least significant r bits of the state vector.

Debug modes/verbosity
**********************
All non-essential and non-error verbosity in MPCD is enabled by enabling "debug mode".
This occurs if the ``# define DBG`` line is uncommented in ``mpcd/headers/definitions.h``.
For production-level simulation, such as on clusters, this should usually be disabled by commenting out that line prior to recompiling. 
However it can be helpful to keep it enabled when testing and or debugging simulation input scripts.

If debug mode is enabled, then the verbosity level of the MPCD simulator is controlled by the integer assigned to the input ``.json`` option called ``"debugOut"``. 
Possible values for this are also defined in ``mpcd/headers/definitions.h``, and generally a given value will also output all messages of lower values.
Commonly used debug levels include:

0. Run --- Outputs nothing but the sim start and sim end.
1. Warn --- Outputs warnings.
2. Init --- Outputs initialisation messages.
3. Steps --- Outputs iteration steps.
4. Title --- Gives a title to some major methods and outputs them when called.

The following debug levels are less frequently used, and are most helpful when debugging the code itself:

5. Wait --- Unused. 
6. Therm --- Outputs thermostat related properties, such as system temperature. Only outputs at this verbosity level.
7. Hist --- Outputs whenever cell quantities are binned for output as a histogram.
8. BC Count --- Outputs number of particles that bounce between many BCs in a single timestep. Only outputs at this verbosity level.
9. MPC BC --- States operations occurred during MPC particle to BC interactions. Only outputs at this verbosity level.
10. BC MPC --- State operations occurred during BC to MPC particle interactions. Only outputs at this verbosity level.
11. BC BC --- States operations occurred during BC to BC interactions. Only outputs at this verbosity level.
12. BC Maximal --- States a summary of all BC collision related operations that occur during a particle-BC collision. Only outputs at this verbosity level.
13. LC Collision --- Outputs intermediate states that occur during the liquid crystal collision algorithm. Only outputs at this verbosity level.
14. BC Ori --- Outputs intermediate states that occur during the BC re-orientation algorithm. Only outputs at this verbosity level.
15. Jeffery's Torque --- Outputs intermediate states that occur during the Jeffery's torque algorithm. Only outputs at this verbosity level.
16. Magnetic --- Outputs intermediary debug information when magnetic fields are applied to liquid crystals and swimmers. Only outputs at this verbosity level.
17. Binary --- Unused.
18. Swimmer --- Outputs intermediary debug information when swimmers are present. Only outputs at this verbosity level.
19. Run Tumble --- Outputs intermediary debug information when swimmers perform the run-tumble algorithm. Only outputs at this verbosity level.
20. Escape --- Legacy: Used to be used for debugging particle escapes from the domain. Currently identical to LC Collision.
21. Swimmer Details --- Outputs high level detail of all intermediary debug information when swimmers are present. Only outputs at this verbosity level.
22. Swimmer Torque --- Outputs intermediary debug information on swimmer torques. Only outputs at this verbosity level.
23. Incompressibility --- Outputs intermediary debug information for the incompressibility algorithm. Only outputs at this verbosity level.

Floating-point error handling
*****************************
.. warning:: 
    This feature only works on Linux (and WSL) systems. It does not function on Mac OS.

.. note:: 
    For normal operation, this should be kept disabled.

Primarily for debugging purposes, it is possible to set MPCD to halt execution upon receipt of a ``SIGFPE`` hardware signal. 
This hardware signal is thrown by the CPU when a floating-point error occurs, such as a division by zero or an overflow.
Typically, such errors will cause the program to be flooded with ``NaN`` values, which will either cause the simulation to crash, or propagate through the system to give garbled output.

To enable this feature, uncomment the ``// #define FPE`` line in ``mpcd/headers/definitions.h``.
This is most useful when running through a debugger, such as ``gdb``, to identify the source of the error.

File flush
***************
.. note:: 
    There are few drawbacks to enabling this, and it is recommended to keep this on. 
    There is a slight slowdown associated with writting to file before the buffer is full. 

By default, codes written in the *C* programming language do not immediatly write to terminal or file.
Instead, they write to a buffer controlled by the operating system, which will periodically flush the buffer. 
It is possible to force the operating system to flush the buffer immediatly, which in MPCD is controlled by the ``# define FFLSH`` line in ``mpcd/headers/definitions.h``.
To disable this option, comment out the definition. 