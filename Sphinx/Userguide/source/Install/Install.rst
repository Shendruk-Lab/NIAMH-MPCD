.. _installAndRun:

Installation & run guide
##########################

This section explains how to install and run MPCD, makefile options, and arguments used.

.. note:: 
    Key information for advanced users:

    * MPCD has no external dependencies
    * Compilation is done with a standard makefile
    * Common UNIX C compilers successfully compile MPCD.

.. _installation:

Installation instructions
*************************
In this section, we describe how to install MPCD, and how to prepare to install it on a new system.
It is assumed that you have basic understanding of a `command-line terminal <https://ubuntu.com/tutorials/command-line-for-beginners#1-overview>`_. 

Depending on your operating system, please follow these instructions to install MPCD:

* If you are running on **Windows**, and do not have Windows Subsystem for Linux (WSL) installed, please follow the instructions in the :ref:`WSL` section.
* If you are running on **Mac OS**, and do not have dev tools such as make or a compiler installed, please follow the instructions in the :ref:`Mac` section.
* If you are running on native **Linux**, or **otherwise**, please follow the instructions in the :ref:`InstallCommon` section.

.. _WSL:

Windows - WSL setup
-------------------
.. warning:: 
    The Shendruk Lab does not support MPCD natively on Windows without WSL. 

    This is due to lacking support for ``make`` on Windows, and a requirement of the MD subsystem not compiling using MinGW or MSVC.
    Ultimately, this is likely to be fixed in a future update.

Running MPCD on Windows must be done through the Windows Subsystem for Linux (WSL) compatability layer.
This is a feature of Windows 10 and 11, allowing for native Linux calls to be made from the Windows kernel, effectively acting as a low-level VM.
More practically, it gives you a UNIX terminal within Windows.

The remainder of this section roughly follows the explanation given by Microsoft in their `official installation guide for WSL <https://learn.microsoft.com/en-us/windows/wsl/install>`_, and `tutorial for setting up a developer environment in WSL <https://learn.microsoft.com/en-us/windows/wsl/setup/environment#set-up-your-linux-username-and-password>`_.

For basic install, open up a command prompt or powershell terminal as administrator, and run

.. code-block:: powershell

    wsl --install

By default, this installs Ubuntu; however, alternative distributions can be installed via terminal commands, or via the Microsoft store app (search for WSL). 
You will be able to then open WSL by opening the start menu, and searching for the distribution you installed.
You can also directly open a WSL instance from the new Windows terminal.

.. This stupid thing is the only way to make this show up side-by-side
.. list-table:: 
    :header-rows: 0
    :widths: 50 50
    :align: center

    * - .. image:: ./WSLStartMenu.png
            :width: 95%
            :align: center
      - .. image:: ./WSLTerminal.png
            :width: 95%
            :align: center

Once open, WSL will ask for you to make a username and password for your Linux user account.
This is separate from your Windows account, and can be anything you like.
Note that the password will be required in the future for ``sudo`` commands, such as to install packages.

Note that the default WSL file system is located within your appdata folder and is not easily accessible from the Windows file system.
Instead, from within WSL, you can open a Windows file explorer window by running

.. code-block:: console

    explorer.exe .

This will open a file explorer window in the current directory of your WSL instance.

An alternative method is to save files in a Windows directory easily accessible from a root drive, such as ``C:\\``. 
These can be accessed within WSL by navigating to ``/mnt/c/`` (or whichever drive you saved the files to).

Before continuing, ensure git and a compiler (such as gcc) is installed on your WSL instance by running

.. code-block:: console

    git --version
    gcc --version

If you recieve a warning about either git or gcc not being installed after running these commands, then you will need to install them via your Linux distribution's package manager.
In Ubuntu, this will be done by running

.. code-block:: console

    sudo apt install git
    sudo apt install gcc

Once installed, you can continue with the :ref:`InstallCommon` section.

.. _Mac:

Mac OS X - Dev tools setup
--------------------------
To run MPCD on a Mac, the only requirements are a C compiler and the GNU make utility.
These can be installed via the XCode command line tools:
Open a terminal (under Applications/Utilities), and run

.. code-block:: console

    xcode-select --install

In the pop-up window, click Install and agree to the terms of service.
This will install both make and a C compiler.

An alternative method is to install using homebrew, a package manager for Mac OS.
To install make and a C compiler using homebrew, run

.. code-block:: console

    brew install make gcc

.. note::
    By default on Mac, gcc is an alias for the clang compiler. 
    This can be verified by running ``gcc -v`` in a terminal and verifying output.

    See :ref:`CompilerOptions` for more information on compilers.

Once done, you can continue with the :ref:`InstallCommon` section.

.. _InstallCommon:

Common installation instructions
--------------------------------

MPCD is distributed via it's `Github repository <https://github.com/Shendruk-Lab/MPCD>`_, from which it can be downloaded by either clicking the green "Code" button and selecting "Download ZIP", or by cloning the repository using git. 
To clone via git, open a UNIX terminal and run

.. code-block:: console

    git clone https://github.com/Shendruk-Lab/MPCD.git

This will create a folder called ``MPCD`` in your current directory, containing the source code for MPCD. 
To compile, navigate within this folder (such that the ``Makefile`` is at the same level as your current working directory) and then call 

.. code-block:: console

    make

**This is all that is required to compile the code**.

Resulting executable
^^^^^^^^^^^^^^^^^^^^

The compiled executable file will be called ``mpcd.out`` and will be located in the same directory as the ``Makefile``.
The program can then be run by calling ``./mpcd.out``, with arguments as detailed in the :ref:`ProgArgs` section.

Makefile options
****************
The functionality of the ``make`` call to compile MPCD is entirely controled by the ``Makefile`` in the root MPCD directory. 
There are two main ways to control the compilation process: By adjustng Makefile variables, or calling make phonies.
In this section we highlight a specific important phony call, before explaining both methods.

Make clean
----------
A clean operation is a particularly important phony call.
This will remove compiled binaries and "object files", which are used by the compiler.
Performing a make clean after any significant code change, or before a new compilation, is highly recommended.

To perform a clean, simply run

.. code-block:: console

    make clean

Makefile phonies 
----------------
Makefile phonies are the recommended way to compile non-standard builds of MPCD.
These are pre-defined rules that can be called by running ``make <phony>``, of which ``make clean`` is the most important example.

Phonies are defiend in the ``Makefile`` by a line of format

.. code-block:: makefile

    .PHONY: <phony>

The most up-to-date list of these will always be present in the ``Makefile``, but a summary of particularly important ones include

* ``clean``: 
    Removes all compiled binaries and object files.
* ``guide``: 
    Builds this user guide. 
    Requires Sphinx to be installed.
* ``docs``: 
    Builds doxygen code documentation for MPCD.
    Requires doxygen to be installed.
* ``debug``: 
    Compiles the code with debug symbols and optimisation disabled.
    Designed for use with any general debugger, and is only used for code development.
* ``debug+``: 
    Compiles the code with debug symbols and optimisation disabled, but with extra information specifically read by the ``gdb`` debugger and the ``valgrind`` tool.
    Only used for code development.
* ``prof``: 
    Compiles the code with debug symbols and profiling logging enabled.
    Designed for easier use with the ``gprof`` profiler, and is used to optimise the simulator during development.


Makefile variables
------------------
Direct editing of the ``Makefile`` is highly discouraged.
Instead, variables can be appended to the ``make`` call to change compilation behaviour.

There are three variables that may be helpful to modify:

* ``cc``: 
    The C compiler to use. By default, this is ``gcc``, however it can be changed to any C compiler installed on your system.
* ``cflags``: 
    The compilation flags to use. By default, this is just ``-Wall``, to ensure all warnings are shown. 
    This can be changed to any flags supported by your C compiler.
* ``opt``: 
    This option is for optimisation flags. 
    Compilers will try to optimise and speed up the binaries they produce by making shortcuts in code, however this can cause issues with the compiled code.
    By default, this is set to ``-O3``, which is the highest level of optimisation, so it may be necessary to lower it to lower levels such as ``-O2``, ``-O1``, or ``-O0``.

An example of a ``make`` call setting all of these variables is

.. code-block:: console

    make cc=clang cflags="-Wall -Wextra -Werror" opt=-O2

.. _CompilerOptions:

Compilers 
---------
MPCD has been tested with a number of compilers, and is known to work with the following

- ``gcc`` --- The GNU C compiler
- ``clang`` --- The LLVM C compiler
- ``icc`` --- The Intel C compiler

Of these, ``gcc`` and ``clang`` are the most common and we endeavour to support both primarily. 
We have qualatative evidence that ``clang`` is more resilient to code "undefined behaviour", and is slighly more performant.

These compilers can be set by setting the ``cc`` parameter of the ``Makefile``: 

.. code-block:: console

    make cc=<compiler>

.. note:: 
    The Intel C Compiler is known to produce very performant code, but only if you are running on Intel CPU hardware. 
    Double check the brand of CPU that your computer or cluster is Intel before using this compiler.

.. _ProgArgs:

Program arguments & input files
*********************************
MPCD is designed to be run from command line, and as such has an arguments interface.
In this section, we will describe the arguments, and also describe the required input files for the two input modes.

Command line arguments
----------------------
.. note:: 
    All arguments are case sensitive.

MPCD arguments are purely programmatic --- There is no GUI, or interacive mode, and all physics is specified by the input files.
There are only two required arguments for MPCD

* Input files, which can either use ``-i`` and point to a ``.json`` file for :ref:`JSONInput`, or use :ref:`LegacyInput`.
* An output directory, which can be set using ``-o``.

Optional arguments include

* ``-h``: 
    Prints a help message, explaining arguments and their usage, then exits.
* ``-v``: 
    Prints a legacy version message.

Some examples of valid calls are

.. code-block:: console

    ./mpcd.out -i ./path/to/input.json -o ./path/to/output
    ./mpcd.out -h 
    ./mpcd.out -v

.. _JSONInput:

JSON input
----------

The recommended way to run MPCD is by using a JSON input file. 
This is a single file that is .json formatted, and contains all of the physics information for the simulation.

A guide to all input parameters is provided on Github `here <https://github.com/Shendruk-Lab/MPCD/blob/master/docs/InputGuide.md>`_.
Furthermore, the :ref:`Usage <usageChapter>` section and the :ref:`Tutorials <tutorials>` in this user guide give explanations on how to set up input files for specific simulations.

.. note:: 
    If you want to use the molecular dynamics (MD) subsystem, you will need to also provide an MD .inp file along with a .json input file. 
    See `the MD input guide <https://github.com/Shendruk-Lab/MPCD/blob/master/docs/MDguide.md>`_ for more information.

.. _LegacyInput:

Legacy input
------------

.. warning:: 
    Legacy input files have been considered depreciated since the introduction of JSON input files in summer 2021.
    They are still supported, but no features implemented since 2020 are supported in them.

    Use of legacy input is strongly discouraged.

In order to use legacy input files, you pass them to the ``mpcd.out`` using the legacy input argument ``-Li`` followed by the path to the input files.
For example

.. code-block:: console

    ./mpcd.out -Li ./path/to/input/files/directory -o ./path/to/output

Legacy input files are a series of 5 ``.inp`` files, which are read in order to set up the simulation. 
These include

* input.inp
* bc.inp
* printcom.inp
* swimmer.inp
* md.inp

These files are read in order, and are all required for the simulation to run.
Furthermore, these files expect parameters input in a particular order to function.
Examples are provided within the ``sampleInputs`` folder of the MPCD repository, and an incomplete guide is provided `on Github <https://github.com/Shendruk-Lab/MPCD/blob/master/docs/legacyInputSummary.txt>`_.
