Analysis
########

In the root MPCD repository, we provide a series of scripts to help visualise simulations output by MPCD. 
These scripts are written in Python, and require Python 3, and required packages are included (``defectHandler.py``, ``shendrukGroupFormat.py``).

In this section, we briefly describe their usage, including their arguments and what they output.
Examples of their usage can be found in :ref:`the tutorials section <tutorials>`.

.. error:: 
    Tim's comment: An example of their output should be provided either here, or in the tutorials section. 

.. note:: 
    These analysis scripts are provided primarily to quickly verify simulation output for first time users, and to give an example of how to parse MPCD's output files.
    Users are expected to develop their own analysis scripts appropriate for their own purposes.

Common
******

.. warning:: 
    This section does not apply to ``localS.py`` due to arguments being handled differently.

Most scripts make use of Python's ``argparse`` library to parse command line arguments.
Arguments have tried to be kept consistent between scripts, despite some differences in their usage.
Here, we briefly describe some common features

Help
====

In general, you can always give the ``--help`` (shorthand ``-h``) flag to any script to see a list of arguments and their usage.

For example:

.. code-block:: console

    $ python3 analysisScripts/orientationField2Danimated.py -h

gives output

.. code-block:: console

    usage: orientationField2Danimated.py [-h] [--qx QX] [--qy QY] [-c LENGTH] [-a MYASPECT] [-k KEEPFRAMES] [-d DEFECTDATA] dataname inputname start finish avdim

    Orientation field rendering script.

    positional arguments:
    dataname              Path to the data (should be directorfield.dat)
    inputname             Path to input .json file
    start                 Starting timestep for averaging
    finish                Finishing timestep for averaging
    avdim                 Dimension to 'slice' over

    optional arguments:
    -h, --help            show this help message and exit
    --qx QX               Only show every qx arrow in x
    --qy QY               Only show every qy arrow in y
    -c LENGTH, --length LENGTH
                          Length of director lines
    -a MYASPECT, --myAspect MYASPECT
                          'auto' or 'equal'
    -k KEEPFRAMES, --keepFrames KEEPFRAMES
                          0=don't keep (delete) frames; 1=keep frames
    -d DEFECTDATA, --defectData DEFECTDATA
                          Path to defect data (if any)

The output may vary depending on script, but will always be the most up-to-date description of the script's arguments.

.. _ConsistentArguments:
Consistent arguments
====================

The following arguments are (generally) consistent between scripts:

- ``dataname`` 
    Path to the input data file.
    Usually a ``.dat`` file.
    The exact file will vary from script to script.
- ``inputname`` 
    Path to the input ``.json`` file corresponding to the given data.
    This is used to read in simulation parameters.
- ``start`` 
    Starting timestep index for plotting.
- ``finish`` 
    Ending timestep index for plotting.
    Set to a very large number to plot until the end of the simulation.
- ``avdim`` 
    For scripts that will "slice" over a dimension, this is the slicing dimension.
    These scripts will average in this dimension, hence the name ``avdim``.
    Accepts ``"x"``, ``"y"``, or ``"z"``.
- ``-a``/ ``--myAspect`` 
    (Optional, default ``"auto"``)
    The aspect ratio to produce plots in. 
    Either ``"auto"`` for automatic aspect ratio, or ``"equal"`` for equal aspect ratio.
- ``-k``/ ``--keepFrames`` 
    (Optional, default ``0``)
    When producing frames for a movie, whether to keep the frames or delete them after the movie is produced.
    Set to ``1`` to keep the frames, and ``0`` to delete frames.

density2Danimated.py 
********************

.. include:: density2Danimated.txt

flowFieldAnimated.py 
********************

.. include:: flowFieldAnimated.txt

localS.py 
*********

.. include:: localS.txt

orientationField2DAnimated.py
*****************************

.. include:: orientationField2Danimated.txt

probDistScalar.py
*****************

.. include:: probDistScalar.txt

probDistVector.py
*****************

.. include:: probDistVector.txt

swimmerRefFrame_flow2D.py
*************************

.. include:: swimmerRefFrame_flow2D.txt

swimmerTraj_fields2D.py
***********************

.. include:: swimmerTraj_fields2D.txt

topochargeField2Danimated.py
****************************

.. include:: topochargeField2Danimated.txt
