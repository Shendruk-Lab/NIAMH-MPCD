MPCD supports over 20 unique collision operators :math:`\vec\Xi`.
These are set in the input ``.json`` file using either the ``"rTech"`` or ``"collOp"`` keywords as follows:

.. code-block:: console

    "collOp": <int value>

Where ``<int value>`` is an integer code corresonding to the desired collision operator.

In this section, we describe the various families of collision operators available in MPCD, state their associated keycodes, and describe the parameters they make use of.

.. _SRD:

Stochastic rotation dynamics
============================
The original derivation of MPCD is known as Stochastic Rotation Dynamics (SRD) **TODO: cite**.
Translating to modern terminology, SRD's collision operator :math:`\vec\Xi` performs
a rotation through an angle :math:`\theta` about a randomly chosen axis such that conservation of energy, isotropy and a Maxwell Boltzmann velocity distribution are met in the continuum limit.
It is this random rotation operation that gives this MPCD variant the name Stochastic Rotation Dynamics (SRD). 

More formally, in SRD the collision operator is

.. math:: 
    :name: eq:SRD

    \vec{\Xi}_i^\mathrm{SRD}(t) = \Omega_\theta \cdot \left( 
        \vec{v}_i(t) - \vec{v}_{cm}(t)    
    \right)

where :math:`\Omega_\theta` is a rotation tensor through an angle :math:`\theta` about a given axis, and :math:`\vec{v}_{cm}` is the center of mass velocity of the particles in the collision cell.

To be specific, the rotation tensor :math:`\Omega_\theta` is defined as

.. math:: 
    :name: eq:SRD_rotation_tensor

    \begin{align}
        \Omega_\theta &= \begin{pmatrix}
            n_x^2+(1-n_x^2)\cos\theta & n_xn_y\left(1 - \cos\theta\right)-n_z \sin\theta & n_xn_z\left(1-\cos\theta\right)+n_y \sin\theta \\
            n_xn_y \left(1 - \cos\theta\right)+n_z \sin\theta & n_y^2 + \left(1-n_y^2\right)\cos\theta & n_yn_z \left(1 - \cos\theta\right)-n_x \sin\theta \\
            n_xn_z\left(1 - \cos\theta\right)-n_y \sin\theta & n_yn_z\left(1 - \cos\theta\right)+n_x \sin\theta  &n_z^2 + \left(1-n_z^2\right)\cos\theta 
        \end{pmatrix}. 
    \end{align}

where the vector :math:`\vec{\hat{n}} = [n_x, n_y, n_z]` defines a unit vector representing the axis of rotation.

**TODO: SRD schematic from Tyler's thesis**

This version of SRD where :math:`\vec{\hat{n}}` is chosen from the unit sphere is implemented with a keycode of ``0``. 
To use it add the following to your input ``.json`` file:

.. code-block:: console

    "collOp": 0

However, it is found that it is sufficient to choose the axis of rotation :math:`\vec{\hat{n}}` randomly from the the cartesian axes. 
ie, :math:`\vec{\hat{n}} = [1, 0, 0]`, :math:`[0, 1, 0]` or :math:`[0, 0, 1]`.
This is quicker, and in bulk is sufficient to reproduce a fluid with the correct hydrodynamic properties **TODO: cite**.

This version of SRD is implemented with a keycode of ``1``.
To use it add the following to your input ``.json`` file:

.. code-block:: console

    "collOp": 1

It is possible to verify that momentum and energy are conserved before and after collision by explicitly computing the them.

.. _Andersen:

Andersen-thermostatted MPCD
===========================
While SRD is sufficient for some applications, it is not without its drawbacks:

- The SRD collision operator does not generally conserve angular momentum, as the positions of particles within the cell are not considered.
- An explicit external thermostat is required to control the system energy when external forces or activity is applied.

Following other mesoscale simulation techniques, such as DPD, where a thermostat is included intrisically as part of the simulation technique, thermostats have been incoroporated into MPCD collision operators.
That is, :math:`\vec\Xi` is both the collision operator, and a local thermostat for the cell, such that no velocity rescaling techniques are necessary.

The primary thermostatted collision operator of use in this simulator is the Andersen thermostatted collision operator **TODO: cite**.
This collision operator is defined as

.. math:: 
    :name: eq:AndersenOp

    \vec{\Xi}_i^\mathrm{A}(t) = 
    \vec\xi_i -
    \frac{\sum_j^{N_C} m_j \vec\xi_j}{\sum_j^{N_C} m_j}

where each component of :math:`\vec\xi_i` is randomly generated from a Gaussian distribution with variance :math:`\sqrt{k_B T/ m}`.

The final term in :ref:`the Andersen collision operator <eq:AndersenOp>` is often referred to simply as the residual, and denoted :math:`\langle \vec\xi_j \rangle_{N_C}`. 
Each fluid particle in a given cell is assigned a new, randomly chosen velocity during each collision.
However, the presence of the residual means that each particle's velocity is the cell’s centre of mass velocity plus its new randomly assigned value, minus the average of all the new random values for that cell.
This ensures that the center of mass velocity of each cell does not change during collision, and thus momentum is conserved.

.. list-table:: 
    :header-rows: 0
    :widths: 50 50
    :align: center

    * - .. image:: ./CollisionOperators/Andersen1.png
            :width: 95%
            :align: center
      - .. image:: ./CollisionOperators/Andersen2.png
            :width: 95%
            :align: center
    * - **(a)** Bin MPCD particles into cells.
      - **(b)** Calculate the center of mass velocity :math:`\vec{v}^\mathrm{cm}` of the cell.
    * - .. image:: ./CollisionOperators/Andersen3.png
            :width: 95%
            :align: center
      - .. image:: ./CollisionOperators/Andersen4.png
            :width: 95%
            :align: center
    * - **(c)** For each MPCD particle, generate a random velocity that sums to :math:`\vec{v}^\mathrm{cm}`.
      - **(d)** Apply to each MPCD particle.

The Andersen thermostatted collision operator is implemented with a keycode of ``2``.
To use it, add the following to your input ``.json`` file:

.. code-block:: console

    "collOp": 2

Note that the :ref:`the basic Andersen collision operator <eq:AndersenOp>` conserves mass and translational momentum, but not angular momentum (or energy, as it is thermostatted).
However, it can be extended to do so: The basic Andersen collision operator introduces a small change in angular momentum every timestep, which we can denote :math:`\delta \vec{L}`. 
This can be cancelled out by applying a small counter-rotation to each particle, such that the total angular momentum change is zero.
If the particles in the cell have a given intertia tensor :math:`I` about the center of mass, then the counter-rotation required will be :math:`\vec \omega = I \cdot \delta \vec{L}`.

A term can be added to the Andersen collision operator that performs this counter-rotation, as such:

.. math:: 
    :name: eq:AndersenOpAngular

    \vec{\Xi}_i^\mathrm{A}(t) = 
    \vec\xi_i -
    \langle \vec\xi_j \rangle_{N_C} +
    \left[
        I \cdot \left(
            \sum_j^{N_C} m_j \left\lbrace
                \vec{x}'_j \times (\vec{v}_j - \vec\xi_j)
            \right\rbrace
        \right)
    \right] \times \vec{x}'_j

.. note:: 
    The angular-conserving Andersen thermostatted collision operator above is the default, and preferred, collision operator in this MPCD simulator for traditional fluids.

The angular-conserving Andersen thermostatted collision operator is implemented with a keycode of ``3``.
To use it, add the following to your input ``.json`` file:

.. code-block:: console

    "collOp": 3

.. note:: 
    The angular-conserving Andersen thermostatted collision operator is the basis used for the :ref:`Nematic MPCD <chapter8>` algorithm.

.. _Vicsek:

Vicsek MPCD
===========
**TODO: Olek has been writing this. Chate Vicsek should be in here too I guess**

.. _ActiveNematic:

Active nematic MPCD
===================
.. warning:: 
    The Active-nematic MPCD algorithm relies on liquid crystal mode being enabled. See the section on :ref:`Nematic MPCD <chapter8>` for more information.

    This can be enabled by adding the following to your input ``.json`` file:
    
    .. code-block:: console

        "lc": 1

The simulator also has the capability to simulate active-nematic fluids **TODO: cite Tim**.
This builds on the :ref:`Nematic MPCD collision operator <eq:NematicCollision>`, henceforth referred to as :math:`\vec\Xi^\mathrm{N}`.

The fundamental requirement for a nematic to be turned into an active-nematic is the presence of local force dipoles. 
In order to induce these, a planar collision operator is applied to the nematic particles in the cell:
A plane normal to the cell's director, :math:`\vec{n}_C`, is placed at the cell center of mass.
Particles forward of the plane get a kick forward, and particles backward of the plane get a kick backward.

.. list-table:: 
    :header-rows: 0
    :widths: 50 50
    :align: center

    * - .. image:: ./CollisionOperators/ActivityDiagram-0.png
            :width: 95%
            :align: center
      - .. image:: ./CollisionOperators/ActivityDiagram-1.png
            :width: 95%
            :align: center
    * - **(a)** Bin MPCD particles into cells.
      - **(b)** Place a plane normal to the cell director at the cell center of mass.
    * - .. image:: ./CollisionOperators/ActivityDiagram-2.png
            :width: 95%
            :align: center
      - .. image:: ./CollisionOperators/ActivityDiagram-3.png
            :width: 95%
            :align: center
    * - **(c)** For each MPCD particle forward of the plane, give a kick forward, and vica-versa.
      - **(d)** Apply kicks to the post-nematic velocity collision.

This is sufficient to reproduce a local force dipole within each cell. 
Furthermore, the positioning of the plane ensures that this operation conserves linear momentum, despite the local injection of energy. 

Mathematically, this collision operation is expressed as:

.. math:: 
    :name: eq:ActiveNematicCollision

    \vec{\Xi}_i^\mathrm{AN}(t) =
    \vec{\Xi}_i^\mathrm{N}(t) +
    \alpha_C \delta t \left( 
        \frac{\kappa_i}{m_i} - \langle \frac{\kappa_j}{m_j} \rangle_C
    \right)

where :math:`\kappa_i=\pm1` represent whether a particle is forward or backward of the plane **TODO: cite**.

The magnitude of the kicks, ie the active energy input :math:`\alpha_C`, has 4 different formulations:

1. **Active-Sum**
2. **Active-Average**
3. **Sigmoidal-Average**
4. **Sigmoidal-Sum**

These are described in the sub-sections below.

.. note:: 
    For optimal density behaviour, we recommend using Sigmoid-Sum in the general case.

.. _ActiveSum:

Active-sum
----------

Active-Sum is the simplest formulation, and is akin to a "particle-carried activity" **TODO: cite**.
The activity of every particle within the cell is summed, and a kick of this magnitude is then applied to every particle. ie:

.. math::

    \alpha_C^\mathrm{Sum} = \sum_j^{N_C} \alpha_j

This is implemented with a keycode of ``16``, and can be used by adding the following to your input ``.json`` file:

.. code-block:: console

    "collOp": 16

.. _ActiveAv:

Active-average
--------------

In contrast to active-sum, with "particle-carried activity", active-average instead has "cell-carried activity" **TODO: cite**.
Mathematicaly, the sum of particle activities within a cell is averaged, which is then applied as an active force magnitude to every particle. 
ie:

.. math:: 

    \alpha_C^\mathrm{Av} = \frac{1}{N_C} \sum_j^{N_C} \alpha_j

This leads to smaller active forces than active-sum, resulting in less severe density fluctuations while still maintaining active nematic scaling laws **TODO: cite**.

This is implemented with a keycode of ``17``, and can be used by adding the following to your input ``.json`` file:

.. code-block:: console

    "collOp": 17

.. _SigmoidalANMPCD:

Sigmoidally modulated active-nematic MPCD
-----------------------------------------

In an effort to reduce the density fluctuations seen in Active-Sum and Active-Average, the code implements a density modulated active-nematic MPCD algorithm **TODO: cite**.
A sigmoidally shaped modulation function is used

.. math:: 

    \mathcal{S}_C(\rho_C; \sigma_p, \sigma_w) = 
    \frac{1}{2}
    \left(
        1 - \tanh \left(
            \frac{\rho_C - \langle\rho_C\rangle_{N_C} \left( 1 + \sigma_p \right)}{\langle\rho_C\rangle_{N_C} \sigma_w}
        \right)
    \right)

This modulation function takes the local cell number density, :math:`\rho_C`, and returns a values in the range :math:`[0, 1]` to modulate it.
This is controlled by two parameters:

- :math:`\sigma_p` is a decimale repesenting where the modulation function is centred relative, in units of mean density.
- :math:`\sigma_w` is a decimal representing how wide the modulation drop-off should be, in units of the mean density.

These are configured with species-based parameters ``"sigPos"`` and ``"sigWidth"`` in the input ``.json`` respectively.

**TODO: show figure of modulation function**

This modulation function is then applied to each cells active force magnitude, giving two new formulations:

- **Sigmoidal-Average**: :math:`\alpha_C = \mathcal{S}_C(\rho_C) \alpha_C^\mathrm{Av}`
- **Sigmoidal-Sum**: :math:`\alpha_C = \mathcal{S}_C(\rho_C) \alpha_C^\mathrm{Sum}`

These are implemented with keycodes ``20`` and ``21`` respectively, and can be used by adding the following to your input ``.json`` file:

.. code-block:: console

    "collOp": 20

.. code-block:: console

    "collOp": 21
