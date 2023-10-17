MPCD supports over 20 unique collision operators :math:`\vec\Xi`.
These are set in the input ``.json`` file using either the ``"rTech"`` or ``"collOp"`` keywords as follows:

.. code-block:: console

    "collOp": <int value>

Where ``<int value>`` is an integer code corresonding to the desired collision operator.

In this section, we describe the various families of collision operators available in MPCD, state their associated keycodes, and describe the parameters they make use of.

.. _SRD:

Stochastic Rotation Dynamics
----------------------------
**TODO**

.. _Andersen:

Andersen Thermostatted MPCD
---------------------------
**TODO**

.. note:: 
    **TODO: Mention that if you use angular conserving andersen, you can enable liquid crystals. Link to section in here too**

.. _Vicsek:

Vicsek MPCD
-----------
**TODO: Olek has been writing this. Chate Vicsek should be in here too I guess**

.. _ActiveNematic:

Active Nematic MPCD
-------------------
.. warning:: 
    **TODO: mention that this assumes liquid crystal is on**

**TODO**
**TODO: mention it is dervied from angular conserving andersen**
