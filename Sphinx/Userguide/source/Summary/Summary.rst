Executive summary
##################


*Multi-Particle Collision Dynamics* (**MPCD**) is a coarse-grained, particle-based, simulation technique for simulating fluids that reproduces thermally fluctuating, hydrodynamics interactions. 
MPCD is ideal for situations involving moderate Peclet numbers---where diffusion and advection are both significant. 

This implementation of MPCD is written in *C* and requires no external libraries. 
It is serial but versatile, and is able to simulate complex boundary conditions, polymer suspensions, nematic fluids, colloids and swimming bacteria. 

In addition to this userguide, a number of reviews have been published including:

* Gompper, Ihle, Kroll & Winkler (2009). `Multi-Particle Collision Dynamics: A Particle-Based Mesoscale Simulation Approach to the Hydrodynamics of Complex Fluids <https://link.springer.com/chapter/10.1007/978-3-540-87706-6_1>`_. In *Advanced Computer Simulation Approaches for Soft Matter Sciences III. Advances in Polymer Science*, vol 221.
* Kapral (2008). `Multiparticle Collision Dynamics: Simulation of Complex Systems on Mesoscales <https://onlinelibrary.wiley.com/doi/10.1002/9780470371572.ch2>`_. In *Advances in Chemical Physics*. 
* Yeomans (2006). `Mesoscale Simulations: Lattice Boltzmann and Particle Algorithms <https://www.sciencedirect.com/science/article/pii/S0378437106004067>`_. *Physica A: Statistical Mechanics and its Applications*, vol 369 (1). 
* Howard, Nikoubashman & Palmer (2019). `Modeling Hydrodynamic Interactions in Soft Materials with Multiparticle Collision Dynamics <https://www.sciencedirect.com/science/article/pii/S2211339819300024>`_. *Current Opinion in Chemical Engineering*, vol 23.
