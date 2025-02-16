# Rheologically-tuned-modes-of-collective-transport-in-active-viscoelastic-films
This repository contains Python scripts used to investigate the spatiotemporal patterns in viscous, viscoelastic, or elastic thin active films. 

<h2>REQUIREMENTS</h2>

Python 3

tested with Python 3.10.12


<h2>USAGE INSTRUCTIIONS</h2>

All scripts are executed via the command line. For simplicity, parameter values are set directly within the scripts.

The scripts starting "*simulation...*" are used for the main calculations. They solve the evolution equations for the polar order parameter field, the velocity field and the displacement field via a pseudo spectral method and perform additional analysis of the dynamics.

*determineStateDiagram.py* <br>
This script determines the boundaries between regions in the state diagram of the system. These are calculated based on the linear stability of the stationary states (isotropic and polar) and the existence of the rotational solutions. The script moves through different values of the displacement relaxation time and determines the values of the active force strength that mark the boundaries. For the calculations, we employ bisection methods, that is, we start from a large interval of possible values of active force strength and successively half the size of this interval until an appropriate accuracy is achieved.

*simulation_determineCorrelations.py* <br>
This script solves the evolution equations for a spatially extended system and calculates the local angle of polar order parameter, velocity and displacement. In the case of local rotations, these angles are shifted with respect to each other. This is used to determine the spatial correlations of the patterns in terms of the phase shift between polar order parameter and velocity.

*simulation_tracerMotion.py* <br>
This script solves the evolution equations for a spatially extended system and simultaneously tracks the motion of tracers advected by the velocity field. The resulting trajectories can then be used to visualize the modes of collective transport. The script allows for periodic changing of the displacement relaxation time. This results in temporary directed transport and intermediate storage times.

*simulation_determineInstabilityStripePatterns.py* <br>
This script finds the value of active force strength needed to destabilize stripe-like patterns for a given value of displacement relaxation time. It employs a bisection method, that is, we start from a large interval of possible values of active force strength and successively half the size of this interval until an appropriate accuracy is achieved. For every sampling value of active force strength, the script solves the evolution equations for a spatially extended system starting from an initially polar state with stripe-like modulations. Then, we check whether the polar order parameter locally starts to rotate to determine if the stripe-patterns become unstable.


