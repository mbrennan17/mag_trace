# Planetary Magnetic Field Line Tracer

The mag_trace function propogates the magnetic field line from user defined starting position(s), using the planetary internal and external field model functions (also provided by the user). Additional optional inputs can be provided to change default settings: 
    altitude for field line termination 
    planet/body radii of triaxial ellipsoid
    max propogation distance
    max radius of field line from planet/body
    min radius for applying external magnetic field model

The vectorized function allows single or multiple input positions. Sets of field line points are output as a cell array, with a field line array for each input position. The current configuration uses the native ode45 integrator, leveraging the dynamic step size propogation for faster computation, as well as the built-in termination functionality. Other integrators are also being considered for easier translation to other platforms (python, IDL, etc). 

This is part of a community code project: Magnetospheres of the Outer Planets Group Community Code

Authors: M. Brennan C. Lawler, and R.J. Wilson
