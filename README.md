# Planetary Magnetic Field Line Tracer

The mag_trace function propogates the magnetic field line from user defined starting position(s), using the planetary internal and external field model functions (also provided by the user). Additional optional inputs can be provided to change default settings: 
    Altitude for field line termination 
    Planet/body radii of triaxial ellipsoid
    Max propogation distance
    Max radius of field line from planet/body
    Min radius for applying external magnetic field model

The vectorized function allows single or multiple input positions. Sets of field line points are output as a cell array, with a field line array for each input position. The current configuration uses the native ode45 integrator, leveraging the variable step size propogation for faster computation, as well as the built-in termination functionality. Other integrators are also being considered for easier translation to other platforms (python, IDL, etc). 

Authors: M. Brennan C. Lawler, and R.J. Wilson

This is part of a community code project: Magnetospheres of the Outer Planets Group Community Code and is formulated to work seamlessly with the Jupiter internal and external magnetic field models as outlined in the January 2023 in Space Science Reviews paper:

Wilson, R.J., Vogt, M.F., Provan, G. et al. Internal and External Jovian Magnetic Fields: Community Code to Serve the Magnetospheres of the Outer Planets Community. Space Sci Rev 219, 15 (2023). https://doi.org/10.1007/s11214-023-00961-3


