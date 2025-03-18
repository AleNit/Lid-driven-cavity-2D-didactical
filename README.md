# Lid-driven-cavity-2D-didactical
Simulate the 2d Lid-Driven-Cavity test by solving incompressible Navier-Stokes equations by means a fractional step method.  
Spatial derivatives are discretized by 2nd-order accurate centered differences on a staggered grid. Time integration is carried out by a 2nd order Adams-Bashforth method. The steady solution at Re=400 is varified against the reference data from Ghia et al. (1982). Boundary conditions are prescribed with a ghost cell approach.
