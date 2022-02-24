# Random seeding of particles at the inflow boundary

Current setup:
At each time step there is a 1/2 chance of a particle being seeded at the inflow boundary (x=0). The initial y-coordinate of the particle is chosen using the numpy random module.


# TODO
-[ ] Include Stokes drag for particle motion
-[ ] Related to above; store additional time step 
-[ ] Related to above; store in solver or particle class?
