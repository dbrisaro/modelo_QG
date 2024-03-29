#............................................................................
# This program solves the barotropic vorticity equation in non-dimensional
# form using finite differences.
#
# The model has incorporated the "partial" slipping boundary conditions.
#............................................................................

im=152           # number of grid points in the zonal direction
jm=152           # number of grid points in the meridional direction
ds=0.05          # grid step
dt=0.05          # time step
Ro=0.005           # Rossby number (measures non-linearity of the flow)
eps=0.00          # non-dimensional coefficient representing bottom friction
Ah=0.001           # non-dimensional coeff. of horizontal Laplacian mixing
Bh=0.0           # non-dimensional coeff. of horizontal bi-harmonic mixing
gamma=0.0        # coeff. of "intermediate slipping" used as boundary cond.
nst=1            # start time step number
nend=50        # end time step number
nlpt=100         # frequency (time steps) for saving output
ncrit=4000       # number of steps allowed to do the relaxation (sub. helm)
pcrit=0.1        # criterium to stop the relaxation
BFP=1            # Beta (BFP=1) or F plane (BFP=0)
GYR=1            # Simple Gyre (GYR=1) or Double Gyre (GYR=2)
HEM=-1           # North Hemisphere Gyre (HEM=1) or South Hemisphere Gyre (HEM=-1)
