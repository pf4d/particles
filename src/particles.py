from pylab import *

# PARTICLES MODULE:
# This set of classes was developed to model interacting spherical particles in
# 3 dimensions. There are presently two types of interaction forces possible. 
# Each is a class
#
# GranularMaterialsForce: This is a force that represents elastic forces between
# particles that are overlapping, and a damping, or dissipative force. It does
# not include the often used rotational degrees of freedom. Also note that the
# constraint forces are taken care of with an extra function call.
# 
# LennardJonesForce: This force describes the weak attraction at large distances
# and strong repulsion experienced at close distances by mono-atomic gases.
#
# The equation of motion are integrated by the Verlet method, which is presently
# the only integration scheme supported.
#
# All particle data (position, velocity, acceleration, radii, and distances) are
# stored in the Particles class. 

# 9/23/08 JVJ and Tim Bocek

class GranularMaterialForce(object):

  def __init__(self, k=1.5, gamma=0.3, g=0.25):
    # parameters in force model
    self.k     = k       # Elastic 'bounce'
    self.gamma = gamma   # Energy dissipation/loss
    self.g     = g       # Gravity
    self.f     = 1.0

  def __call__(self, p):
    # Find position differences
    d, dx, dy, dz = p.distanceMatrix(p.x, p.y, p.z)

    # Compute overlap
    dr = d - p.sumOfRadii
    fill_diagonal(dr, 0)

    # No forces arising in no overlap cases
    dr[dr > 0] = 0

    # Compute elastic particle/particle forces
    mag_r = self.k * dr

    # Velocity differences
    dv, dvx, dvy, dvz = p.distanceMatrix(p.vx, p.vy, p.vz)
    da, dax, day, daz = p.distanceMatrix(p.ax, p.ay, p.az)
    
    # Damping terms
    vijDotrij             = dvx*dx + dvy*dy + dvz*dz
    vijDotrij[dr==0]      = 0 

    # Damping is subtracted from force
    mag_r += self.gamma * vijDotrij / d

    # floor components of acceleration :
    crx, cry, crz, ctx, cty, ctz = self.floorConstraint(p)

    # Project onto components, sum all forces on each particle
    p.ax = sum(mag_r * dx/d * p.ratioOfRadii, axis=1) + crx
    p.ay = sum(mag_r * dy/d * p.ratioOfRadii, axis=1) + cry - self.g 
    p.az = sum(mag_r * dz/d * p.ratioOfRadii, axis=1) + crz
    
    omegax = tile(p.omegax, (p.N, 1))
    omegay = tile(p.omegay, (p.N, 1))
    omegaz = tile(p.omegaz, (p.N, 1))
    
    omegax = omegax + omegax.T 
    omegay = omegax + omegax.T 
    omegaz = omegax + omegax.T 
    
    omegax[dr==0] = 0
    omegay[dr==0] = 0
    omegaz[dr==0] = 0

    # projection of a onto the tangent plane to r (tangential acceleration) :
    atx = dax - (dax * dx) / d**2 * dx
    aty = day - (day * dy) / d**2 * dy
    atz = daz - (daz * dz) / d**2 * dz
    
    # projection of v onto the tangent plane to r (tangential velocity) : 
    vtx = dvx - (dvx * dx) / d**2 * dx
    vty = dvy - (dvy * dy) / d**2 * dy
    vtz = dvz - (dvz * dz) / d**2 * dz
    
    # calculate the radius vector to the point of torque :
    radx = dx/d * p.r
    rady = dy/d * p.r
    radz = dz/d * p.r
   
    # no angular forces where particles do not touch :
    atx[dr==0] = 0
    aty[dr==0] = 0
    atz[dr==0] = 0
    vtx[dr==0] = 0
    vty[dr==0] = 0
    vtz[dr==0] = 0

    # calculate torque (r x F) :
    taux = rady*atz - aty*radz
    tauy = radx*atz - atx*radz
    tauz = radx*aty - atx*rady

    # calculate tangential velocity parallel to torque (r x vt) :
    epix = rady*vtz - vty*radz
    epiy = radx*vtz - vtx*radz
    epiz = radx*vty - vtx*rady

    # linear velocity frictional contribution to torque coefficient :
    f = 0.0

    # angular velocity damping coefficient from medium (air) :
    g = 0.1

    # moment of inertia for a sphere :
    I = 0.4*p.r**2
    
    # angular momentum exchange coefficient :
    kappa = 0.9*p.r
   
    # project onto components, sum all angular forces on each particle
    p.alphax = + sum(-(taux + f*epix) / I - kappa*omegax, axis=1) \
               - g*p.omegax \
               + ctx
    p.alphay = + sum(-(tauy + f*epiy) / I - kappa*omegay, axis=1) \
               - g*p.omegay \
               + cty
    p.alphaz = + sum(-(tauz + f*epiz) / I - kappa*omegaz, axis=1) \
               - g*p.omegaz \
               + ctz

  def floorConstraint(self, p):
    """ 
    This is a highly specific function for a floor that responds (elasticity 
    and damping) the same way a particle does. Presently, if constraints are 
    to change, one would have to rewrite the function.
    """
    r          = 3.0 # This is how 'hard' the floor is
    fd         = p.y + p.L/2 - p.r
    fd[fd > 0] = 0
    
    floorForce_r     = -self.k * fd 
    floorDamping_r   = -self.gamma * p.vy * fd
    floorForce_r     = floorForce_r - floorDamping_r
    crx = 0
    cry = floorForce_r * r / p.r
    crz = 0
    
    floorDamping_tx  = p.omegax.copy()
    floorDamping_ty  = p.omegay.copy()
    floorDamping_tz  = p.omegaz.copy()
    floorDamping_tx[fd == 0] = 0 
    floorDamping_ty[fd == 0] = 0 
    floorDamping_tz[fd == 0] = 0 

    ctx = -self.f*floorDamping_tx
    cty = -self.f*floorDamping_ty
    ctz = -self.f*floorDamping_tz

    return crx, cry, crz, ctx, cty, ctz


class VerletIntegrator(object):

  def __init__(self, dt=0.01):
    # Time step
    self.dt = dt

  def __call__(self, force, p):
    dt = self.dt

    # Position update
    p.x = p.x + p.vx*dt + 0.5*p.ax*dt**2
    p.y = p.y + p.vy*dt + 0.5*p.ay*dt**2
    p.z = p.z + p.vz*dt + 0.5*p.az*dt**2

    # angular roation update :
    p.thetax = p.thetax + p.omegax*dt + 0.5*p.alphax*dt**2
    p.thetay = p.thetay + p.omegay*dt + 0.5*p.alphay*dt**2
    p.thetaz = p.thetaz + p.omegaz*dt + 0.5*p.alphaz*dt**2

    # Update periodic BC
    p.pbcUpdate()

    # Store accelerations for averaging that is done 
    ax = p.ax
    ay = p.ay
    az = p.az
    
    # store angular accelerations for crank-nicolson scheme :
    alphax = p.alphax
    alphay = p.alphay
    alphaz = p.alphaz

    force(p) # Force update with new positions

    # Velocity updates
    p.vx = p.vx + 0.5*(ax + p.ax)*dt
    p.vy = p.vy + 0.5*(ay + p.ay)*dt
    p.vz = p.vz + 0.5*(az + p.az)*dt

    # angular velocity updates :
    p.omegax = p.omegax + 0.5*(alphax + p.alphax)*dt
    p.omegay = p.omegay + 0.5*(alphay + p.alphay)*dt
    p.omegaz = p.omegaz + 0.5*(alphaz + p.alphaz)*dt


class Particles(object):

  def __init__(self, L, force, periodicX=1, periodicY=1, periodicZ=1):
    # Container size
    self.L = L
    # Total Number of particles
    self.N = 0
    # type
    self.type = 'float32'
    # Positions
    self.x  = array([],dtype=self.type)
    self.y  = array([],dtype=self.type)
    self.z  = array([],dtype=self.type)
    # Velocities
    self.vx = array([],dtype=self.type)
    self.vy = array([],dtype=self.type)
    self.vz = array([],dtype=self.type)
    # Forces
    self.ax = array([],dtype=self.type)
    self.ay = array([],dtype=self.type)
    self.az = array([],dtype=self.type)
    # Radii
    self.r  = array([],dtype=self.type)
    # angular rotation :
    self.thetax = array([],dtype=self.type)
    self.thetay = array([],dtype=self.type)
    self.thetaz = array([],dtype=self.type)
    # angular velocity :
    self.omegax = array([],dtype=self.type)
    self.omegay = array([],dtype=self.type)
    self.omegaz = array([],dtype=self.type)
    # angular acceleration :
    self.alphax = array([],dtype=self.type)
    self.alphay = array([],dtype=self.type)
    self.alphaz = array([],dtype=self.type)
    # Periodic on?
    self.periodicX = periodicX 
    self.periodicY = periodicY 
    self.periodicZ = periodicZ 
    # Force function
    self.f = force
     
  def addParticle(self, x, y, z, vx, vy, vz, r,
                  thetax, thetay, thetaz, 
                  omegax, omegay, omegaz): 
    self.x  = hstack((self.x,x))
    self.y  = hstack((self.y,y))
    self.z  = hstack((self.z,z))
    self.vx = hstack((self.vx,vx))
    self.vy = hstack((self.vy,vy))
    self.vz = hstack((self.vz,vz))
    self.ax = hstack((self.ax,0))
    self.ay = hstack((self.ay,0))
    self.az = hstack((self.az,0))
    self.r  = hstack((self.r,r))
    self.thetax = hstack((self.thetax,thetax))
    self.thetay = hstack((self.thetay,thetay))
    self.thetaz = hstack((self.thetaz,thetaz))
    self.omegax = hstack((self.omegax,omegax))
    self.omegay = hstack((self.omegay,omegay))
    self.omegaz = hstack((self.omegaz,omegaz))
    self.alphax = hstack((self.alphax,0))
    self.alphay = hstack((self.alphay,0))
    self.alphaz = hstack((self.alphaz,0))
    self.N  = self.N+1
    temp    = tile(self.r,(self.N,1))
    self.sumOfRadii   = temp + temp.T
    self.ratioOfRadii = temp / temp.T
    self.f(self)

  def pbcUpdate(self):
    """
    Moves paricles across periodic boundary
    """
    if self.periodicX:
      self.x[self.x >  self.L/2] = self.x[self.x >  self.L/2] - self.L
      self.x[self.x < -self.L/2] = self.x[self.x < -self.L/2] + self.L
    if self.periodicY:
      self.y[self.y >  self.L/2] = self.y[self.y >  self.L/2] - self.L
      self.y[self.y < -self.L/2] = self.y[self.y < -self.L/2] + self.L
    if self.periodicZ:
      self.z[self.z >  self.L/2] = self.z[self.z >  self.L/2] - self.L
      self.z[self.z < -self.L/2] = self.z[self.z < -self.L/2] + self.L

  def distanceMatrix(self, x, y, z):
    """
    Computes distances between all particles and places the result in a 
    matrix such that the ij th matrix entry corresponds to the distance 
    between particle i and j
    """ 
    xtemp = tile(x, (self.N,1))
    ytemp = tile(y, (self.N,1))
    ztemp = tile(z, (self.N,1))
    dx    = xtemp - xtemp.T
    dy    = ytemp - ytemp.T
    dz    = ztemp - ztemp.T
  
    # Particles 'feel' each other across the periodic boundaries
    if self.periodicX:
      dx[dx >  self.L/2] = dx[dx >  self.L/2] - self.L
      dx[dx < -self.L/2] = dx[dx < -self.L/2] + self.L
    if self.periodicY:
      dy[dy >  self.L/2] = dy[dy >  self.L/2] - self.L
      dy[dy < -self.L/2] = dy[dy < -self.L/2] + self.L
    if self.periodicZ:
      dz[dz >  self.L/2] = dz[dz >  self.L/2] - self.L
      dz[dz < -self.L/2] = dz[dz < -self.L/2] + self.L

    # Total Distances
    d = sqrt(dx**2 + dy**2 + dz**2)

    # Mark zero entries with negative 1 to avoid divergences
    d[d==0] = -1

    return d, dx, dy, dz



