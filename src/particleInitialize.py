from pylab import *
def particleInitialize(p,initialization,L):
  vel=1.0
  dist = L/10.
  if initialization == 'one':
    p.addParticle(0,dist,0,0,0,0,1)

  elif initialization == 'two':
    p.addParticle(-dist,0.,0.,vel,0.,0.,1.)
    p.addParticle(dist,0.,0,-vel,0.,0.,1.)

  elif initialization == 'three':
    p.addParticle(0.,dist*sqrt(3)/2,0.,0.,-vel,0.,1.)
    p.addParticle(-dist,0.,0.,vel,0.,0.,1.0)
    p.addParticle(dist,0.,0,-vel,0.,0.,1.)

  elif initialization == 'four':
    p.addParticle(-dist,0.,0.,vel,0.,0.,1.)
    p.addParticle(dist,0.,0,-vel,0.,0.,1.)
    p.addParticle(0.,dist,0,0.,-vel,0.,1.)
    p.addParticle(0.,-dist,0,0.,vel,0.,1.)

  elif initialization == 'six':
    p.addParticle(0.,dist,0.,0.,-vel,0.,1.)
    p.addParticle(0.,-dist,0.,0.,vel,0.,1.)
    p.addParticle(dist/sqrt(2),dist/sqrt(2),0,-vel/sqrt(2),-vel/sqrt(2),0.,1.)
    p.addParticle(-dist/sqrt(2),dist/sqrt(2),0,vel/sqrt(2),-vel/sqrt(2),0.,1.)
    p.addParticle(-dist/sqrt(2),-dist/sqrt(2),0,vel/sqrt(2),vel/sqrt(2),0.,1.)
    p.addParticle(dist/sqrt(2),-dist/sqrt(2),0,-vel/sqrt(2),vel/sqrt(2),0.,1.)       
  elif initialization == 'eight':
    p.addParticle(-dist,0.,0.,vel,0.,0.,1.)
    p.addParticle(dist,0.,0,-vel,0.,0.,1.)
    p.addParticle(0.,dist,0,0.,-vel,0.,1.)
    p.addParticle(0.,-dist,0,0.,vel,0.,1.)

    p.addParticle(dist/sqrt(2),dist/sqrt(2),0,-vel/sqrt(2),-vel/sqrt(2),0.,1.)
    p.addParticle(-dist/sqrt(2),dist/sqrt(2),0,vel/sqrt(2),-vel/sqrt(2),0.,1.)
    p.addParticle(-dist/sqrt(2),-dist/sqrt(2),0,vel/sqrt(2),vel/sqrt(2),0.,1.)
    p.addParticle(dist/sqrt(2),-dist/sqrt(2),0,-vel/sqrt(2),vel/sqrt(2),0.,1.)       
     
  else:
    p.addParticle(-dist,0.,0.,vel,0.,0.,1.)
    p.addParticle(dist,0.,0,-vel,0.,0.,1.)

