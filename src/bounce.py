from particles          import *
from particleInitialize import *
from pylab              import *
from OpenGL.GL          import *
from OpenGL.GLUT        import *
from OpenGL.GLE         import *
from OpenGL.GLU         import *
import sys, time


# this program is a 'driver' for a simple simulation of partilces in a box with
# periodic boundary conditions. Your objective will be to complete the code here
# so that you can 'see' the particles with OpenGL.

tStart = t0 = time.time()

rotx = 0
roty = 0
rotz = 0

dt        = 0.01    # time step taken by the time integration routine.
L         = 10.0   # size of the box.
t         = 0      # initial time
vy        = 0      # vertical velocity
vx        = 0      # horizontal velocity
vz        = 0      # depth velocity

k         = 30.0   # elastic 'bounce'
gamma     = 0.1    # energy dissipation/loss
k         = 1.5    # elastic 'bounce'
gamma     = 0.1    # energy dissipation/loss
g         = 0.25   # downward acceleration

on        = True   # start / stop adding particles
trans     = False  # transparency enable
partInt   = 6/dt   # how often to add a new particle
radiusDiv = 1      # radius divisor
massive   = False  # the big ball.

# particle update data:
COUNT         = 1  # number of time steps computed
UPDATE_FRAMES = 2  # how often to redraw screen

# how resolved are the spheres?
STACKS = 30
SLICES = 30

# instantiate the forces function between particles
f = GranularMaterialForce(k=k, g=g, gamma=gamma)
# create some particles and a box
p = Particles(L, f, periodicY=0, periodicZ=1, periodicX=1)
#  def addParticle(self, x, y, z, vx, vy, vz, r,
#                  thetax, thetay, thetaz, 
#                  omegax, omegay, omegaz): 
p.addParticle(0,2*L,0,0,0,0,1,0,0,0,0,0,0)
# instantiate Integrator
integrate = VerletIntegrator(dt)

def init():
  glEnable(GL_COLOR_MATERIAL)
  glEnable(GL_BLEND)
  glShadeModel(GL_SMOOTH)
  #glShadeModel(GL_FLAT)
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
  glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE)

  glEnable(GL_POLYGON_OFFSET_FILL) # Prevents hidden line problems when drawing
  glPolygonOffset(1.0, 1.0)        # a wireframe on top of filled polygons.

  glEnable(GL_CULL_FACE)
  glEnable(GL_DEPTH_TEST)
  
  glEnable(GL_LIGHTING)
  glEnable(GL_LIGHT0)
  glEnable(GL_LIGHT1)
  
def display():
  glColor(0.5, 0.8, 0.0)
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

  glPushMatrix() 
  glMaterial(GL_FRONT, GL_EMISSION,  [0.0, 0.0, 0.0, 0.0])
  glMaterial(GL_FRONT, GL_SPECULAR,  [0.5, 0.5, 0.5, 0.0])
  glMaterial(GL_FRONT, GL_SHININESS, 100.0)
  
  for i in range(p.N):
    if trans:
      mag = sqrt(p.vx[i]**2 + p.vy[i]**2 + p.vz[i]**2) + 0.02
    else:
      mag = 1.0
  
    glLoadIdentity()
    gluLookAt(0,0,10,   # Camera Position
              0,0,0,    # Point the Camera looks at
              0,1,0)    # the Up-Vector
    
    glRotate(rotx,1,0,0)
    glRotate(roty,0,1,0)
    glRotate(rotz,0,0,1)
    if p.ax[i] >= 0.5 or p.ax[i] <= -0.5:
      glColor(p.r[i]/1.5, p.r[i]/1.5, p.r[i]/1.5, mag)
    elif p.vx[i] >= 0.5:
      glColor(0.8, 0.4, 0.0, mag)
    elif p.vx[i] <= -0.5:
      glColor(0.0, 0.5, 0.8, mag)
    else:
      glColor(p.r[i]/2, p.r[i]/2, 0.0, mag)
    glTranslate(p.x[i], p.y[i], p.z[i])
    glRotate(p.thetax[i],1,0,0)
    glRotate(p.thetay[i],0,1,0)
    glRotate(p.thetaz[i],0,0,1)
    glutSolidSphere(p.r[i]/radiusDiv, SLICES, STACKS)
    glColor(0.0,0.0,0.0,1.0)
    glutWireSphere(p.r[i]/radiusDiv*1.01, SLICES/3, STACKS/3)
    glColor(0.5, 0.5, 0.5)
  print p.alphax[0], p.alphay[0], p.alphaz[0]
  glPopMatrix()

  glLoadIdentity()
  gluLookAt(0,0,10,   # Camera Position
            0,0,0,    # Point the Camera looks at
            0,1,0)    # the Up-Vector
  
  glRotatef(rotx,1,0,0)
  glRotatef(roty,0,1,0)
  glRotatef(rotz,0,0,1)
  
  glPushMatrix()
  
  lx1 = 0.0
  ly1 = 2*L + 2
  lz1 = 0.0
  
  lx2 = L
  ly2 = 2*L + 2
  lz2 = L
 
  glMaterial(GL_FRONT, GL_EMISSION,  [1.0, 1.0, 1.0, 0.0])
  glTranslate(lx1, ly1, lz1)
  glutSolidSphere(0.5, SLICES, STACKS)
  glTranslate(lx2, ly2, lz2)
  glutSolidSphere(0.5, SLICES, STACKS)
  glPopMatrix()

  # spot light :
  glLight(GL_LIGHT0, GL_AMBIENT,               [0.2,  0.2, 0.2, 1.0])
  glLight(GL_LIGHT0, GL_SPECULAR,              [1.0,  1.0, 1.0, 1.0])
  glLight(GL_LIGHT0, GL_DIFFUSE,               [1.0,  1.0, 1.0, 1.0])
  glLight(GL_LIGHT0, GL_POSITION,              [lx1,  ly1, lz1, 1.0])
  glLight(GL_LIGHT0, GL_SPOT_DIRECTION,        [0.0, -1.0, 0.0])
  glLight(GL_LIGHT0, GL_SPOT_CUTOFF,           10.0)
  glLight(GL_LIGHT0, GL_SPOT_EXPONENT,          3.0)
  glLight(GL_LIGHT0, GL_CONSTANT_ATTENUATION,   1.2)
  glLight(GL_LIGHT0, GL_LINEAR_ATTENUATION,     0.0)
  glLight(GL_LIGHT0, GL_QUADRATIC_ATTENUATION,  0.002)
  
  # regular light
  glLight(GL_LIGHT1, GL_AMBIENT,               [0.2,  0.2, 0.2, 1.0])
  glLight(GL_LIGHT1, GL_SPECULAR,              [1.0,  1.0, 1.0, 1.0])
  glLight(GL_LIGHT1, GL_DIFFUSE,               [1.0,  1.0, 1.0, 1.0])
  glLight(GL_LIGHT1, GL_POSITION,              [lx2,  ly2, lz2, 1.0])
  glLight(GL_LIGHT1, GL_CONSTANT_ATTENUATION,   2.0)
  glLight(GL_LIGHT1, GL_LINEAR_ATTENUATION,     0.0)
  glLight(GL_LIGHT1, GL_QUADRATIC_ATTENUATION,  0.0)
  
  glutSwapBuffers()
  glFlush()

def reshape(width, height):
  glViewport(0, 0, width, height)
  glMatrixMode(GL_PROJECTION)
  glLoadIdentity()
  glOrtho(-L, L, -L, L, -4*L, 4*L)
  glMatrixMode(GL_MODELVIEW)
  glLoadIdentity()

def idle():
  global COUNT, vy, vx, vz, massive
  for i in range(UPDATE_FRAMES):
    integrate(f,p) # Move the system forward in time
    COUNT = COUNT + 1 
    if mod(COUNT,partInt) == 0:
      # syntax is addParticle(x,y,z,vx,vy,vz,radius)
      # note y is into page.
      if massive:
        r = L/4
      else:
        r = 0.3*randn() + 1.
      if on:
        px = 0.25*randn()
        py = 2*L
        pz = 0
        p.addParticle(px, py, pz, vx, vy, vz, r,0,0,0,0,0,0)
  glutPostRedisplay()

def key(k, x, y):
  global trans, on, radiusDiv, massive

  if k == 'q':
    print "'q' was pressed"
    exit(0)
  
  if k == 't':
    if trans == True:
      trans = False
    else:
      trans = True
    print "'t' was pressed: Trans =", trans
    
  if k == 'o':
    if on == True:
      on = False
    else:
      on = True
    print "'o' was pressed: On =", on
    
  if k == '=':
    radiusDiv += 0.05
    print "'+' was pressed: radiusDiv =", radiusDiv
  
  if k == '-':
    if radiusDiv-1 == 0:
      print "RADIUS DIVISOR AT MIN VALUE"
    else:
      radiusDiv -= 0.05
      print "'-' was pressed: radiusDiv =", radiusDiv
  
  if k == 'm':
    if massive == True:
      massive = False
    else:
      massive = True
    print "'m' was pressed: massive =", massive
  
  if k == 'n':
    print "'n' was pressed: n =", p.N
    

def special(k, x, y):
  global vy, vx, vz, partInt
  
  if k == GLUT_KEY_UP:
    vy += .5
    print 'UP    key was pressed: vy =', vy
  
  if k == GLUT_KEY_DOWN:
    vy -= .5
    print 'DOWN  key was pressed: vy =', vy
  
  if k == GLUT_KEY_RIGHT:
    vx += .5
    print 'RIGHT key was pressed: vx =', vx
  
  if k == GLUT_KEY_LEFT:
    vx -= .5
    print 'LEFT  key was pressed: vx =', vx
  
  if k == GLUT_KEY_PAGE_UP:
    vz += .5
    print 'PGUP  key was pressed: vz =', vz
  
  if k == GLUT_KEY_PAGE_DOWN:
    vz -= .5
    print 'PGDWN key was pressed: vz =', vz
  
  if k == GLUT_KEY_HOME:
    partInt += 1/dt
    print 'HOME  key was pressed: partInt =', partInt
  
  if k == GLUT_KEY_INSERT:
    if partInt - 1/dt <= 0:
      print "ADD PARTICLE INTERVAL DIVISOR AT MIN VALUE"
    else:
      partInt -= 1/dt
      print 'INS   key was pressed: partInt =', partInt

def mouse(button,state,x,y):
  global beginx,beginy,rotate
  if button == GLUT_LEFT_BUTTON and state == GLUT_DOWN:
    print "Mouseclick: ",x,"x> ",y,"yv"
    rotate = 1
    beginx = x
    beginy = y
  if button == GLUT_LEFT_BUTTON and state == GLUT_UP:
    rotate = 0

def motion(x,y):
  global rotx,roty,beginx,beginy,rotate
  if rotate:
    rotx = rotx + (y - beginy)
    roty = roty + (x - beginx)
    beginx = x
    beginy = y
    glutPostRedisplay()
  

if __name__ == '__main__':

    i      = 80
    width  = i*int(L)
    height = i*int(L)
    
    sx = 600
    sy = 250

    # open a window
    glutInit(sys.argv)
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH)
    glutInitWindowPosition(sx, sy)
    glutInitWindowSize(width, height)
    glutCreateWindow("bounce")
    glutDisplayFunc(display)
    glutMouseFunc(mouse)
    glutMotionFunc(motion)
    glutReshapeFunc(reshape)
    glutIdleFunc(idle)
    glutKeyboardFunc(key)
    glutSpecialFunc(special)
  
    
    # initialize
    init()

    # hand off control to event loop
    glutMainLoop()


