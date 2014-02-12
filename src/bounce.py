from particles          import *
from particleInitialize import *
from pylab              import *
from OpenGL.GL          import *
from OpenGL.GLUT        import *
import sys, time


# this program is a 'driver' for a simple simulation of partilces in a box with
# periodic boundary conditions. Your objective will be to complete the code here
# so that you can 'see' the particles with OpenGL.

tStart = t0 = time.time()

dt        = 0.1    # time step taken by the time integration routine.
L         = 10.0   # size of the box.
t         = 0      # initial time
vy        = 0      # vertical velocity
vx        = 0      # horizontal velocity
vz        = 0      # depth velocity

on        = True   # start / stop adding particles
trans     = False  # transparency enable
partInt   = 400    # how often to add a new particle
radiusDiv = 1      # radius divisor
massive   = False  # the big ball.

# particle update data:
COUNT         = 1  # number of time steps computed
UPDATE_FRAMES = 2  # how often to redraw screen

# how resolved are the spheres?
STACKS = 10
SLICES = 25

# instantiate the forces function between particles
f = GranularMaterialForce(g=0.1)
# create some particles and a box
p = Particles(L, f, periodicY=0, periodicZ=1, periodicX=1)
p.addParticle(0,L,L/2,0,0,0,1)
# instantiate Integrator
integrate = VerletIntegrator(dt)

def init():
  glutReshapeFunc(reshape)
  glutIdleFunc(idle)
  glutKeyboardFunc(key)
  glutSpecialFunc(special)
  glutDisplayFunc(draw)
  
  glEnable(GL_LIGHTING)
  glEnable(GL_LIGHT0)
  glEnable(GL_COLOR_MATERIAL)
  glEnable(GL_BLEND)
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
  glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE)

  glEnable(GL_POLYGON_OFFSET_FILL)  # Prevents hidden line problems when drawing
  glPolygonOffset(1.0, 1.0)         # a wireframe on top of filled polygons.

  glEnable(GL_DEPTH_TEST)
  
def draw():
  glColor(0.5, 0.8, 0.0)
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
  for i in range(p.N):
    if trans:
      mag = sqrt(p.vx[i]**2 + p.vy[i]**2 + p.vz[i]**2) + 0.02
    else:
      mag = 1.0
    
    glLoadIdentity()
    if p.ax[i] >= 0.5 or p.ax[i] <= -0.5:
      glColor(p.r[i]/1.5, p.r[i]/1.5, p.r[i]/1.5, mag)
    elif p.vx[i] >= 0.5:
      glColor(0.8, 0.4, 0.0, mag)
    elif p.vx[i] <= -0.5:
      glColor(0.0, 0.5, 0.8, mag)
    else:
      glColor(p.r[i]/2, p.r[i]/2, 0.0, mag)
    glTranslate(p.x[i], p.y[i], p.z[i])
    glutSolidSphere(p.r[i]/radiusDiv, SLICES, STACKS)
    glColor(0.5, 0.5, 0.5)
  
  glutSwapBuffers()

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
        p.addParticle(.25*randn(), L, L/2, vx, vy, vz, r)
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
    partInt += 10
    print 'HOME  key was pressed: partInt =', partInt
  
  if k == GLUT_KEY_INSERT:
    if partInt-10 == 0:
      print "ADD PARTICLE INTERVAL DIVISOR AT MIN VALUE"
    else:
      partInt -= 10
      print 'INS   key was pressed: partInt =', partInt

def reshape(width, height):
  glViewport(0,0,width,height)
  glMatrixMode(GL_PROJECTION)
  glLoadIdentity()
  glOrtho(-L/1.5, L/1.5, -L/1.5, L/1.5, -L, L)
  glMatrixMode(GL_MODELVIEW)
  glShadeModel(GL_SMOOTH)
  #glLight(GL_LIGHT0, GL_AMBIENT, [0.2, 0.2, 0.2 , 1.0])
  #glLight(GL_LIGHT0, GL_DIFFUSE, [0.8, 0.8, 0.8 , 1.0])
  #glLight(GL_LIGHT0, GL_SPECULAR, [1.0, 1.0, 1.0 , 1.0])
  glLight(GL_LIGHT0, GL_POSITION, [0.7, 1.0, -1.0, 0.0])
  

def visible(vis):
  None

if __name__ == '__main__':

    width  = 500
    height = 500
    
    sx = 600
    sy = 250

    # open a window
    glutInit(sys.argv)
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE)
    glutInitWindowPosition(sx, sy)
    glutInitWindowSize(width, height)
    glutCreateWindow("bounce")
    
    # initialize
    init()

    # hand off control to event loop
    glutMainLoop()


