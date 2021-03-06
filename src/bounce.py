from particles          import *
from particleInitialize import *
from pylab              import *
from OpenGL.GL          import *
from OpenGL.GLUT        import *
from OpenGL.GLE         import *
from OpenGL.GLU         import *
from FTGL               import *
from time               import time
import sys

rotx      = 0      # camera x rotation
roty      = 0      # camera y rotation
rotz      = 0      # camera z rotation

frames    = 0      # for spf calculation
lastTime  = time() # current time
fps       = 1.0    # current frames per second

dt        = 0.10   # time step taken by the time integration routine.
L         = 10.0   # size of the box.
t         = 0      # initial time
vy        = 0      # vertical velocity
vx        = 0      # horizontal velocity
vz        = 0      # depth velocity

k         = 30.0   # elastic 'bounce'
gamma     = 1.0    # energy dissipation/loss
#k         = 1.0    # elastic 'bounce'
#gamma     = 0.5    # energy dissipation/loss
g         = 0.25   # downward acceleration

on        = False  # start / stop adding particles
trans     = False  # transparency enable
partInt   = 6/dt   # how often to add a new particle
radiusDiv = 1      # radius divisor
massive   = False  # the big ball.

# particle update data:
COUNT         = 1  # number of time steps computed
UPDATE_FRAMES = 1  # how often to redraw screen

# how resolved are the spheres?
STACKS = 30
SLICES = 30

# instantiate the forces function between particles
f = GranularMaterialForce(k=k, g=g, gamma=gamma)
# create some particles and a box
p = Particles(L, f, periodicY=0, periodicZ=1, periodicX=1)
#  addParticle(x, y, z, vx, vy, vz, r,
#              thetax, thetay, thetaz, 
#              omegax, omegay, omegaz): 
p.addParticle(0, L,0,0,0,0,1.0/2,0,0,0,0,0,0)
#p.addParticle(4,  -L/2+1/2.0,0,0,0,0,1.0/2,0,0,0,1,1,0)
#p.addParticle(2,  -L/2+1/2.0,0,0,0,0,1.0/2,0,0,0,1,0,0)
#p.addParticle(0,  -L/2+1/2.0,0,0,0,0,1.0/2,0,0,0,0,1,0)
#p.addParticle(-2, -L/2+1/2.0,0,0,0,0,1.0/2,0,0,0,0,0,1)

# instantiate Integrator
integrate = VerletIntegrator(dt)

def init():
  glEnable(GL_COLOR_MATERIAL)
  glEnable(GL_BLEND)
  glShadeModel(GL_SMOOTH)
  glShadeModel(GL_FLAT)
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
  glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE)

  glEnable(GL_POLYGON_OFFSET_FILL) # Prevents hidden line problems when drawing
  glPolygonOffset(1.0, 1.0)        # a wireframe on top of filled polygons.

  glEnable(GL_CULL_FACE)
  glEnable(GL_DEPTH_TEST)
  
  glEnable(GL_LIGHTING)
  glEnable(GL_LIGHT0)
  glEnable(GL_LIGHT1)

def print_p0_stats(dx, dy):
  # disable lighting :
  glDisable(GL_LIGHTING)
  
  # print statistics :
  glPushMatrix()
  glLoadIdentity()
  
  glColor(1.0,1.0,1.0,1.0) 
  font = BitmapFont('ProggySquareSZ.ttf')
  font.FaceSize(16)
  
  t1 = 'red particle statistics (x,y,z) :'
  t2 = 'theta        : %.2E, %.2E, %.2E' % (p.thetax[0],p.thetay[0],p.thetaz[0])
  t3 = 'omega        : %.2E, %.2E, %.2E' % (p.omegax[0],p.omegay[0],p.omegaz[0])
  t4 = 'alpha        : %.2E, %.2E, %.2E' % (p.alphax[0],p.alphay[0],p.alphaz[0])
  t5 = 'position     : %.2E, %.2E, %.2E' % (p.x[0],  p.y[0],  p.z[0])
  t6 = 'velocity     : %.2E, %.2E, %.2E' % (p.vx[0], p.vy[0], p.vz[0])
  t7 = 'acceleration : %.2E, %.2E, %.2E' % (p.ax[0], p.ay[0], p.az[0])
  glRasterPos2f(-L+dy,-(L-1)+dy*3)
  font.Render(t1)
  glRasterPos2f(-L+dy,-(L-1)+dy*2.5)
  font.Render(t2)
  glRasterPos2f(-L+dy,-(L-1)+dy*2)
  font.Render(t3)
  glRasterPos2f(-L+dy,-(L-1)+dy*1.5)
  font.Render(t4)
  glRasterPos2f(-L+dy,-(L-1)+dy*1.0)
  font.Render(t5)
  glRasterPos2f(-L+dy,-(L-1)+dy*0.5)
  font.Render(t6)
  glRasterPos2f(-L+dy,-(L-1)+dy*0.0)
  font.Render(t7)
  glPopMatrix()
  
  # re-enable lighting :
  glEnable(GL_LIGHTING)
  
def print_stats(dx, dy):
  # save projection matrix
  glMatrixMode(GL_PROJECTION)
  glPushMatrix()

  # switch to orthographic projection :
  glLoadIdentity()
  glOrtho(-L, L, -L, L, -4*L, 4*L)

  # back to the modelview matrix mode, so that we can translate/scale text :
  glMatrixMode(GL_MODELVIEW)

  glDisable(GL_LIGHTING)   # disable lighting
  
  # print statistics :
  glPushMatrix()
  glLoadIdentity()
  
  glColor(1.0,1.0,1.0,1.0) 
  glRasterPos2f(L-dx, L-dy)
  font = BitmapFont('ProggySquareSZ.ttf')
  font.FaceSize(16)
  #font = TextureFont('ProggySquareSZ.ttf')
  #font.FaceSize(13)
  #glScale(0.05, 0.05, 0.05)
  font.Render("n = %i" % p.N)
  glRasterPos2f(-L+dy, L-dy)
  font.Render("%i FPS" % fps)
  glPopMatrix()
  
  # re-enable lighting :
  glEnable(GL_LIGHTING)

  # get back to old perspective matrix :
  glMatrixMode(GL_PROJECTION)
  glPopMatrix()
  glMatrixMode(GL_MODELVIEW)

def draw_velocity_vectors():
  """
  draw velocity vectors.
  """
  glLineWidth(1.0)
  glPopMatrix()  
  glColor4f(1.0,1.0,1.0,1.0)
  glDisable(GL_LIGHTING)
  glBegin(GL_LINES)
  for i in range(p.N):
    v_mag = sqrt(p.vx[i]**2 + p.vy[i]**2 + p.vz[i]**2) + 1e-16
    xyz1 = array([p.x[i],  p.y[i],  p.z[i]])
    vxyz = array([p.vx[i], p.vy[i], p.vz[i]])
    vxyz = vxyz / v_mag * 2.0*p.r[i]
    xyz2 = xyz1 + vxyz
    glVertex3fv(xyz1)
    glVertex3fv(xyz2)
  glEnd()
  glEnable(GL_LIGHTING)
  glPushMatrix()

def draw_acceleration_vectors():
  """
  draw acceleration vectors.
  """
  glLineWidth(1.0)
  glPopMatrix()  
  glColor4f(1.0,0.0,0.0,1.0)
  glDisable(GL_LIGHTING)
  glBegin(GL_LINES)
  for i in range(p.N):
    a_mag = sqrt(p.ax[i]**2 + p.ay[i]**2 + p.az[i]**2) + 1e-16
    xyz1 = array([p.x[i],  p.y[i],  p.z[i]])
    axyz = array([p.ax[i], p.ay[i], p.az[i]])
    axyz = axyz / a_mag * 2.0*p.r[i]
    xyz2 = xyz1 + axyz
    glVertex3fv(xyz1)
    glVertex3fv(xyz2)
  glEnd()
  glEnable(GL_LIGHTING)
  glPushMatrix()

def draw_rotation_vectors():
  """
  draw rotation vectors.
  """
  glLineWidth(1.0)
  glPopMatrix()  
  glColor4f(0.0,1.0,1.0,1.0)
  glDisable(GL_LIGHTING)
  glBegin(GL_LINES)
  for i in range(p.N):
    xyz1  = array([p.x[i],  p.y[i],  p.z[i]])
    fxyz  = array([0, 0, 1])            # initial forward vector
    #uxyz  = array([0, 1, 0])            # initial up vector
    uxyz  = array([p.omegax[i], p.omegay[i], p.omegaz[i]])
    txyz  = array([p.thetax[i], p.thetay[i], p.thetaz[i]])
    fxyz  = rotate_vector(fxyz, txyz)
    uxyz  = rotate_vector(uxyz, txyz)
    fxyz  = fxyz / norm(fxyz) * 2.0*p.r[i]
    uxyz  = uxyz / norm(uxyz) * 2.0*p.r[i]
    fxyz2 = xyz1 + fxyz
    uxyz2 = xyz1 + uxyz
    glVertex3fv(xyz1)
    glVertex3fv(fxyz2)
    glVertex3fv(xyz1)
    glVertex3fv(uxyz2)
  glEnd()
  glEnable(GL_LIGHTING)
  glPushMatrix()

def draw_angular_velocity_vectors():  
  """
  draw angular velocity vectors.
  """
  glLineWidth(1.0)
  glPopMatrix()
  glColor4f(0.0,1.0,0.0,1.0)
  glDisable(GL_LIGHTING)
  glBegin(GL_LINES)
  for i in range(p.N):
    omega_mag = sqrt(p.omegax[i]**2 + p.omegay[i]**2 + p.omegaz[i]**2) + 1e-16
    xyz1 = array([p.x[i],      p.y[i],      p.z[i]])
    vxyz = -array([p.omegax[i], p.omegay[i], p.omegaz[i]]) 
    vxyz = vxyz / omega_mag * (p.r[i]+0.5)
    xyz2 = xyz1 + vxyz
    glVertex3fv(xyz1)
    glVertex3fv(xyz2)
  glEnd()
  glEnable(GL_LIGHTING)
  glPushMatrix()
  
def draw_angular_acceleration_vectors():
  """
  draw angular acceleration vectors.
  """ 
  glLineWidth(1.0)
  glPopMatrix()
  glColor4f(1.0,1.0,0.0,1.0)
  glDisable(GL_LIGHTING)
  glBegin(GL_LINES)
  for i in range(p.N):
    alpha_mag = sqrt(p.alphax[i]**2 + p.alphay[i]**2 + p.alphaz[i]**2) + 1e-16
    xyz1 = array([p.x[i],      p.y[i],      p.z[i]])
    vxyz = -array([p.alphax[i], p.alphay[i], p.alphaz[i]])
    vxyz = vxyz / alpha_mag * (p.r[i]+0.5)
    xyz2 = xyz1 + vxyz
    glVertex3fv(xyz1)
    glVertex3fv(xyz2)
  glEnd()
  glEnable(GL_LIGHTING)
  glPushMatrix()

def rotate_vector(v, r):
  """
  rotate vector <v> about the x, y, and z axes by angles provided in <r> array.
  """
  rx = r[0]
  ry = r[1]
  rz = r[2]
  c  = cos(rx)
  s  = sin(rx)
  Rx = array([[1, 0,  0],
              [0, c, -s],
              [0, s,  c]])
  c  = cos(ry)
  s  = sin(ry)
  Ry = array([[ c, 0, s],
              [ 0, 1, 0],
              [-s, 0, c]])
  c  = cos(rz)
  s  = sin(rz)
  Rz = array([[c, -s, 0],
              [s,  c, 0],
              [0,  0, 1]])
  R  = dot(Rx, dot(Ry, Rz))
  vn = dot(R, v)
  return vn

def Rx(theta):
  c  = cos(theta)
  s  = sin(theta)
  Rx = array([[1, 0,  0],
              [0, c, -s],
              [0, s,  c]])
  return Rx

def Ry(theta):
  c  = cos(theta)
  s  = sin(theta)
  Ry = array([[ c, 0, s],
              [ 0, 1, 0],
              [-s, 0, c]])
  return Ry

def Rz(theta):
  c  = cos(theta)
  s  = sin(theta)
  Rz = array([[c, -s, 0],
              [s,  c, 0],
              [0,  0, 1]])
  return Rz


def display():
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
 
  # camera viewpoint :
  glLoadIdentity()
  gluLookAt(0,0,10,   # Camera Position
            0,0,0,    # Point the Camera looks at
            0,1,0)    # the Up-Vector
  
  glRotate(rotx,1,0,0)
  glRotate(roty,0,1,0)
  glRotate(rotz,0,0,1)
  
  # interval to space text : 
  dx = 0.2 * L
  dy = 0.1 * L

  # print statistics :
  print_p0_stats(dx,dy)
  print_stats(dx,dy)
  
  # draw the spheres :
  glPushMatrix() 
  glColor(0.6, 0.1, 0.1)
  glMaterial(GL_FRONT, GL_EMISSION,  [0.0, 0.0, 0.0, 0.0])
  glMaterial(GL_FRONT, GL_SPECULAR,  [0.5, 0.5, 0.5, 0.0])
  glMaterial(GL_FRONT, GL_SHININESS, 100.0)
  for i in range(p.N):
    if trans:
      mag = sqrt(p.vx[i]**2 + p.vy[i]**2 + p.vz[i]**2) + 0.02
    else:
      mag = 1.0
    
    #if abs(p.alphax[i]) > 0.01 or \
    #   abs(p.alphay[i]) > 0.01 or \
    #   abs(p.alphaz[i]) > 0.01:
    #  glColor(p.r[i]/2.0, p.r[i]/4.0, p.r[i]/2.0, mag)
    if (p.ax[i] > 0.5 or p.ax[i] < -0.5) and i != 0:
      #glColor(p.r[i]/2.0, p.r[i]/2.0, p.r[i]/2.0, mag)
      glColor(1/2.0, 1/2.0, 1/2.0, mag)
    #elif (p.vx[i] > 0.5 or p.vy[i] > 0.5 or p.vz[i] > 0.5) and i != 0:
    #  glColor(0.8, 0.4, 0.0, mag)
    #elif p.vx[i] > 0.5:
    #  glColor(0.8, 0.4, 0.0, mag)
    #elif p.vx[i] < -0.5:
    #  glColor(0.0, 0.5, 0.8, mag)
    elif i != 0:
      #glColor(p.r[i]/1, p.r[i]/2, 0.0, mag)
      glColor(1, 1/2.0, 0.0, mag)
    
    glPushMatrix()
    glTranslate(p.x[i], p.y[i], p.z[i])
    
    # rotation :
    mvm = glGetFloatv(GL_MODELVIEW_MATRIX)
    M   = dot(p.theta[i], mvm[:3,:3])
    mvm[:3,:3] = M
    glLoadMatrixf(mvm)

    # draw spehere :
    glMaterial(GL_FRONT, GL_SPECULAR,  [0.5, 0.5, 0.5, 0.0])
    glMaterial(GL_FRONT, GL_SHININESS, 100.0)
    glutSolidSphere(p.r[i]/radiusDiv, SLICES, STACKS)
    glColor(0.0,0.0,0.0,1.0)
    glMaterial(GL_FRONT, GL_SPECULAR,  [0.0, 0.0, 0.0, 0.0])
    glMaterial(GL_FRONT, GL_SHININESS, 0.0)
    glutWireSphere(p.r[i]/radiusDiv*1.01, SLICES/6, STACKS/6)
    glPopMatrix()

  # draw vectors on particles :
  #draw_velocity_vectors()
  #draw_acceleration_vectors()
  #draw_rotation_vectors()
  draw_angular_velocity_vectors()
  #draw_angular_acceleration_vectors()   
   
  # draw the lights : 
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
  global COUNT, vy, vx, vz, massive, frames, lastTime, fps
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
        r = 1.0/2
      if on:
        px = 0.25*randn()
        py = 2*L
        pz = 0
        p.addParticle(px, py, pz, vx, vy, vz, r,0,0,0,0,0,0)
  glutPostRedisplay()

  # calculate fps :
  currentTime = time()
  frames     += 1
  if (currentTime - lastTime) >= 1.0:
    fps       = frames
    frames    = 0
    lastTime += 1.0
    #print fps

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

    i      = 70
    width  = i*int(L)
    height = i*int(L)
    
    sx = 600# + 1920
    sy = 300# + 100

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


