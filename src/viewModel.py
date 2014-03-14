import sys, pygame
from pygame.locals import *
from pygame.constants import *
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
from particles import *
import numpy
from numpy import *
from random import randint

#import obj loader
from objLoader import *

pygame.init()
dt = .01
radiusDiv = 1
SLICES = 30
STACKS = 30
viewport = (800,600)
hx = viewport[0]/2
hy = viewport[1]/2
srf = pygame.display.set_mode(viewport, OPENGL | DOUBLEBUF)
L = 10.0
px = 0
py = 2*L
pz = 0
pvx = 0
pvy = 0
pvz = 0
pr = 1
f = GranularMaterialForce(k=1.5, g=.000001, gamma=.1)
p = Particles(L, f, periodicY=0, periodicZ=1, periodicX=1)
integrate = VerletIntegrator(dt)
#p.addParticle(px,py,pz,pvx,pvy,pvz,pr,0,0,0,0,0,0)

def initParticles():
    #particle = p.addParticle(0,py,0,pvx,pvy,pvz,pr,0,0,0,0,0,0)
    for i in range(0,10):
        spx = randint(1,40)
        #spy = randint(1,40)
        #spz = randint(1,40)
        #p.addParticle(spx,py,spz,pvx,pvy,pvz,pr,0,0,0,0,0,0)

def drawParticles():
    glTranslate(p.x[0],p.y[0],p.z[0])
    glutSolidSphere(p.r[0]/radiusDiv, SLICES, STACKS)
    for i in range(0,10):
        glTranslate(randint(1,40),py,randint(1,40))
        glutSolidSphere(p.r[0]/radiusDiv, SLICES, STACKS)



#  def addParticle(self, x, y, z, vx, vy, vz, r,
#                  thetax, thetay, thetaz, 
#                  omegax, omegay, omegaz): 


#particle = p.addParticle(px,py,pz,pvx,pvy,pvz,pr,0,0,0,0,0,0)

glLightfv(GL_LIGHT0, GL_POSITION, (-40, 200, 100, 0.0))
glLightfv(GL_LIGHT0, GL_AMBIENT, (0.5, 0.5, 0.5, 1.0))
glLightfv(GL_LIGHT0, GL_DIFFUSE, (0.5, 0.5, 0.5, 1.0))
glEnable(GL_LIGHT0)
glEnable(GL_LIGHTING)
glEnable(GL_COLOR_MATERIAL)
glEnable(GL_DEPTH_TEST)
glShadeModel(GL_SMOOTH)
glutInit()

obj = OBJ(sys.argv[1], swapyz=True)

clock = pygame.time.Clock()

glMatrixMode(GL_PROJECTION)
glLoadIdentity()
width, height = viewport
gluPerspective(90.0, width/float(height), 1, 100.0)
glEnable(GL_DEPTH_TEST)
glMatrixMode(GL_MODELVIEW)

rx, ry = (0,0)
tx, ty = (0,0)
zpos = 5
rotate = move = False

p.addParticle(px,py,pz,pvx,pvy,pvz,pr,0,0,0,0,0,0)

#initParticles()
while 1:
    clock.tick(30)
    for e in pygame.event.get():
        if e.type == QUIT:
            sys.exit()
        elif e.type == KEYDOWN and e.key == K_ESCAPE:
            sys.exit()
        elif e.type == KEYDOWN and e.key == K_UP:
            #paz -= .2
            p.az[0] -= 10
        elif e.type == KEYDOWN and e.key == K_DOWN:
            #paz += .2
            p.az[0] += 10
        elif e.type == KEYDOWN and e.key == K_a:
            pvx -= .2
            p.ax[0] -= 10
        elif e.type == KEYDOWN and e.key == K_d:
            pvx += .2
            p.ax[0] += 10
        elif e.type == KEYDOWN and e.key == K_w:
            pvy += .2
            p.ay[0] += 10
        elif e.type == KEYDOWN and e.key == K_s:
            pvy -= .2
            p.ay[0] -= 10
        elif e.type == MOUSEBUTTONDOWN:
            if e.button == 4: zpos = max(1, zpos-1)
            elif e.button == 5: zpos += 1
            elif e.button == 1: rotate = True
            elif e.button == 3: move = True
        elif e.type == MOUSEBUTTONUP:
            if e.button == 1: rotate = False
            elif e.button == 3: move = False
        elif e.type == MOUSEMOTION:
            i, j = e.rel
            if rotate:
                rx += i
                ry += j
            if move:
                tx += i
                ty -= j

    integrate(f,p)
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
    glLoadIdentity()

    #px += pvx
    #py += pvy
    #pz += pvz

    #p.x[0] += pvx
    #p.y[0] += pvy
    #p.z[0] += pvz
    particlePosition = numpy.array([p.x[0],p.y[0],p.z[0]])
    cameraDist = particlePosition + numpy.array([20*cos(pvz),10,20*sin(pvz)])
    cameraDist = particlePosition + numpy.array([0,3,-10])
    cameraTarget = particlePosition
    cameraUp = numpy.array([0,1,0])

    print("x: %f  y: %f  z: %f" % (p.x[0],p.y[0],p.z[0]))

    # render object
    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()
    #gluLookAt(cameraDist, cameraTarget, cameraUp)
    gluLookAt(cameraDist[0],cameraDist[1],cameraDist[2],cameraTarget[0],cameraTarget[1],cameraTarget[2],cameraUp[0],cameraUp[1],cameraUp[2]) 

    glTranslate(p.x[0],p.y[0],p.z[0])
    glRotate(ry, 1, 0, 0)
    glRotate(rx, 0, 1, 0)
    glCallList(obj.gl_list)

    pygame.display.flip()


