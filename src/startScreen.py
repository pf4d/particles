from OpenGL.GL          import *
from OpenGL.GLUT        import *
from OpenGL.GLE         import *
from OpenGL.GLU         import *

buttonX1 = 300
buttonX2 = 700

startY1 = 900
startY2 = 700

helpY1 = 600
helpY2 = 400

exitY1 = 300
exitY2 = 100

def display():
    glLoadIdentity();
    glClear(GL_COLOR_BUFFER_BIT)
    glColor(1, 0, 0)
    drawButton(buttonX1, buttonX2, startY1, startY2)
    printChars(buttonX1+10, startY2+40, GLUT_STROKE_ROMAN, 'START', 1, 1, 1)
    glColor(1, 0, 0)
    drawButton(buttonX1, buttonX2, helpY1, helpY2)
    printChars(buttonX1+10, helpY2+40, GLUT_STROKE_ROMAN, 'HELP', 1, 1, 1)
    glColor(1, 0, 0)
    drawButton(buttonX1, buttonX2, exitY1, exitY2)
    printChars(buttonX1+10, exitY2+40, GLUT_STROKE_ROMAN, 'EXIT', 1, 1, 1)
    glutSwapBuffers()

def keyboard(key, x, y):
    if key == 'q':
        exit(0)

def setProjection():
    glMatrixMode(GL_PROJECTION)
    glOrtho(0, 1000, 0, 1000, -1, 1)
    glMatrixMode(GL_MODELVIEW)

def drawButton(x1, x2, y1, y2):
    glBegin(GL_TRIANGLES)
    glVertex(x1, y1, 0)
    glVertex(x2, y1, 0)
    glVertex(x2, y2, 0)
    glVertex(x2, y2, 0)
    glVertex(x1, y2, 0)
    glVertex(x1, y1, 0)
    glEnd()

def printChars(x, y, font, text, r, g, b):
    glColor(r, g, b)
    glPushMatrix()
    #glRasterPos(x, y)
    glTranslate(x, y, 0)
    for char in text:
        #glutBitmapCharacter(font, ctypes.c_int( ord(char)))
        glutStrokeCharacter(font, ctypes.c_int(ord(char)))
    glPopMatrix()

if __name__ == '__main__':
    glutInit()
    glutInitWindowSize(1000, 1000)
    glutCreateWindow("hello")
    glutDisplayFunc(display)
    glutKeyboardFunc(keyboard)
    setProjection()
    glutMainLoop()
