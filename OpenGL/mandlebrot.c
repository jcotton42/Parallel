#include <stdio.h>
#include <GL/glut.h>
#include <GL/freeglut_ext.h>
#include <complex.h>
#define SX 800
#define SY 600
#define CXMIN -2.0
#define CXMAX 2.0
#define CYMIN -1.5
#define CYMAX 1.5

int cache[SX][SY];
int max_iteration = 100;

int cart2screenX(double cx) {
    return (int)((cx - CXMIN) / (CXMAX - CXMIN) * SX);
}

int cart2screenY(double cy) {
    return (int)((cy - CYMIN) / (CYMAX - CYMIN) * SY);
}

double screen2cartX(int sx) {
    return ((double)sx / SX) * (CXMAX - CXMIN) + CXMIN;
}

double screen2cartY(int sy) {
    return ((double)sy / SY) * (CYMAX - CYMIN) + CYMIN;
}

void initcache(void);
void mandlebrot(void);
void mousefunc(int button, int state, int x, int y);
void motionfunc(int x, int y);

int main(int argc, char* argv[]) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(SX, SY);
    glutInitWindowPosition(100, 50);
    glutCreateWindow("");
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glShadeModel(GL_SMOOTH);

    glViewport(0, 0, (GLsizei)SX, (GLsizei)SY);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0, 1.0 * SX, 0.0, 1.0 * SY);
    glMatrixMode(GL_MODELVIEW);

    initcache();

    glutDisplayFunc(mandlebrot);
    glutMouseFunc(mousefunc);
    glutMotionFunc(NULL);
    glutIdleFunc(NULL);
    glutReshapeFunc(NULL);
    glutCloseFunc(NULL);
    glutMainLoop();
    return 0;
}

void initcache(void) {
    int iteration = 0;
    double complex z, c;
    glClear(GL_COLOR_BUFFER_BIT);
    for(int xp = 0; xp < SX; xp++) {
        for(int yp = 0; yp < SY; yp++) {
            c = screen2cartX(xp) + screen2cartY(yp) * 1i;
            z = 0 + 0i;
            while(cabs(z) < 2 && iteration < max_iteration) {
                z = cpow(z, 2) + c;
                iteration++;
            }
            cache[xp][yp] = iteration;
            iteration = 0;
        }
    }
}

void mandlebrot(void) {
    for(int xp = 0; xp < SX; xp++) {
        for(int yp = 0; yp < SY; yp++) {
            if(cache[xp][yp] == max_iteration) {
                glColor3d(0.0, 0.0, 0.0);
            } else {
                double ratio = (double)cache[xp][yp] / max_iteration;
                glColor3d(ratio, ratio, ratio);
            }
            glBegin(GL_POINTS);
            glVertex2d(xp, yp);
            glEnd();
        }
    }
    glutSwapBuffers();
}

void mousefunc(int button, int state, int x, int y) {
    glutPostRedisplay();
    y = SY - y;
    if(state != GLUT_DOWN) return;
    double _Complex z, c;
    c = screen2cartX(x) + screen2cartY(y) * 1i;
    z = 0 + 0i;
    int iteration = 0;
    int max_iterations = 100;
    glColor3d(1.0, 1.0, 0.0);
    glBegin(GL_LINE_STRIP);
    glVertex2d(SX / 2, SY / 2);
    while(cabs(z) < 2 && iteration < max_iterations) {
        z = cpow(z, 2.0) + c;
        glVertex2d(cart2screenX(creal(z)), cart2screenY(cimag(z)));
    }
    glEnd();
    glutSwapBuffers();
}

void motionfunc(int x, int y) {

}
