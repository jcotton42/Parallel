#include <iostream>
#include <complex>
#include <GL/glut.h>
#include <GL/freeglut_ext.h>
#include <mpi.h>

using std::complex;
using std::pair;
using std::cout;
using std::endl;

int screenWidth = 800;
int screenHeight = 600;
double CXMIN = -2.0;
double CXMAX = 2.0;
double CYMIN = -1.5;
double CYMAX = 1.5;

int **cache = nullptr;
int max_iteration = 100;
int clickedX;
int clickedY;
int dragX;
int dragY;
int drawMode;
bool rebuildCache;
enum {
    DRAW_NONE,
    DRAW_LINE,
    DRAW_RECT
};

enum {
    KILL,
    BUILD_CACHE,
    RESHAPED,
    ZOOMED
};

int size;
int rank;
MPI_Status status;

complex<double> screen2cart(int sx, int sy);
int cart2screenX(complex<double> c);
int cart2screenY(complex<double> c);

void display();
void mousefunc(int button, int state, int x, int y);
void motionfunc(int x, int y);
void drawPointPath();
void worker();
void manager(int* argc, char* argv[]);
void reshape(int width, int height);
void close();

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0) manager(&argc, argv);
    else worker();
    MPI_Finalize();
}

void worker() {
    int iteration = 0;
    int* result = new int[screenWidth];
    complex<double> z(0, 0);
    complex<double> c;

    while(true) {
        MPI_Recv(nullptr , 0, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        switch(status.MPI_TAG) {
            case BUILD_CACHE:
                int row;
                MPI_Recv(&row, 1, MPI_INT, 0, BUILD_CACHE, MPI_COMM_WORLD, &status);;
                for(int col = 0; col < screenWidth; col++) {
                    z = 0;
                    c = screen2cart(col, row);
                    while(std::abs(z) < 2 && iteration < max_iteration) {
                        z = std::pow(z, 2) + c;
                        iteration++;
                    }
                    result[col] = iteration == max_iteration ? 0 : iteration;
                    iteration = 0;
                }
                MPI_Send(&row, 1, MPI_INT, 0, BUILD_CACHE, MPI_COMM_WORLD);
                MPI_Send(result, screenWidth, MPI_INT, 0, BUILD_CACHE, MPI_COMM_WORLD);
                break;
            case RESHAPED:
                MPI_Recv(&screenWidth, 1, MPI_INT, 0, RESHAPED, MPI_COMM_WORLD, &status);
                MPI_Recv(&screenHeight, 1, MPI_INT, 0, RESHAPED, MPI_COMM_WORLD, &status);
                break;
            case ZOOMED:
                MPI_Recv(&CXMIN, 1, MPI_DOUBLE, 0, ZOOMED, MPI_COMM_WORLD, &status);
                MPI_Recv(&CXMAX, 1, MPI_DOUBLE, 0, ZOOMED, MPI_COMM_WORLD, &status);
                MPI_Recv(&CYMIN, 1, MPI_DOUBLE, 0, ZOOMED, MPI_COMM_WORLD, &status);
                MPI_Recv(&CYMAX, 1, MPI_DOUBLE, 0, ZOOMED, MPI_COMM_WORLD, &status);
                break;
            case KILL:
                return;
        }
    }
}

void manager(int* argc, char* argv[]) {
    glutInit(argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(screenWidth, screenHeight);
    glutInitWindowPosition(100, 50);
    glutCreateWindow("");
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glShadeModel(GL_SMOOTH);

    glViewport(0, 0, screenWidth, screenHeight);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0, 1.0 * screenWidth, 0.0, 1.0 * screenHeight);
    glMatrixMode(GL_MODELVIEW);

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutCloseFunc(close);
    glutMainLoop();
}

void display() {
    int row;
    if(rebuildCache) {
        int nextRank = 1;
        for(row = 0; row < screenHeight; row++) {
            if(nextRank == size) nextRank = 1;
            // this line may cause slowdowns...
            MPI_Send(nullptr, 0, MPI_INT, nextRank, BUILD_CACHE, MPI_COMM_WORLD);
            MPI_Send(&row, 1, MPI_INT, nextRank, BUILD_CACHE, MPI_COMM_WORLD);
        }
    }
    // TODO if there is no "progress" behavior, move glBegin/glEnd inside the for loop
    double ratio;
    glBegin(GL_POINTS);
    for(int i = 0; i < screenHeight; i++) {
        if(rebuildCache) {
            MPI_Recv(&row, 1, MPI_INT, MPI_ANY_SOURCE, BUILD_CACHE, MPI_COMM_WORLD, &status);
            MPI_Recv(cache[row], screenWidth, MPI_INT, status.MPI_SOURCE, BUILD_CACHE, MPI_COMM_WORLD, &status);
        } else {
            row = i;
        }
        for(int col = 0; col < screenWidth; col++) {
            ratio = cache[row][col] * 1.0 / max_iteration;
            glColor3d(ratio, ratio, ratio);
            glVertex2i(col, row);
        }
    }
    glEnd();
    glutSwapBuffers();
}

void reshape(int width, int height) {
    if(cache)
        for(int i = 0; i < screenHeight; i++)
            delete[] cache[i];
    delete[] cache;
    cache = new int*[height];
    for(int i = 0; i < height; i++)
        cache[i] = new int[width];
    screenWidth = width;
    screenHeight = height;
    for(int i = 1; i < size; i++) {
        MPI_Send(nullptr, 0, MPI_INT, i, RESHAPED, MPI_COMM_WORLD);
        MPI_Send(&width, 1, MPI_INT, i, RESHAPED, MPI_COMM_WORLD);
        MPI_Send(&height, 1, MPI_INT, i, RESHAPED, MPI_COMM_WORLD);
    }
    rebuildCache = true;
    glutPostRedisplay();
}

void close() {
    for(int i = 1; i < size; i++)
        MPI_Send(nullptr, 0, MPI_INT, i, KILL, MPI_COMM_WORLD);
    MPI_Finalize();
}

complex<double> screen2cart(int sx, int sy) {
    return complex<double>(
        ((double)sx / screenWidth) * (CXMAX - CXMIN) + CXMIN,
        ((double)sy / screenHeight) * (CYMAX - CYMIN) + CYMIN
        );
}

int cart2screenX(complex<double> c) {
    return (int)((c.real() - CXMIN) / (CXMAX - CXMIN) * screenWidth);
}

int cart2screenY(complex<double> c) {
    return (int)((c.imag() - CYMIN) / (CYMAX - CYMIN) * screenHeight);
}
