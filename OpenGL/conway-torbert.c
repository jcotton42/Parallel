// 
// Torbert, 1.5.2015
// 
// Conway's Game of Life in Serial
// 
// Reads from a file, loads Gosper's glider gun.
// 
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <GL/freeglut.h>
#include <mpi.h>

#define WID 800
#define HEI 600
#define DIM 5

int w = WID, h = HEI;
int n = WID / DIM, m = HEI / DIM;
int dn = DIM, dm = DIM;

//int c[WID / DIM][HEI / DIM]; // 36 9
int** c;
int xo = 12, yo = 6;

int flag = 0;

int size;
int rank;
MPI_Status status;

void display(void);
void onestep(void);
void idle(void);
void mouse(int button, int state, int xscr, int yscr);
void keyfunc(unsigned char key, int xscr, int yscr);
void reshape(int wscr, int hscr);
void manager(int argc, char* argv[]);
void worker(void);
int getNeighbor(int col, int row);

enum {
    KILL,
    RESIZE,
    INFO,
    STEP,
    INIT
};

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0) manager(argc, argv);
    else worker();
    MPI_Finalize();
    return 0;
}

void manager(int argc, char* argv[]) {
    FILE* fin;
    int x, y;
    char ch;

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(w, h);
    glutInitWindowPosition(100, 50);
    glutCreateWindow("Conway's Game of Life");

    glClearColor(0.0, 0.0, 0.0, 0.0);
    glShadeModel(GL_SMOOTH);

    for(x = 0; x < n; x++)
        for(y = 0; y < m; y++)
            c[x][y] = 0;

    fin = fopen("glider_gun.txt", "r");
    for(y = 0; y < 9; y++) {
        for(x = 0; x < 36; x++) {
            fread(&ch, sizeof(char), 1, fin);
            if(ch != '.')
                c[xo + x][yo + y] = 1;
        }
        fread(&ch, sizeof(char), 1, fin); // newline
    }
    fclose(fin);

    glutDisplayFunc(display);
    glutIdleFunc(idle);
    glutMouseFunc(mouse);
    glutKeyboardFunc(keyfunc);
    glutReshapeFunc(reshape);

    glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_GLUTMAINLOOP_RETURNS);
    glutMainLoop();
}

void worker(void) {
    MPI_Recv(&w, 1, MPI_INT, 0, INIT, MPI_COMM_WORLD, &status);
    MPI_Recv(&h, 1, MPI_INT, 0, INIT, MPI_COMM_WORLD, &status);

    c = malloc(sizeof(int*) * w);
    for(int i = 0; i < w; i++) {
        c[i] = malloc(sizeof(int) * h);
        MPI_Recv(c[i], h, MPI_INT, 0, INIT, MPI_COMM_WORLD, &status);
    }

    while(true) {
        MPI_Recv(NULL, 0, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        switch(status.MPI_TAG) {
            case STEP:
                for(int col = 0; col < w; col++) {
                    for(int row = 0; row < h; row++) {
                        int numLiveNeighbors =
                            getNeighbor(col - 1, row) +
                            getNeighbor(col + 1, row) +
                            getNeighbor(col, row - 1) +
                            getNeighbor(col, row + 1);
                        if(numLiveNeighbors < 2 || numLiveNeighbors > 3) c[col][row] = 0;
                        else if(numLiveNeighbors == 3) c[col][row] = 1;
                    }
                }
                MPI_Send(NULL, 0, MPI_INT, 0, STEP, MPI_COMM_WORLD);
                for(int col = 0; col < w; col++)
                    MPI_Send(c[col], h, MPI_INT, 0, STEP, MPI_COMM_WORLD);
                break;
            case KILL:
                return;
            case RESIZE:
                MPI_Recv(&w, 1, MPI_INT, status.MPI_SOURCE, RESIZE, MPI_COMM_WORLD, &status);
                MPI_Recv(&h, 1, MPI_INT, status.MPI_SOURCE, RESIZE, MPI_COMM_WORLD, &status);
                // resize and remap buffer here
                break;
            case INFO:
                int pts[2];
                MPI_Recv(pts, 2, MPI_INT, status.MPI_SOURCE, INFO, MPI_COMM_WORLD, &status);
                int col = pts[0], row = pts[1];
                if(col < 0) col += w;
                MPI_Send(c[col][row], 1, MPI_INT, status.MPI_SOURCE, INFO, MPI_COMM_WORLD);
                break;
            default:
                fprintf(stderr, "WTF happened?");
                break;
        }
    }
}

int getNeighbor(int col, int row) {
    if(row < 0) row += h;
    else if(row >= h) row -= h;
    if(col >= 0 && col < w) return c[col][row];
    int dest;
    if(col < 0) { 
        dest = rank == 1 ? size - 1 : rank - 1;
    } else if(col >= w) {
        col -= w;
        dest = rank == size - 1 ? 1 : rank + 1;
    }
    int pts[2] = {col, row};
    int val;
    MPI_Send(NULL, 0, MPI_INT, dest, INFO, MPI_COMM_WORLD);
    MPI_Send(pts, 2, MPI_INT, dest, INFO, MPI_COMM_WORLD);
    MPI_Recv(&val, 1, MPI_INT, dest, INFO, MPI_COMM_WORLD, &status);
    return val;
}

void display(void) {
    int x, y;
    int j, k;
    glClear(GL_COLOR_BUFFER_BIT);

    glColor3f(1.0, 1.0, 1.0);
    for(x = 0; x < n; x++) {
        for(y = 0; y < m; y++) {
            if(c[x][y]) {
                glBegin(GL_POINTS);
                for(j = 0; j <= dn; j++)
                    for(k = 0; k <= dm; k++)
                        glVertex2f(x*dn + j, y*dm + k);
                glEnd();
            }
        }
    }
    //
    glutSwapBuffers();
}

void onestep() {
    div_t info = div(w, size - 1);
    int overflow = info.rem;
    int workerWidth = info.quot;
    for(int i = 1; i < size; i++)
        MPI_Send(NULL, 0, MPI_INT, i, STEP, MPI_COMM_WORLD);
    for(int i = 0; i < size; i++) {
        MPI_Recv(NULL, 0, MPI_INT, MPI_ANY_SOURCE, STEP, MPI_COMM_WORLD, &status);
        int col = workerWidth * (status.MPI_SOURCE - 1);
        int end = col + workerWidth;
        // last worker handles anything that doesn't divide cleanly
        if(status.MPI_SOURCE == size - 1) end += overflow;
        while(col < end) {
            MPI_Recv(c[col], h, MPI_INT, status.MPI_SOURCE, STEP, MPI_COMM_WORLD, &status);
            col++;
        }
    }
}

void idle(void) {
    if(flag) {
        onestep();

        glutPostRedisplay();
    }
}

void mouse(int button, int state, int xscr, int yscr) {
    if(button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
        flag = (flag + 1) % 2;
}

void keyfunc(unsigned char key, int xscr, int yscr) {
    if(key == ' ') {
        onestep();

        glutPostRedisplay();
    } else if(key == 'q') {
        exit(0);
    }
}

void reshape(int wscr, int hscr) {
    w = wscr; h = hscr;
    glViewport(0, 0, (GLsizei)w, (GLsizei)h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    gluOrtho2D(0, w - 1, 0, h - 1);
    glMatrixMode(GL_MODELVIEW);
    // infrom workers
}
