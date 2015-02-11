#include <stdio.h>
#include <stdlib.h>
#include <GL/freeglut.h>
#include <mpi.h>

enum {
    KILL,
    DRAW,
    NEW_DATA,
    NEED_DATA
};

int size;
int rank;
int width; // varies by worker
int height = 600;
int** cells;

void manager(int* argc, char** argv);
void worker(void);

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0) manager(&argc, argv);
    else worker();
    MPI_Finalize();
}

void manager(int* argc, char** argv) {
    FILE* f;
    width = 800;
    cells = calloc(sizeof(int*), width);
    for(int i = 0; i < width; i++)
        cells[i] = calloc(sizeof(int), height);

    glutInit(argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(width, height);
    glutInitWindowPosition(100, 50);
    glutCreateWindow("");
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glShadeModel(GL_SMOOTH);
    glViewport(0, 0, width, height);
    glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_GLUTMAINLOOP_RETURNS);

    glutMainLoop();

    for(int i = 1; i < size; i++)
        MPI_Send(NULL, 0, MPI_INT, i, KILL, MPI_COMM_WORLD);
}
