#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <immintrin.h>
#include <stdbool.h>

#define W 640
#define H 480

typedef struct {
    double* origin;
    double* direction;
} Ray;

typedef struct {
    int* rgb;
    bool(*rayIntersect)(Ray*, struct Object*, double*);
} Object;

typedef struct {
    int* rgb;
    bool(*rayIntersect)(Ray*, struct Sphere*, double*);
    // it's 4 for a good reason, don't ask :)
    double origin[4];
    double radius;
} Sphere;

typedef struct {
    int* rgb;
    bool(*rayIntersect)(Ray*, struct Plane*, double*);
    double normal[4];
    double point[4];
} Plane;

double light[4] = {0.500000, 0.500000, -1.000000, 0.000000};
double eye[4] = {0.500000, 0.500000, -1.000000, 0.000000};

void setPixel(int x, int y, Ray* ray, int*** pixels);
bool rayIntersectSphere(Ray* ray, Sphere* sphere, double* t);
bool rayIntersectPlane(Ray* ray, Plane* plane, double* t);
void writeppm(int*** data);
void vectorNorm(double a[4], double out[4]);
double vectorDot(double a[4], double b[4]);
void vectorSub(double a[4], double b[4], double out[4]);
void vectorAdd(double a[4], double b[4], double out[4]);
double distance(double a[4], double b[4]);

#define NUM_OBJECTS 4

Object* objects[NUM_OBJECTS] = {
    &(Sphere){
        (int[]){0, 0, 255},
        rayIntersectSphere,
        {0.500000, 0.500000, 0.166667, 0.000000},
        0.166667
    },
    &(Sphere){
        (int[]){0, 255, 0},
        rayIntersectSphere,
        {0.833333, 0.500000, 0.500000, 0.000000},
        0.166667
    },
    &(Sphere){
        (int[]){255, 0, 0},
        rayIntersectSphere,
        {0.333333, 0.666667, 0.666667, 0.000000},
        0.333333
    },
    &(Plane){
        (int[]){255, 255, 255},
        rayIntersectPlane,
        {0, 1, 0, 0},
        {1, 0.333333, 1, 0}
        }
};

int main(void) {
    int*** data = malloc(sizeof(int**) * H);
    Ray ray = {eye, _aligned_malloc(sizeof(double) * 4, 32)};
    for(int yp = 0; yp < H; yp++) {
        double y = 1.0 * yp / H;
        data[H - yp - 1] = malloc(sizeof(int*) * W);
        for(int xp = 0; xp < W; xp++) {
            data[H - yp - 1][xp] = calloc(3, sizeof(int));
            double x = 1.0 * xp / W;
            ray.direction = (double[]) { x, y, 0, 0 };
            vectorSub(ray.direction, ray.origin, ray.direction);
            vectorNorm(ray.direction, ray.direction);

            setPixel(xp, yp, &ray, data);
        }
    }
    writeppm(data);
    return 0;
}

void setPixel(int x, int y, Ray* ray, int*** pixels) {
    double t = INFINITY;
    double tmp;
    int* color = NULL;

    for(int i = 0; i < NUM_OBJECTS; i++) {
        if(objects[i]->rayIntersect(ray, objects[i], &tmp)) {
            if(tmp < t) {
                t = tmp;
                color = objects[i]->rgb;
            }
        }
    }

    if(color == NULL) {
        pixels[H - y - 1][x][0] = 0;
        pixels[H - y - 1][x][1] = 0;
        pixels[H - y - 1][x][2] = 0;
    } else {
        pixels[H - y - 1][x][0] = color[0];
        pixels[H - y - 1][x][1] = color[1];
        pixels[H - y - 1][x][2] = color[2];
    }
}

void vectorNorm(double a[4], double out[4]) {
    double length = sqrt(vectorDot(a, a));
    out[0] = a[0] / length;
    out[1] = a[1] / length;
    out[2] = a[2] / length;
}

bool rayIntersectSphere(Ray* ray, Sphere* sphere, double* t) {
    // d.dt^2 + 2d.(p0 - c)t + (p0 - c).(p0 - c) - r^2 = 0
    double A, B, C;
    A = vectorDot(ray->direction, ray->direction);

    double diff[4]; // p0 - c
    vectorSub(ray->origin, sphere->origin, diff);
    B = 2 * vectorDot(ray->direction, diff);

    C = vectorDot(diff, diff) - pow(sphere->radius, 2);

    double discr; // discriminant (sp?)
    discr = pow(B, 2) - 4 * A * C;
    if(discr < 0) {
        *t = 0;
        return false;
    } else {
        // (-b +- sqrt(discr))/(4ac)
        double minus = -B - sqrt(discr) / (4 * A * C);
        double plus = -B + sqrt(discr) / (4 * A * C);
        if(minus < 0) minus = INFINITY;
        if(plus < 0) plus = INFINITY;
        *t = fmin(minus, plus);
        return true;
    }
}

bool rayIntersectPlane(Ray* ray, Plane* plane, double* t) {
    // d = ((p0 - l0).n)/(l.n)
    // l.n == 0: line and plane are parallel
    //      (p0 - l0).n == 0: line is inside the plane
    //      (p0 - l0).n != 0: no intersection
    // l.n != 0: single point of intersection
    //      t = distance(ray->origin, dl + l0)
    // n: normal vector
    // p0: point on plane
    // l: direction of ray
    // l0: point on ray
    double numerator;
    double denominator;
    double tmp[4];
    denominator = vectorDot(ray->direction, plane->normal);
    vectorSub(plane->point, ray->origin, tmp);
    numerator = vectorDot(tmp, plane->normal);
    if(denominator == 0) {
        *t = INFINITY;
        return numerator == 0;
    } else {
        double d = numerator / denominator;
        tmp[0] = ray->direction[0] * d;
        tmp[1] = ray->direction[1] * d;
        tmp[2] = ray->direction[2] * d;
        tmp[3] = 0;
        vectorAdd(tmp, ray->origin, tmp);
        *t = distance(ray->origin, tmp);
        return true;
    }
}

double vectorDot(double a[4], double b[4]) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

void vectorSub(double a[4], double b[4], double out[4]) {
    out[0] = a[0] - b[0];
    out[1] = a[1] - b[1];
    out[2] = a[2] - b[2];
}

void vectorAdd(double a[4], double b[4], double out[4]) {
    out[0] = a[0] + b[0];
    out[1] = a[1] + b[1];
    out[2] = a[2] + b[2];
}

double distance(double a[4], double b[4]) {
    return sqrt(
        pow(b[0] - a[0], 2) +
        pow(b[1] - a[1], 2) +
        pow(b[2] - a[2], 2)
        );
}

void writeppm(int*** data) {
    FILE* ppm = fopen("out.ppm", "w");

    fprintf(ppm, "P3\n");
    fprintf(ppm, "%d %d\n", W, H);
    fprintf(ppm, "255\n");

    for(int y = 0; y < H; y++) {
        for(int x = 0; x < W; x++) {
            fprintf(
                ppm,
                "%d %d %d\n",
                data[y][x][0],
                data[y][x][1],
                data[y][x][2]
                );
        }
    }

    fclose(ppm);
}
