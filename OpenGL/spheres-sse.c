#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <immintrin.h>
#include <stdbool.h>

#define W 640
#define H 480

typedef struct {
    // it's 4 for a good reason, don't ask :)
    __m256d origin;
    double radius;
} Sphere;

typedef struct {
    __m256d origin;
    __m256d direction;
} Ray;

Sphere sphereB = {
    {0.500000, 0.500000, 0.166667, 0.000000},
    0.166667
};
Sphere sphereG = {
    {0.833333, 0.500000, 0.500000, 0.000000},
    0.166667
};
Sphere sphereR = {
    {0.333333, 0.666667, 0.666667, 0.000000},
    0.333333
};
double light[4] = {0.500000, 0.500000, -1.000000, 0.000000};
__declspec(align(32)) double eye[4] = {0.500000, 0.500000, -1.000000, 0.000000};

bool rayIntersectsSphere(Ray* ray, Sphere* sphere);
void writeppm(int*** data);
void vectorNorm(double a[4], double out[4]);
double vectorDot(double a[4], double b[4]);

int main(void) {
    int*** data = malloc(sizeof(int**) * H);
    Ray ray;
    ray.origin = _mm256_load_pd(eye);
    for(int yp = 0; yp < H; yp++) {
        double y = 1.0 * yp / H;
        data[H - yp - 1] = malloc(sizeof(int*) * W);
        for(int xp = 0; xp < W; xp++) {
            data[H - yp - 1][xp] = calloc(3, sizeof(int));
            double x = 1.0 * xp / W;
            __declspec(align(32)) double tmp[] = {x, y, 0, 0};
            ray.direction = _mm256_load_pd(tmp);
            ray.direction = _mm256_sub_pd(ray.direction, ray.origin);
            vectorNorm(ray.direction, &ray.direction);
            if(
                rayIntersectsSphere(&ray, &sphereB) ||
                rayIntersectsSphere(&ray, &sphereG) ||
                rayIntersectsSphere(&ray, &sphereR)
                ) {
                data[H - yp - 1][xp][0] = data[H - yp - 1][xp][1] = data[H - yp - 1][xp][2] = 255;
            }
        }
    }
    writeppm(data);
    return 0;
}

void vectorNorm(__m256d a, __m256d* out) {
    double length = sqrt(vectorDot(a, a));
    out[0] = a[0] / length;
    out[1] = a[1] / length;
    out[2] = a[2] / length;
}

bool rayIntersectsSphere(Ray* ray, Sphere* sphere) {
    // d.dt^2 + 2d.(p0 - c)t + (p0 - c).(p0 - c) - r^2 = 0
    double A, B, C;
    A = vectorDot(ray->direction, ray->direction);

    double diff[4]; // p0 - c
    vectorSub(ray->origin, sphere->origin, diff);
    B = 2 * vectorDot(ray->direction, diff);

    C = vectorDot(diff, diff) - pow(sphere->radius, 2);

    double discr; // discriminant (sp?)
    discr = pow(B, 2) - 4 * A * C;
    return discr >= 0;
}

double vectorDot(__m256d a, __m256d b) {
    __m256d c = _mm256_mul_pd(a, b);
    c = _mm256_hadd_pd(_mm256_hadd_pd(c, c));
    c = _mm256_shuffle_pd()
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
