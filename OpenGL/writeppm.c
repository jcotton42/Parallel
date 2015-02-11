//
// Torbert, 3 February 2014
//
#include <stdio.h>
#include <stdlib.h>

#define WIDTH 640
#define HIGHT 480

int main(void) {
    int*** rgb = malloc(sizeof(int**) * HIGHT); // red-green-blue for each pixel

    int y, x;

    FILE* ppm;

    for(y = 0; y < HIGHT; y++) {
        rgb[y] = malloc(sizeof(int*) * WIDTH);
        for(x = 0; x < WIDTH; x++) {
            rgb[y][x] = malloc(sizeof(int) * 3);
            rgb[y][x][0] = 0; // red
            rgb[y][x][1] = 255; // green
            rgb[y][x][2] = 0; // blue
        }
    }

    ppm = fopen("out.ppm", "w");

    fprintf(ppm, "P3\n");
    fprintf(ppm, "%d %d\n", WIDTH, HIGHT);
    fprintf(ppm, "255\n");

    for(y = 0; y < HIGHT; y++) {
        for(x = 0; x < WIDTH; x++) {
            fprintf(
                ppm,
                "%d %d %d\n",
                rgb[y][x][0],
                rgb[y][x][1],
                rgb[y][x][2]);
        }
    }
    fclose(ppm);

    return 0;
}
