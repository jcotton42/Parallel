#pragma warning(disable:4756)
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <functional>

const int W = 640;
const int H = 480;

class Vector {
    double x, y, z;
public:
    Vector(double x, double y, double z) {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    Vector& operator+=(const Vector& rhs) {
        this->x += rhs.x;
        this->y += rhs.y;
        this->z += rhs.z;
        return *this;
    }

    Vector& operator-=(const Vector& rhs) {
        this->x -= rhs.x;
        this->y -= rhs.y;
        this->z -= rhs.z;
        return *this;
    }

    friend Vector operator+(Vector lhs, const Vector& rhs) {
        return lhs += rhs;
    }

    friend Vector operator-(Vector lhs, const Vector& rhs) {
        return lhs -= rhs;
    }

    double dot(const Vector& other) {
        return this->x * other.x + this->y * other.y + this->z * other.z;
    }

    Vector cross(const Vector& other) {
        return Vector(
            this->y * other.z - this->z * other.y,
            this->z * other.x - this->x * other.z,
            this->x * other.y - this->y * other.x
            );
    }

    double distance(const Vector& other) {
        return (other - *this).magnitude();
    }

    double magnitude() {
        return sqrt(
            pow(this->x, 2) +
            pow(this->y, 2) +
            pow(this->z, 2)
            );
    }

    Vector scale(double d) {
        return Vector(this->x * d, this->y * d, this->z * d);
    }

    Vector normal() {
        double len = this->magnitude();
        return Vector(this->x / len, this->y / len, this->z / len);
    }

    double angle(Vector& other) {
        return acos(
            this->dot(other) /
            (this->magnitude() * other.magnitude())
            );
    }
};

class Ray {
public:
    Vector origin;
    Vector direction;
    Ray(const Vector& origin, const Vector& direction)
        : origin(origin), direction(direction) {}
};

class Light {
public:
    Vector origin;
    Vector direction;
    std::function<int*()> color;

    Light(const Vector& origin, const Vector& direction, const std::function<int*()>& color)
        : origin(origin), direction(direction), color(color) {}
};

class Object {
public:
    std::function<int*(const Vector&, const Object&)> color;
    std::function<bool(const Ray&, const Object&, const Vector&)> rayIntersect;
    
    Object(std::function<int*(const Vector&, const Object&)> color,
        std::function<bool(const Ray&, const Object&, const Vector&)> rayIntersect)
        : color(color), rayIntersect(rayIntersect) {}
};

class Sphere : public Object {
public:
    const Vector origin;
    const double radius;

    Sphere(std::function<int*(const Vector&, const Object&)> color,
        std::function<bool(const Ray&, const Object&, const Vector&)> rayIntersect,
        const Vector& origin, const double& radius)
        : Object(color, rayIntersect), origin(origin), radius(radius) {}
};

class Plane : public Object {
public:
    const Vector normal;
    const Vector point;

    Plane(std::function<int*(const Vector&, const Object&)> color,
        std::function<bool(const Ray&, const Object&, const Vector&)> rayIntersect,
        const Vector& normal, const Vector& point)
        : Object(color, rayIntersect), normal(normal), point(point) {}
};

bool rayIntersectSphere(const Ray& ray, const Sphere& sphere, double* t);
bool rayIntersectPlane(const Ray& ray, const Plane& plane, double* t);

static Object* objects[] = {
    &(Sphere) {
        (int[]) {0, 0, 255},
        rayIntersectSphere,
        {0.500000, 0.500000, 0.166667, 0.000000},
        0.166667
    },
    &(Sphere) {
        (int[]) {0, 255, 0},
        rayIntersectSphere,
        {0.833333, 0.500000, 0.500000, 0.000000},
        0.166667
    },
    &(Sphere) {
        (int[]) {255, 0, 0},
        rayIntersectSphere,
        {0.333333, 0.666667, 0.666667, 0.000000},
        0.333333
    },
    &(Plane) {
        (int[]) {255, 255, 255},
        rayIntersectPlane,
        {0, 0.333333, 0, 0},
        {0, 0.333333, 0, 0}
    }
};

void setPixel(int x, int y, Ray* eye, Light* light, int*** pixels);
void writeppm(int*** data);

#define NUM_OBJECTS sizeof(objects) / sizeof(Object*)

int main(void) {
    int*** data = malloc(sizeof(int**) * H);
    Ray eye = {
        (double[]) {0.500000, 0.500000, -1.000000, 0.000000},
        malloc(sizeof(double) * 4)
    };
    Light light = {
        (double[]){0.000000, 1.000000, -0.500000, 0.000000},
        (double[]){0.000000, 0.000000, 0.000000, 0.000000},
        (int[]){255, 255, 255}
    };
    for(int yp = 0; yp < H; yp++) {
        double y = 1.0 * yp / H;
        data[H - yp - 1] = malloc(sizeof(int*) * W);
        for(int xp = 0; xp < W; xp++) {
            data[H - yp - 1][xp] = calloc(3, sizeof(int));
            double x = 1.0 * xp / W;
            eye.direction = (double[]) { x, y, 0, 0 };
            vectorSub(eye.direction, eye.origin, eye.direction);
            vectorNorm(eye.direction, eye.direction);

            setPixel(xp, yp, &eye, &light, data);
        }
    }

    writeppm(data);
    system("imdisplay.exe out.ppm");
    return 0;
}

void setPixel(int x, int y, Ray* eye, Light* light, int*** pixels) {
    double t = INFINITY;
    double tmp;
    int ti = 0;
    Ray ray = {
        (double[]){0,0,0,0},
        (double[]){0,0,0,0}
    };
    int* color = NULL;

    for(int i = 0; i < NUM_OBJECTS; i++) {
        if(objects[i]->rayIntersect(eye, objects[i], &tmp)) {
            if(tmp < t && tmp > 0) {
                t = tmp;
                if(color == NULL) color = alloca(sizeof(int) * 3);
                color[0] = objects[i]->rgb[0];
                color[1] = objects[i]->rgb[1];
                color[2] = objects[i]->rgb[2];
                ti = i;
            }
        }
    }

    if(light != NULL && !isinf(t)) {
        vectorScale(eye->direction, t, ray.origin);
        vectorAdd(ray.origin, eye->origin, ray.origin);
        vectorSub(light->origin, ray.origin, ray.direction);
        vectorNorm(ray.direction, ray.direction);
        t = INFINITY;
        for(int i = 0; i < NUM_OBJECTS; i++)
            if(i != ti && objects[i]->rayIntersect(&ray, objects[i], &tmp))
                if(tmp < t && tmp > 0)
                    t = tmp;
        if(!isinf(t)) {
            vectorScale(eye->direction, t, ray.origin);
            vectorAdd(ray.origin, eye->origin, ray.origin);
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

bool rayIntersectSphere(const Ray& ray, const Sphere& sphere, double* t) {
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
        double minus = (-B - sqrt(discr)) / (2 * A);
        double plus = (-B + sqrt(discr)) / (2 * A);

        if(minus < 0) minus = INFINITY;
        if(plus < 0) plus = INFINITY;
        *t = fmin(minus, plus);
        return true;
    }
}

bool rayIntersectPlane(const Ray& ray, const Plane& plane, double* t) {
    // http://www.cs.princeton.edu/courses/archive/fall00/cs426/lectures/raycast/sld017.htm
    // t = -(P0.N + d) / (V.N)
    // P0: ray origin
    // V: ray direction
    // N: normal to plane
    // t: distance from ray origin to plane
    // if bottom is zero
    //      if top is zero
    //          intersects everywhere
    //      else
    //          intersects nowhere
    // else
    //      intersects at one point

    double d = -vectorDot(plane->normal, plane->point);
    double numerator = -(vectorDot(ray->origin, plane->normal) + d);
    double denominator = vectorDot(ray->direction, plane->normal);
    if(denominator == 0.0) {
        if(numerator == 0.0) {
            *t = 0.0;
            return true;
        } else {
            *t = 0.0;
            return false;
        }
    } else {
        *t = numerator / denominator;
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

void vectorScale(double a[4], double b, double out[4]) {
    out[0] = a[0] * b;
    out[1] = a[1] * b;
    out[2] = a[2] * b;
    out[3] = a[3] * b;
}

double vectorAngle(double a[4], double b[4]) {
    // acos((a.b)/(|a||b|))
    return acos(
        vectorDot(a, b) /
        sqrt(vectorDot(a, a) * vectorDot(b, b))
        );
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
