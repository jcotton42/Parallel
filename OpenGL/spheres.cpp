#pragma warning(disable:4756)
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <functional>
#include <vector>
#include <array>

const int W = 640;
const int H = 480;

class Vector {
public:
    double x, y, z;
    Vector() : x(0), y(0), z(0) {}

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

    double dot(const Vector& other) const {
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
    Ray(Vector origin, Vector direction)
        : origin(origin), direction(direction) {}
};

class Light {
public:
    const Vector origin;
    const Vector direction;
    std::function<std::array<int ,3>()> color;

    Light(const Vector& origin, const Vector& direction, const std::function<int*()>& color)
        : origin(origin), direction(direction), color(color) {}
};

class Object {
protected:
    std::function<std::array<int, 3>(const Vector&, const Object&)> _color;
    std::function<bool(const Ray&, const Object&, double* const)> _rayIntersect;
public:
    Object(std::function<std::array<int, 3>(const Vector&, const Object&)> color,
        std::function<bool(const Ray&, const Object&, double* const)> rayIntersect)
        : _color(color), _rayIntersect(rayIntersect) {}

    std::array<int, 3> color(const Vector& position) {
        return this->_color(position, *this);
    }

    bool rayIntersect(const Ray& ray, double* const t) {
        return this->_rayIntersect(ray, *this, t);
    }
};

class Sphere : public Object {
public:
    const Vector origin;
    const double radius;

    Sphere(std::function<std::array<int, 3>(const Vector&, const Object&)> color,
        std::function<bool(const Ray&, const Object&, double* const)> rayIntersect,
        const Vector& origin, const double& radius)
        : Object(color, rayIntersect), origin(origin), radius(radius) {}
};

class Plane : public Object {
public:
    const Vector normal;
    const Vector point;

    Plane(std::function<std::array<int, 3>(const Vector&, const Object&)> color,
        std::function<bool(const Ray&, const Object&, double* const)> rayIntersect,
        const Vector& normal, const Vector& point)
        : Object(color, rayIntersect), normal(normal), point(point) {}
};

bool rayIntersectSphere(const Ray& ray, const Sphere& sphere, double* const t);
bool rayIntersectPlane(const Ray& ray, const Plane& plane, double* const t);

static std::vector<Object*> objects;

void setPixel(int x, int y, Ray& eye, Light& light, int*** pixels);
void writeppm(int*** data);

void initObjects();

int main(void) {
    initObjects();
    int*** data = new int**[H];
    Ray eye(
        Vector(0.500000, 0.500000, -1.000000),
        Vector()
        );

    Light light(
        Vector(0.000000, 1.000000, -0.500000),
        Vector(0.000000, 0.000000, 0.000000),
        []()->std::array<int, 3> { 
            return {255, 255, 255};
        }
    );
    for(int yp = 0; yp < H; yp++) {
        double y = 1.0 * yp / H;
        data[H - yp - 1] = new int*[W];
        for(int xp = 0; xp < W; xp++) {
            data[H - yp - 1][xp] = new int[3]();
            double x = 1.0 * xp / W;
            eye.direction.x = x;
            eye.direction.y = y;
            eye.direction -= eye.origin;
            eye.direction = eye.direction.normal();

            setPixel(xp, yp, eye, light, data);
        }
    }

    writeppm(data);
    system("imdisplay.exe out.ppm");
    return 0;
}

void initObjects() {
    objects.push_back(new Sphere(
        []()->std::array<int, 3> { return{0, 0, 255}; },
        rayIntersectSphere,
        Vector(0.500000, 0.500000, 0.166667),
        0.166667
        ));
    objects.push_back(new Sphere(
        []()->std::array<int, 3> { return{0, 255, 0}; },
        rayIntersectSphere,
        Vector(0.833333, 0.500000, 0.500000),
        0.166667
        ));
    objects.push_back(new Sphere(
        []()->std::array<int, 3> { return{255, 0, 0}; },
        rayIntersectSphere,
        Vector(0.333333, 0.666667, 0.666667),
        0.333333
        ));
    objects.push_back(new Plane(
        []()->std::array<int, 3> { return{255, 255, 255}; },
        rayIntersectPlane,
        Vector(0, 0.333333, 0),
        Vector(0, 0.333333, 0)
        ));
}

void setPixel(int x, int y, Ray& eye, Light& light, int*** pixels) {
    double t = INFINITY;
    double tmp;
    int ti = 0;
    Ray ray{
        Vector(),
        Vector()
    };
    int* color = nullptr;

    for(int i = 0; i < objects.size(); i++) {
        if(objects[i]->rayIntersect(eye, &tmp)) {
            if(tmp < t && tmp > 0) {
                t = tmp;
                ray.origin = eye.direction.scale(t) + eye.origin;
                color = objects[i]->color(ray.origin).data();
                ti = i;
            }
        }
    }

    if(!isinf(t)) {
        ray.direction = (light.origin - ray.origin).normal();

        t = INFINITY;
        for(int i = 0; i < objects.size(); i++)
            if(i != ti && objects[i]->rayIntersect(ray, &tmp))
                if(tmp < t && tmp > 0)
                    t = tmp;
        if(!isinf(t)) {
            ray.origin = eye.direction.scale(t) + eye.origin;
        }
    }

    if(!color) {
        pixels[H - y - 1][x][0] = 0;
        pixels[H - y - 1][x][1] = 0;
        pixels[H - y - 1][x][2] = 0;
    } else {
        pixels[H - y - 1][x][0] = color[0];
        pixels[H - y - 1][x][1] = color[1];
        pixels[H - y - 1][x][2] = color[2];
    }
}

bool rayIntersectSphere(const Ray& ray, const Sphere& sphere, double* const t) {
    // d.dt^2 + 2d.(p0 - c)t + (p0 - c).(p0 - c) - r^2 = 0
    double A, B, C;
    A = ray.direction.dot(ray.direction);

    Vector diff; // p0 - c
    diff = ray.origin - sphere.origin;

    B = 2 * ray.direction.dot(diff);

    C = diff.dot(diff) - pow(sphere.radius, 2);

    double discr; // discriminant (sp?)
    discr = pow(B, 2) - 4 * A * C;
    if(discr < 0) {
        *t = 0;
        return false;
    }
    // (-b +- sqrt(discr))/(4ac)

    double minus = (-B - sqrt(discr)) / (2 * A);
    double plus = (-B + sqrt(discr)) / (2 * A);

    if(minus < 0) minus = INFINITY;
    if(plus < 0) plus = INFINITY;
    *t = fmin(minus, plus);
    return true;
}

bool rayIntersectPlane(const Ray& ray, const Plane& plane, double* const t) {
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
    
    double d = -plane.normal.dot(plane.point);
    double numerator = -ray.origin.dot(plane.normal) - d;
    double denominator = ray.direction.dot(plane.normal);

    if(denominator == 0.0) {
        if(numerator == 0.0) {
            *t = 0.0;
            return true;
        }
        *t = 0.0;
        return false;
    }
    *t = numerator / denominator;
    return true;
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
