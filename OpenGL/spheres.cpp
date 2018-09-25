#pragma warning(disable:4756)
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <functional>
#include <vector>
#include <array>

class Object;
const int W = 640;
const int H = 480;

class Vector {
public:
    double x, y, z;

    Vector() : x(0), y(0), z(0) {}

    Vector(double x, double y, double z) : x(x), y(y), z(z) {}

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

    Vector cross(const Vector& other) const {
        return Vector(
            this->y * other.z - this->z * other.y,
            this->z * other.x - this->x * other.z,
            this->x * other.y - this->y * other.x
            );
    }

    double distance(const Vector& other) const {
        return (other - *this).magnitude();
    }

    double magnitude() const {
        return sqrt(
            pow(this->x, 2) +
            pow(this->y, 2) +
            pow(this->z, 2)
            );
    }

    Vector scale(double d) const {
        return Vector(this->x * d, this->y * d, this->z * d);
    }

    Vector normal() const {
        double len = this->magnitude();
        return Vector(this->x / len, this->y / len, this->z / len);
    }

    Vector rotate(Vector axis, double theta) const {
        //http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/
        axis = axis.normal();
        double u = axis.x;
        double v = axis.y;
        double w = axis.z;
        double ct = cos(theta);
        double st = sin(theta);
        double piece = u * x + v * y + w * z;

        return Vector(
            u * piece * (1 - ct) + x * ct + (-w * y + v * z) * st,
            v * piece * (1 - ct) + y * ct + (w * x - u * z) * st,
            w * piece * (1 - ct) + z * ct + (-v * x + u * y) * st
            );
    }

    double angle(const Vector& other) const {
        return acos(
            this->dot(other) /
            (this->magnitude() * other.magnitude())
            );
    }
};

typedef std::function<std::array<int, 3>()> lightColorFunc;

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
    lightColorFunc color;

    Light(const Vector& origin, const Vector& direction, const lightColorFunc color)
        : origin(origin), direction(direction), color(color) {}
};

typedef std::function<std::array<int, 3>(const Vector&, const Object&)> objectColorFunc;

class Object {
protected:
    objectColorFunc _color;
    explicit Object() {}
public:
    explicit Object(objectColorFunc color)
        : _color(color) {}

    virtual std::array<int, 3> color(const Vector& position) const {
        return this->_color(position, *this);
    }

    virtual bool rayIntersect(const Ray& ray, double* const t) const = 0;
};

class Sphere : public Object {
public:
    const Vector origin;
    const double radius;

    virtual bool rayIntersect(const Ray& ray, double* const t) const override {
        // d.dt^2 + 2d.(p0 - c)t + (p0 - c).(p0 - c) - r^2 = 0
        double A, B, C;
        A = ray.direction.dot(ray.direction);

        Vector diff; // p0 - c
        diff = ray.origin - this->origin;

        B = 2 * ray.direction.dot(diff);

        C = diff.dot(diff) - pow(this->radius, 2);

        double discr; // discriminant (sp?)
        discr = pow(B, 2) - 4.0 * A * C;
        if(discr < 0) {
            *t = 0;
            return false;
        }
        // (-b +- sqrt(discr))/(2a)

        double minus = (-B - sqrt(discr)) / (2.0 * A);
        double plus = (-B + sqrt(discr)) / (2.0 * A);

        if(minus < 0) minus = INFINITY;
        if(plus < 0) plus = INFINITY;
        *t = fmin(minus, plus);
        return true;
    }

    Sphere(objectColorFunc color,
        const Vector& origin, const double& radius)
        : Object(color), origin(origin), radius(radius) {}
};

class Plane : public Object {
public:
    const Vector normal;
    const Vector point;
    
    virtual bool rayIntersect(const Ray& ray, double* const t) const override {
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

        double d = -this->normal.dot(this->point);
        double numerator = -ray.origin.dot(this->normal) - d;
        double denominator = ray.direction.dot(this->normal);

        if(denominator == 0.0) {
            *t = 0.0;
            return numerator == 0.0;
        }
        *t = numerator / denominator;
        return true;
    }

    Plane(objectColorFunc color,
        const Vector& normal, const Vector& point)
        : Object(color), normal(normal), point(point) {}
};

class CompoundObject : Object {
    std::vector<Sphere> spheres;
    std::vector<Plane> sides;
public:
    explicit CompoundObject(const std::vector<Sphere>& spheres)
        : spheres(spheres), sides(std::vector<Plane>(6)) {
        double xmin = INFINITY, ymin = INFINITY, zmin = INFINITY;
        double xmax = -INFINITY, ymax = -INFINITY, zmax = -INFINITY;
        for(auto&& s : spheres) {
            xmin = std::min(xmin, s.origin.x - s.radius);
            ymin = std::min(ymin, s.origin.y - s.radius);
            zmin = std::min(zmin, s.origin.z - s.radius);

            xmax = std::max(xmax, s.origin.x + s.radius);
            ymax = std::max(ymax, s.origin.y + s.radius);
            zmax = std::max(zmax, s.origin.z + s.radius);
        }
        Vector p1(xmin, ymin, zmin), p2(xmax, ymax, zmax);
        sides.push_back(Plane(nullptr, ))
    }

    virtual bool rayIntersect(const Ray& ray, double* const t) const override {
        for(auto&& side : sides)
            if(side.rayIntersect(ray, t))
                goto inBox;
        return false;
    inBox:
        double tmp;
        for(auto&& s : spheres)
            if(s.rayIntersect(ray, &tmp))
                if(tmp < *t && tmp > 0)
                    *t = tmp;
        return true;
    }

    virtual std::array<int, 3> color(const Vector& position) const override {
        return std::array<int, 3> { { 0, 0, 255 }};
    }
};

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
        []()->std::array<int, 3> { return std::array<int, 3>{{255, 255, 255}}; }
    );
    for(int yp = 0; yp < H; yp++) {
        double y = 1.0 * yp / H;
        data[H - yp - 1] = new int*[W];
        for(int xp = 0; xp < W; xp++) {
            data[H - yp - 1][xp] = new int[3]();
            double x = 1.0 * xp / W;
            eye.direction.x = x;
            eye.direction.y = y;
            eye.direction.z = 0;
            eye.direction -= eye.origin;
            eye.direction = eye.direction.normal();

            setPixel(xp, yp, eye, light, data);
        }
        printf("set %d\n", yp);
    }

    writeppm(data);
    system("imdisplay.exe out.ppm");
    return 0;
}

/*void loadFile(const char* file) {
    double x, y, z, r;
    FILE* f = fopen(file, "r");
    auto func = [](const Vector& a, const Object& b)->std::array<int, 3> { return std::array<int, 3>{ { 0, 0, 255 }}; };
    std::vector<Object*> objs;
    while(fscanf(f, "%lf %lf %lf %lf\n", &x, &y, &z, &r) != EOF)
        objs.push_back(new Sphere(
        func,
        Vector(x, y, z),
        r
        ));
    fclose(f);
    objects.push_back(new CompoundObject(objs, ))
    printf("loaded: %llu\n", objects.size());
}*/

void initObjects() {
    objects.push_back(new Sphere(
        [](const Vector& a, const Object& b)->std::array<int, 3> { return std::array<int, 3>{{0, 0, 255}}; },
        Vector(0.500000, 0.500000, 0.166667),
        0.166667
        ));
    objects.push_back(new Sphere(
        [](const Vector& a, const Object& b)->std::array<int, 3> { return std::array<int, 3>{{0, 255, 0}}; },
        Vector(0.833333, 0.500000, 0.500000),
        0.166667
        ));
    objects.push_back(new Sphere(
        [](const Vector& a, const Object& b)->std::array<int, 3> { return std::array<int, 3>{{255, 0, 0}}; },
        Vector(0.333333, 0.666667, 0.666667),
        0.333333
        ));
    objects.push_back(new Plane(
        [](const Vector& a, const Object& b)->std::array<int, 3> { return std::array<int, 3>{{255, 255, 255}}; },
        Vector(0, 0.333333, 0),
        Vector(0, 0.333333, 0)
        ));
//    loadFile("helix.txt");
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
            color[0] /= 2;
            color[1] /= 2;
            color[2] /= 2;
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
