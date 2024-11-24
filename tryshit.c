#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define BOX_SIZE 1003.0
#define NUM_OBJECTS 10
#define TRACE_DEPTH 5

static enum Shape { SPHERE, CUBE };

float random() {
    return (rand() / (RAND_MAX + 1.0));
}

typedef struct {
    float x, y, z;
} Vector3;

// Sum of vectors
Vector3 add(Vector3* vec1, Vector3* vec2) {
    Vector3 result = {vec1->x + vec2->x, vec1->y + vec2->y, vec1->z + vec2->z};
    return result;
}

// Hadamand product
Vector3 hadprod(Vector3* vec1, Vector3* vec2) {
    Vector3 result = {vec1->x * vec2->x, vec1->y * vec2->y, vec1->z * vec2->z};
    return result;
}

// Dot product
float dotprod(Vector3* vec1, Vector3* vec2) {
    return vec1->x * vec2->x + vec1->y * vec2->y + vec1->z * vec2->z;
}

// Make vector longer by f
Vector3 scale(Vector3* vec, float f) {
    Vector3 result = {vec->x * f, vec->y * f, vec->z * f};
    return result;
}

// Normalize the vector
Vector3 normalize(Vector3* vec) {
    Vector3 result = scale(vec, 1/sqrt(vec->x*vec->x + vec->y*vec->y + vec->z*vec->z));
    return result;
}

Vector3 crossprod(Vector3* vec1, Vector3* vec2) {
    Vector3 result = { vec1->y * vec2->z - vec1->z * vec2->y, vec1->z * vec2->x - vec1->x * vec2->z, vec1->x * vec2->y - vec1->y * vec2->x};
    return result;
}

typedef struct {
    // vector of origin
    Vector3 origin;
    // normalized vector of direction
    Vector3 direction;
} Ray;

// Get a point along the ray over distance of t
Vector3 point_at(Ray* r, float t) {
    Vector3 s = scale(&r->direction, t);
    return add(&s, &r->origin);
}



typedef struct {
    Vector3 center;
    float radius;
} Sphere;

int collide_sphere(Sphere* s, Ray* r, float* t) {
    Vector3 ray_reverse = scale(&r->origin, -1);
    Vector3 from_center = add(&ray_reverse, &s->center);
    float b = dotprod(&r->direction, &from_center);
    float c = dotprod(&from_center, &from_center) - s->radius;
    float D = b * b - c;
    if (D < 0) return 0;
    float solve = b - sqrt(D);
    if (solve > 0 && solve < *t) {
        *t = solve;
        return 1;
    }
    return 0;
}

typedef struct {
    Vector3 min_corner;
    Vector3 max_corner;
} Cube;

// Collision detection for AABB (Axis-Aligned Bounding Box)
int collide_cube(Cube* c, Ray* q, float* t) {
    float t_min = 0, t_max = *t;

    for (int i = 0; i < 3; ++i) {
        float invD = 1.0f / (i == 0 ? q->direction.x : (i == 1 ? q->direction.y : q->direction.z));
        float t0 = ((i == 0 ? c->min_corner.x : (i == 1 ? c->min_corner.y : c->min_corner.z)) - (i == 0 ? q->origin.x : (i == 1 ? q->origin.y : q->origin.z))) * invD;
        float t1 = ((i == 0 ? c->max_corner.x : (i == 1 ? c->max_corner.y : c->max_corner.z)) - (i == 0 ? q->origin.x : (i == 1 ? q->origin.y : q->origin.z))) * invD;

        if (invD < 0) {
            // Swap t0 and t1 when the ray direction is negative
            float temp = t0;
            t0 = t1;
            t1 = temp;
        }

        // Update t_min and t_max
        t_min = t0 > t_min ? t0 : t_min;
        t_max = t1 < t_max ? t1 : t_max;

        // If no intersection on this axis, return 0
        if (t_max <= t_min) return 0;
    }

    // Update t with the closest intersection distance
    *t = t_min;
    return 1; // Intersection detected
}

int to_color(float c) {
    return pow(c < 0 ? 0 : c > 1 ? 1 : c, .45) * 255 + .5;
}

Ray create_random(float x, float y) {
    Vector3 e = { x, -y, 1 };
    e = normalize(&e);
    e = scale(&e, 4);
    float a = 6 * random();
    float c = .2 * sqrt(random());
    float b = sin(a) * c;
    a = cos(a) * c;
    e.x -= a;
    e.y -= b;
    Vector3 position = { a, b, 0 };
    e = normalize(&e);
    Ray result = { position, e };
    return result;
}

Sphere box[NUM_OBJECTS] = {
        {
            {0, -2, 5}, 1
        },
        {
            {0, -BOX_SIZE, 0}, 1e6
        },
        {
            {0, BOX_SIZE, 0}, 1e6
        },
        {
            {BOX_SIZE, 0, 0}, 1e6
        },
        {
            {-BOX_SIZE, 0, 0}, 1e6
        },
        {
            {0, 0, -BOX_SIZE}, 1e6
        },
        {
            {0, 0, BOX_SIZE + 3}, 1e6
        },
        {
            {-2, -2, 4}, 2
        },
        {
            {2, -3, 4}, 1
        },
        {
            {2, -1, 4}, 1
        }
};

Cube cubes[1] = {
    {
        {0, 0, 0}
    }
};

int find_intersection(Ray* a, float* t) {
    int n = -1;
    for (int m = 0; m < NUM_OBJECTS; m++) {
        if (collide_sphere(&box[m], a, t)) n = m;
    }
    return n;
}

Vector3 trace_ray(Ray* a, int b) {
    float t = 1e6;
    int n = find_intersection(a, &t);
    if (b > TRACE_DEPTH || n < 0) {
        Vector3 result = {0, 0, 0};
        return result;
    }
     if (!n) {
        Vector3 result = {.9, .5, .1};
        return result;
     }

    Vector3 P = point_at(a, t);

    Vector3 sphere_ray_normal = box[n].center;
    sphere_ray_normal = scale(&sphere_ray_normal, -1);
    sphere_ray_normal = add(&sphere_ray_normal, &P);
    sphere_ray_normal = normalize(&sphere_ray_normal);

    if (n > 6) {
        Vector3 inverse_direction = scale(&sphere_ray_normal, -2 * dotprod(&sphere_ray_normal, &a->direction));
        Vector3 reflection = add(&a->direction, &inverse_direction);
        Ray r = { P, reflection };
        Vector3 next_vec = trace_ray(&r, b+1);
        return scale(&next_vec, (n - 6.5)/2);
    }

    float O = 6 * random();
    float A = sqrt(random());
    Vector3 up = scale(&a->direction, -1);
    up = crossprod(&up, &sphere_ray_normal);
    up = normalize(&up);
    Vector3 T = crossprod(&up, &sphere_ray_normal);

    Vector3 direction = scale(&T, cos(O) * A);
    Vector3 k1 = scale(&up, sin(O) * A);
    direction = add(&direction, &k1);
    Vector3 k2 = scale(&sphere_ray_normal, sqrt(1 - A * A));
    direction = add(&direction, &k2);
    direction = normalize(&direction);

    Ray result = {P, direction};

    Vector3 j = {1, 1, 1};
    if (n == 3)
        j.x = j.z = 0;
    if (n == 4) j.y = j.z = 0;
    Vector3 next_vec = trace_ray(&result, b+1);
    next_vec = hadprod(&next_vec, &j);
    return next_vec;
}

int main() {
        FILE *file = fopen("cornell_box.ppm", "w");
        fprintf(file, "P3\n512 512\n255\n");
        for (int m = 0; m < 512; m++) {
            for (int n = 0; n < 512; n++) {
                Vector3 q = {0, 0, 0};
                for (int k = 0; k < 100; k++) {
                    Ray j = create_random(n / 256.0 - 1, m / 256.0 - 1);
                    Vector3 traced = trace_ray(&j, 0);
                    traced = scale(&traced, .02);
                    q = add(&q, &traced);
                }
                fprintf(file, "%d %d %d ", to_color(q.x), to_color(q.y), to_color(q.z));
            }
        }
}




