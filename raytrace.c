#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>

#define BOX_SIZE 1003.0
#define NUM_OBJECTS 10
#define TRACE_DEPTH 7
#define NUM_ITERATIONS 100
#define COLLISION_THRESHOLD 1e-7

#ifndef INF
#define INF 1e6
#endif

#ifndef M_PI
#define M_PI 3.141592653589
#endif

typedef enum { SPHERE, CUBE } Shape;

double random() {
    return (rand() / (RAND_MAX + 1.0));
}

typedef struct {
    double x, y, z;
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
double dotprod(Vector3* vec1, Vector3* vec2) {
    return vec1->x * vec2->x + vec1->y * vec2->y + vec1->z * vec2->z;
}

// Make vector longer by f
Vector3 scale(Vector3* vec, double f) {
    return (Vector3){vec->x * f, vec->y * f, vec->z * f};
}

// Normalize the vector
Vector3 normalize(Vector3* vec) {
    Vector3 result = scale(vec, pow(sqrt(vec->x*vec->x + vec->y*vec->y + vec->z*vec->z), -1));
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
Vector3 point_at(Ray* r, double t) {
    Vector3 s = scale(&r->direction, t);
    return add(&s, &r->origin);
}



typedef struct {
    Vector3 center;
    double radius;
} Sphere;

int collide_sphere(Sphere* s, Ray* r, double* t) {
    Vector3 ray_reverse = scale(&r->origin, -1);
    Vector3 from_center = add(&ray_reverse, &s->center);
    double b = dotprod(&r->direction, &from_center);
    double c = dotprod(&from_center, &from_center) - s->radius;
    double D = b * b - c;
    if (D < 0) return 0;
    double solve = b - sqrt(D);
    if (solve > 0 && solve < *t) {
        *t = solve;
        return 1;
    }
    return 0;
}

typedef struct {
    Vector3 min_corner;
    Vector3 max_corner;
    double angle_y;
} Cube;

Vector3 rotate_y(Vector3* vec, double angle, Vector3* center) {
    double sin_a = sin(angle);
    double cos_a = cos(angle);

    // Translate to origin
    Vector3 translated = {
        vec->x - center->x,
        vec->y - center->y,
        vec->z - center->z
    };

    Vector3 rotated = {
        translated.x * cos_a - translated.z * sin_a,
        translated.y,
        translated.x * sin_a + translated.z * cos_a
    };

    return (Vector3){
        rotated.x + center->x,
        rotated.y + center->y,
        rotated.z + center->z
    };
}

Vector3 inverse_rotate_y(Vector3* vec, double angle, Vector3* center) {
    double sin_a = sin(-angle);
    double cos_a = cos(-angle);

    Vector3 rotated = {
        vec->x * cos_a - vec->z * sin_a,
        vec->y,
        vec->x * sin_a + vec->z * cos_a
    };
    return rotated;
}



// Collision detection for AABB (Axis-Aligned Bounding Box)
int collide_cube(Cube* c, Ray* q, double* t) {
    double t_min = 0, t_max = *t;

    for (int i = 0; i < 3; ++i) {
        double invD = 1.0f / (i == 0 ? q->direction.x : (i == 1 ? q->direction.y : q->direction.z));
        double t0 = ((i == 0 ? c->min_corner.x : (i == 1 ? c->min_corner.y : c->min_corner.z)) - (i == 0 ? q->origin.x : (i == 1 ? q->origin.y : q->origin.z))) * invD;
        double t1 = ((i == 0 ? c->max_corner.x : (i == 1 ? c->max_corner.y : c->max_corner.z)) - (i == 0 ? q->origin.x : (i == 1 ? q->origin.y : q->origin.z))) * invD;

        if (invD < 0) {
            double temp = t0;
            t0 = t1;
            t1 = temp;
        }

        t_min = t0 > t_min ? t0 : t_min;
        t_max = t1 < t_max ? t1 : t_max;

        if (t_max <= t_min) return 0;
    }

    *t = t_min;
    return 1; // Intersection detected
}

int collide_rotated_cube(Cube* cube, Ray* ray, double* t, double angle) {
    //Vector3 center = {0, 0, 0};
    Vector3 center = {
        (cube->min_corner.x + cube->max_corner.x) / 2.0,
        (cube->min_corner.y + cube->max_corner.y) / 2.0,
        (cube->min_corner.z + cube->max_corner.z) / 2.0
    };

    // Rotate the ray
    Ray rotated_ray = {
        rotate_y(&ray->origin, -angle, &center),
        inverse_rotate_y(&ray->direction, angle, &center)
    };

    // AABB collision
    return collide_cube(cube, &rotated_ray, t);
}




int to_color(double c) {
    return pow(c < 0 ? 0 : c > 1 ? 1 : c, .45) * 255 + .5;
}

Ray create_random(double x, double y) {
    Vector3 e = { x, -y, 1 };
    e = normalize(&e);
    e = scale(&e, 4);
    Vector3 position = {0, 0, 0};
    e = normalize(&e);
    Ray result = { position, e };
    return result;
}

Sphere box[NUM_OBJECTS] = {
        {
            {0, 2, 4}, 1
        },
        {
            {0, -BOX_SIZE, 0}, INF
        },
        {
            {0, BOX_SIZE, 0}, INF
        },
        {
            {BOX_SIZE, 0, 0}, INF
        },
        {
            {-BOX_SIZE, 0, 0}, INF
        },
        {
            {0, 0, -BOX_SIZE}, INF
        },
        {
            {0, 0, BOX_SIZE + 3}, INF
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
        {-1, -3, 3}, {1, -1, 5}, M_PI / 6
    }
};


int find_intersection(Ray* a, Shape* type, double* t) {
    int n = -1;
    for (int m = 0; m < NUM_OBJECTS; m++) {
        if (collide_sphere(&box[m], a, t)) {
            n = m;
            *type = SPHERE;
        }
    }

    for (int m = 0; m < 1; m++) {
        if (collide_rotated_cube(&cubes[m], a, t, cubes[m].angle_y)) {
            n = m;
            *type = CUBE;
        }
    }
    return n;
}

Vector3 trace_ray(Ray* a, int b) {
    double t = INF;
    Shape sh = -1;
    int n = find_intersection(a, &sh, &t);
    if (b > TRACE_DEPTH || n < 0) {
        Vector3 result = {0, 0, 0};
        return result;
    }
     if (!n && sh == SPHERE) {
        Vector3 result = {.9, .5, .1};
        return result;
     }

    Vector3 P = point_at(a, t);

    Vector3 normal;
    if (sh == SPHERE) {
        normal = box[n].center;
        normal = scale(&normal, -1);
        normal = add(&normal, &P);
        normal = normalize(&normal);
    }
    else {
        Cube* c = &cubes[n];
        Vector3 center = {
            (c->min_corner.x + c->max_corner.x) / 2.0,
            (c->min_corner.y + c->max_corner.y) / 2.0,
            (c->min_corner.z + c->max_corner.z) / 2.0
        };
        normal = center;
        normal = scale(&normal, -1);
        normal = add(&normal, &P);
        normal = normalize(&normal);
        /*
        if (fabs(P.x - c->min_corner.x) < COLLISION_THRESHOLD) normal = (Vector3){-1, 0, 0};
        else if (fabs(P.x - c->max_corner.x) < COLLISION_THRESHOLD) normal = (Vector3){1, 0, 0};
        else if (fabs(P.y - c->min_corner.y) < COLLISION_THRESHOLD) normal = (Vector3){0, -1, 0};
        else if (fabs(P.y - c->max_corner.y) < COLLISION_THRESHOLD) normal = (Vector3){0, 1, 0};
        else if (fabs(P.z - c->min_corner.z) < COLLISION_THRESHOLD) normal = (Vector3){0, 0, -1};
        else if (fabs(P.z - c->max_corner.z) < COLLISION_THRESHOLD) normal = (Vector3){0, 0, 1};
        normal = rotate_y(&normal, c->angle_y, &center);
        normal = normalize(&normal);
        */
    }

        if (n > 6) {
            Vector3 inverse_direction = scale(&normal, -2 * dotprod(&normal, &a->direction));
            Vector3 reflection = add(&a->direction, &inverse_direction);
            Ray r = { P, reflection };
            Vector3 next_vec = trace_ray(&r, b+1);
            return scale(&next_vec, (n - 6.5)/2);
        }

        double O = 2 * M_PI * random();
        double A = sqrt(random());
        Vector3 up = scale(&a->direction, -1);
        up = crossprod(&up, &normal);
        up = normalize(&up);
        Vector3 T = crossprod(&up, &normal);

        Vector3 direction = scale(&T, cos(O) * A);
        Vector3 k1 = scale(&up, sin(O) * A);
        direction = add(&direction, &k1);
        Vector3 k2 = scale(&normal, sqrt(1 - A * A));
        direction = add(&direction, &k2);
        direction = normalize(&direction);

        Ray result = (Ray){ P, direction};

    Vector3 j = {1, 1, 1};
    if (n == 3)
        j.x = j.z = 0;
    if (n == 4) j.y = j.z = 0;
    if (n == 0 && sh == CUBE) {
        j = (Vector3) { 0.2, 0.4, 0.6 };
    }
    Vector3 next_vec = trace_ray(&result, b+1);
    next_vec = hadprod(&next_vec, &j);
    return next_vec;
}

int main() {
        srand(time(0));
        FILE *file = fopen("cornell_box.ppm", "w");
        fprintf(file, "P3\n512 512\n255\n");
        #pragma omp parallel for collapse(2) schedule(dynamic)
        for (int m = 0; m < 512; m++) {
            for (int n = 0; n < 512; n++) {
                Vector3 q = {0, 0, 0};
                for (int k = 0; k < NUM_ITERATIONS; k++) {
                    Ray j = create_random(n / 256.0 - 1, m / 256.0 - 1);
                    Vector3 traced = trace_ray(&j, 0);
                    traced = scale(&traced, .02);
                    q = add(&q, &traced);
                }
                q = scale(&q, 0.9);
                fprintf(file, "%d %d %d ", to_color(q.x), to_color(q.y), to_color(q.z));
            }
        }
}




