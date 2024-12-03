#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>

#define BOX_SIZE 1003
#define NUM_SPHERES 10
#define NUM_CUBES 1
#define TRACE_DEPTH 10
#define NUM_ITERATIONS 1
#define LIGHTS_CNT 1

#define COLLISION_THRESHOLD 1e-6
#define EPSILON 1e-8


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
    Vector3 color;
    double reflect;
    double refract;
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
    Vector3 color;
    double angle_y;
    double reflect;
    double refract;
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
        double delim = (i == 0 ? q->direction.x : (i == 1 ? q->direction.y : q->direction.z));
        double invD = 1.0f / ((fabs(delim) > EPSILON) ? delim : EPSILON);

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

Ray create_ray(double x, double y) {
    Vector3 e = { x, -y , 1 };
    e = normalize(&e);
    //e = scale(&e, 4);
    Vector3 position = {0, 0, 0};;
    Ray result = { position, e };
    return result;
}

typedef struct {
    Vector3 position;
    Vector3 color;
} LightSource;

LightSource lights[LIGHTS_CNT] = {
        {
            {0, 2, 4}, {.9, .5, .1}
        },
        // {
        //     {1, -2, 4}, {.2, .5, .7}
        // },
};

Sphere spheres[NUM_SPHERES] = {
        {
            {0, -BOX_SIZE, 0}, INF, {1, 0, 1}, 0, 0
        },
        {
            {0, BOX_SIZE, 0}, INF, {1, 1, 1}, 0, 0
        },
        {
            {BOX_SIZE, 0, 0}, INF, {0, 1, 0}, 0, 0
        },
        {
            {-BOX_SIZE, 0, 0}, INF, {1, 0, 0}, 0, 0
        },
        {
            {0, 0, -BOX_SIZE}, INF, {1, 1, 1}, 0, 0
        },
        {
            {0, 0, BOX_SIZE + 3}, INF, {1, 1, 1}, 0, 0
        },
        {
            {-1.4, -2, 3}, 0.2, {.3, .5, .7}, 0, 1.5
        },
        {
            {0, -0.2, 4}, 0.2, {.1, .2, .3}, 1, 0
        },
        {
            {2, -1, 4}, 1, {1, 1, 1}, .1, 0
        },
        {
            {0, -0.5, 3}, 1, {1, 0, 1}, .4, 0
        }
};

Cube cubes[NUM_CUBES] = {
    // {
    //     {-3, -1, 3}, {-3.9, 1, 5}, {1, 0, 1}, 0, 0, 0
    // }
    {
        {-1, -3, 3}, {1, -1, 5}, {1, .9, 1}, M_PI / 6, 0, 0
    }
};

// Cube border_cubes[6] = {
//     {{-BOX_SIZE, -BOX_SIZE, -BOX_SIZE}, {BOX_SIZE, BOX_SIZE, -BOX_SIZE}, 0},  // Front
//     {{-BOX_SIZE, -BOX_SIZE, BOX_SIZE}, {BOX_SIZE, BOX_SIZE, BOX_SIZE}, 0},   // Back
//     {{-BOX_SIZE, -BOX_SIZE, -BOX_SIZE}, {-BOX_SIZE, BOX_SIZE, BOX_SIZE}, M_PI / 2}, // Left
//     {{BOX_SIZE, -BOX_SIZE, -BOX_SIZE}, {BOX_SIZE, BOX_SIZE, BOX_SIZE}, M_PI / 2},  // Right
//     {{-BOX_SIZE, -BOX_SIZE, -BOX_SIZE}, {BOX_SIZE, -BOX_SIZE, BOX_SIZE}, M_PI}, // Bottom
//     {{-BOX_SIZE, BOX_SIZE, -BOX_SIZE}, {BOX_SIZE, BOX_SIZE, BOX_SIZE}, M_PI}  // Top
// };




int find_intersection(Ray* a, Shape* type, double* t) {
    int n = -1;
    for (int m = 0; m < NUM_SPHERES; m++) {
        if (collide_sphere(&spheres[m], a, t)) {
            n = m;
            *type = SPHERE;
        }
    }

    for (int m = 0; m < NUM_CUBES; m++) {
        if (collide_rotated_cube(&cubes[m], a, t, cubes[m].angle_y)) {
            n = m;
            *type = CUBE;
        }
    }

    return n;
}

Vector3 refract(Vector3* direction, Vector3* normal, double eta) {
    double cos_i = -dotprod(direction, normal);

    Vector3 local_normal = *normal;

    if (cos_i < 0) {
        cos_i *= -1.0;
        local_normal = scale(normal, -1);
        eta = 1/eta;
    }

    double sin2_t = eta * eta * (1 - cos_i * cos_i); // Snell's Law to compute sin^2(theta_t)

    double cos_t = sqrt(1.0 - sin2_t); // Cosine of the refraction angle

    if (cos_t < 0) {
        return (Vector3){0, 0, 0};
    }

    // Calculate the refracted direction using Snell's Law
    Vector3 refracted = scale(direction, eta);
    Vector3 scaled_normal = scale(&local_normal, eta * cos_i - cos_t);
    return add(&refracted, &scaled_normal);
}

Vector3 reflect(Vector3* direction, Vector3* normal) {
    Vector3 inverse_direction = scale(normal, -2 * dotprod(normal, direction));
    Vector3 reflection = add(direction, &inverse_direction);
    return normalize(&reflection);
}

int is_point_in_shadow(Vector3* point, LightSource* light) {
    Vector3 to_light = {light->position.x - point->x,
                        light->position.y - point->y,
                        light->position.z - point->z};
    Ray shadow_ray = {*point, normalize(&to_light)};

    double t = sqrt(dotprod(&to_light, &to_light));
    Shape type;
    int hit = find_intersection(&shadow_ray, &type, &t);

    // If t is smaller than the distance to the light, the point is shadowed.
    return hit && t < sqrt(dotprod(&to_light, &to_light));
}


Vector3 trace_ray(Ray* a, int b) {
    if (b > TRACE_DEPTH) return (Vector3){0, 0, 0};
    Vector3 next_vec = {0, 0, 0};
    double t = INF;
    Shape sh;
    int n = find_intersection(a, &sh, &t);

    if (n < 0) {
        Vector3 result = {0, 0, 0};
        return result;
    }
    //  if (sh == SPHERE && dotprod(&spheres[n].emitted, &e) != 0) {
    //     Vector3 result = spheres[n].emitted;
    //     return result;
    //  }

    Vector3 P = point_at(a, t);

    Vector3 normal;
    Vector3 view_dir = scale(&a->direction, -1);

    if (sh == SPHERE) {
        normal = scale(&spheres[n].center, -1);
        normal = add(&normal, &P);
        normal = normalize(&normal);
    } else {
        Cube* c = &cubes[n];
        if (fabs(P.x - c->min_corner.x) < COLLISION_THRESHOLD) {
            normal = (Vector3){-1.0, 0.0, 0.0};
        } else if (fabs(P.x - c->max_corner.x) < COLLISION_THRESHOLD) {
            normal = (Vector3){1.0, 0.0, 0.0};
        } else if (fabs(P.y - c->min_corner.y) < COLLISION_THRESHOLD) {
            normal = (Vector3){0.0, -1.0, 0.0};
        } else if (fabs(P.y - c->max_corner.y) < COLLISION_THRESHOLD) {
            normal = (Vector3){0.0, 1.0, 0.0};
        } else if (fabs(P.z - c->min_corner.z) < COLLISION_THRESHOLD) {
            normal = (Vector3){0.0, 0.0, -1.0};
        } else {
            normal = (Vector3){0.0, 0.0, 1.0};
        }
    }

    double refr, refl;
    Vector3 color;

    if (sh == SPHERE) {
        refr = spheres[n].refract;
        refl = spheres[n].reflect;
        color = spheres[n].color;
    }
    else {
        refr = cubes[n].refract;
        refl = cubes[n].reflect;
        color = cubes[n].color;
    }
    // Vector3 ambient = scale(&color, 0.1);
    // next_vec = add(&next_vec, &ambient);


    if (refr > 0) {
        Vector3 refracted_direction = refract(&a->direction, &normal, 1.5);
        Ray refracted_ray = {P, refracted_direction};
        Vector3 refraction = trace_ray(&refracted_ray, b+1);
        refraction = scale(&refraction, refr);
        next_vec = add(&next_vec, &refraction);
    }

    if (refl > 0) {
        Vector3 reflected_direction = reflect(&a->direction, &normal);
        Ray reflected_ray = {P, reflected_direction};
        Vector3 reflection = trace_ray(&reflected_ray, b+1);
        reflection = scale(&reflection, refl);
        next_vec = add(&next_vec, &reflection);
    }

    Vector3 diffuse = {0, 0, 0};

    for (int i = 0; i < LIGHTS_CNT; i++) {
        if (is_point_in_shadow(&P, &lights[i]))
            continue;

        Vector3 light_dir = {
            lights[i].position.x - P.x,
            lights[i].position.y - P.y,
            lights[i].position.z - P.z
        };
        light_dir = normalize(&light_dir);

        double dot = fmax(0.0, dotprod(&normal, &light_dir));
        Vector3 diff_intensity = scale(&lights[i].color, dot);
        next_vec = add(&next_vec, &diff_intensity);

        Vector3 reflect_dir = reflect(&light_dir, &normal);
        double spec_intensity = pow(fmax(0, dotprod(&reflect_dir, &view_dir)), 64);
        Vector3 specular = scale(&lights[i].color, spec_intensity);
        next_vec = add(&next_vec, &specular);
    }

    next_vec = add(&next_vec, &diffuse);

    Vector3 ambient = scale(&color, 0.2);
    next_vec = add(&next_vec, &ambient);

    next_vec = hadprod(&next_vec, &color);

    return next_vec;
}

int main() {
        srand(time(0));
        FILE *file = fopen("cornell_box.ppm", "w");
        fprintf(file, "P3\n512 512\n255\n");
        //#pragma omp parallel for collapse(2) schedule(dynamic, 16)
        for (int m = 0; m < 512; m++) {
            for (int n = 0; n < 512; n++) {
                Vector3 q = {0, 0, 0};
                for (int k = 0; k < NUM_ITERATIONS; k++) {
                    Ray j = create_ray(n / 256.0 - 1, m / 256.0 - 1);
                    Vector3 traced = trace_ray(&j, 1);
                    q = add(&q, &traced);
                }
                double exposure = 1;
                q.x *= exposure / NUM_ITERATIONS;
                q.y *= exposure / NUM_ITERATIONS;
                q.z *= exposure / NUM_ITERATIONS;

                fprintf(file, "%d %d %d ", to_color(q.x), to_color(q.y), to_color(q.z));
            }
        }
}




