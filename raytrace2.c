#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>

#define SPHERE_ALIGNMENT 1003
#define NUM_SPHERES 9
#define NUM_CUBES 1
#define TRACE_DEPTH 10
#define NUM_ITERATIONS 1
#define LIGHTS_CNT 2

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
    double c = dotprod(&from_center, &from_center);
    if (s->radius == INF) {
        c -= s->radius;
    }
    else {
        c -= s->radius * s->radius;
    }
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
    Vector3 center;
    double edge_size;
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
    Vector3 half_size = {c->edge_size / 2, c->edge_size / 2, c->edge_size / 2};
    Vector3 min_corner = {c->center.x - half_size.x, c->center.y - half_size.y, c->center.z - half_size.z};
    Vector3 max_corner = {c->center.x + half_size.x, c->center.y + half_size.y, c->center.z + half_size.z};

    double tmin = (min_corner.x - q->origin.x) / q->direction.x;
    double tmax = (max_corner.x - q->origin.x) / q->direction.x;
    if (tmin > tmax) { double temp = tmin; tmin = tmax; tmax = temp; }

    double tymin = (min_corner.y - q->origin.y) / q->direction.y;
    double tymax = (max_corner.y - q->origin.y) / q->direction.y;
    if (tymin > tymax) { double temp = tymin; tymin = tymax; tymax = temp; }

    if ((tmin > tymax) || (tymin > tmax)) return 0;
    if (tymin > tmin) tmin = tymin;
    if (tymax < tmax) tmax = tymax;

    double tzmin = (min_corner.z - q->origin.z) / q->direction.z;
    double tzmax = (max_corner.z - q->origin.z) / q->direction.z;
    if (tzmin > tzmax) { double temp = tzmin; tzmin = tzmax; tzmax = temp; }

    if ((tmin > tzmax) || (tzmin > tmax)) return 0;
    if (tzmin > tmin) tmin = tzmin;
    if (tzmax < tmax) tmax = tzmax;

    if (tmin > 0 && tmin < *t) {
        *t = tmin;
        return 1;
    }
    return 0;
}

// Collision detection for rotated cubes
int collide_rotated_cube(Cube* cube, Ray* ray, double* t, double angle) {
    // Transform the ray to cube's local coordinate system
    Vector3 rotated_origin = inverse_rotate_y(&ray->origin, angle, &cube->center);
    Vector3 rotated_direction = inverse_rotate_y(&ray->direction, angle, &cube->center);
    Ray local_ray = {rotated_origin, rotated_direction};

    // Perform AABB collision
    int result = collide_cube(cube, &local_ray, t);

    return result;
}

void calculate_cube_normal(Cube* c, Vector3* hit_point, Vector3* normal) {
    Vector3 local_hit = inverse_rotate_y(hit_point, c->angle_y, &c->center);

    double dx = fabs(local_hit.x - c->center.x);
    double dy = fabs(local_hit.y - c->center.y);
    double dz = fabs(local_hit.z - c->center.z);

    if (dx > dy && dx > dz) {
        normal->x = (local_hit.x > c->center.x) ? 1 : -1;
        normal->y = 0;
        normal->z = 0;
    } else if (dy > dx && dy > dz) {
        normal->x = 0;
        normal->y = (local_hit.y > c->center.y) ? 1 : -1;
        normal->z = 0;
    } else {
        normal->x = 0;
        normal->y = 0;
        normal->z = (local_hit.z > c->center.z) ? 1 : -1;
    }

    *normal = rotate_y(normal, c->angle_y, &c->center);
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
            {-2, 2, 3.5}, {.9, .5, .1}
        },
        {
            {1, -2, 4}, {.2, .5, .7}
        },
};

Sphere spheres[NUM_SPHERES] = {
        {
            {0, -SPHERE_ALIGNMENT, 0}, INF, {1, 0, 1}, 0, 0
        },
        {
            {0, SPHERE_ALIGNMENT, 0}, INF, {1, 1, 1}, 0, 0
        },
        {
            {SPHERE_ALIGNMENT, 0, 0}, INF, {0, 1, 0}, 1, 0
        },
        {
            {-SPHERE_ALIGNMENT, 0, 0}, INF, {1, 0, 0}, 0, 0
        },
        {
            {0, 0, -SPHERE_ALIGNMENT}, INF, {1, 1, 1}, 0, 0
        },
        {
            {0, 0, SPHERE_ALIGNMENT + 3}, INF, {1, 1, 1}, 0, 0
        },
        {
            {0, -0.2, 4}, 0.2, {.1, .2, .3}, 1, 0
        },
        {
            {2, -1, 4}, 1, {1, 1, 1}, .1, 0
        },
        {
            {0, -0.5, 3}, 1, {1, 0, 1}, .4, 1
        }
};

Cube cubes[NUM_CUBES] = {
    {{1, -1, 3}, 1.0, {0.8, 0.2, 0.2}, M_PI / 4, 0, 0.5}
};


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


Vector3 reflect(Vector3* direction, Vector3* normal) {
    Vector3 inverse_direction = scale(normal, -2 * dotprod(normal, direction));
    Vector3 reflection = add(direction, &inverse_direction);
    return normalize(&reflection);
}

Vector3 refract(Vector3* direction, Vector3* normal, double eta) {
    double cos_i = -dotprod(direction, normal);
    double sin2_t = eta * eta * (1.0 - cos_i * cos_i);

    // Total internal reflection
    if (sin2_t > 1.0) {
        return (Vector3){0, 0, 0}; // No refraction in this case
    }

    double cos_t = sqrt(1.0 - sin2_t);
    Vector3 refracted_parallel = scale(direction, eta);
    Vector3 refracted_perpendicular = scale(normal, eta * cos_i - cos_t);
    Vector3 refracted = add(&refracted_parallel, &refracted_perpendicular);

    return normalize(&refracted);
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
    } else if (sh == CUBE) {
        calculate_cube_normal(&cubes[n], &P, &normal);
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

    Vector3 eps_normal = scale(&normal, -EPSILON);
    P = add(&P, &eps_normal);

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

    Vector3 ambient = (Vector3){.001, .001, .001};
    next_vec = add(&next_vec, &ambient);

    next_vec = hadprod(&next_vec, &color);

    return next_vec;
}


int main()
{
    srand(time(0));

    unsigned char image[512][512][3];
    #pragma omp parallel for collapse(2) schedule(dynamic)
    for (int m = 0; m < 512; m++)
    {
        for (int n = 0; n < 512; n++)
        {
            Vector3 q = {0, 0, 0};
            for (int k = 0; k < NUM_ITERATIONS; k++)
            {
                Ray j = create_ray(n / 256.0 - 1, m / 256.0 - 1);
                Vector3 traced = trace_ray(&j, 1);
                q = add(&q, &traced);
            }
            double exposure = 2;
            q.x *= exposure / NUM_ITERATIONS;
            q.y *= exposure / NUM_ITERATIONS;
            q.z *= exposure / NUM_ITERATIONS;
            image[m][n][0] = to_color(q.x);
            image[m][n][1] = to_color(q.y);
            image[m][n][2] = to_color(q.z);
        }
    }

    FILE *file = fopen("cornell_box2.ppm", "wb");
    fprintf(file, "P3\n512 512\n255\n");

    for (int m = 0; m < 512; m++)
    {
        for (int n = 0; n < 512; n++)
        {
            fprintf(file, "%d %d %d ", image[m][n][0], image[m][n][1], image[m][n][2]);
        }
        fprintf(file, "\n");
    }

    fclose(file);
}





