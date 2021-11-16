#ifndef MYSPHERE_H
#define MYSPHERE_H

#include "vec3.h"
#include "ray.h"

double hit_sphere(const point3& center, double radius, const ray& r) {
    vec3 oc = r.origin() - center;
    auto a = r.direction().length_squared();
    auto half_b = dot(oc, r.direction());
    auto c = oc.length_squared() - radius*radius;
    auto discriminant = (half_b*half_b) - (a*c);

    if (discriminant < 0) {
        return -1.0;
    } else {
        return (-half_b - sqrt(discriminant) ) / a;
    }
}

color ray_color(const ray& r) {

    auto t = hit_sphere(point3(0, 0, -1), 0.5, r);

    if (t > 0.0) {
        vec3 N = unit_vector(r.at(t) - vec3(0, 0, -1));
        return 0.5 * color(N.x() + 1, N.y() + 1, N.z() + 1);
    }
    
    vec3 unit_direction = unit_vector(r.direction());
    
    t = 0.5 * (unit_direction.y() + 1.0);

    return (1.0 - t) * color(1.0,1.0,1.0) + t*color(0.5,0.7,1.0);
}

#endif