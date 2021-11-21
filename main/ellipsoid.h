#ifndef ELLIPSOID_H
#define ELLIPSOID_H

#include "hittable.h"
#include "vec3.h"

class ellipsoid : public hittable {
    public:
        ellipsoid() {}
        ellipsoid(point3 cen, double aa, double bb, double cc) : center(cen), aa(aa), bb(bb), cc(cc){};

        virtual bool hit(
            const ray& r, double t_min, double t_max, hit_record& rec) const override;

    public:
        point3 center;
        double aa;
        double bb;
        double cc;

};

bool ellipsoid::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    vec3 oc = r.origin() - center;


    auto a = (r.direction().x()/aa)*(r.direction().x()/aa) + (r.direction().y()/bb)*(r.direction().y()/bb) + (r.direction().z()/cc)*(r.direction().z()/cc);
    auto half_b = (oc.x()*r.direction().x()/(aa*aa)) + (oc.y()*r.direction().y()/(bb*bb)) + (oc.z()*r.direction().z()/(cc*cc));
    auto c = (oc.x()/aa)*(oc.x()/aa) + (oc.y()/bb)*(oc.y()/bb) + (oc.z()/cc)*(oc.z()/cc) - 1;


//    auto a = r.direction().length_squared();
//    auto half_b = dot(oc, r.direction());
//    auto c = oc.length_squared() - radius*radius;

//
//    std::cerr << "r.direction().x():" << std::to_string(r.direction().x()) << "\n";
//    std::cerr << "oc.x():" << std::to_string(oc.x()) << "\n";
//    std::cerr << "oc.y():" << std::to_string(oc.y()) << "\n";
//    std::cerr << "oc.z():" << std::to_string(oc.z()) << "\n";
//    std::cerr << "a:" << std::to_string(a) << "\n";
//    std::cerr << "a_meu:" << std::to_string(a_meu) << "\n";
//    std::cerr << "half_b:" << std::to_string(half_b) << "\n";
//    std::cerr << "half_b_meu:" << std::to_string(half_b_meu) << "\n";
//    std::cerr << "c:" << std::to_string(c) << "\n";
//    std::cerr << "c_meu:" << std::to_string(c_meu) << "\n";

    auto discriminant = (half_b * half_b) - a*c;
    if (discriminant < 0) return false;
    auto sqrtd = sqrt(discriminant);

    auto root = (-half_b - sqrtd) / a;
    if (root < t_min || t_max < root) {
        root = (-half_b + sqrtd) / a;
        if (root < t_min || t_max < root)
            return false;
    }

    rec.t = root;
    rec.p = r.at(rec.t);
    vec3 outward_normal = (rec.p - center) / aa;
    rec.set_face_normal(r, outward_normal);
    // rec.normal = (rec.p - center) / radius;

    return true;

}

#endif