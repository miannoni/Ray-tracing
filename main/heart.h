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

1/320 (

bool ellipsoid::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    vec3 oc = r.origin() - center;

    double Ax =
    double Ay =
    double Az =
    double dx =
    double dy =
    double dz =
    double t6 = -5*std::pow((4*std::pow(dx,2)+9*std::pow(dy,2)+4*std::pow(dz,2)),3)
    double t5 = -1920*Ax*std::pow(dx,5)-21870*Ay*std::pow(dy,5)-1920*Az*std::pow(dz,5)-9720*Ax*dx*std::pow(dy,4)-1920*Ax*dx*std::pow(dz,4)-4320*Ay*dy*std::pow(dz,4)-19440*Ay*std::pow(dx,2)*std::pow(dy,3)-3840*Az*std::pow(dx,2)*std::pow(dz,3)+320*std::pow(dx,2)*std::pow(dz,3)-8640*Az*std::pow(dy,2)*std::pow(dz,3)+36*std::pow(dy,2)*std::pow(dz,3)-8640*Ax*std::pow(dx,3)*std::pow(dy,2)-3840*Ax*std::pow(dx,3)*std::pow(dz,2)-19440*Ay*std::pow(dy,3)*std::pow(dz,2)-8640*Ax*dx*std::pow(dy,2)*std::pow(dz,2)-8640*Ay*std::pow(dx,2)*dy*std::pow(dz,2)-4320*Ay*std::pow(dx,4)*dy-1920*Az*std::pow(dx,4)*dz-9720*Az*std::pow(dy,4)*dz-8640*Az*std::pow(dx,2)*std::pow(dy,2)*dz
    double t4 = -4800*std::pow(Ax,2)*std::pow(dx,4)-2160*std::pow(Ay,2)*std::pow(dx,4)-960*std::pow(Az,2)*std::pow(dx,4)+960*std::pow(dx,4)-4860*std::pow(Ax,2)*std::pow(dy,4)-54675*std::pow(Ay,2)*std::pow(dy,4)-4860*std::pow(Az,2)*std::pow(dy,4)+4860*std::pow(dy,4)-960*std::pow(Ax,2)*std::pow(dz,4)-2160*std::pow(Ay,2)*std::pow(dz,4)-4800*std::pow(Az,2)*std::pow(dz,4)+960*std::pow(dz,4)-38880*Ax*Ay*dx*std::pow(dy,3)+640*Ax*dx*std::pow(dz,3)-7680*Ax*Az*dx*std::pow(dz,3)+72*Ay*dy*std::pow(dz,3)-17280*Ay*Az*dy*std::pow(dz,3)-12960*std::pow(Ax,2)*std::pow(dx,2)*std::pow(dy,2)-29160*std::pow(Ay,2)*std::pow(dx,2)*std::pow(dy,2)-4320*std::pow(Az,2)*std::pow(dx,2)*std::pow(dy,2)+4320*std::pow(dx,2)*std::pow(dy,2)-5760*std::pow(Ax,2)*std::pow(dx,2)*std::pow(dz,2)-4320*std::pow(Ay,2)*std::pow(dx,2)*std::pow(dz,2)-5760*std::pow(Az,2)*std::pow(dx,2)*std::pow(dz,2)+960*Az*std::pow(dx,2)*std::pow(dz,2)+1920*std::pow(dx,2)*std::pow(dz,2)-4320*std::pow(Ax,2)*std::pow(dy,2)*std::pow(dz,2)-29160*std::pow(Ay,2)*std::pow(dy,2)*std::pow(dz,2)-12960*std::pow(Az,2)*std::pow(dy,2)*std::pow(dz,2)+108*Az*std::pow(dy,2)*std::pow(dz,2)+4320*std::pow(dy,2)*std::pow(dz,2)-17280*Ax*Ay*dx*dy*std::pow(dz,2)-17280*Ax*Ay*std::pow(dx,3)*dy-7680*Ax*Az*std::pow(dx,3)*dz-38880*Ay*Az*std::pow(dy,3)*dz-17280*Ax*Az*dx*std::pow(dy,2)*dz-17280*Ay*Az*std::pow(dx,2)*dy*dz
    double t3 = -6400*std::pow(Ax,3)*std::pow(dx,3)-8640*Ax*std::pow(Ay,2)*std::pow(dx,3)-3840*Ax*std::pow(Az,2)*std::pow(dx,3)+3840*Ax*std::pow(dx,3)-72900*std::pow(Ay,3)*std::pow(dy,3)-19440*Ay*std::pow(Az,2)*std::pow(dy,3)-19440*std::pow(Ax,2)*Ay*std::pow(dy,3)+19440*Ay*std::pow(dy,3)-6400*std::pow(Az,3)*std::pow(dz,3)+320*std::pow(Ax,2)*std::pow(dz,3)+36*std::pow(Ay,2)*std::pow(dz,3)-3840*std::pow(Ax,2)*Az*std::pow(dz,3)-8640*std::pow(Ay,2)*Az*std::pow(dz,3)+3840*Az*std::pow(dz,3)-8640*std::pow(Ax,3)*dx*std::pow(dy,2)-58320*Ax*std::pow(Ay,2)*dx*std::pow(dy,2)-8640*Ax*std::pow(Az,2)*dx*std::pow(dy,2)+8640*Ax*dx*std::pow(dy,2)-3840*std::pow(Ax,3)*dx*std::pow(dz,2)-8640*Ax*std::pow(Ay,2)*dx*std::pow(dz,2)-11520*Ax*std::pow(Az,2)*dx*std::pow(dz,2)+3840*Ax*dx*std::pow(dz,2)+1920*Ax*Az*dx*std::pow(dz,2)-19440*std::pow(Ay,3)*dy*std::pow(dz,2)-25920*Ay*std::pow(Az,2)*dy*std::pow(dz,2)-8640*std::pow(Ax,2)*Ay*dy*std::pow(dz,2)+8640*Ay*dy*std::pow(dz,2)+216*Ay*Az*dy*std::pow(dz,2)-19440*std::pow(Ay,3)*std::pow(dx,2)*dy-8640*Ay*std::pow(Az,2)*std::pow(dx,2)*dy-25920*std::pow(Ax,2)*Ay*std::pow(dx,2)*dy+8640*Ay*std::pow(dx,2)*dy-3840*std::pow(Az,3)*std::pow(dx,2)*dz+960*std::pow(Az,2)*std::pow(dx,2)*dz-11520*std::pow(Ax,2)*Az*std::pow(dx,2)*dz-8640*std::pow(Ay,2)*Az*std::pow(dx,2)*dz+3840*Az*std::pow(dx,2)*dz-8640*std::pow(Az,3)*std::pow(dy,2)*dz+108*std::pow(Az,2)*std::pow(dy,2)*dz-8640*std::pow(Ax,2)*Az*std::pow(dy,2)*dz-58320*std::pow(Ay,2)*Az*std::pow(dy,2)*dz+8640*Az*std::pow(dy,2)*dz-34560*Ax*Ay*Az*dx*dy*dz
    double t2 = -4800*std::pow(Ax,4)*std::pow(dx,2)-4860*std::pow(Ay,4)*std::pow(dx,2)-960*std::pow(Az,4)*std::pow(dx,2)+320*std::pow(Az,3)*std::pow(dx,2)+5760*std::pow(Ax,2)*std::pow(dx,2)-12960*std::pow(Ax,2)*std::pow(Ay,2)*std::pow(dx,2)+4320*std::pow(Ay,2)*std::pow(dx,2)-5760*std::pow(Ax,2)*std::pow(Az,2)*std::pow(dx,2)-4320*std::pow(Ay,2)*std::pow(Az,2)*std::pow(dx,2)+1920*std::pow(Az,2)*std::pow(dx,2)-960*std::pow(dx,2)-2160*std::pow(Ax,4)*std::pow(dy,2)-54675*std::pow(Ay,4)*std::pow(dy,2)-2160*std::pow(Az,4)*std::pow(dy,2)+36*std::pow(Az,3)*std::pow(dy,2)+4320*std::pow(Ax,2)*std::pow(dy,2)-29160*std::pow(Ax,2)*std::pow(Ay,2)*std::pow(dy,2)+29160*std::pow(Ay,2)*std::pow(dy,2)-4320*std::pow(Ax,2)*std::pow(Az,2)*std::pow(dy,2)-29160*std::pow(Ay,2)*std::pow(Az,2)*std::pow(dy,2)+4320*std::pow(Az,2)*std::pow(dy,2)-2160*std::pow(dy,2)-960*std::pow(Ax,4)*std::pow(dz,2)-4860*std::pow(Ay,4)*std::pow(dz,2)-4800*std::pow(Az,4)*std::pow(dz,2)+1920*std::pow(Ax,2)*std::pow(dz,2)-4320*std::pow(Ax,2)*std::pow(Ay,2)*std::pow(dz,2)+4320*std::pow(Ay,2)*std::pow(dz,2)-5760*std::pow(Ax,2)*std::pow(Az,2)*std::pow(dz,2)-12960*std::pow(Ay,2)*std::pow(Az,2)*std::pow(dz,2)+5760*std::pow(Az,2)*std::pow(dz,2)+960*std::pow(Ax,2)*Az*std::pow(dz,2)+108*std::pow(Ay,2)*Az*std::pow(dz,2)-960*std::pow(dz,2)-38880*Ax*std::pow(Ay,3)*dx*dy-17280*Ax*Ay*std::pow(Az,2)*dx*dy-17280*std::pow(Ax,3)*Ay*dx*dy+17280*Ax*Ay*dx*dy-7680*Ax*std::pow(Az,3)*dx*dz+1920*Ax*std::pow(Az,2)*dx*dz-7680*std::pow(Ax,3)*Az*dx*dz-17280*Ax*std::pow(Ay,2)*Az*dx*dz+7680*Ax*Az*dx*dz-17280*Ay*std::pow(Az,3)*dy*dz+216*Ay*std::pow(Az,2)*dy*dz-38880*std::pow(Ay,3)*Az*dy*dz-17280*std::pow(Ax,2)*Ay*Az*dy*dz+17280*Ay*Az*dy*dz
    double t1 = -1920*Ax^5*dx-9720*Ax*sd::pow(Ay,4)*dx-1920*Ax*sd::pow(Az,4)*dx+3840*sd::pow(Ax,3)*dx+640*Ax*sd::pow(Az,3)*dx-8640*sd::pow(Ax,3)*sd::pow(Ay,2)*dx+8640*Ax*sd::pow(Ay,2)*dx-3840*sd::pow(Ax,3)*sd::pow(Az,2)*dx-8640*Ax*sd::pow(Ay,2)*sd::pow(Az,2)*dx+3840*Ax*sd::pow(Az,2)*dx-1920*Ax*dx-21870*Ay^5*dy-4320*Ay*sd::pow(Az,4)*dy-19440*sd::pow(Ax,2)*sd::pow(Ay,3)*dy+19440*sd::pow(Ay,3)*dy+72*Ay*sd::pow(Az,3)*dy-19440*sd::pow(Ay,3)*sd::pow(Az,2)*dy-8640*sd::pow(Ax,2)*Ay*sd::pow(Az,2)*dy+8640*Ay*sd::pow(Az,2)*dy-4320*sd::pow(Ax,4)*Ay*dy+8640*sd::pow(Ax,2)*Ay*dy-4320*Ay*dy-1920*Az^5*dz-3840*sd::pow(Ax,2)*sd::pow(Az,3)*dz-8640*sd::pow(Ay,2)*sd::pow(Az,3)*dz+3840*sd::pow(Az,3)*dz+960*sd::pow(Ax,2)*sd::pow(Az,2)*dz+108*sd::pow(Ay,2)*sd::pow(Az,2)*dz-1920*sd::pow(Ax,4)*Az*dz-9720*sd::pow(Ay,4)*Az*dz+3840*sd::pow(Ax,2)*Az*dz-8640*sd::pow(Ax,2)*sd::pow(Ay,2)*Az*dz+8640*sd::pow(Ay,2)*Az*dz-1920*Az*dz
    double t0 = -320*Ax^6-3645*Ay^6-320*Az^6+960*std::pow(Ax,4)-4860*std::pow(Ax,2)*std::pow(Ay,4)+4860*std::pow(Ay,4)-960*std::pow(Ax,2)*std::pow(Az,4)-2160*std::pow(Ay,2)*std::pow(Az,4)+960*std::pow(Az,4)+320*std::pow(Ax,2)*std::pow(Az,3)+36*std::pow(Ay,2)*std::pow(Az,3)-960*std::pow(Ax,2)-2160*std::pow(Ax,4)*std::pow(Ay,2)+4320*std::pow(Ax,2)*std::pow(Ay,2)-2160*std::pow(Ay,2)-960*std::pow(Ax,4)*std::pow(Az,2)-4860*std::pow(Ay,4)*std::pow(Az,2)+1920*std::pow(Ax,2)*std::pow(Az,2)-4320*std::pow(Ax,2)*std::pow(Ay,2)*std::pow(Az,2)+4320*std::pow(Ay,2)*std::pow(Az,2)-960*std::pow(Az,2)+320
    t6 = t6/320
    t5 = t5/320
    t4 = t4/320
    t3 = t3/320
    t2 = t2/320
    t1 = t1/320
    t0 = t0/320
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