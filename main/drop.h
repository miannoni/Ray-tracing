#ifndef DROP_H
#define DROP_H


#include "hittable.h"
#include "vec3.h"
#include <limits>
#include <math.h>
#include "../aux/poly34.h"

#define PI 3.14159265
// poly34.cpp : solution of cubic and quartic equation
// (c) Khashin S.I. http://math.ivanovo.ac.ru/dalgebra/Khashin/index.html
// khash2 (at) gmail.com
// Thanks to Alexandr Rakhmanin <rakhmanin (at) gmail.com>
// public domain
//
//#include <math.h>

//#include "../aux/poly34.h"		   // solution of cubic and quartic equation
#define	TwoPi  6.28318530717958648
const double eps=1e-14;

//=============================================================================
// _root3, root3 from http://prografix.narod.ru
//=============================================================================
static double _root3 ( double x )
{
    double s = 1.;
    while ( x < 1. )
    {
        x *= 8.;
        s *= 0.5;
    }
    while ( x > 8. )
    {
        x *= 0.125;
        s *= 2.;
    }
    double r = 1.5;
    r -= 1./3. * ( r - x / ( r * r ) );
    r -= 1./3. * ( r - x / ( r * r ) );
    r -= 1./3. * ( r - x / ( r * r ) );
    r -= 1./3. * ( r - x / ( r * r ) );
    r -= 1./3. * ( r - x / ( r * r ) );
    r -= 1./3. * ( r - x / ( r * r ) );
    return r * s;
}

double root3 ( double x )
{
    if ( x > 0 ) return _root3 ( x ); else
    if ( x < 0 ) return-_root3 (-x ); else
    return 0.;
}


// x - array of size 2
// return 2: 2 real roots x[0], x[1]
// return 0: pair of complex roots: x[0]�i*x[1]
int   SolveP2(double *x, double a, double b) {			// solve equation x^2 + a*x + b = 0
	double D = 0.25*a*a - b;
	if (D >= 0) {
		D = sqrt(D);
		x[0] = -0.5*a + D;
		x[1] = -0.5*a - D;
		return 2;
	}
	x[0] = -0.5*a;
	x[1] = sqrt(-D);
	return 0;
}
//---------------------------------------------------------------------------
// x - array of size 3
// In case 3 real roots: => x[0], x[1], x[2], return 3
//         2 real roots: x[0], x[1],          return 2
//         1 real root : x[0], x[1] � i*x[2], return 1
int SolveP3(double *x,double a,double b,double c) {	// solve cubic equation x^3 + a*x^2 + b*x + c = 0
	double a2 = a*a;
    double q  = (a2 - 3*b)/9;
	double r  = (a*(2*a2-9*b) + 27*c)/54;
	// equation x^3 + q*x + r = 0
    double r2 = r*r;
	double q3 = q*q*q;
	double A,B;
	if (r2 <= (q3 + eps)) {//<<-- FIXED!
		double t=r/sqrt(q3);
		if( t<-1) t=-1;
		if( t> 1) t= 1;
        t=acos(t);
        a/=3; q=-2*sqrt(q);
        x[0]=q*cos(t/3)-a;
        x[1]=q*cos((t+TwoPi)/3)-a;
        x[2]=q*cos((t-TwoPi)/3)-a;
        return(3);
    } else {
        //A =-pow(fabs(r)+sqrt(r2-q3),1./3);
        A =-root3(fabs(r)+sqrt(r2-q3));
		if( r<0 ) A=-A;
		B = (A==0? 0 : B=q/A);

		a/=3;
		x[0] =(A+B)-a;
        x[1] =-0.5*(A+B)-a;
        x[2] = 0.5*sqrt(3.)*(A-B);
		if(fabs(x[2])<eps) { x[2]=x[1]; return(2); }
        return(1);
    }
}// SolveP3(double *x,double a,double b,double c) {
//---------------------------------------------------------------------------
// a>=0!
void  CSqrt( double x, double y, double &a, double &b) // returns:  a+i*s = sqrt(x+i*y)
{
	double r  = sqrt(x*x+y*y);
	if( y==0 ) {
		r = sqrt(r);
		if(x>=0) { a=r; b=0; } else { a=0; b=r; }
	} else {		// y != 0
		a = sqrt(0.5*(x+r));
		b = 0.5*y/a;
	}
}
//---------------------------------------------------------------------------
int   SolveP4Bi(double *x, double b, double d)	// solve equation x^4 + b*x^2 + d = 0
{
	double D = b*b-4*d;
	if( D>=0 )
	{
		double sD = sqrt(D);
		double x1 = (-b+sD)/2;
		double x2 = (-b-sD)/2;	// x2 <= x1
		if( x2>=0 )				// 0 <= x2 <= x1, 4 real roots
		{
			double sx1 = sqrt(x1);
			double sx2 = sqrt(x2);
			x[0] = -sx1;
			x[1] =  sx1;
			x[2] = -sx2;
			x[3] =  sx2;
			return 4;
		}
		if( x1 < 0 )				// x2 <= x1 < 0, two pair of imaginary roots
		{
			double sx1 = sqrt(-x1);
			double sx2 = sqrt(-x2);
			x[0] =    0;
			x[1] =  sx1;
			x[2] =    0;
			x[3] =  sx2;
			return 0;
		}
		// now x2 < 0 <= x1 , two real roots and one pair of imginary root
			double sx1 = sqrt( x1);
			double sx2 = sqrt(-x2);
			x[0] = -sx1;
			x[1] =  sx1;
			x[2] =    0;
			x[3] =  sx2;
			return 2;
	} else { // if( D < 0 ), two pair of compex roots
		double sD2 = 0.5*sqrt(-D);
		CSqrt(-0.5*b, sD2, x[0],x[1]);
		CSqrt(-0.5*b,-sD2, x[2],x[3]);
		return 0;
	} // if( D>=0 )
} // SolveP4Bi(double *x, double b, double d)	// solve equation x^4 + b*x^2 d
//---------------------------------------------------------------------------
#define SWAP(a,b) { t=b; b=a; a=t; }
static void  dblSort3( double &a, double &b, double &c) // make: a <= b <= c
{
	double t;
	if( a>b ) SWAP(a,b);	// now a<=b
	if( c<b ) {
		SWAP(b,c);			// now a<=b, b<=c
		if( a>b ) SWAP(a,b);// now a<=b
	}
}
//---------------------------------------------------------------------------
int   SolveP4De(double *x, double b, double c, double d)	// solve equation x^4 + b*x^2 + c*x + d
{
	//if( c==0 ) return SolveP4Bi(x,b,d); // After that, c!=0
	if( fabs(c)<1e-14*(fabs(b)+fabs(d)) ) return SolveP4Bi(x,b,d); // After that, c!=0

	int res3 = SolveP3( x, 2*b, b*b-4*d, -c*c);	// solve resolvent
	// by Viet theorem:  x1*x2*x3=-c*c not equals to 0, so x1!=0, x2!=0, x3!=0
	if( res3>1 )	// 3 real roots,
	{
		dblSort3(x[0], x[1], x[2]);	// sort roots to x[0] <= x[1] <= x[2]
		// Note: x[0]*x[1]*x[2]= c*c > 0
		if( x[0] > 0) // all roots are positive
		{
			double sz1 = sqrt(x[0]);
			double sz2 = sqrt(x[1]);
			double sz3 = sqrt(x[2]);
			// Note: sz1*sz2*sz3= -c (and not equal to 0)
			if( c>0 )
			{
				x[0] = (-sz1 -sz2 -sz3)/2;
				x[1] = (-sz1 +sz2 +sz3)/2;
				x[2] = (+sz1 -sz2 +sz3)/2;
				x[3] = (+sz1 +sz2 -sz3)/2;
				return 4;
			}
			// now: c<0
			x[0] = (-sz1 -sz2 +sz3)/2;
			x[1] = (-sz1 +sz2 -sz3)/2;
			x[2] = (+sz1 -sz2 -sz3)/2;
			x[3] = (+sz1 +sz2 +sz3)/2;
			return 4;
		} // if( x[0] > 0) // all roots are positive
		// now x[0] <= x[1] < 0, x[2] > 0
		// two pair of comlex roots
		double sz1 = sqrt(-x[0]);
		double sz2 = sqrt(-x[1]);
		double sz3 = sqrt( x[2]);

		if( c>0 )	// sign = -1
		{
			x[0] = -sz3/2;
			x[1] = ( sz1 -sz2)/2;		// x[0]�i*x[1]
			x[2] =  sz3/2;
			x[3] = (-sz1 -sz2)/2;		// x[2]�i*x[3]
			return 0;
		}
		// now: c<0 , sign = +1
		x[0] =   sz3/2;
		x[1] = (-sz1 +sz2)/2;
		x[2] =  -sz3/2;
		x[3] = ( sz1 +sz2)/2;
		return 0;
	} // if( res3>1 )	// 3 real roots,
	// now resoventa have 1 real and pair of compex roots
	// x[0] - real root, and x[0]>0,
	// x[1]�i*x[2] - complex roots,
	// x[0] must be >=0. But one times x[0]=~ 1e-17, so:
	if (x[0] < 0) x[0] = 0;
	double sz1 = sqrt(x[0]);
	double szr, szi;
	CSqrt(x[1], x[2], szr, szi);  // (szr+i*szi)^2 = x[1]+i*x[2]
	if( c>0 )	// sign = -1
	{
		x[0] = -sz1/2-szr;			// 1st real root
		x[1] = -sz1/2+szr;			// 2nd real root
		x[2] = sz1/2;
		x[3] = szi;
		return 2;
	}
	// now: c<0 , sign = +1
	x[0] = sz1/2-szr;			// 1st real root
	x[1] = sz1/2+szr;			// 2nd real root
	x[2] = -sz1/2;
	x[3] = szi;
	return 2;
} // SolveP4De(double *x, double b, double c, double d)	// solve equation x^4 + b*x^2 + c*x + d
//-----------------------------------------------------------------------------
double N4Step(double x, double a,double b,double c,double d)	// one Newton step for x^4 + a*x^3 + b*x^2 + c*x + d
{
	double fxs= ((4*x+3*a)*x+2*b)*x+c;	// f'(x)
	if (fxs == 0) return x;	//return 1e99; <<-- FIXED!
	double fx = (((x+a)*x+b)*x+c)*x+d;	// f(x)
	return x - fx/fxs;
}
//-----------------------------------------------------------------------------
// x - array of size 4
// return 4: 4 real roots x[0], x[1], x[2], x[3], possible multiple roots
// return 2: 2 real roots x[0], x[1] and complex x[2]�i*x[3],
// return 0: two pair of complex roots: x[0]�i*x[1],  x[2]�i*x[3],
int   SolveP4(double *x,double a,double b,double c,double d) {	// solve equation x^4 + a*x^3 + b*x^2 + c*x + d by Dekart-Euler method
	// move to a=0:
	double d1 = d + 0.25*a*( 0.25*b*a - 3./64*a*a*a - c);
	double c1 = c + 0.5*a*(0.25*a*a - b);
	double b1 = b - 0.375*a*a;
	int res = SolveP4De( x, b1, c1, d1);
	if( res==4) { x[0]-= a/4; x[1]-= a/4; x[2]-= a/4; x[3]-= a/4; }
	else if (res==2) { x[0]-= a/4; x[1]-= a/4; x[2]-= a/4; }
	else             { x[0]-= a/4; x[2]-= a/4; }
	// one Newton step for each real root:
	if( res>0 )
	{
		x[0] = N4Step(x[0], a,b,c,d);
		x[1] = N4Step(x[1], a,b,c,d);
	}
	if( res>2 )
	{
		x[2] = N4Step(x[2], a,b,c,d);
		x[3] = N4Step(x[3], a,b,c,d);
	}
	return res;
}
//-----------------------------------------------------------------------------
#define F5(t) (((((t+a)*t+b)*t+c)*t+d)*t+e)
//-----------------------------------------------------------------------------
double SolveP5_1(double a,double b,double c,double d,double e)	// return real root of x^5 + a*x^4 + b*x^3 + c*x^2 + d*x + e = 0
{
	int cnt;
	if( fabs(e)<eps ) return 0;

	double brd =  fabs(a);			// brd - border of real roots
	if( fabs(b)>brd ) brd = fabs(b);
	if( fabs(c)>brd ) brd = fabs(c);
	if( fabs(d)>brd ) brd = fabs(d);
	if( fabs(e)>brd ) brd = fabs(e);
	brd++;							// brd - border of real roots

	double x0, f0;					// less than root
	double x1, f1;					// greater than root
	double x2, f2, f2s;				// next values, f(x2), f'(x2)
	double dx;

	if( e<0 ) { x0 =   0; x1 = brd; f0=e; f1=F5(x1); x2 = 0.01*brd; }	// positive root
	else	  { x0 =-brd; x1 =   0; f0=F5(x0); f1=e; x2 =-0.01*brd; }	// negative root

	if( fabs(f0)<eps ) return x0;
	if( fabs(f1)<eps ) return x1;

	// now x0<x1, f(x0)<0, f(x1)>0
	// Firstly 10 bisections
	for( cnt=0; cnt<10; cnt++)
	{
		x2 = (x0 + x1) / 2;					// next point
		//x2 = x0 - f0*(x1 - x0) / (f1 - f0);		// next point
		f2 = F5(x2);				// f(x2)
		if( fabs(f2)<eps ) return x2;
		if( f2>0 ) { x1=x2; f1=f2; }
		else       { x0=x2; f0=f2; }
	}

	// At each step:
	// x0<x1, f(x0)<0, f(x1)>0.
	// x2 - next value
	// we hope that x0 < x2 < x1, but not necessarily
	do {
		if(cnt++>50) break;
		if( x2<=x0 || x2>= x1 ) x2 = (x0 + x1)/2;	// now  x0 < x2 < x1
		f2 = F5(x2);								// f(x2)
		if( fabs(f2)<eps ) return x2;
		if( f2>0 ) { x1=x2; f1=f2; }
		else       { x0=x2; f0=f2; }
		f2s= (((5*x2+4*a)*x2+3*b)*x2+2*c)*x2+d;		// f'(x2)
		if( fabs(f2s)<eps ) { x2=1e99; continue; }
		dx = f2/f2s;
		x2 -= dx;
	} while(fabs(dx)>eps);
	return x2;
} // SolveP5_1(double a,double b,double c,double d,double e)	// return real root of x^5 + a*x^4 + b*x^3 + c*x^2 + d*x + e = 0
//-----------------------------------------------------------------------------
int   SolveP5(double *x,double a,double b,double c,double d,double e)	// solve equation x^5 + a*x^4 + b*x^3 + c*x^2 + d*x + e = 0
{
	double r = x[0] = SolveP5_1(a,b,c,d,e);
	double a1 = a+r, b1=b+r*a1, c1=c+r*b1, d1=d+r*c1;
	return 1+SolveP4(x+1, a1,b1,c1,d1);
} // SolveP5(double *x,double a,double b,double c,double d,double e)	// solve equation x^5 + a*x^4 + b*x^3 + c*x^2 + d*x + e = 0
//-----------------------------------------------------------------------------











class drop : public hittable {
    public:
        drop() {}
        drop(point3 cen, double size_x, double size_y, double size_z) : center(cen), size_x(size_x), size_y(size_y), size_z(size_z){};

        virtual bool hit(
            const ray& r, double t_min, double t_max, hit_record& rec) const override;



    public:
        point3 center;
        double size_x;
        double size_y;
        double size_z;

};

bool drop::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    vec3 oc = r.origin() - center;

    double dz = r.direction().x()*size_x;
//    std::cerr << "dx: " << std::to_string(dx) << "\n";
    double dy = r.direction().z()*size_z;
    double dx = -r.direction().y()*size_y;
    double ez = oc.x()*size_x;
    double ex = -oc.y()*size_y;
    double ey = oc.z()*size_z;


    double t4 = std::pow(dx, 4);
    double a = 4*ex*std::pow(dx, 3) - std::pow(dx, 3);
    double b = 6*std::pow(ex, 2)*std::pow(dx, 2) - 3*ex*std::pow(dx, 2) + std::pow(dy, 2) + std::pow(dz, 2);
    double c = 4*std::pow(ex, 3)*dx - 3*std::pow(ex, 2)*dx + 2*ey*dy + 2*ez*dz;
    double d = std::pow(ex, 4) - std::pow(ex, 3) + std::pow(ey, 2) + std::pow(ez, 2);

    a = a/t4;
    b = b/t4;
    c = c/t4;
    d = d/t4;


//α^4 - α^3 + β^2 + γ^2 +
//ρ^4 t^4 +
//4 α ρ^3 t^3 - ρ^3 t^3 +
//6 α^2 ρ^2 t^2 - 3 α ρ^2 t^2 + σ^2 t^2 + t^2 τ^2 +
//4 α^3 ρ t - 3 α^2 ρ t + 2 β σ t + 2 γ t τ = 0

    double min_root = std::numeric_limits<int>::max();
    double roots[4];
    int real_roots = SolveP4( roots, a, b, c, d);
    if (real_roots == 4){
        int found = 0;
        for (int i = 0; i < 4; i++){
            if (roots[i] > 0){
                found = 1;
                if (min_root > roots[i]){
                    min_root = roots[i];
                }
            }
        }
        if (!found){
            return false;
        }
    } else if (real_roots == 2){
        int found = 0;
        for (int i = 0; i < 2; i++){
            if (roots[i] > 0){
                found = 1;
                if (min_root > roots[i]){
                    min_root = roots[i];
                }
            }
        }
        if (!found){
            return false;
        }
    } else {
        return false;
    }




    rec.t = min_root;
    rec.p = r.at(rec.t);
    double u = asin ((rec.p - center).y()/(rec.p - center).length()) * 180.0 / PI;
    double v = asin ((rec.p - center).z()/(rec.p - center).length()) * 180.0 / PI;
    double tu_x = 0.5*sin(u*PI/180)*sin(v*PI/180)*(cos(v*PI/180)-1);
    double tu_y = -0.5*cos(u*PI/180)*sin(v*PI/180)*(cos(v*PI/180)-1);
    double tu_z = 0;
    double tv_x = cos(u*PI/180)*(0.5*std::pow(sin(v*PI/180), 2)-0.5*(cos(v*PI/180) - 1)*cos(v*PI/180));
    double tv_y = sin(u*PI/180)*(0.5*std::pow(sin(v*PI/180), 2)-0.5*(cos(v*PI/180) - 1)*cos(v*PI/180));
    double tv_z = -sin(v*PI/180);
    vec3 tu(tu_x,tu_y,tu_z);
    vec3 tv(tv_x,tv_y,tv_z);

//    tu.x() = 0.5*sin(u*PI/180)*sin(v*PI/180)*(cos(v*PI/180)-1);
//    tu.y() = -0.5*cos(u*PI/180)*sin(v*PI/180)*(cos(v*PI/180)-1);
//    tu.z() = 0;
//    tv.x() = cos(u*PI/180)*(0.5*std::pow(sin(v*PI/180), 2)-0.5*(cos(v*PI/180) - 1)*cos(v*PI/180));
//    tv.y() = sin(u*PI/180)*(0.5*std::pow(sin(v*PI/180), 2)-0.5*(cos(v*PI/180) - 1)*cos(v*PI/180));
//    tv.z() = -sin(v*PI/180);







//    std::cerr << "a: " << std::to_string(a) << "\n";
//    std::cerr << "b: " << std::to_string(b) << "\n";
//    std::cerr << "c: " << std::to_string(c) << "\n";
//    std::cerr << "d: " << std::to_string(d) << "\n";
//    std::cerr << "x[0]: " << std::to_string(x[0]) << "\n";
//    std::cerr << "x[1]: " << std::to_string(x[1]) << "\n";
//    std::cerr << "x[2]: " << std::to_string(x[2]) << "\n";
//    std::cerr << "x[3]: " << std::to_string(x[3]) << "\n";

//    auto a = (r.direction().x()/aa)*(r.direction().x()/aa) + (r.direction().y()/bb)*(r.direction().y()/bb) + (r.direction().z()/cc)*(r.direction().z()/cc);
//    auto half_b = (oc.x()*r.direction().x()/(aa*aa)) + (oc.y()*r.direction().y()/(bb*bb)) + (oc.z()*r.direction().z()/(cc*cc));
//    auto c = (oc.x()/aa)*(oc.x()/aa) + (oc.y()/bb)*(oc.y()/bb) + (oc.z()/cc)*(oc.z()/cc) - 1;


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

//    auto discriminant = (half_b * half_b) - a*c;

//    if (discriminant < 0) return false;
//    auto sqrtd = sqrt(discriminant);
//
//    auto root = (-half_b - sqrtd) / a;
//    if (root < t_min || t_max < root) {
//        root = (-half_b + sqrtd) / a;
//        if (root < t_min || t_max < root)
//            return false;
//    }


//    vec3 outward_normal = (rec.p - center) / a;
    vec3 outward_normal = cross(tu, tv);
    rec.set_face_normal(r, outward_normal);
//    rec.normal = (rec.p - center);
//    rec.normal = rec.normal / rec.normal.length();
    rec.normal = outward_normal / outward_normal.length();

//    for (int i = 0; i < 3; i++){
//        std::cerr << "tu: " << std::to_string(tu[i]) << "\n";
//        std::cerr << "tv: " << std::to_string(tv[i]) << "\n";
////        std::cerr << "cross: " << std::to_string(crozzz[i]) << "\n";
//    }

    return true;
//    return false;
}

#endif