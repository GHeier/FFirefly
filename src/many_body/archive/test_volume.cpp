#include <iostream>
#include <cmath>

struct Vec3 {
    double x, y, z;
    Vec3(double x_=0, double y_=0, double z_=0) : x(x_), y(y_), z(z_) {}
    Vec3 operator-(const Vec3& v) const { return Vec3(x-v.x, y-v.y, z-v.z); }
    double dot(const Vec3& v) const { return x*v.x + y*v.y + z*v.z; }
    Vec3 cross(const Vec3& v) const {
        return Vec3(y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x);
    }
};

double volume(Vec3 k[4]) {
    Vec3 v1 = k[1] - k[0];
    Vec3 v2 = k[2] - k[0];
    Vec3 v3 = k[3] - k[0];
    return std::abs(v1.dot(v2.cross(v3))) / 6.0;
}

int main() {
    // Triangle in 2D (degenerate tetrahedron)
    Vec3 k[4];
    k[0] = Vec3(0, 0, 0);
    k[1] = Vec3(1, 0, 0);
    k[2] = Vec3(0, 1, 0);
    k[3] = Vec3(0, 0, 0); // duplicate
    
    std::cout << "Triangle volume (should be 0.5): " << volume(k) << "\n";
    std::cout << "Expected area = 0.5 * 1 * 1 = 0.5\n";
    
    return 0;
}
