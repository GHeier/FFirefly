
float integrate_susceptibility(Vec q, float T, float fermi_energy, float w, int num_points) {
    //return imaginary_integration(q, T, fermi_energy, w, num_points, 0.000);
    if (q.norm() < 0.0001) {
        vector<Vec> FS = get_FS(fermi_energy);
        float sum = 0; for (auto x : FS) sum += x.area / vp(x.n, x);
        return sum / pow(2*k_max,dim);
    }
    if (q.norm() < 0.0001) q = Vec(0.01,0.01,0.01);
    auto func = [q, w, T](Vec k) -> float {
        return ratio(k, q, w, T);
    };
    auto e_diff = [q](Vec k) -> float {
        return epsilon(q.n, k+q) - epsilon(q.n, k);
    };
    auto denom = [q, w](Vec k) -> float { 
        return w - (epsilon(q.n, k+q) - epsilon(q.n, k)); 
    };
    auto denom_diff = [q, w](Vec k) -> float { 
        return vp_diff(k.n, k, q); 
    };

    float a, b; get_surface_transformed_bounds(b, a, e_diff);
    vector<float> spacing; get_spacing_vec(spacing, w, a, b, num_points);
    //float c = tetrahedron_sum_continuous(denominator, denominator_diff, q, spacing, w, T);
    float c = surface_transform_integral(func, denom, denom_diff, spacing);
    return c / pow(2*k_max,dim);
}

complex<float> complex_susceptibility_integration(Vec q, float T, float fermi_energy, complex<float> w, int num_points) {
    auto func = [T, q, w, fermi_energy](float x, float y, float z) -> complex<float> {
        Vec k(x,y,z);
        float e_k = epsilon(k.n, k) - fermi_energy;
        float e_kq = epsilon(k.n, k+q) - fermi_energy;
        float f_kq = fermi_dirac(e_kq, T);
        float f_k = fermi_dirac(e_k, T);
        if (fabs(e_kq - e_k) < 0.0001 and fabs(w.imag()) < 0.0001 and fabs(w.real()) < 0.0001) {
            if (T == 0 or exp(e_k/T) > 1e6) return e_k < 0;
            return 1/T * exp(e_k/T) / pow( exp(e_k/T) + 1,2);
        }
        return (f_kq - f_k) / (w - (e_kq - e_k));
    };
    auto func_r = [T, q, w, fermi_energy](Vec k) -> float {
        float e_k = epsilon(k.n, k) - fermi_energy;
        float e_kq = epsilon(k.n, k+q) - fermi_energy;
        float f_kq = fermi_dirac(e_kq, T);
        float f_k = fermi_dirac(e_k, T);
        if (fabs(e_kq - e_k) < 0.0001 and fabs(w.imag()) < 0.0001 and fabs(w.real()) < 0.0001) {
            if (T == 0 or exp(e_k/T) > 1e6) return e_k < 0;
            return 1/T * exp(e_k/T) / pow( exp(e_k/T) + 1,2);
        }
        return (f_kq - f_k)*(w.real() - (e_kq - e_k)) / 
            ((w.real() - (e_kq - e_k))*(w.real() - (e_kq - e_k)) + w.imag()*w.imag());
    };
    auto func_i = [T, q, w, fermi_energy](Vec k) -> float {
        float e_k = epsilon(k.n, k) - fermi_energy;
        float e_kq = epsilon(k.n, k+q) - fermi_energy;
        float f_kq = fermi_dirac(e_kq, T);
        float f_k = fermi_dirac(e_k, T);
        if (fabs(e_kq - e_k) < 0.0001 and fabs(w.imag()) < 0.0001 and fabs(w.real()) < 0.0001) {
            if (T == 0 or exp(e_k/T) > 1e6) return e_k < 0;
            return 1/T * exp(e_k/T) / pow( exp(e_k/T) + 1,2);
        }
        return -(f_kq - f_k) * w.imag() /
            ((w.real() - (e_kq - e_k))*(w.real() - (e_kq - e_k)) + w.imag()*w.imag());
    };
    auto denom = [q, w](Vec k) -> float { 
        return w.real() - (epsilon(k.n, k+q) - epsilon(k.n, k)); 
    };
    auto denom_diff = [q, w](Vec k) -> float { 
        return vp_diff(k.n, k, q); 
    };
    float d = pow(2*k_max,dim);
    float a, b; get_surface_transformed_bounds(b, a, denom);
    vector<float> spacing; get_spacing_vec(spacing, w.real(), a, b, num_points);
    float c_real = surface_transform_integral(func_r, denom, denom_diff, spacing);
    float c_imag = surface_transform_integral(func_i, denom, denom_diff, spacing);
    complex<float> c(c_real, c_imag);
    return c / d;

    return complex_trapezoidal_integration(func, -k_max, k_max, -k_max, k_max, 
            -k_max, k_max, num_points) / d;
}

vector<vector<vector<float>>> chi_cube(float T, float fermi_energy, float w, string message) {
    int num_points = (dim == 3) ? 50 : 1000; // Number of integral surfaces
    int m_z = m*(dim%2) + 3*((dim+1)%2);
    vector<vector<vector<float>>> cube(m, vector<vector<float>> (m, vector<float> (m_z)));
    unordered_map<string, float> map;
    float empty_val = -98214214.0;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            for (int k = 0; k < m_z; k++) {
                Vec q((2.0*k_max*i)/(1.0*m-1), (2.0*k_max*j)/(1.0*m-1), (2.0*k_max*k)/(1.0*m_z-1));
                Vec q2 = to_IBZ_2(q);
                if (map.find(vec_to_string(q2)) == map.end())
                    map[vec_to_string(q2)] = empty_val;
            }
        }
    }
    
    //cout << "Taking " << map.size() << " integrals in " << dim << " dimensions.\n";
    //#pragma omp parallel for
    for(unsigned int i = 0; i < map.size(); i++) {
        auto datIt = map.begin();
        advance(datIt, i);
        string key = datIt->first;
        map[key] = integrate_susceptibility(string_to_vec(key), T, fermi_energy, w, num_points);
        progress_bar(1.0 * i / (map.size()-1), message);
    }
    cout << endl;

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            for (int k = 0; k < m_z; k++) {
                Vec q((2*k_max*i)/(m-1), (2*k_max*j)/(m-1), (2*k_max*k)/(m_z-1));
                Vec q2 = to_IBZ_2(q);
                cube[i][j][k] = map[vec_to_string(q2)];
            }
        }
    }

    return cube;
}

float calculate_chi_from_cube(const vector<vector<vector<float>>> &chi_cube, Vec q) {
    Vec v = to_IBZ_2(q);
    float d = 2*k_max/(m-1);

    float x = v(0), y = v(1), z = v(2);
    if (dim == 2) z = 0;

    int i = floor(x / d);
    int j = floor(y / d);
    int k = floor(z / d);

    float x1 = i * d; 
    float y1 = j * d; 
    float z1 = k * d; 

    float x2 = x1 + d; 
    float y2 = y1 + d; 
    float z2 = z1 + d; 

    float dx = 0, dy = 0, dz = 0, wx = 0, wy = 0, wz = 0, w0 = 0;

    // Make sure there's no issue with indexing
    //cout << q << q(2) << endl;
    //int s = chi_cube.size()-1; 
    //assert( i < s and j < s and k < chi_cube[0][0].size()-1);

    float f1 = chi_cube[i][j][k], f2 = chi_cube[i+1][j][k];
    float f3 = chi_cube[i+1][j+1][k], f4 = chi_cube[i][j+1][k];
    float f5 = chi_cube[i][j][k+1], f6 = chi_cube[i+1][j][k+1];
    float f7 = chi_cube[i+1][j+1][k+1], f8 = chi_cube[i][j+1][k+1];

    if (x - x1 <= z2 - z and x - x1 >= y - y1) {// blue @ 1
        w0 = f1;
        wx = (f2 - f1) / d; dx = x - x1;
        wy = (f3 - f2) / d; dy = y - y1;
        wz = (f5 - f1) / d; dz = z - z1;
    }

    else if (y + z <= z1 + y2 and x - x1 <= y - y1) {// orange @ 1
        w0 = f1;
        wx = (f3 - f4) / d; dx = x - x1;
        wy = (f4 - f1) / d; dy = y - y1;
        wz = (f5 - f1) / d; dz = z - z1;
    }

    else if (x + z >= z1 + x2 and y + z <= z1 + y2) {// red @ 2
        w0 = f2;
        wx = (f6 - f5) / d; dx = x - x2;
        wy = (f3 - f2) / d; dy = y - y1;
        wz = (f6 - f2) / d; dz = z - z1;
    }

    else if (y + z >= z1 + y2 and x + z <= z1 + x2) {// purple @ 4
        w0 = f4;
        wx = (f3 - f4) / d; dx = x - x1;
        wy = (f8 - f5) / d; dy = y - y2;
        wz = (f8 - f4) / d; dz = z - z1;
    }

    else if (x - x1 >= y - y1 and y + z >= z1 + y2) {// teal @ 7
        w0 = f7;
        wx = (f6 - f5) / d; dx = x - x2;
        wy = (f7 - f6) / d; dy = y - y2;
        wz = (f7 - f3) / d; dz = z - z2;
    }

    else if (x - x1 <= y - y1 and y + z >= z1 + y2) {// green @ 7
        w0 = f7;
        wx = (f7 - f8) / d; dx = x - x2;
        wy = (f8 - f5) / d; dy = y - y2;
        wz = (f7 - f3) / d; dz = z - z2;
    }
    else return f1 + (f2-f1)/d*x + (f4-f1)/d*y + (f5-f1)/d*z;

    return w0 + wx*dx + wy*dy + wz*dz;
}


float trap_cube(auto &f, float x0, float x1, float y0, float y1, float z0, float z1) {
    return (x1-x0)*(y1-y0)*(z1-z0) / 8 * (f(x0,y0,z0) + f(x0,y0,z1) + f(x0,y1,z0) + f(x0,y1,z1) + f(x1,y0,z0) + f(x1,y0,z1) + f(x1,y1,z0) + f(x1,y1,z1));
}

float trap_8_cubes(auto &f, float x0, float x1, float y0, float y1, float z0, float z1) {
    return trap_cube(f, x0, (x0+x1)/2, y0, (y0+y1)/2, z0, (z0+z1)/2) 
        + trap_cube(f, x0, (x0+x1)/2, y0, (y0+y1)/2, (z0+z1)/2, z1) 
        + trap_cube(f, x0, (x0+x1)/2, (y0+y1)/2, y1, z0, (z0+z1)/2) 
        + trap_cube(f, x0, (x0+x1)/2, (y0+y1)/2, y1, (z0+z1)/2, z1) 
        + trap_cube(f, (x0+x1)/2, x1, y0, (y0+y1)/2, z0, (z0+z1)/2) 
        + trap_cube(f, (x0+x1)/2, x1, y0, (y0+y1)/2, (z0+z1)/2, z1) 
        + trap_cube(f, (x0+x1)/2, x1, (y0+y1)/2, y1, z0, (z0+z1)/2) 
        + trap_cube(f, (x0+x1)/2, x1, (y0+y1)/2, y1, (z0+z1)/2, z1);
}

float adaptive_trapezoidal(auto &f, float x0, float x1, float y0, float y1, float z0, float z1, int xdivs, int ydivs, int zdivs, float error_relative) {
    float sum = 0;

    float dx = (x1 - x0) / (xdivs);
    float dy = (y1 - y0) / (ydivs);
    float dz = (z1 - z0) / (zdivs);
    if (dim == 2) dz = 1;

    for (int i = 0; i < xdivs; i++) {
        float x = x0 + i*dx;
        for (int j = 0; j < ydivs; j++) {
            float y = y0 + j*dy;
            for (float k = 0; k < zdivs; k++) {
                float z = z0 + k*dz;

                float t1 = trap_cube(f, x, x+dx, y, y+dy, z, z+dz);
                float t2 = trap_8_cubes(f, x, x+dx, y, y+dy, z, z+dz);

                if (fabs(t1 - t2) < error_relative * fabs(t2) or fabs(t1 - t2) < 0.0001) {
                    sum += t2;
                }
                else {
//                    cout << t1 << " " << t2 << " " << fabs(t1 - t2) << " " << error_relative * fabs(t2) << endl;
                    float new_zdiv = 2 * (dim % 2) + 1 * ((dim+1)%2);
                    sum += adaptive_trapezoidal(f, x, x+dx, y, y+dy, z, z+dz, 2, 2, new_zdiv, error_relative);
                }

            }
        }
    }
    return sum;
}

float iteratively_splitting_cubes(auto &f, float x0, float x1, float y0, float y1, float z0, float z1, float error_total, float error_relative) {

    float total_sum = 0;
    float dx = (x1 - x0) / 2;
    bool no_errors = false;

    for (int i = 1; not no_errors; i++) {
        total_sum = 0;
        no_errors = true;
        int iters = pow(2,i);
        //#pragma omp parallel for reduction(+:total_sum)
        for (int j = 0; j < iters; j++) {
            for (int k = 0; k < iters; k++) {
                for (int l = 0; l < iters; l++) {
                    float t1 = trap_cube(f, x0+j*dx, x0+(j+1)*dx, y0+k*dx, y0+(k+1)*dx, z0+l*dx, z0+(l+1)*dx);
                    float t2 = trap_8_cubes(f, x0+j*dx, x0+(j+1)*dx, y0+k*dx, y0+(k+1)*dx, z0+l*dx, z0+(l+1)*dx); 
                    float err = fabs(t2 - t1);
                    if ( err > fabs(error_relative*t2) and err > error_total / pow(2,3*i) )
                        no_errors = false;

                    total_sum += t2;
                }
            }
        }
        //cout << "Splits: " << i << endl;
        dx /= 2;
    }
    return total_sum;
}

// Note: Should only be used for 2D, fails in 3D
vector<Vec> sort_by_adjacent(vector<Vec> &FS) {
    vector<Vec> sorted_FS;
    Vec empty;
    Vec k0 = FS[0];
    sorted_FS.push_back(k0);
    FS[0] = empty;
    while (sorted_FS.size() < FS.size()) {
        float min_dist = 100000;
        Vec min_vec;
        for (Vec k : FS) {
            if (k == empty) continue;
            float dist = (k - k0).norm();
            if (dist < min_dist) {
                min_dist = dist;
                min_vec = k;
            }
        }
        sorted_FS.push_back(min_vec);
        cout << sorted_FS.size() << " ";
        k0 = min_vec;
        for (unsigned int i = 0; i < FS.size(); i++) {
            if (FS[i] == min_vec) {
                FS[i] = empty;
                break;
            }
        }
    }
    return sorted_FS;
}

vector<VecAndEnergy> FS_to_q_shifted(vector<Vec> &FS, Vec q) {
    vector<VecAndEnergy> q_shifted;
    for (Vec k : FS) {
        float energy = epsilon(k+q);
        q_shifted.push_back({k, energy});
    }
    return q_shifted;
}

// NOTE: Changes delta value as well
int get_num_points_from_delta(float &delta) {
    if (delta == 0) {
        delta = 0.0001;
    }
    int pts = 10*k_max/delta + 1;
    if (delta > 0.1) {
        delta = 0;
    }
    return pts;
}

void sanitize_I_vals(float &V1, float &V2, float &V3, float &V4) {
    if (fabs(V1 - V2) < 1e-3) V2 = V1;
    if (fabs(V1 - V3) < 1e-3) V3 = V1;
    if (fabs(V1 - V4) < 1e-3) V4 = V1;
    if (fabs(V2 - V3) < 1e-3) V3 = V2;
    if (fabs(V2 - V4) < 1e-3) V4 = V2;
    if (fabs(V3 - V4) < 1e-3) V4 = V3;
}

vector<float> getUnique(float a, float b, float c, float d) {
    // Use a set to find unique values 
    set<float, greater<float>> uniqueValues = {a, b, c, d};
    // Copy the sorted unique values to a vector
    vector<float> result(uniqueValues.begin(), uniqueValues.end());
    return result;
}

bool check_two_equal(float V1, float V2, float V3, float V4) {
    return (V1 == V2 and V3 == V4) or (V1 == V3 and V2 == V4) or (V1 == V4 and V2 == V3);
}

// Integral value of the tetrahedron method when interpolated linearly across each small cube
float get_I(float D1, float D2, float D3, float V1, float V2, float V3, float V4) {
    sanitize_I_vals(V1, V2, V3, V4);
    vector<float> V = getUnique(V1, V2, V3, V4);
    if (find(V.begin(), V.end(), 0) != V.end()) {
        printf("Option 0\n");
        return 0;
    }
    if (V.size() == 1) {
        printf("Option 1\n");
        float r = 1 / V[0];
        return r;
    }
    if (V.size() == 2 and check_two_equal(V1, V2, V3, V4)) {
        printf("Option 2\n");
        float t1 = 2 * V[0] * V[1] / pow(V[0] - V[1], 3) * log(fabs(V[1] / V[0]));
        float t2 = (V[0] + V[1]) / (pow(V[0] - V[1], 2));
        float r = 3 * (t1 + t2);
        return 3 * (t1 + t2);
    }
    else if (V.size() == 2) {
        printf("Option 3\n");
        float t1 = V[1]*V[1] / (pow(V[0] - V[1],3)) * log(fabs(V[0]/V[1]));
        float t2 = (1.5 * V[1]*V[1] + 0.5 * V[0]*V[0] - 2 * V[0]*V[1]) / pow(V[0] - V[1],3);
        float r = 3 * (t1 + t2);
        return r;
    }
    if (V.size() == 3) {
        printf("Option 4\n");
        float t1 = V[1] * V[1] / (pow(V[1] - V[0], 2) * (V[1] - V[2])) * log(fabs(V[1] / V[0]));
        float t2 = V[2] * V[2] / (pow(V[2] - V[0], 2) * (V[2] - V[1])) * log(fabs(V[2] / V[0]));
        float t3 = V[0] / ((V[1] - V[0]) * (V[2] - V[0]));
        float r = 3 * (t1 + t2 + t3);
        return r;
    }
    printf("Option 5\n");
    float t1 = (V1*V1/D1*log(fabs(V1/V4)));
    float t2 = (V2*V2/D2*log(fabs(V2/V4)));
    float t3 = (V3*V3/D3*log(fabs(V3/V4)));
    float r = 3 * (t1 + t2 + t3);
    return r;
}

bool scattering_available(Vec q, vector<Vec> p) {
    for (Vec v : p) {
        if (epsilon(v) < mu and epsilon(v+q) < mu) return false;
        if (epsilon(v) > mu and epsilon(v+q) > mu) return false;
    }
    return true;
}

// Computes the sum analytically, which should be quite a bit faster
// Done at zero temperature is the only caveat
float analytic_tetrahedron_linear_energy_method(Vec q, float w, int num_pts) {
    vector<vector<float>> tetrahedrons {
        {1, 2, 3, 5}, 
        {1, 3, 4, 5},
        {2, 5, 6, 3},
        {4, 5, 8, 3},
        {5, 8, 7, 3},
        {5, 6, 7, 3}
    };

    float upper = 1.1;
    float lower = 1.0;
    float width = (upper - lower) / 2;
    double sum = 0;
    double Omega = pow(2*width,3) / (6*num_pts*num_pts*num_pts);
    float d3x = 1 / pow(num_pts, 3);
    if (dim == 2) Omega = pow(2*width,2) / (2*n*n);
    //#pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < num_pts; i++) {
        for (int j = 0; j < num_pts; j++) {
            for (int k = 0; k < num_pts * (dim%2) + 1 * ((dim+1)%2); k++) {
                //Vec temp(2*width * i / (num_pts-1) - lower, 2*width * j / (num_pts-1) - lower, 2*width * k / (num_pts-1) - lower);
                Vec temp(i * (upper - lower) / (num_pts-1) + lower, j * (upper - lower) / (num_pts-1) + lower, k * (upper - lower) / (num_pts-1) + lower);
                vector<Vec> points = points_from_indices(epsilon, i, j, k, num_pts);
                //sum += points[0](0) * Omega * 6;

                for (int c = 0; c < 6; c++) {

                    vector<Vec> ep_points(4);
                    for (int p = 0; p < 4; p++) {
                        ep_points[p] = points[tetrahedrons[c][p]-1];
                    }
                    //if (not scattering_available(q, ep_points)) continue;
                    // Sort by e(k)
                    sort(ep_points.begin(), ep_points.end(), [](Vec a, Vec b) { 
                            return a.w < b.w; 
                    });
                    float Theta_kq = epsilon(ep_points[3]+q) < mu ? 1 : 0;
                    float Theta_k = ep_points[0].w < mu ? 1 : 0;

                    float V1 = e_diff(ep_points[3],q) - w;
                    float V2 = e_diff(ep_points[2],q) - w;
                    float V3 = e_diff(ep_points[1],q) - w;
                    float V4 = e_diff(ep_points[0],q) - w;

                    float D1 = (V1 - V4) * (V1 - V3) * (V1 - V2);
                    float D2 = (V2 - V4) * (V2 - V3) * (V2 - V1);
                    float D3 = (V3 - V4) * (V3 - V2) * (V3 - V1);

                    float I = get_I(D1, D2, D3, V1, V2, V3, V4);
                    //sum += (Theta_k - Theta_kq) * Omega * I;
                    sum += Omega * I;
                    //sum += 1 / (q.norm() * q.norm() + 2*q*temp) * Omega;
                    //sum += Omega;
                    //sum += temp(0) * d3x;
                }
            }
        }
    }
    assert(not isnan(sum));
    return (float)sum;
}

