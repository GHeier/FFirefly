#include <complex>
#include <vector>

#include "src/config/load/cpp_config.hpp"
#include "src/objects/CMField/fields.hpp"
#include "src/objects/vec.hpp"
#include "vertex.hpp"

using namespace std;


//void call_flex() {
//    string filename = outdir + prefix + "_chi." + filetype;
//    printf("Reading chi from %s\n", filename.c_str());
//    Field_C chi(filename);
//    float U = U0;
//    float nx = q_mesh[0], ny = q_mesh[1], nz = q_mesh[2];
//    int chidim = chi.cmf.data.dimension;
//    if (chidim == 2) nz = 1;
//
//    vector<Vec> points;
//    vector<float> wpts = chi.cmf.data.w_points;
//    if (wpts.size() == 0) {
//        wpts.push_back(0.0);
//    }
//    printv("wpts size: %d\n", wpts.size());
//    vector<complex<Vec>> values;
//    vector<vector<vector<float>>> vec_values(1);
//
//    printf("Computing vertex\n");
//    for (int i = 0; i < nx; i++) {
//        for (int j = 0; j < ny; j++) {
//            for (int k = 0; k < nz; k++) {
//                Vec q = brillouin_zone * Vec(i / nx - 0.5, j / ny - 0.5, k / nz - 0.5);
//                q.dimension = chidim;
//                for (int l = 0; l < wpts.size(); l++) {
//                    float w = wpts[l];
//                    if (chidim == 3) q.w = w;
//                    else q.z = w;
//                    points.push_back(q);
//                    complex<float> X = chi(q, w);
//                    complex<float> val = (U * U * X) / complex<float>(1.0f - U * X) + (U * U * U * X * X) / complex<float>(1.0f - U * U * X * X);
//                    if (filetype == "dat" || filetype == "txt")
//                        values.push_back(complex<Vec>(Vec(val), Vec(val.imag())));
//                    else if (filetype == "h5" || filetype == "hdf5")
//                        vec_values[0].push_back({val, val.imag()});
//                    if (abs(U * X) >= 1) {
//                        printf("Geometric series not convergent: U*X = %f\n", U * X);
//                        exit(1);
//                    }
//                }
//            }
//        }
//    }
//    printf("Saving Vertex\n");
//    string file = outdir + prefix + "_vertex." + filetype;
//    //if (filetype == "dat" || filetype == "txt")
//    //    save_to_file(file, points, values, chi.cmf.data.dimension, chi.cmf.data.with_w, chi.cmf.data.with_n, chi.cmf.data.is_complex, chi.cmf.data.is_vector);
//    if (filetype == "hdf5" || filetype == "h5") {
//        vector<vector<float>> BZ = brillouin_zone;
//        BZ.resize(dimension);
//        for (auto &row : BZ)
//            row.resize(dimension);
//        Vec first = BZ * Vec(-0.5, -0.5, -0.5);
//        vector<int> mesh = q_mesh;
//        if (chi.cmf.data.dimension == 2)
//            mesh = {q_mesh[0], q_mesh[1]};
//        //save_to_field(file, vec_values, BZ, mesh, chi.cmf.data.w_points, chi.cmf.data.is_complex, chi.cmf.data.is_vector);
//    }
//    cout << "Saved to " << outdir + prefix + "_vertex." + filetype << endl;
//}

void call_flex() {
    string filename = outdir + prefix + "_chi." + filetype;
    printf("Reading chi from %s\n", filename.c_str());

    Field_R chi(filename);
    float U = U0;
    int chidim = chi.get_data()->dimension;

    vector<float> wpts = chi.get_data()->w_points;
    if (wpts.size() == 0) {
        wpts.push_back(0.0);
    }
    vector<float> vals;
    vector<float> singlet_vals;

    if (!chi.get_data()->as_mesh) {
        printf("Using stored point data\n");
        // Loop over stored points
        for (const auto& point : chi.get_data()->points) {
            Vec q(point.data(), chidim);
            q.dimension = chidim;
            for (size_t l = 0; l < wpts.size(); l++) {
                float w = wpts[l];

                float X = chi(q, w);
                //float val = (U * U * X) / float(1.0f - U * X) + (U * U * U * X * X) / float(1.0f - U * U * X * X);
                float val = 1.5 * (U * U * X) / float(1.0f - U * X) + 0.5 * U * U * X / (1 + U * X) - U * U * X;
                vals.push_back(val);

                float singlet_val = 1.5 * (U * U * X) / float(1.0f - U * X) - 0.5 * U * U * X / (1 + U * X);
                singlet_vals.push_back(singlet_val);

                if (abs(U * X) >= 1) {
                    printf("Geometric series not convergent: U*X = %f\n", U * X);
                    exit(1);
                }
            }
        }
    } else {
        printf("Using mesh data\n");
        // Loop over mesh
        float nx = q_mesh[0], ny = q_mesh[1], nz = q_mesh[2];
        if (chidim == 2) nz = 1;

        // w loop outermost to generate data in w-k order (matching CMF_search expectations)
        for (size_t l = 0; l < wpts.size(); l++) {
            float w = wpts[l];
            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < ny; j++) {
                    for (int k = 0; k < nz; k++) {
                        Vec q = brillouin_zone * Vec(i / nx - 0.5, j / ny - 0.5, k / nz - 0.5);
                        q.dimension = chidim;

                        float X = chi(q, w);
                        //float val = (U * U * X) / float(1.0f - U * X) + (U * U * U * X * X) / float(1.0f - U * U * X * X);
                        float val = 1.5 * (U * U * X) / float(1.0f - U * X) + 0.5 * U * U * X / (1 + U * X) - U * U * X;
                        vals.push_back(val);

                        float singlet_val = 1.5 * (U * U * X) / float(1.0f - U * X) - 0.5 * U * U * X / (1 + U * X);
                        singlet_vals.push_back(singlet_val);

                        if (abs(U * X) >= 1.0) {
                            printf("Geometric series not convergent: U*X = %f\n", U * X);
                            exit(1);
                        }
                    }
                }
            }
        }
    }

    float max_val = 0.0f;
    for (const auto& v : vals) {
        if (abs(v) > max_val) {
            max_val = abs(v);
        }
    }
    printf("Computed %zu vertex values\n", vals.size());
    printf("Max vertex magnitude: %f\n", max_val);

    printf("Saving Vertex\n");
    string file = outdir + prefix + "_vertex." + filetype;
    string singlet_file = outdir + prefix + "_vertex_singlet." + filetype;
    if (filetype == "hdf5" || filetype == "h5") {
        if (chi.get_data()->as_mesh) {
            vector<int> mesh_vec(q_mesh.begin(), q_mesh.begin() + chidim);
            save_data(file, vals, vector<int>{}, mesh_vec, brillouin_zone, wpts);
            save_data(singlet_file, singlet_vals, vector<int>{}, mesh_vec, brillouin_zone, wpts);
        } else {
            save_data(file, vals, vector<int>{}, vector<int>{}, vector<vector<float>>{{}}, wpts, chi.get_data()->points);
            save_data(singlet_file, singlet_vals, vector<int>{}, vector<int>{}, vector<vector<float>>{{}}, wpts, chi.get_data()->points);
        }
    }
    cout << "Saved to " << file << endl;
    cout << "Saved to " << singlet_file << endl;
}
