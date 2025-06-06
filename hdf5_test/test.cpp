#include "H5Cpp.h"

int main() {
    const std::string filename = "example.h5";
    const hsize_t dims[1] = {5};
    int data[5] = {10, 20, 30, 40, 50};

    try {
        H5::H5File file(filename, H5F_ACC_TRUNC);

        // Create dataspace and dataset
        H5::DataSpace dataspace(1, dims);
        H5::DataSet dataset = file.createDataSet("my_dataset", H5::PredType::NATIVE_INT, dataspace);

        // Write the array to the dataset
        dataset.write(data, H5::PredType::NATIVE_INT);

    } catch (H5::Exception& error) {
        error.printErrorStack();
        return -1;
    }

    return 0;
}

