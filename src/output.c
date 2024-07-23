#include "output.h"
#define NPARAMS 10
#define STR_SIZE 50

void open_file(input_t* input){

    herr_t status;
    hid_t attribute_id;

    hsize_t dim_pos_dims[] = {input->pos_dims};
    hsize_t dim_pos_params[] = {input->pos_dims, NPARAMS};
    hsize_t dim_mom_dims[] = {input->mom_dims};
    hsize_t dim_mom_params[] = {input->mom_dims, NPARAMS};

    hsize_t str_dim[] = {STR_SIZE};
    hsize_t str_dim_pos_ndims[] = {input->pos_dims, STR_SIZE};
    hsize_t str_dim_mom_ndims[] = {input->mom_dims, STR_SIZE};
    hsize_t str_dim_pos_params[] = {input->pos_dims, NPARAMS, STR_SIZE};
    hsize_t str_dim_mom_params[] = {input->mom_dims, NPARAMS, STR_SIZE};
    hsize_t str_dim_fields[] = {input->n_fields, STR_SIZE};

    hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    hid_t file_id = H5Fcreate(input->filename, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);

    /*

    hid_t scalar_dataspace_id = H5Screate(H5S_SCALAR);
    hid_t pos_dims_dataspace_id = H5Screate_simple(1, dim_pos_dims, NULL);
    hid_t mom_dims_dataspace_id = H5Screate_simple(1, dim_mom_dims, NULL);
    hid_t pos_params_dataspace_id = H5Screate_simple(2, dim_pos_params, NULL);
    hid_t mom_params_dataspace_id = H5Screate_simple(2, dim_mom_params, NULL);

    hid_t string_scalar_dataspace_id = H5Screate_simple(1, str_dim, NULL);
    hid_t string_pos_dims_dataspace_id = H5Screate_simple(2, str_dim_pos_ndims, NULL);
    hid_t string_mom_dims_dataspace_id = H5Screate_simple(2, str_dim_mom_ndims, NULL);
    hid_t string_pos_params_dataspace_id = H5Screate_simple(2, str_dim_pos_params, NULL);
    hid_t string_mom_params_dataspace_id = H5Screate_simple(2, str_dim_mom_params, NULL);
    hid_t string_fields_dataspace_id = H5Screate_simple(2, str_dim_fields, NULL);
    
    attribute_id = H5Acreate (file_id, "Time Delta", H5T_NATIVE_REAL, scalar_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attribute_id, H5T_NATIVE_REAL, &(input->deltaT));

    attribute_id = H5Acreate (file_id, "Number of Position Dimensions", H5T_NATIVE_INT, scalar_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attribute_id, H5T_NATIVE_INT, &(input->pos_dims));

    attribute_id = H5Acreate (file_id, "Number of Momentum Dimensions", H5T_NATIVE_INT, scalar_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attribute_id, H5T_NATIVE_INT, &(input->mom_dims));

    attribute_id = H5Acreate (file_id, "Number of Timesteps", H5T_NATIVE_INT, scalar_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attribute_id, H5T_NATIVE_INT, &(input->n_timesteps));

    attribute_id = H5Acreate (file_id, "Diagnostic Frequency", H5T_NATIVE_INT, scalar_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attribute_id, H5T_NATIVE_INT, &(input->diag_freq));

    attribute_id = H5Acreate (file_id, "Number of Fields", H5T_NATIVE_INT, scalar_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attribute_id, H5T_NATIVE_INT, &(input->n_fields));

    attribute_id = H5Acreate (file_id, "Number of Species", H5T_NATIVE_INT, scalar_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attribute_id, H5T_NATIVE_INT, &(input->n_species));

    attribute_id = H5Acreate (file_id, "Position Min.", H5T_NATIVE_REAL, pos_dims_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attribute_id, H5T_NATIVE_REAL, input->pos_min);

    attribute_id = H5Acreate (file_id, "Momentum Min.", H5T_NATIVE_REAL, mom_dims_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attribute_id, H5T_NATIVE_REAL, input->mom_min);

    attribute_id = H5Acreate (file_id, "Position Max.", H5T_NATIVE_REAL, pos_dims_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attribute_id, H5T_NATIVE_REAL, input->pos_max);

    attribute_id = H5Acreate (file_id, "Momentum Max.", H5T_NATIVE_REAL, mom_dims_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attribute_id, H5T_NATIVE_REAL, input->mom_max);

    attribute_id = H5Acreate (file_id, "Position Delta", H5T_NATIVE_REAL, pos_dims_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attribute_id, H5T_NATIVE_REAL, input->pos_delta);

    attribute_id = H5Acreate (file_id, "Momentum Delta", H5T_NATIVE_REAL, mom_dims_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attribute_id, H5T_NATIVE_REAL, input->mom_delta);

    attribute_id = H5Acreate (file_id, "Number of Position Points", H5T_NATIVE_INT, pos_dims_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attribute_id, H5T_NATIVE_INT, input->pos_points);

    attribute_id = H5Acreate (file_id, "Number of Momentum Points", H5T_NATIVE_INT, mom_dims_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attribute_id, H5T_NATIVE_INT, input->mom_points);

    attribute_id = H5Acreate (file_id, "Position Initial Condition Parameters", H5T_NATIVE_REAL, pos_params_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attribute_id, H5T_NATIVE_REAL, input->pos_init_params[0]);

    attribute_id = H5Acreate (file_id, "Momentum Initial Condition Parameters", H5T_NATIVE_REAL, pos_params_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attribute_id, H5T_NATIVE_REAL, input->mom_init_params[0]);

    attribute_id = H5Acreate (file_id, "Dispersion Relation", H5T_NATIVE_CHAR, string_scalar_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attribute_id, H5T_NATIVE_CHAR, input->dispersion_names[0]);

    attribute_id = H5Acreate (file_id, "Kernels", H5T_NATIVE_CHAR, string_scalar_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attribute_id, H5T_NATIVE_CHAR, input->kernel_names[0]);

    attribute_id = H5Acreate (file_id, "Time-advancing algorithm", H5T_NATIVE_CHAR, string_scalar_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attribute_id, H5T_NATIVE_CHAR, input->pusher);

    attribute_id = H5Acreate (file_id, "Space-discretization algorithm", H5T_NATIVE_CHAR, string_scalar_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attribute_id, H5T_NATIVE_CHAR, input->operator);

    attribute_id = H5Acreate (file_id, "Initial Condition Functions", H5T_NATIVE_CHAR, string_fields_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attribute_id, H5T_NATIVE_CHAR, input->pos_init_names[0]);
    
    status = H5Aclose(attribute_id);
    status = H5Sclose(scalar_dataspace_id);
    status = H5Sclose(pos_dims_dataspace_id);
    status = H5Sclose(mom_dims_dataspace_id);
    status = H5Sclose(pos_params_dataspace_id);
    status = H5Sclose(mom_params_dataspace_id);

    status = H5Sclose(string_scalar_dataspace_id);
    status = H5Sclose(string_pos_dims_dataspace_id);
    status = H5Sclose(string_mom_dims_dataspace_id);
    status = H5Sclose(string_pos_params_dataspace_id);
    status = H5Sclose(string_mom_params_dataspace_id);
    status = H5Sclose(string_fields_dataspace_id);

    */

    status = H5Fclose(file_id);
}

void write_solution(input_t* input, REAL** results, int timestep){

    hid_t file_id, dataspace_id, dataset_id;
    herr_t status;

    file_id = H5Fopen(input->filename, H5F_ACC_RDWR, H5P_DEFAULT);
    
    char* name = malloc(50);
    sprintf(name, "/Species%d", timestep);

    hsize_t dims[] = {input->n_species, input->pos_total, input->mom_total};
    dataspace_id = H5Screate_simple(3, dims, NULL);
    dataset_id = H5Dcreate2(file_id, name, H5T_NATIVE_REAL, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_REAL, H5S_ALL, H5S_ALL, H5P_DEFAULT, results[0]);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    status = H5Fclose(file_id);
    free(name);
}

void write_fields(input_t* input, COMPLEX** fields, int timestep){

    hid_t complex_id = H5Tcreate (H5T_COMPOUND, sizeof (hdf5_complex_t));
    H5Tinsert (complex_id, "real", HOFFSET(hdf5_complex_t, re), H5T_NATIVE_REAL);
    H5Tinsert (complex_id, "imaginary", HOFFSET(hdf5_complex_t, im), H5T_NATIVE_REAL);

    hid_t file_id, dataspace_id, dataset_id;
    herr_t status;

    file_id = H5Fopen(input->filename, H5F_ACC_RDWR, H5P_DEFAULT);
    
    char* name = malloc(50);
    sprintf(name, "/Fields%d", timestep);

    hsize_t dims[] = {input->n_fields, input->pos_total};
    dataspace_id = H5Screate_simple(2, dims, NULL);
    dataset_id = H5Dcreate2(file_id, name, complex_id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, complex_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, fields[0]);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    status = H5Fclose(file_id);
    free(name);
}


void write_sources(input_t* input, COMPLEX** sources, COMPLEX** aux, int* indices, int timestep){

    MPI_Barrier(MPI_COMM_WORLD);

    // Eliminate ghost cells
    for(int k = 0; k < input->n_species; ++k){
        boundary_shift(input, sources[k], aux[k], indices);
    }

    // Create complex IO datatype
    hid_t complex_id = H5Tcreate (H5T_COMPOUND, sizeof (hdf5_complex_t));
    H5Tinsert (complex_id, "real", HOFFSET(hdf5_complex_t, re), H5T_NATIVE_REAL);
    H5Tinsert (complex_id, "imaginary", HOFFSET(hdf5_complex_t, im), H5T_NATIVE_REAL);

    // Declare ID and error variables
    hid_t file_id, dataspace_id, dataset_id, fapl_id, xf_id;
    herr_t status;

    // Open File
    fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    file_id = H5Fopen(input->filename, H5F_ACC_RDWR, fapl_id);
    
    /*// Create parallel access to file and dataset
    char* name = malloc(50);
    sprintf(name, "/Sources%d", timestep);
    hsize_t dims[] = {input->n_species, input->pos_total};
    dataspace_id = H5Screate_simple(2, dims, NULL);
    dataset_id = H5Dcreate2(file_id, name, complex_id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    xf_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(xf_id, H5FD_MPIO_INDEPENDENT);

    // Select parallel subset
    hsize_t offset[2];
    hsize_t stride[2];
    hsize_t count[2];
    hsize_t block[2];
    offset[0] = 0;
    stride[0] = 1;
    count[0] = 1;
    block[0] = input->n_species;
    offset[1] = input->rank * (input->pos_points[input->pos_dims-1]-2*input->padding);
    stride[1] = (input->size-1)*offset[1];
    count[1] = 1;
    block[1] = (input->pos_points[input->pos_dims-1]-2*input->padding);
    for(int i = 0; i < input->pos_dims-1; ++i){
        count[1] *= (input->pos_points[i]-2*input->padding);
    }
    
    hid_t memspace_id = H5Screate_simple (2, block, NULL);
    status = H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset, stride, count, block);

    // Write to file
    status = H5Dwrite(dataset_id, complex_id, memspace_id, dataspace_id, xf_id, aux[0]);
    
    // Close all spaces
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
    status = H5Sclose(memspace_id);
    status = H5Pclose(fapl_id);
    status = H5Pclose(xf_id);*/
    status = H5Fclose(file_id);

    for(int k = 0; k < input->n_species; ++k){
        inverse_boundary_shift(input, aux[k], sources[k], indices);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    
    //free(name);
}