#include "umat.hpp"

/** @brief Boolean indicating whether the current increment is the first call to the umat function. */
bool first_call = true;

/** @brief Unique pointer to store the memory address of the chosen model. */
std::unique_ptr<Model> model;

/** @brief Eigen vector with six indices. */
typedef Eigen::Vector<double, 6> Vector6d;

/** @brief Eigen vector with six indices. */
typedef Eigen::Matrix<double, 3, 3> Matrix3d;

/** @brief Eigen vector with six indices. */
typedef Eigen::Matrix<double, 6, 6> Matrix6d;

extern "C" void umat(
    double *stress,
    double *statev,
    double *ddsdde,
    double *sse,
    double *spd,
    double *scd,
    double *rpl,
    double *ddsddt,
    double *drplde,
    double *drpldt,
    double *stran,
    double *dstran,
    double *time,
    double *dtime,
    double *temp,
    double *dtemp,
    double *predef,
    double *dpred,
    char *cmname,
    int *ndi,
    int *nshr,
    int *ntens,
    int *nstatv,
    double *props,
    int *nprops,
    double *coords,
    double *drot,
    double *pnewdt,
    double *celent,
    double *dfgrd0,
    double *dfgrd1,
    int *noel,
    int *npt,
    int *layer,
    int *kspt,
    int *kstep,
    int *kinc) {   
    // On first call, initialise logger and instantiate model.
    if (first_call) {
        // Initialise logger.
        char log_filename[] = "umat.log";
        std::remove(log_filename);
        plog::init(plog::debug, log_filename);

        // Instantiate model.
        PLOG_DEBUG << *ndi << "D problem defined with " << *ntens << " stress variables.";
        PLOG_DEBUG << "Attempting to instantiate " << cmname << " model.";
        if (strcmp(cmname, "MCC") == 0) {
            model.reset(new MCC);   
        } else if (strcmp(cmname, "SMCC") == 0) {
            model.reset(new SMCC);    
        } else {
            PLOG_FATAL << "Model name given not implemented. Check name given in CAE / input file.";
        }
        first_call = false;
    } 

    // Create maps to data.
    Eigen::Map<Vector6d> map_to_stress(stress);
    Eigen::Map<Matrix3d> map_to_jacobian(ddsdde);
    std::vector<double> state(statev, statev+*nstatv);
    std::vector<double> parameters(props, props+*nprops);

    // Create a native Eigen type using the map as initialisation.
    Vector6d Eigen_sigma = map_to_stress;
    Matrix3d Eigen_jacobian = map_to_jacobian;

    // Set variables within model.
    model->set_parameters(parameters);
    model->set_sigma(Eigen_sigma);
    model->set_jacobian(Eigen_jacobian);
    model->set_state_variables(state);
    model->compute_invariants();

    // Do some work with it... (i.e. stress integration).
    
    // Perform stress integration.

    // Equate map to updated variable in order to map back to input variable.
    map_to_stress = model->get_sigma();
    map_to_jacobian = model->get_jacobian();
    std::vector<double> new_state = model->get_state_variables();
    std::cout << new_state[0] << "\n";

    std::cout << model->get_elastic_matrix() << "\n";
    
    // Verify that original stress variable from Abaqus has been updated.
    for (int i=0; i < *ntens; ++i) {
        std::cout << stress[i] << " ";
    }
    std::cout << "\n"; 

    // Verify that original stress variable from Abaqus has been updated.
    for (int i=0; i < 9; ++i) {
        std::cout << ddsdde[i] << " ";
    }
    std::cout << "\n";  
}