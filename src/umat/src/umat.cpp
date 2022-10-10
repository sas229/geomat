#include "umat.hpp"

/** @brief Boolean indicating whether the current increment is the first call to the umat function. */
bool first_call = true;

/** @brief Unique pointer to store the memory address of the chosen model. */
std::unique_ptr<Model> model;

/** @brief Eigen vector with six indices. */
typedef Eigen::Vector<double, 6> Tensor3D;

/** @brief Eigen vector with three indices. */
typedef Eigen::Vector<double, 3> Tensor2D;

Eigen::Map<Tensor2D> _map_to_2D_stress_tensor(NULL);
Eigen::Map<Tensor3D> _map_to_3D_stress_tensor(NULL);

Tensor2D _stress_2D;
Tensor3D _stress_3D;

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
    if (*ntens < 6) {
        // 2D problem.

        // Create maps to pointers.
        new (&_map_to_2D_stress_tensor) Eigen::Map<Tensor2D>(stress); 

        // Create a native Eigen type using the map as initialisation.
        _stress_2D = _map_to_2D_stress_tensor; 

        // Set variables within model.
        model->set_stress(_stress_2D);

        // Equate map to updated variable in order to map back to input variable.
        _map_to_2D_stress_tensor = model->get_stress();
    } else {
        // 3D problem.

        // Create maps to pointers.
        new (&_map_to_3D_stress_tensor) Eigen::Map<Tensor3D>(stress);  

        // Create a native Eigen type using the map as initialisation.
        _stress_3D = _map_to_3D_stress_tensor;

        // Set variables within model.
        model->set_stress(_stress_3D);

        // Equate map to updated variable in order to map back to input variable.
        _map_to_3D_stress_tensor = model->get_stress();
    }

    // Do some work with it... (i.e. stress integration).
    
    // Perform stress integration.
    
    // Verify that original stress variable from Abaqus has been updated.
    for (int i=0; i < *ntens; ++i) {
        std::cout << stress[i] << " ";
    }
    std::cout << "\n";  
}