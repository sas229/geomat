#include "umat.hpp"

// Booleans to control model instantiation.
bool first_call = true;
bool model_instantiated = true;

// Unique pointer to store a pointer to the chosen model.
std::unique_ptr<Model> model;

// Debug output to file.
std::ofstream debug("umat.log");

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
    // Instantiate model if first call to umat function.
    if (first_call) {
        if (strcmp(cmname, "MCC") == 0) {
            model.reset(new MCC);   
        } else if (strcmp(cmname, "SMCC") == 0) {
            model.reset(new SMCC);    
        } else {
            model_instantiated = false;
        }
        // Log debug message.
        if (model_instantiated) {
            debug << model->get_name() << " model instantiated with " << model->get_nparams() << " parameters.\n";
        } else {
            debug << "Error: Model name given not implemented. Check name given in input file.\n";
        }
        first_call = false;
    } 
    // If model instantiated perform stress integration.
    if (model_instantiated) {
        // Create map to data.
        Eigen::Map<Vector6d> _map_to_stress(stress);

        // Create a native Eigen type using the map as initialisation.
        Vector6d _stress{_map_to_stress};

        // Set variables within model.
        model->set_stress(_stress);

        // Do some work with it... (i.e. stress integration).
        
        // Perform stress integration.
        
        // Equate map to updated variable in order to map back to input variable.
        _map_to_stress = model->get_stress();

        // Verify that original stress variable from Abaqus has been updated.
        for (int i=0; i < *ntens; ++i) {
            std::cout << stress[i] << " ";
        }
        std::cout << "\n";  
    } 
}

