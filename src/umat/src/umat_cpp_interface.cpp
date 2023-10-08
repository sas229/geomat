#include "umat_cpp_interface.hpp"

/** 
 * @brief Boolean indicating whether the current increment is the first call to the umat function. 
 */
bool first_call = true;

/** 
 * @brief Unique pointer to store the memory address of the chosen model. 
 */
std::unique_ptr<Model> model;

void umat_cpp_interface(
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
    
    // Create maps to data.
    PLOG_INFO << "Creating maps to data...";
    Eigen::Map<Voigt> map_to_stress(stress);
    Eigen::Map<Voigt> map_to_strain_increment(dstran);
    Eigen::Map<State> map_to_state(statev, *nstatv, 1);
    Eigen::Map<Parameters> map_to_parameters(props, *nprops, 1);

    // Create a native Eigen type using the map as initialisation.
    Voigt Eigen_sigma = map_to_stress;
    Voigt Eigen_dstran = map_to_strain_increment;
    State state = map_to_state;
    Parameters parameters = map_to_parameters;

    // On first call, initialise logger.
    if (first_call) {
        // Initialise logger.
        char log_filename[] = "umat.log";
        std::remove(log_filename);
        plog::init(plog::debug, log_filename);
        first_call = false;
    }

    // Instantiate model.
    PLOG_INFO << *ndi << "D problem defined with " << *ntens << " stress variables.";
    PLOG_INFO << "Attempting to instantiate " << cmname << " model.";
    if (strcmp(cmname, "LinearElastic") == 0) {
        model.reset(new LinearElastic(parameters, state, "verbose"));    
    } else if (strcmp(cmname, "C2MC") == 0) {
        model.reset(new C2MC(parameters, state, "verbose"));  
    } else if (strcmp(cmname, "EMC") == 0) {
        model.reset(new EMC(parameters, state, "verbose"));   
    } else if (strcmp(cmname, "MCC") == 0) {
        model.reset(new MCC(parameters, state, "verbose"));   
    } else if (strcmp(cmname, "SMCC") == 0) {
        model.reset(new SMCC(parameters, state, "verbose"));    
    } else {
        PLOG_FATAL << "Model name given not implemented. Check name given in CAE / input file.";
        assert(true);
    }
        
    // Set initial stress state.
    model->set_sigma_prime_tilde(-Eigen_sigma);

    // Set variables within model.    
    model->set_Delta_epsilon_tilde(-Eigen_dstran);
    model->solve();

    // // Equate map to updated variable in order to map back to input variable.
    map_to_stress = -model->get_sigma_prime_tilde();
    map_to_state = model->get_state_variables();

    // // The map to the jacobian below somehow overwrites the stress state in the LinearElastic model class!!!
    // Jacobian jacobian = model->get_jacobian();
    // ddsdde = jacobian.data();
    
    // statev = state.data();
    // std::cout << "Jacobian:\n" << jacobian << "\n";
    // std::cout << "Updated stress after sign change:\n" << map_to_stress << "\n";
    // for (auto i=0; i<6; i++) {
    //     std::cout << stress[i] << " ";
    // }
    // std::cout << "\n";
    // if (*nstatv > 0) {
    //     std::cout << "statev[0]: " << statev[0] << "\n";
    // }   
}