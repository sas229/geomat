#include "umat.hpp"

/** 
 * @brief Boolean indicating whether the current increment is the first call to the umat function. 
 */
bool first_call = true;

/** 
 * @brief Unique pointer to store the memory address of the chosen model. 
 */
std::unique_ptr<Model> model;

void umat(
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
    Eigen::Map<Voigt> map_to_stress(stress);
    Eigen::Map<Voigt> map_to_strain_increment(dstran);
    Eigen::Map<State> map_to_state(statev, *nstatv, 1);
    Eigen::Map<Parameters> map_to_parameters(props, *nprops, 1);

    // Create a native Eigen type using the map as initialisation.
    Voigt Eigen_sigma = map_to_stress;
    Voigt Eigen_dstran = map_to_strain_increment;
    State state = map_to_state;
    Parameters parameters = map_to_parameters;

    // On first call, initialise logger and instantiate model.
    if (first_call) {
        // Initialise logger.
        char log_filename[] = "umat.log";
        std::remove(log_filename);
        plog::init(plog::debug, log_filename);

        // Instantiate model.
        PLOG_INFO << *ndi << "D problem defined with " << *ntens << " stress variables.";
        PLOG_INFO << "Attempting to instantiate " << cmname << " model.";
        if (strcmp(cmname, "LinearElastic") == 0) {
            model.reset(new LinearElastic(parameters, state, "verbose"));    
        } else if (strcmp(cmname, "MCC") == 0) {
            model.reset(new MCC(parameters, state, "verbose"));   
        } else if (strcmp(cmname, "SMCC") == 0) {
            model.reset(new SMCC(parameters, state, "verbose"));    
        } else {
            PLOG_FATAL << "Model name given not implemented. Check name given in CAE / input file.";
            assert(true);
        }
        first_call = false;
    } 

    // Set variables within model.    
    model->set_sigma_prime_tilde(-Eigen_sigma);
    model->set_Delta_epsilon_tilde(-Eigen_dstran);

    // model->update_stress_invariants();
    // model->update_principal_stresses();
    // double c, s, s_bar;
    // model->compute_lode(850.0, 9000.0, c, s, s_bar);
    // model->update_lode();

    // Do some work with it... (i.e. stress integration).
    // int i = 0;
    // while (i<1000) {
    //     model->solve();
    //     model->set_Delta_epsilon_tilde(Eigen_dstran);
    //     i += 1;
    // }
    model->solve();

    // // Equate map to updated variable in order to map back to input variable.
    // std::cout << "Stress prior to update:\n" << map_to_stress << "\n";
    // map_to_stress = model->get_sigma_prime();
    // std::cout << "Updated stress after remapping:\n" << map_to_stress << "\n";

    // // The map to the jacobian below somehow overwrites the stress state in the LinearElastic model class!!!
    // Jacobian jacobian = model->get_jacobian();
    // ddsdde = jacobian.data();
    
    // statev = state.data();
    // std::cout << "Jacobian:\n" << jacobian << "\n";
    // std::cout << "Updated stress after sign change:\n" << map_to_stress << "\n";
    // std::cout << "From C array:\n";
    // for (auto i=0; i<6; i++) {
    //     std::cout << stress[i] << "\n";
    // }
    // if (*nstatv > 0) {
    //     std::cout << "statev[0]: " << statev[0] << "\n";
    // }   
}