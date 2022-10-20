#include "umat.hpp"

/** 
 * @brief Boolean indicating whether the current increment is the first call to the umat function. 
 */
bool first_call = true;

/** 
 * @brief Unique pointer to store the memory address of the chosen model. 
 */
std::unique_ptr<Model> model;

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
    // Create maps to data.
    Eigen::Map<Vector6d> map_to_stress(stress);
    Eigen::Map<Vector6d> map_to_strain_increment(dstran);
    Eigen::Map<Eigen::Matrix3d> map_to_jacobian(ddsdde);
    std::vector<double> state(statev, statev+*nstatv);
    std::vector<double> parameters(props, props+*nprops);

    // Create a native Eigen type using the map as initialisation.
    Voigt Eigen_sigma = map_to_stress;
    Voigt Eigen_dstran = map_to_strain_increment;
    Eigen::Matrix3d Eigen_jacobian = map_to_jacobian;

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
            model.reset(new LinearElastic(parameters, state));    
        } else if (strcmp(cmname, "MCC") == 0) {
            model.reset(new MCC(parameters, state));   
        } else {
            PLOG_FATAL << "Model name given not implemented. Check name given in CAE / input file.";
            assert(true);
        }
        first_call = false;
    } 

    // Set variables within model.    
    model->set_sigma_prime(Eigen_sigma);
    model->set_strain_increment(Eigen_dstran);
    model->set_jacobian(Eigen_jacobian);

    model->update_stress_invariants();
    model->update_principal_stresses();
    double c, s, s_bar;
    model->compute_lode(850.0, 9000.0, c, s, s_bar);
    model->update_lode();

    // Do some work with it... (i.e. stress integration).
    model->solve();

    // Equate map to updated variable in order to map back to input variable.
    map_to_stress = model->get_sigma_prime();
    map_to_jacobian = model->get_jacobian();
    statev = state.data();
    std::cout << "Updated stress after sign change:\n" << map_to_stress << "\n";
    std::cout << "From C array:\n";
    for (auto i=0; i<6; i++) {
        std::cout << stress[i] << "\n";
    }
    if (*nstatv > 0) {
        std::cout << "statev[0]: " << statev[0] << "\n";
    }   
}