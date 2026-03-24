#include "umat_cpp_interface.hpp"

/** @brief Thread-safe one-time logger initialisation. */
static std::once_flag log_init_flag;

static void initialise_umat_log() {
    static char log_filename[] = "umat.log";
    std::remove(log_filename);
    plog::init(plog::debug, log_filename);
}

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

    // Thread-safe one-time logger initialisation.
    std::call_once(log_init_flag, initialise_umat_log);

    // Create maps to data.
    PLOG_INFO << "Creating maps to data...";
    Eigen::Map<Voigt> map_to_stress(stress);
    Eigen::Map<Voigt> map_to_strain_increment(dstran);
    Eigen::Map<State> map_to_state(statev, *nstatv, 1);
    Eigen::Map<Parameters> map_to_parameters(props, *nprops, 1);

    // Create native Eigen types using the maps.
    Voigt Eigen_sigma = map_to_stress;
    Voigt Eigen_dstran = map_to_strain_increment;
    State state = map_to_state;
    Parameters parameters = map_to_parameters;

    // Instantiate model locally (thread-safe — no global state).
    PLOG_INFO << *ndi << "D problem defined with " << *ntens << " stress variables.";
    PLOG_INFO << "Attempting to instantiate " << cmname << " model.";
    std::unique_ptr<Model> model = create_model(std::string(cmname), parameters, state, "verbose");

    // Set initial stress state (negate for Abaqus tension-positive convention).
    model->set_sigma_prime_tilde(-Eigen_sigma);

    // Set strain increment and solve (negate for sign convention).
    model->set_Delta_epsilon_tilde(-Eigen_dstran);
    model->solve();

    // Map results back to Abaqus arrays (negate stress back to tension-positive).
    map_to_stress = -model->get_sigma_prime_tilde();
    map_to_state = model->get_state_variables();
}