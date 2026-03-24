#include "udsm_cpp_interface.hpp"

/** @brief Thread-safe one-time logger initialisation. */
static std::once_flag udsm_log_init_flag;

static void initialise_udsm_log() {
    static char log_filename[] = "udsm.log";
    std::remove(log_filename);
    plog::init(plog::debug, log_filename);
}

void udsm_cpp_interface(
    int *IDTask, int *iMod, int *IsUndr, int *iStep, int *iTer,
    int *iEl, int *Int, double *X, double *Y, double *Z,
    double *Time0, double *dTime, double *Props, double *Sig0,
    double *Swp0, double *StVar0, double *dEps, double *D,
    double *Bulk_W, double *Sig, double *Swp, double *StVar,
    int *ipl, int *nStat, int *NonSym, int *iStrsDep,
    int *iTimeDep, int *iTang, char *iPrjDir, int *iPrjLen,
    int *iAbort)
{
    // Thread-safe one-time logger initialisation.
    std::call_once(udsm_log_init_flag, initialise_udsm_log);
    *iAbort = 0;

    switch (*IDTask) {
        case 1: {
            // Initialise state variables.
            PLOG_INFO << "IDTask 1: Initialising state variables for model index " << *iMod << ".";
            auto model = create_model_for_metadata(*iMod);
            int n_state = model->get_state_required();
            for (int i = 0; i < n_state; i++) {
                StVar0[i] = 0.0;
            }
            break;
        }
        case 2: {
            // Calculate constitutive stresses.
            PLOG_INFO << "IDTask 2: Calculating constitutive stresses for model index " << *iMod << ".";
            auto meta = create_model_for_metadata(*iMod);
            int n_params = meta->get_parameters_required();
            int n_state = meta->get_state_required();

            Parameters parameters = Eigen::Map<Parameters>(Props, n_params, 1);
            State state = (n_state > 0) ? State(Eigen::Map<State>(StVar0, n_state, 1)) : State::Zero(0);

            auto model = create_model(get_model_names()[*iMod - 1], parameters, state, "none");

            // Set stress (negate: Plaxis compression-negative -> geomat compression-positive).
            Voigt sigma_voigt;
            for (int i = 0; i < 6; i++) sigma_voigt[i] = -Sig0[i];
            model->set_sigma_prime_tilde(sigma_voigt);

            // Set strain increment (negate for sign convention).
            Voigt dstran;
            for (int i = 0; i < 6; i++) dstran[i] = -dEps[i];
            model->set_Delta_epsilon_tilde(dstran);

            // Solve.
            model->solve();

            // Write results (negate stress back to Plaxis convention).
            Voigt result_sigma = model->get_sigma_prime_tilde();
            for (int i = 0; i < 6; i++) Sig[i] = -result_sigma[i];

            // Write state variables.
            if (n_state > 0) {
                State result_state = model->get_state_variables();
                for (int i = 0; i < n_state; i++) StVar[i] = result_state[i];
            }

            // Set plasticity indicator (0 = elastic).
            *ipl = 0;
            break;
        }
        case 3: {
            // Create effective material stiffness matrix.
            PLOG_INFO << "IDTask 3: Creating effective stiffness matrix for model index " << *iMod << ".";
            auto meta = create_model_for_metadata(*iMod);
            int n_params = meta->get_parameters_required();
            int n_state = meta->get_state_required();

            Parameters parameters = Eigen::Map<Parameters>(Props, n_params, 1);
            State state = (n_state > 0) ? State(Eigen::Map<State>(StVar0, n_state, 1)) : State::Zero(0);

            auto model = create_model(get_model_names()[*iMod - 1], parameters, state, "none");

            // Set stress state for stress-dependent stiffness.
            Voigt sigma_voigt;
            for (int i = 0; i < 6; i++) sigma_voigt[i] = -Sig0[i];
            model->set_sigma_prime_tilde(sigma_voigt);

            // Compute stiffness matrix via a zero-strain solve to populate the Jacobian.
            Voigt zero_strain = Voigt::Zero();
            model->set_Delta_epsilon_tilde(zero_strain);
            model->solve();
            Constitutive D_e = model->get_jacobian();

            // Write to D array (column-major, Fortran order).
            for (int j = 0; j < 6; j++) {
                for (int i = 0; i < 6; i++) {
                    D[j * 6 + i] = D_e(i, j);
                }
            }
            break;
        }
        case 4: {
            // Return number of state variables.
            PLOG_INFO << "IDTask 4: Returning number of state variables for model index " << *iMod << ".";
            auto model = create_model_for_metadata(*iMod);
            *nStat = model->get_state_required();
            break;
        }
        case 5: {
            // Return matrix attributes.
            PLOG_INFO << "IDTask 5: Returning matrix attributes.";
            *NonSym = 0;
            *iStrsDep = 1;
            *iTimeDep = 0;
            *iTang = 0;
            break;
        }
        case 6: {
            // Create elastic material stiffness matrix.
            PLOG_INFO << "IDTask 6: Creating elastic stiffness matrix for model index " << *iMod << ".";
            auto meta = create_model_for_metadata(*iMod);
            int n_params = meta->get_parameters_required();
            int n_state = meta->get_state_required();

            Parameters parameters = Eigen::Map<Parameters>(Props, n_params, 1);
            State state = (n_state > 0) ? State(Eigen::Map<State>(StVar0, n_state, 1)) : State::Zero(0);

            auto model = create_model(get_model_names()[*iMod - 1], parameters, state, "none");

            // Set stress state for stress-dependent stiffness.
            Voigt sigma_voigt;
            for (int i = 0; i < 6; i++) sigma_voigt[i] = -Sig0[i];
            model->set_sigma_prime_tilde(sigma_voigt);

            // Compute elastic stiffness matrix via a zero-strain solve.
            Voigt zero_strain = Voigt::Zero();
            model->set_Delta_epsilon_tilde(zero_strain);
            model->solve();
            Constitutive D_e = model->get_jacobian();

            // Write to D array (column-major, Fortran order).
            for (int j = 0; j < 6; j++) {
                for (int i = 0; i < 6; i++) {
                    D[j * 6 + i] = D_e(i, j);
                }
            }
            break;
        }
        default: {
            PLOG_ERROR << "Unknown IDTask: " << *IDTask << ".";
            *iAbort = 1;
            break;
        }
    }
}
