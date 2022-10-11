#include "SMCC.hpp"

/**
 * @brief Soft Modified Cam Clay (SMCC) model constructor.
 * */
SMCC::SMCC() {
    set_name("SMCC");
    set_n_parameters(9);
    set_n_state_variables(3);
    PLOG_DEBUG << name << " model instantiated with " << n_parameters << " parameters and " << n_state_variables << " state variables.";
}