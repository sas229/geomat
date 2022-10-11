#include "MCC.hpp"

MCC::MCC() {
    set_name("MCC");
    set_n_parameters(6);
    set_n_state_variables(2);
    PLOG_DEBUG << name << " model instantiated with " << n_parameters << " parameters and " << n_state_variables << " state variables.";
}