#include "SMCC.hpp"

/**
 * @brief Soft Modified Cam Clay (SMCC) model constructor.
 * */
SMCC::SMCC() {
    set_name("SMCC");
    set_nparams(9);
    set_nstatev(3);
    PLOG_DEBUG << _name << " model instantiated with " << _nparams << " parameters and " << _nstatev << " state variables.";
}