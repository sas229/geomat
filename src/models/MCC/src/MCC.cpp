#include "MCC.hpp"

MCC::MCC() {
    set_name("MCC");
    set_nparams(6);
    set_nstatev(2);
    PLOG_DEBUG << _name << " model instantiated with " << _nparams << " parameters and " << _nstatev << " state variables.";
}