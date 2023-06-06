#include "Types.hpp"

Voigt to_voigt(Cauchy cauchy) {
    Voigt voigt = Voigt::Zero();
    voigt(0) = cauchy(0,0);
    voigt(1) = cauchy(1,1);
    voigt(2) = cauchy(2,2);
    voigt(3) = cauchy(0,1);
    voigt(4) = cauchy(0,2);
    voigt(5) = cauchy(1,2);
    return voigt;
}

Cauchy to_cauchy(Voigt voigt) {
    Cauchy cauchy = Cauchy::Zero();
    cauchy(0,0) = voigt(0);
    cauchy(1,1) = voigt(1);
    cauchy(2,2) = voigt(2);
    cauchy(0,1) = cauchy(1,0) = voigt(3);
    cauchy(0,2) = cauchy(2,0) = voigt(4);
    cauchy(1,2) = cauchy(2,1) = voigt(5);
    return cauchy;
}