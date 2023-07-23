#include "Tensor.hpp"

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

Cauchy dev(Cauchy cauchy) {
    Cauchy eye = Cauchy::Identity();
    return cauchy - 1.0/3.0*eye*cauchy.trace();
}

double double_dot_product(Cauchy cauchy) {
    return cauchy.cwisePow(2.0).sum();
}

double double_dot_product(Cauchy a, Cauchy b) {
    return a.cwiseProduct(b).sum();
}

double tr(Cauchy cauchy) {
    return cauchy.trace();
}