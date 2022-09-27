#ifndef UTILITIES_H
#define UTILITIES_H

#include <memory>
#include <string>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string.h>

#include "Model.hpp"
#include "MCC.hpp"
#include "SMCC.hpp"

namespace utilities {

    bool instantiate_model(char* cmname, std::unique_ptr<Model>& model, std::ofstream& debug);

}

#endif