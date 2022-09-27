#include "utilities.hpp"

bool utilities::instantiate_model(char* cmname, std::unique_ptr<Model>& model, std::ofstream& debug) {
    // Compare model name supplied with available model names. Instantiate model if found, otherwise report error.
    if (strcmp(cmname, "MCC") == 0) {
        model.reset(new MCC);   
    } else if (strcmp(cmname, "SMCC") == 0) {
        model.reset(new SMCC);    
    } else {
        debug << "Error: Model name given not implemented. Check name given in input file.\n";
        return false;
    }
    debug << model->get_name() << " model instantiated with " << model->get_nparams() << " parameters.\n";
    return true;
}
