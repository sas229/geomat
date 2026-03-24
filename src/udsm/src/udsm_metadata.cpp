#include "udsm_metadata.hpp"
#include "ModelLookup.hpp"
#include "Model.hpp"
#include <cstring>
#include <string>

/** @brief Helper to write a string into a Plaxis shortstring buffer (256 bytes).
 *
 * The first byte is the length, followed by the characters.
 */
static void write_shortstring(char* dest, const std::string& src) {
    int len = static_cast<int>(src.size());
    if (len > 255) len = 255;
    dest[0] = static_cast<char>(len);
    std::memcpy(dest + 1, src.c_str(), len);
}

void PLAXIS_CALL GetModelCount(int* count) {
    *count = get_model_count();
}

void PLAXIS_CALL GetModelName(int* iModel, char* name) {
    auto model = create_model_for_metadata(*iModel);
    write_shortstring(name, model->get_model_name());
}

void PLAXIS_CALL GetParamCount(int* iModel, int* count) {
    auto model = create_model_for_metadata(*iModel);
    *count = model->get_parameters_required();
}

void PLAXIS_CALL GetParamName(int* iModel, int* iParam, char* name) {
    auto model = create_model_for_metadata(*iModel);
    auto names = model->get_param_names();
    int idx = *iParam - 1;
    if (idx >= 0 && idx < static_cast<int>(names.size())) {
        write_shortstring(name, names[idx]);
    }
}

void PLAXIS_CALL GetParamUnit(int* iModel, int* iParam, char* units) {
    auto model = create_model_for_metadata(*iModel);
    auto unit_list = model->get_param_units();
    int idx = *iParam - 1;
    if (idx >= 0 && idx < static_cast<int>(unit_list.size())) {
        write_shortstring(units, unit_list[idx]);
    }
}

void PLAXIS_CALL GetStateVarCount(int* iModel, int* count) {
    auto model = create_model_for_metadata(*iModel);
    *count = model->get_state_required();
}

void PLAXIS_CALL GetStateVarName(int* iModel, int* iParam, char* name) {
    auto model = create_model_for_metadata(*iModel);
    auto names = model->get_state_var_names();
    int idx = *iParam - 1;
    if (idx >= 0 && idx < static_cast<int>(names.size())) {
        write_shortstring(name, names[idx]);
    }
}

void PLAXIS_CALL GetStateVarUnit(int* iModel, int* iParam, char* units) {
    auto model = create_model_for_metadata(*iModel);
    auto unit_list = model->get_state_var_units();
    int idx = *iParam - 1;
    if (idx >= 0 && idx < static_cast<int>(unit_list.size())) {
        write_shortstring(units, unit_list[idx]);
    }
}
