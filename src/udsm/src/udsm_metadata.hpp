#ifndef UDSM_METADATA_H
#define UDSM_METADATA_H

#ifdef _WIN32
  #define PLAXIS_CALL __stdcall
#else
  #define PLAXIS_CALL
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Return the number of available constitutive models in the DLL.
 * @param[out] count Pointer to an integer that receives the model count.
 */
void PLAXIS_CALL GetModelCount(int* count);

/** @brief Return the name of a model by its 1-based index.
 * @param[in] iModel Pointer to the 1-based model index.
 * @param[out] name Character buffer that receives the model name.
 */
void PLAXIS_CALL GetModelName(int* iModel, char* name);

/** @brief Return the number of parameters for a given model.
 * @param[in] iModel Pointer to the 1-based model index.
 * @param[out] count Pointer to an integer that receives the parameter count.
 */
void PLAXIS_CALL GetParamCount(int* iModel, int* count);

/** @brief Return the name of a parameter for a given model.
 * @param[in] iModel Pointer to the 1-based model index.
 * @param[in] iParam Pointer to the 1-based parameter index.
 * @param[out] name Character buffer that receives the parameter name.
 */
void PLAXIS_CALL GetParamName(int* iModel, int* iParam, char* name);

/** @brief Return the unit string of a parameter for a given model.
 * @param[in] iModel Pointer to the 1-based model index.
 * @param[in] iParam Pointer to the 1-based parameter index.
 * @param[out] units Character buffer that receives the unit string.
 */
void PLAXIS_CALL GetParamUnit(int* iModel, int* iParam, char* units);

/** @brief Return the number of state variables for a given model.
 * @param[in] iModel Pointer to the 1-based model index.
 * @param[out] count Pointer to an integer that receives the state variable count.
 */
void PLAXIS_CALL GetStateVarCount(int* iModel, int* count);

/** @brief Return the name of a state variable for a given model.
 * @param[in] iModel Pointer to the 1-based model index.
 * @param[in] iParam Pointer to the 1-based state variable index.
 * @param[out] name Character buffer that receives the state variable name.
 */
void PLAXIS_CALL GetStateVarName(int* iModel, int* iParam, char* name);

/** @brief Return the unit string of a state variable for a given model.
 * @param[in] iModel Pointer to the 1-based model index.
 * @param[in] iParam Pointer to the 1-based state variable index.
 * @param[out] units Character buffer that receives the unit string.
 */
void PLAXIS_CALL GetStateVarUnit(int* iModel, int* iParam, char* units);

#ifdef __cplusplus
}
#endif

#endif
