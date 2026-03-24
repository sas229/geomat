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

void PLAXIS_CALL GetModelCount(int* count);
void PLAXIS_CALL GetModelName(int* iModel, char* name);
void PLAXIS_CALL GetParamCount(int* iModel, int* count);
void PLAXIS_CALL GetParamName(int* iModel, int* iParam, char* name);
void PLAXIS_CALL GetParamUnit(int* iModel, int* iParam, char* units);
void PLAXIS_CALL GetStateVarCount(int* iModel, int* count);
void PLAXIS_CALL GetStateVarName(int* iModel, int* iParam, char* name);
void PLAXIS_CALL GetStateVarUnit(int* iModel, int* iParam, char* units);

#ifdef __cplusplus
}
#endif

#endif
