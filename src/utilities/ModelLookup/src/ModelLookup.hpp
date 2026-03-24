#ifndef MODELLOOKUP_H
#define MODELLOOKUP_H

#include <memory>
#include <string>
#include <vector>
#include "Types.hpp"

class Model;

/** @brief Create a model instance by name.
 *
 * @param name Model name (e.g. "C2MC", "MCC").
 * @param parameters Model parameters.
 * @param state Initial state variables.
 * @param log_severity Log severity level.
 * @return std::unique_ptr<Model>
 */
std::unique_ptr<Model> create_model(const std::string& name, Parameters parameters, State state, std::string log_severity="none");

/** @brief Get the number of available models.
 * @return int
 */
int get_model_count();

/** @brief Get the list of available model names.
 * @return std::vector<std::string>
 */
std::vector<std::string> get_model_names();

/** @brief Create a model instance for metadata queries (zero-filled parameters and state).
 *
 * @param name Model name.
 * @return std::unique_ptr<Model>
 */
std::unique_ptr<Model> create_model_for_metadata(const std::string& name);

/** @brief Create a model instance for metadata queries by index (1-based, for Plaxis iMod).
 *
 * @param index 1-based model index.
 * @return std::unique_ptr<Model>
 */
std::unique_ptr<Model> create_model_for_metadata(int index);

#endif
