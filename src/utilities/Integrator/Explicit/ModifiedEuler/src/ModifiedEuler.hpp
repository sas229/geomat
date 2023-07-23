#ifndef MODIFIEDEULER_H
#define MODIFIEDEULER_H

#include <plog/Log.h>
#include "Tensor.hpp"
#include "Types.hpp"
#include "Explicit.hpp"


class ModifiedEuler : public Explicit {

    public:

        /**
         * @brief ModifiedEuler integrator constructor.
         * 
         * @param settings Settings.
         * @param mf Model functions.
         */
        ModifiedEuler(Settings *settings, ModelFunctions *mf);

        /**
         * @brief ModifiedEuler integrator destructor.
         */
        ~ModifiedEuler() {}

    private:

        /**
         * @brief Overriden method to compute an initial estimate for the current increment.
         * 
         */
        void compute_initial_estimate(void) override;

};

#endif
