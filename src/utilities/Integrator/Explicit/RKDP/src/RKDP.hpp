#ifndef RKDP_H
#define RKDP_H

#include <plog/Log.h>
#include "Types.hpp"
#include "Explicit.hpp"


class RKDP : public Explicit {

    public:

        /**
         * @brief RKDP integrator constructor.
         * 
         * @param settings Settings.
         * @param mf Model functions.
         */
        RKDP(Settings *settings, ModelFunctions *mf);

        /**
         * @brief RKDP integrator destructor.
         */
        ~RKDP() {}

    private:

        /**
         * @brief Overriden method to compute an initial estimate for the current increment.
         * 
         */
        void compute_initial_estimate(void) override;

};

#endif
