#ifndef MODEL_H
#define MODEL_H

#include <string>
#include <Eigen/Eigen>

typedef Eigen::Vector<double, 6> Vector6d;

class Model {
    public:
    
        // Setters.
        void set_name(std::string name);
        void set_nparams(int nparams);
        void set_stress(Vector6d stress);
        
        // Getters.
        std::string get_name(void);
        int get_nparams(void);
        Vector6d get_stress(void);
        
    private:
        std::string _name;
        Vector6d _stress;
        int _nparams;
};

#endif