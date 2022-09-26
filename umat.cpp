#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Eigen>
#include <vector>
#include <memory>
#include <typeinfo>

class Model {
    public:
        Model(std::string name) {
            _name = name;
        }
        std::string get_name(void) {return _name;}
        virtual double get_lambda(void) {return _lambda;}
    private:
        std::string _name;
        double _lambda;
};

class MCC : public Model {
    public: 
       MCC(std::string name = "MCC") : Model(name) {}
       double get_lambda(void) {return 2.541;} 
    private:
        int _nparams = 6;
};

class SMCC : public Model {
    public: 
       SMCC(std::string name = "SMCC") : Model(name) {}
       double get_lambda(void) {return 0.442;} 
    private:
        int _nparams = 9;
};

// Global booleans to control model instantiation.
bool first_call = true;
bool model_instantiated = false;

// Global unique pointer to store a pointer to the chosen model.
std::unique_ptr<Model> model;

// Debug output to file.
std::ofstream debug("umat.log");

void instantiate_model(char* cmname) {
    // Set first_call boolean to false so that model instantiation is not repeated on subsequent calls to the umat function.
    first_call = false;

    // Compare model name supplied with available model names. Instantiate model if found, otherwise report error.
    if (strcmp(cmname, "MCC") == 0) {
        model.reset(new MCC);   
    } else if (strcmp(cmname, "SMCC") == 0) {
        model.reset(new SMCC);    
    } else {
        debug << "Error: Model name given not implemented. Check name given in input file.\n";
        return;
    }
    model_instantiated = true;
    debug << cmname << " model instantiated.\n";
}

/** Abaqus umat interface.
 * 
 * This function allows Abaqus to call the library of user-defined materials defined within this library. 
 * Fortran passes variables by reference by default, hence pointers are passed for all variables. In the 
 * stress and strain arrays and in the matrices `ddsdde', `ddsddt', and `drplde', direct components are 
 * stored first, followed by shear components. There are `ndi' direct and `nshr' engineering shear components.
 * At a minimum the function must define `ddsdde', `stress', `statev'. 
 * 
 * @brief Interface to allow Abaqus to call umat material models via a C++ implementation.
 * @param stress Stress tensor at the beginning of the increment to be updated.
 * @param statev An array of solution dependent state variables to be updated.
 * @param ddsdde Jacobian matrix to be updated.
 * @param sse Specific elastic strain energy.
 * @param spd Specific plastic dissipation.
 * @param scd Specific creep dissipation.
 * @param rpl Volumetric heat generation per unit time.
 * @param ddsddt Variation of the stress increments with respect to the temperature.
 * @param drplde Variation of rpl with respect to the strain increments.
 * @param drpldt Variation of rpl with respect to the temperature.
 * @param stran The total strains at the beginning of the increment.
 * @param dstran Array of strain increments.
 * @param time Value of step time (1) and total time (2) at the beginning of the increment.
 * @param dtime Time increment.
 * @param temp Temperature at the start of the increment.
 * @param predef Array of predefined field variables at this point at the start of the increment, interpolated from nodes.
 * @param dpred Array of increments of predefined field variables.
 * @param cmname User-defined material name (avoid ABQ_ prefix).
 * @param ndi Number of direct stress components at this point.
 * @param nshr Number of engineering shear stress components at this point.
 * @param ntens Size of the stress or strain component array (`ndi' + `nshr').
 * @param nstatv Number of solution-dependent state variables that are associated with this material type.
 * @param props User-specified array of material constants associated with this user material.
 * @param nprops Number of user-defined material constants associated with this user material.
 * @param coords An array containing the coordinates of this point.
 * @param drot Rotation increment matrix.
 * @param celent Characteristic element length, which is a typical length of a line across an element for a first-order element; it is half of the same typical length for a second-order element.
 * @param dfgrd0 Array containing the deformation gradient at the beginning of the increment.
 * @param dfgrd1 Array containing the deformation gradient at the end of the increment.
 * @param noel Element number.
 * @param npt Integration point number.
 * @param layer Layer number.
 * @param kspt Section point number within the current layer.
 * @param kstep Procedure type key; nlgeom boolean; linear perturbation boolean.
 * @param kinc Increment number.
 * */
extern "C" void umat(
    double *stress,
    double *statev,
    double *ddsdde,
    double *sse,
    double *spd,
    double *scd,
    double *rpl,
    double *ddsddt,
    double *drplde,
    double *drpldt,
    double *stran,
    double *dstran,
    double *time,
    double *dtime,
    double *temp,
    double *dtemp,
    double *predef,
    double *dpred,
    char *cmname,
    int *ndi,
    int *nshr,
    int *ntens,
    int *nstatv,
    double *props,
    int *nprops,
    double *coords,
    double *drot,
    double *pnewdt,
    double *celent,
    double *dfgrd0,
    double *dfgrd1,
    int *noel,
    int *npt,
    int *layer,
    int *kspt,
    int *kstep,
    int *kinc
    )
{   
    if (first_call) {
        // Instantiate model if first call to umat function.
        instantiate_model(cmname);

        // Set parameters.

        // Set state variables.

        // Test mapping from C arrays to Eigen structures.
        typedef Eigen::Vector<double, 6> Vector6d;

        // Map to Eigen structure.
        Eigen::Map<Vector6d> map(stress);
        std::cout << map.transpose() << "\n";
        
        // Modify Eigen structure.
        map(2) = 100.0;
        std::cout << map.transpose() << "\n";

        // Check this is reflected in the C array.
        for(int i = 0; i < 6; ++i) std::cout << stress[i] << " ";
        std::cout << "\n";
        if (stress[2] == 100.0) {
            std::cout << "Mapping performed successfully.\n";
        } else {
            std::cout << "Mapping failed\n";
        }

        // Make Eigen matrix and initialise with Map.
        Vector6d s(map);

        // Manipulate Eigen matrix;
        s(3) = 200.0;
        
        // Update Eigen Map.
        map = s;
        std::cout << map.transpose() << "\n";
        
        // Verify C array has been modified.
        for(int i = 0; i < 6; ++i) std::cout << stress[i] << " ";
        std::cout << "\n";
        if (stress[3] == 200.0) {
            std::cout << "Mapping performed successfully.\n";
        } else {
            std::cout << "Mapping failed\n";
        }
        
    } 
    if (model_instantiated) {
        // Perform stress integration.
        std::string name;
        name = model->get_name();
        debug << "Model name: " << name << "\n";

        double lambda;
        lambda = model->get_lambda();
        debug << "Lambda = " << lambda << "\n";

        // Do stress integration with existing model.
        stress[0] = 99;
        debug << "ddsdde: " << ddsdde[0];
        for (int i = 1; i < *ntens; ++i) {
            debug  << ", " << ddsdde[i];
        }
        debug << "\n";   
    }
    
}

