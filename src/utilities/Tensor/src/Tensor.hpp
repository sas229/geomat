
#ifndef TENSOR_H
#define TENSOR_H

#include <Eigen/Core>
#include <iostream>
 
/** 
 * @brief Eigen vector with six indices. 
 */
typedef Eigen::Vector<double, 6> Vector6d;

/** 
 * @brief Eigen matrix with six indices. 
 */
typedef Eigen::Matrix<double, 6, 6> Matrix6d;

/** 
 * @brief Custom Tensor type derived from Eigen::Matrix3d type. 
 * Method added to return the tensor contents in Voigt notation form, 
 */
class Tensor : public Eigen::Matrix3d {

public:

    Tensor(void):Eigen::Matrix3d() {}
 
    // This constructor allows you to construct Tensor from Eigen expressions.
    template<typename OtherDerived>
    Tensor(const Eigen::MatrixBase<OtherDerived>& other) : Eigen::Matrix3d(other){}
 
    // This method allows you to assign Eigen expressions to Tensor.
    template<typename OtherDerived>
    Tensor& operator=(const Eigen::MatrixBase <OtherDerived>& other) {
        this->Eigen::Matrix3d::operator=(other);
        return *this;
    }

    /** @brief Return tensor in Voigt notation form. */
    Eigen::Vector<double, 6> voigt(void) {
        Eigen::Vector<double, 6> voigt;
        voigt(0) = this->operator()(0,0);
        voigt(1) = this->operator()(1,1);
        voigt(2) = this->operator()(2,2);
        voigt(3) = this->operator()(0,1);
        voigt(4) = this->operator()(0,2);
        voigt(5) = this->operator()(1,2);
        return voigt;    
    }

};

#endif