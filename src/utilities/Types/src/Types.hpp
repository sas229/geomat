
#ifndef TYPES_H
#define TYPES_H

#include <Eigen/Core>
#include <iostream>
 
/** 
 * @brief Eigen vector with six indices. 
 */
typedef Eigen::Vector<double, 6> Vector6d;

/** 
 * @brief Eigen square matrix with six indices for storing the constitutive matrix. 
 */
typedef Eigen::Matrix<double, 6, 6> Constitutive;

/** 
 * @brief Eigen square matrix with six indices for storing the Jacobian matrix. 
 */
typedef Eigen::Matrix<double, 6, 6> Jacobian;

/** 
 * @brief Custom Cauchy type derived from Eigen::Matrix3d type. 
 * Method added to return the tensor in Voigt notation form. 
 */
class Cauchy : public Eigen::Matrix3d {

public:

    Cauchy(void):Eigen::Matrix3d() {}
 
    // This constructor allows you to construct Cauchy from Eigen expressions.
    template<typename OtherDerived>
    Cauchy(const Eigen::MatrixBase<OtherDerived>& other) : Eigen::Matrix3d(other){}
 
    // This method allows you to assign Eigen expressions to Cauchy.
    template<typename OtherDerived>
    Cauchy& operator=(const Eigen::MatrixBase <OtherDerived>& other) {
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

/** 
 * @brief Custom Voigt type derived from Eigen::Vector<double, 6> type. 
 * Method added to return the tensor in Cauchy notation form. 
 */
class Voigt : public Vector6d {

public:

    Voigt(void):Vector6d() {}
 
    // This constructor allows you to construct Voigt from Eigen expressions.
    template<typename OtherDerived>
    Voigt(const Eigen::MatrixBase<OtherDerived>& other) : Vector6d(other){}
 
    // This method allows you to assign Eigen expressions to Voigt.
    template<typename OtherDerived>
    Voigt& operator=(const Eigen::MatrixBase <OtherDerived>& other) {
        this->Vector6d::operator=(other);
        return *this;
    }

    /** @brief Return tensor in Cauchy notation form. */
    Cauchy cauchy(void) {
        Cauchy cauchy = Cauchy::Zero();
        cauchy(0,0) = this->operator()(0);
        cauchy(1,1) = this->operator()(1);
        cauchy(2,2) = this->operator()(2);
        cauchy(0,1) = cauchy(1,0) = this->operator()(3);
        cauchy(0,2) = cauchy(2,0) = this->operator()(4);
        cauchy(1,2) = cauchy(2,1) = this->operator()(5);
        return cauchy;    
    }

};

#endif