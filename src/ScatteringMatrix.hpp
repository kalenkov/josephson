// Copyright (c) 2026 Mikhail S. Kalenkov <mikhail.kalenkov@gmail.com>
// Licensed under the MIT License
#ifndef SCATTERING_MATRIX_H
#define SCATTERING_MATRIX_H

#include <eigen3/Eigen/Core>

class ScatteringMatrix {
public:
    ScatteringMatrix(double reflection_probability_up, double reflection_probability_down, double spin_mixing_angle);

    const double& reflection_probability_up() const;
    const double& reflection_probability_down() const;
    const double& transmission_probability_up() const;
    const double& transmission_probability_down() const;
    const double& reflection_amplitude_up() const;
    const double& reflection_amplitude_down() const;
    const double& transmission_amplitude_up() const;
    const double& transmission_amplitude_down() const;
    const double& spin_mixing_angle() const;
    const Eigen::Matrix<std::complex<double>, 8, 8>& Sn() const;
    const Eigen::Matrix<std::complex<double>, 8, 8>& W() const;
    const Eigen::Matrix<std::complex<double>, 8, 8>& SnO3P3Sn_O3P3() const;

private:
    void calculate_matrices();

    Eigen::Matrix<std::complex<double>, 8, 8> Sn_;
    Eigen::Matrix<std::complex<double>, 8, 8> W_;
    Eigen::Matrix<std::complex<double>, 8, 8> SnO3P3Sn_O3P3_;
    double reflection_probability_up_;
    double reflection_probability_down_;
    double transmission_probability_up_;
    double transmission_probability_down_;
    double reflection_amplitude_up_;
    double reflection_amplitude_down_;
    double transmission_amplitude_up_;
    double transmission_amplitude_down_;
    double spin_mixing_angle_;
};

#endif // SCATTERING_MATRIX_H
