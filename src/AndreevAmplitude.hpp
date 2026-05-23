// Copyright (c) 2026 Mikhail S. Kalenkov <mikhail.kalenkov@gmail.com>
// Licensed under the MIT License
#ifndef ANDREEV_AMPLITUDE_H
#define ANDREEV_AMPLITUDE_H

#include "EnergyGrid.hpp"
#include <complex>
#include <eigen3/Eigen/Core>

double sigmoid(double x);

double sigmoid_derivative(double x);

std::complex<double> ae(std::complex<double> e, double Delta);

class AndreevAmplitude {
public:
    AndreevAmplitude(const EnergyGrid& E, double energy_decay, double order_parameter, double order_parameter_cutoff_smoothing);

    const double& energy_decay() const;
    const double& order_parameter() const;
    double order_parameter(int index) const;
    double order_parameter_at(int index) const;
    const double& order_parameter_cutoff_smoothing() const;
    const std::complex<double>& operator[](int index) const;
    std::complex<double> at(int index) const;

private:
    double energy_decay_;
    double order_parameter_;
    double order_parameter_cutoff_smoothing_;
    Eigen::VectorXcd a_;
    Eigen::VectorXd order_parameter_energy_grid_;
};

#endif // ANDREEV_AMPLITUDE_H
