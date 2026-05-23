// Copyright (c) 2026 Mikhail S. Kalenkov <mikhail.kalenkov@gmail.com>
// Licensed under the MIT License
#include "AndreevAmplitude.hpp"
#include "EnergyGrid.hpp"

using namespace std::complex_literals;

double sigmoid(double x)
{
    return (3.0 * x * x - 2.0 * x * x * x);
}

double sigmoid_derivative(double x)
{
    return 6.0 * x * (1.0 - x);
}

std::complex<double> ae(std::complex<double> e, double Delta)
{
    if (std::abs(e) < Delta) {
        return (-e + std::sqrt(e - Delta) * std::sqrt(e + Delta)) / Delta;
    } else {
        return Delta / (-e - std::sqrt(e - Delta) * std::sqrt(e + Delta));
    }
}

AndreevAmplitude::AndreevAmplitude(const EnergyGrid& E, double energy_decay, double order_parameter, double order_parameter_cutoff_smoothing)
    : energy_decay_(energy_decay)
    , order_parameter_(order_parameter)
    , order_parameter_cutoff_smoothing_(order_parameter_cutoff_smoothing)
    , order_parameter_energy_grid_(Eigen::VectorXd::Constant(E.size(), order_parameter))
    , a_(E.size())
{
    if (energy_decay_ <= 0) {
        throw std::invalid_argument("EnergyGrid: energy_decay must be positive.");
    }
    if (order_parameter_ < 0) {
        throw std::invalid_argument("EnergyGrid: order_parameter must be positive.");
    }
    if (order_parameter_cutoff_smoothing_ <= 0) {
        throw std::invalid_argument("EnergyGrid: order_parameter_cutoff_smoothing must be positive.");
    }
    if (order_parameter_cutoff_smoothing_ >= E.energy_cutoff()) {
        throw std::invalid_argument("EnergyGrid: order_parameter_cutoff_smoothing must be less than energy_cutoff.");
    }
    size_t i = 0;
    while (E[i] < -E.energy_cutoff()) {
        order_parameter_energy_grid_[i] = 0;
        order_parameter_energy_grid_[E.size() - 1 - i] = 0;
        i++;
    }
    while (E[i] < -E.energy_cutoff() + order_parameter_cutoff_smoothing) {
        order_parameter_energy_grid_[i] = order_parameter_ * sigmoid((E[i] + E.energy_cutoff()) / order_parameter_cutoff_smoothing_);
        order_parameter_energy_grid_[E.size() - 1 - i] = order_parameter_energy_grid_[i];
        i++;
    }
    for (int i = 0; i < E.size(); i++) {
        a_[i] = ae(E[i] + 1i * energy_decay_, order_parameter_energy_grid_[i]);
    }
}

const double& AndreevAmplitude::energy_decay() const
{
    return energy_decay_;
}
const double& AndreevAmplitude::order_parameter() const
{
    return order_parameter_;
}

double AndreevAmplitude::order_parameter(int index) const
{
    return order_parameter_energy_grid_[index];
}
double AndreevAmplitude::order_parameter_at(int index) const
{
    if (index < 0 or index >= a_.size()) {
        return 0;
    } else {
        return order_parameter_energy_grid_[index];
    }
}
const double& AndreevAmplitude::order_parameter_cutoff_smoothing() const
{
    return order_parameter_cutoff_smoothing_;
}

const std::complex<double>& AndreevAmplitude::operator[](int index) const
{
    return a_[index];
}

std::complex<double> AndreevAmplitude::at(int index) const
{
    if (index < 0 or index >= a_.size()) {
        return 0;
    } else {
        return a_[index];
    }
}
