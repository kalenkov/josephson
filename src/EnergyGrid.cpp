// Copyright (c) 2026 Mikhail S. Kalenkov <mikhail.kalenkov@gmail.com>
// Licensed under the MIT License
#include "EnergyGrid.hpp"
#include <algorithm>
#include <cmath>

EnergyGrid::EnergyGrid(double energy_cutoff, double voltage, double step_size)
    : energy_cutoff_(energy_cutoff)
    , voltage_(voltage)
    , step_size_(step_size)
{
    if (energy_cutoff_ <= 0) {
        throw std::invalid_argument("EnergyGrid: energy_cutoff must be positive.");
    }
    if (voltage_ <= 0) {
        throw std::invalid_argument("EnergyGrid: voltage must be positive.");
    }
    if (step_size_ <= 0) {
        throw std::invalid_argument("EnergyGrid: step_size must be positive.");
    }
    steps_per_half_voltage_ = std::max(1, static_cast<int>(std::ceil(0.5 * voltage_ / step_size_)));
    steps_per_voltage_ = 2 * steps_per_half_voltage_;
    step_size_ = voltage_ / steps_per_voltage_;

    num_voltage_steps_ = 2 * (1 + static_cast<int>(std::ceil(energy_cutoff_ / voltage_)));
    grid_max_energy_ = num_voltage_steps_ * voltage_ / 2.0;
    size_ = num_voltage_steps_ * steps_per_voltage_ + 1;
    grid_ = Eigen::VectorXd::LinSpaced(size_, -grid_max_energy_, grid_max_energy_);
}

const double& EnergyGrid::grid_max_energy() const { return grid_max_energy_; }

const double& EnergyGrid::energy_cutoff() const { return energy_cutoff_; }

const double& EnergyGrid::step_size() const { return step_size_; }

const double& EnergyGrid::voltage() const { return voltage_; }

const int& EnergyGrid::steps_per_voltage() const { return steps_per_voltage_; }

const int& EnergyGrid::steps_per_half_voltage() const { return steps_per_half_voltage_; }

const int& EnergyGrid::num_voltage_steps() const { return num_voltage_steps_; }

const int& EnergyGrid::size() const { return size_; }

const double& EnergyGrid::operator[](int index) const
{
    return grid_(index);
}

double EnergyGrid::at(int index) const
{
    return step_size_ * index - grid_max_energy_;
}
