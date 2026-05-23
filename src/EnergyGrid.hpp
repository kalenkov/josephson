// Copyright (c) 2026 Mikhail S. Kalenkov <mikhail.kalenkov@gmail.com>
// Licensed under the MIT License
#ifndef ENERGY_GRID_H
#define ENERGY_GRID_H

#include <eigen3/Eigen/Core>

class EnergyGrid {
public:
    EnergyGrid(double energy_cutoff, double voltage, double step_size);

    const double& grid_max_energy() const;
    const double& energy_cutoff() const;
    const double& step_size() const;
    const double& voltage() const;
    const int& steps_per_voltage() const;
    const int& steps_per_half_voltage() const;
    const int& num_voltage_steps() const;
    const int& size() const;
    const double& operator[](int index) const;
    double at(int index) const;

private:
    double grid_max_energy_;
    double energy_cutoff_;
    double voltage_;
    double step_size_;
    int steps_per_voltage_;
    int steps_per_half_voltage_;
    int num_voltage_steps_;
    int size_;
    Eigen::VectorXd grid_;
};

#endif // ENERGY_GRID_H
