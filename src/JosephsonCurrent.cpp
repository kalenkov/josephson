// Copyright (c) 2026 Mikhail S. Kalenkov <mikhail.kalenkov@gmail.com>
// Licensed under the MIT License
#include "JosephsonCurrent.hpp"
#include "AlphaCoefficients.hpp"
#include "AndreevAmplitude.hpp"
#include "Config.hpp"
#include "EnergyGrid.hpp"
#include "ScatteringMatrix.hpp"
#include "RuntimeInfo.hpp"
#include "misc.hpp"

JosephsonCurrent::JosephsonCurrent(Config config) : config_(config),
                                                              E(config_.energy_cutoff(), config_.eV(), config_.step_size()),
                                                              S(config_.R_up(), config_.R_down(), config_.spin_mixing_angle()),
                                                              a(E, config_.energy_broadening(), 1.0, config_.order_parameter_cutoff_smoothing()),
                                                              alpha(E, S, a)
{
}

double JosephsonCurrent::evaluate(double phi) const
{
    return evaluateJphi(alpha, E, S, a, config_, {phi})[0];
}

std::vector<double> JosephsonCurrent::evaluate(const std::vector<double> &phi) const
{
    return evaluateJphi(alpha, E, S, a, config_, phi);
}

RuntimeInfo JosephsonCurrent::runtime_info() const
{
    RuntimeInfo info;
    info.step_size = E.step_size();
    info.grid_max_energy = E.grid_max_energy();
    info.grid_size = E.size();
    return info;
}

