// Copyright (c) 2026 Mikhail S. Kalenkov <mikhail.kalenkov@gmail.com>
// Licensed under the MIT License
#ifndef MISC_H
#define MISC_H

#include "AlphaCoefficients.hpp"
#include "Config.hpp"
#include "EnergyGrid.hpp"
#include "ScatteringMatrix.hpp"
#include "RuntimeInfo.hpp"

void print_config(std::ostream& out, const Config& config, const RuntimeInfo& rt_info, const double& phi_min, const double& phi_max, const int& n_phi, const std::string& output_filename, const std::string& cmd_line);

void printJphi(std::ostream& out, const std::vector<double>& phi, const std::vector<double>& J);

std::vector<double> evaluateJphi(const AlphaCoefficients& alpha, const EnergyGrid& E, const ScatteringMatrix& S, const AndreevAmplitude& a, const Config& config, std::vector<double> phi);

#endif // MISC_H
