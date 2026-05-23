// Copyright (c) 2026 Mikhail S. Kalenkov <mikhail.kalenkov@gmail.com>
// Licensed under the MIT License
#ifndef JOSEPHSON_CURRENT_HPP
#define JOSEPHSON_CURRENT_HPP

#include "AlphaCoefficients.hpp"
#include "AndreevAmplitude.hpp"
#include "Config.hpp"
#include "EnergyGrid.hpp"
#include "ScatteringMatrix.hpp"
#include "RuntimeInfo.hpp"
#include <vector>

class JosephsonCurrent
{
public:
    JosephsonCurrent(Config config);
    double evaluate(double phi) const;
    std::vector<double> evaluate(const std::vector<double>& phi) const;
    RuntimeInfo runtime_info() const;

private:
    Config config_;
    EnergyGrid E;
    ScatteringMatrix S;
    AndreevAmplitude a;
    AlphaCoefficients alpha;
};

#endif // JOSEPHSON_CURRENT_HPP