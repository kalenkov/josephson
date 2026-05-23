// Copyright (c) 2026 Mikhail S. Kalenkov <mikhail.kalenkov@gmail.com>
// Licensed under the MIT License
#include "ScatteringMatrix.hpp"
#include "MatrixConstants.hpp"
#include <cmath>
#include <complex>

using namespace std::complex_literals;

ScatteringMatrix::ScatteringMatrix(double reflection_probability_up, double reflection_probability_down, double spin_mixing_angle)
    : reflection_probability_up_(reflection_probability_up)
    , reflection_probability_down_(reflection_probability_down)
    , transmission_probability_up_(1.0 - reflection_probability_up)
    , transmission_probability_down_(1.0 - reflection_probability_down)
    , reflection_amplitude_up_(std::sqrt(reflection_probability_up))
    , reflection_amplitude_down_(std::sqrt(reflection_probability_down))
    , transmission_amplitude_up_(std::sqrt(1.0 - reflection_probability_up))
    , transmission_amplitude_down_(std::sqrt(1.0 - reflection_probability_down))
    , spin_mixing_angle_(spin_mixing_angle)
{
    double r1 = reflection_amplitude_up_;
    double r2 = reflection_amplitude_down_;
    double d1 = transmission_amplitude_up_;
    double d2 = transmission_amplitude_down_;
    std::complex<double> ph = std::exp(-1i * spin_mixing_angle / 2.0);
    Sn_ = (Eigen::Matrix<std::complex<double>, 8, 8>() << // clang-format off
              r1 * ph, -1i * d1 * ph,            0,            0,             0,             0,            0,            0,
        -1i * d1 * ph,       r1 * ph,            0,            0,             0,             0,            0,            0,
        0            ,             0,      r2 * ph, 1i * d2 * ph,             0,             0,            0,            0,
        0            ,             0, 1i * d2 * ph,      r2 * ph,             0,             0,            0,            0,
        0            ,             0,            0,            0,       r2 / ph, -1i * d2 / ph,            0,            0,
        0            ,             0,            0,            0, -1i * d2 / ph,       r2 / ph,            0,            0,
        0            ,             0,            0,            0,             0,             0,      r1 / ph, 1i * d1 / ph,
        0            ,             0,            0,            0,             0,             0, 1i * d1 / ph,      r1 / ph)
    .finished();
    // clang-format on
    calculate_matrices();
}

const double& ScatteringMatrix::reflection_probability_up() const
{
    return reflection_probability_up_;
}

const double& ScatteringMatrix::reflection_probability_down() const
{
    return reflection_probability_down_;
}

const double& ScatteringMatrix::transmission_probability_up() const
{
    return transmission_probability_up_;
}

const double& ScatteringMatrix::transmission_probability_down() const
{
    return transmission_probability_down_;
}

const double& ScatteringMatrix::reflection_amplitude_up() const
{
    return reflection_amplitude_up_;
}

const double& ScatteringMatrix::reflection_amplitude_down() const
{
    return reflection_amplitude_down_;
}

const double& ScatteringMatrix::transmission_amplitude_up() const
{
    return transmission_amplitude_up_;
}

const double& ScatteringMatrix::transmission_amplitude_down() const
{
    return transmission_amplitude_down_;
}

const double& ScatteringMatrix::spin_mixing_angle() const
{
    return spin_mixing_angle_;
}

const Eigen::Matrix<std::complex<double>, 8, 8>& ScatteringMatrix::Sn() const
{
    return Sn_;
}

const Eigen::Matrix<std::complex<double>, 8, 8>& ScatteringMatrix::W() const
{
    return W_;
}

const Eigen::Matrix<std::complex<double>, 8, 8>& ScatteringMatrix::SnO3P3Sn_O3P3() const
{
    return SnO3P3Sn_O3P3_;
}

void ScatteringMatrix::calculate_matrices()
{
    W_ = Sn_ * MatrixConstants::T1;
    SnO3P3Sn_O3P3_ = Sn_.adjoint() * MatrixConstants::O3T3 * Sn_ - MatrixConstants::O3T3;
}
