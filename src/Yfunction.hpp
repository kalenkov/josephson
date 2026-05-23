// Copyright (c) 2026 Mikhail S. Kalenkov <mikhail.kalenkov@gmail.com>
// Licensed under the MIT License
#ifndef Y_FUNCTION_H
#define Y_FUNCTION_H

#include "AlphaCoefficients.hpp"
#include "AndreevAmplitude.hpp"
#include "EnergyGrid.hpp"
#include "MatrixConstants.hpp"
#include "ScatteringMatrix.hpp"
#include <complex>
#include <eigen3/Eigen/Core>
#include <format>
#include <iostream>
#include <tuple>
#include <vector>

using std::complex_literals::operator""i;

class Yfunction {
public:
    Yfunction(const AlphaCoefficients& alpha, const EnergyGrid& E, const ScatteringMatrix& S, const AndreevAmplitude& a, double phi)
        : Y_(E.size())
        , Q_(E.size())
        , E_(E)
        , S_(S)
        , a_(a)
        , alpha_(alpha)
        , phi_(phi)
    {
    }

    void evaluateY(std::atomic<std::size_t>& progress)
    {
        int local_progress;
        int skip = 1000;
        std::vector<Eigen::Matrix<std::complex<double>, 8, 8>> Z_p(E_.size(), MatrixConstants::I);
        std::vector<Eigen::Matrix<std::complex<double>, 8, 8>> Z_m(E_.size(), MatrixConstants::I);
        for (int i = 0; i < E_.size(); i++) {
            Y_[i] = alpha_.z(i);
            local_progress++;
            if (local_progress % skip == 0) {
                local_progress = 0;
                progress++;
            }
        }

        for (int n = 1; n < E_.num_voltage_steps(); n++) {
            for (int i = n * E_.steps_per_voltage(); i < E_.size(); i++) {
                Z_p[i] = alpha_.p(i - (2 * n - 1) * E_.steps_per_half_voltage()) * Z_p[i] * std::exp(1.0i * phi_ / 2.0);
                Y_[i] = Y_[i] + alpha_.z(i - n * E_.steps_per_voltage()) * Z_p[i];
                local_progress++;
                if (local_progress % skip == 0) {
                    local_progress = 0;
                    progress++;
                }
            }
        }

        for (int n = 1; n < E_.num_voltage_steps(); n++) {
            for (int i = E_.size() - 1 - n * E_.steps_per_voltage(); i >= 0; i--) {
                Z_m[i] = alpha_.m(i + (2 * n - 1) * E_.steps_per_half_voltage()) * Z_m[i] * std::exp(-1.0i * phi_ / 2.0);
                Y_[i] = Y_[i] + alpha_.z(i + n * E_.steps_per_voltage()) * Z_m[i];
                local_progress++;
                if (local_progress % skip == 0) {
                    local_progress = 0;
                    progress++;
                }
            }
        }
    }

    void evaluateQ(double Temperature)
    {
        for (int i = E_.steps_per_half_voltage(); i < E_.size() - E_.steps_per_half_voltage(); i++) {
            Q_[i] = MatrixConstants::Pp * (1.0 - std::abs(a_[i + E_.steps_per_half_voltage()] * a_[i + E_.steps_per_half_voltage()])) * std::tanh(E_[i + E_.steps_per_half_voltage()] / (2.0 * Temperature)) + MatrixConstants::Pm * (1.0 - std::abs(a_[i - E_.steps_per_half_voltage()] * a_[i - E_.steps_per_half_voltage()])) * std::tanh(E_[i - E_.steps_per_half_voltage()] / (2.0 * Temperature));
        }
    }

    double ElectricCurrent()
    {
        double current = 0;

        for (int i = E_.steps_per_half_voltage(); i < E_.size() - E_.steps_per_half_voltage(); i++) {
            current += (1.0/8.0)*(Y_[i] * Q_[i] * Y_[i].adjoint() * (S_.W().adjoint() * MatrixConstants::O3T3 * S_.W() + MatrixConstants::O3T3)).trace().real() * E_.step_size() / (2.0 * M_PI);
        }
        return current;
    }

    const Eigen::Matrix<std::complex<double>, 8, 8>& operator[](size_t index) const
    {
        return Y_[index];
    }
    const Eigen::Matrix<std::complex<double>, 8, 8>& at(int index) const
    {
        if (index < 0 or index >= Y_.size()) {
            return Y_[0];
        } else {
            return Y_[index];
        }
    }

private:
    std::vector<Eigen::Matrix<std::complex<double>, 8, 8>> Y_;
    std::vector<Eigen::Matrix<std::complex<double>, 8, 8>> Q_;
    const EnergyGrid& E_;
    const ScatteringMatrix& S_;
    const AndreevAmplitude& a_;
    const AlphaCoefficients& alpha_;
    double phi_;
};

#endif // Y_FUNCTION_H
