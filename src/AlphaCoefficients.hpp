// Copyright (c) 2026 Mikhail S. Kalenkov <mikhail.kalenkov@gmail.com>
// Licensed under the MIT License
#ifndef ALPHA_COEFFICIENTS_H
#define ALPHA_COEFFICIENTS_H

#include "AndreevAmplitude.hpp"
#include "EnergyGrid.hpp"
#include "MatrixConstants.hpp"
#include "ScatteringMatrix.hpp"
#include <eigen3/Eigen/Core>
#include <complex>
#include <vector>

using std::complex_literals::operator""i;

class AlphaCoefficients {
public:
    AlphaCoefficients(const EnergyGrid& E, const ScatteringMatrix& S, const AndreevAmplitude& a)
        : alpha_p_(E.size())
        , alpha_m_(E.size())
        , alpha_0_(E.size())
        , E_(E)
        , S_(S)
        , a_(a)
    {
        for (int i = 0; a.order_parameter(i + E.steps_per_half_voltage()) == 0; i++) {
            alpha_0_[i] = S.W().inverse();
            alpha_0_[E.size() - 1 - i] = alpha_0_[i];
        }
#pragma omp parallel sections
        {
#pragma omp section
            {
                for (int i = E_.size() - 1 - E_.steps_per_voltage(); i >= 0; i--)
                {
                    alpha_p_[i] = -a_[i] * MatrixConstants::Pp * (S_.W() + a_[i + E_.steps_per_voltage()] * alpha_p_[i + E_.steps_per_voltage()] * MatrixConstants::Pm).inverse();
                }
            }
#pragma omp section
            {
                for (int i = E_.steps_per_voltage(); i < E_.size(); i++)
                {
                    alpha_m_[i] = -a_[i] * MatrixConstants::Pm * (S_.W() + a_[i - E_.steps_per_voltage()] * alpha_m_[i - E_.steps_per_voltage()] * MatrixConstants::Pp).inverse();
                }
            }
        }

#pragma omp parallel for
        for (int i = E_.steps_per_half_voltage(); i < E_.size() - E_.steps_per_half_voltage(); i++) {
            alpha_0_[i] = (S_.W() + a_[i - E_.steps_per_half_voltage()] * alpha_m_[i - E_.steps_per_half_voltage()] * MatrixConstants::Pp + a_[i + E_.steps_per_half_voltage()] * alpha_p_[i + E_.steps_per_half_voltage()] * MatrixConstants::Pm).inverse();
        }
    }
    const Eigen::Matrix<std::complex<double>, 8, 8>& p(size_t index) const
    {
        return alpha_p_[index];
    }
    const Eigen::Matrix<std::complex<double>, 8, 8>& m(size_t index) const
    {
        return alpha_m_[index];
    }
    const Eigen::Matrix<std::complex<double>, 8, 8>& z(size_t index) const
    {
        return alpha_0_[index];
    }
    const Eigen::Matrix<std::complex<double>, 8, 8>& p_at(int index) const
    {
        if (index < 0 or index >= alpha_p_.size()) {
            return alpha_p_[0];
        } else {
            return alpha_p_[index];
        }
    }
    const Eigen::Matrix<std::complex<double>, 8, 8>& m_at(int index) const
    {
        if (index < 0 or index >= alpha_m_.size()) {
            return alpha_m_[0];
        } else {
            return alpha_m_[index];
        }
    }
    const Eigen::Matrix<std::complex<double>, 8, 8>& z_at(int index) const
    {
        if (index < 0 or index >= alpha_0_.size()) {
            return alpha_0_[0];
        } else {
            return alpha_0_[index];
        }
    }

protected:
    std::vector<Eigen::Matrix<std::complex<double>, 8, 8>> alpha_p_;
    std::vector<Eigen::Matrix<std::complex<double>, 8, 8>> alpha_m_;
    std::vector<Eigen::Matrix<std::complex<double>, 8, 8>> alpha_0_;
    const EnergyGrid& E_;
    const AndreevAmplitude& a_;
    const ScatteringMatrix& S_;
};

#endif // ALPHA_COEFFICIENTS_H
