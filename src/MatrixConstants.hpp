// Copyright (c) 2026 Mikhail S. Kalenkov <mikhail.kalenkov@gmail.com>
// Licensed under the MIT License
#ifndef MATRIX_CONSTANTS_H
#define MATRIX_CONSTANTS_H

#include <complex>
#include <eigen3/Eigen/Dense>

namespace MatrixConstants {

inline const Eigen::Matrix<std::complex<double>, 8, 8> I = (Eigen::Matrix<std::complex<double>, 8, 8>() << // clang-format off
    1,  0,  0, 0, 0,  0,  0, 0,
    0,  1,  0, 0, 0,  0,  0, 0,
    0,  0,  1, 0, 0,  0,  0, 0,
    0,  0,  0, 1, 0,  0,  0, 0,
    0,  0,  0, 0, 1,  0,  0, 0,
    0,  0,  0, 0, 0,  1,  0, 0,
    0,  0,  0, 0, 0,  0,  1, 0,
    0,  0,  0, 0, 0,  0,  0, 1)
    .finished();
// clang-format on

inline const Eigen::Matrix<std::complex<double>, 8, 8> O3T3 = (Eigen::Matrix<std::complex<double>, 8, 8>() << // clang-format off
    1,  0,  0, 0, 0,  0,  0, 0,
    0, -1,  0, 0, 0,  0,  0, 0,
    0,  0, -1, 0, 0,  0,  0, 0,
    0,  0,  0, 1, 0,  0,  0, 0,
    0,  0,  0, 0, 1,  0,  0, 0,
    0,  0,  0, 0, 0, -1,  0, 0,
    0,  0,  0, 0, 0,  0, -1, 0,
    0,  0,  0, 0, 0,  0,  0, 1)
    .finished();
// clang-format on

inline const Eigen::Matrix<std::complex<double>, 8, 8> Pp = (Eigen::Matrix<std::complex<double>, 8, 8>() << // clang-format off
    1,  0,  0, 0, 0,  0,  0, 0,
    0,  0,  0, 0, 0,  0,  0, 0,
    0,  0,  0, 0, 0,  0,  0, 0,
    0,  0,  0, 1, 0,  0,  0, 0,
    0,  0,  0, 0, 1,  0,  0, 0,
    0,  0,  0, 0, 0,  0,  0, 0,
    0,  0,  0, 0, 0,  0,  0, 0,
    0,  0,  0, 0, 0,  0,  0, 1)
    .finished();
// clang-format on

inline const Eigen::Matrix<std::complex<double>, 8, 8> Pm = (Eigen::Matrix<std::complex<double>, 8, 8>() << // clang-format off
    0,  0,  0, 0, 0,  0,  0, 0,
    0,  1,  0, 0, 0,  0,  0, 0,
    0,  0,  1, 0, 0,  0,  0, 0,
    0,  0,  0, 0, 0,  0,  0, 0,
    0,  0,  0, 0, 0,  0,  0, 0,
    0,  0,  0, 0, 0,  1,  0, 0,
    0,  0,  0, 0, 0,  0,  1, 0,
    0,  0,  0, 0, 0,  0,  0, 0)
    .finished();
// clang-format on

inline const Eigen::Matrix<std::complex<double>, 8, 8> O1 = (Eigen::Matrix<std::complex<double>, 8, 8>() << // clang-format off
    0, 1, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 1, 0)
    .finished();
// clang-format on

inline const Eigen::Matrix<std::complex<double>, 8, 8> S1 = (Eigen::Matrix<std::complex<double>, 8, 8>() << // clang-format off
    0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 1,
    1, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0)
    .finished();
// clang-format on

inline const Eigen::Matrix<std::complex<double>, 8, 8> S3 = (Eigen::Matrix<std::complex<double>, 8, 8>() << // clang-format off
    1, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, -1, 0, 0, 0,
    0, 0, 0, 0, 0, -1, 0, 0,
    0, 0, 0, 0, 0, 0, -1, 0,
    0, 0, 0, 0, 0, 0, 0, -1)
    .finished();
// clang-format on

inline const Eigen::Matrix<std::complex<double>, 8, 8> T1 = (Eigen::Matrix<std::complex<double>, 8, 8>() << // clang-format off
    0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0)
    .finished();
// clang-format on

inline const Eigen::Matrix<std::complex<double>, 8, 8> A = (Eigen::Matrix<std::complex<double>, 8, 8>() << // clang-format off
    1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 1, 0, 0)
    .finished();
// clang-format on

}

#endif // MATRIX_CONSTANTS_H
