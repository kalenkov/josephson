// Copyright (c) 2026 Mikhail S. Kalenkov <mikhail.kalenkov@gmail.com>
// Licensed under the MIT License
#include "Config.hpp"

Config::Config()
{
    R_up_ = 0.2;
    R_down_ = 0.5;
    spin_mixing_angle_ = 2.1;
    energy_cutoff_ = 3.0;
    step_size_ = 0.1;
    eV_ = 0.3;
    T_ = 0.1;
    energy_broadening_ = 0.1;
    order_parameter_cutoff_smoothing_ = 1.0;
    show_progress_ = false;
}



