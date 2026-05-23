// Copyright (c) 2026 Mikhail S. Kalenkov <mikhail.kalenkov@gmail.com>
// Licensed under the MIT License
#ifndef CONFIG_H
#define CONFIG_H

class Config {
public:
    Config(); 
    double R_up() const { return R_up_; }
    void R_up(double v) { R_up_ = v; }
    double R_down() const { return R_down_; }
    void R_down(double v) { R_down_ = v; }
    double spin_mixing_angle() const { return spin_mixing_angle_; }
    void spin_mixing_angle(double v) { spin_mixing_angle_ = v; }
    double energy_cutoff() const { return energy_cutoff_; }
    void energy_cutoff(double v) { energy_cutoff_ = v; }
    double step_size() const { return step_size_; }
    void step_size(double v) { step_size_ = v; }
    double eV() const { return eV_; }
    void eV(double v) { eV_ = v; }
    double T() const { return T_; }
    void T(double v) { T_ = v; }
    double energy_broadening() const { return energy_broadening_; }
    void energy_broadening(double v) { energy_broadening_ = v; }
    double order_parameter_cutoff_smoothing() const { return order_parameter_cutoff_smoothing_; }
    void order_parameter_cutoff_smoothing(double v) { order_parameter_cutoff_smoothing_ = v; }
    bool show_progress() const { return show_progress_; }
    void show_progress(bool v) { show_progress_ = v; }

private:
    double R_up_;
    double R_down_;
    double spin_mixing_angle_;
    double energy_cutoff_; 
    double step_size_;
    double eV_;
    double T_;
    double energy_broadening_;
    double order_parameter_cutoff_smoothing_;
    bool show_progress_;
};

#endif // CONFIG_H
