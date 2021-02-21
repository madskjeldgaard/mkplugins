#include "lpg.hpp"
#include "math.h"
#include "mkutils.hpp"
#include <iostream>

/*
 * The code for the LowpassGate emulation here is based on the work of Julian
 * Parker and Stefano D'angelo and their paper "A Digital Model of the Buchla
 * Lowpass-Gate" from DAFX-13, 2013 and the accompanying Max code.
 *
 * See the paper here for more information:
 * https://www.dafx.de/paper-archive/2013/papers/44.dafx2013_submission_56.pdf
 */

namespace lpg {

// -------------------------------------------
// Low pass gate emulation
// -------------------------------------------
void LowPassGate::setup(float samplerate, float initOffset, float initScale,
                        float initVCA, float initResonance,
                        LpgFilter::Linearity initLinearity) {
  m_samplerate = samplerate;

  // Setup vactrol
  vactrol.setup(samplerate);

  // Setup control circuit
  controlcircuit.setup(initOffset, initScale);

  // Setup filter
  float initResistance = 2500.;
  filter.setup(samplerate, initResistance, initVCA, initResonance,
               initLinearity);
}

inline void LowPassGate::setControlOffset(float offset) {
  controlcircuit.setOffset(offset);
}

inline void LowPassGate::setControlScale(float scale) {
  controlcircuit.setScale(scale);
}

inline void LowPassGate::setVCAAmount(float amount) {
  filter.setVCAAmount(amount);
}

inline void LowPassGate::setResonance(double resonance) {
  filter.setResonance(resonance);
}

inline void LowPassGate::setLowpassMode(LpgFilter::LowpassProcessing mode) {
  filter.setLowpassMode(mode);
}

inline void LowPassGate::setLinearity(LpgFilter::Linearity linearity) {
  filter.setProcessingMode(linearity);
}

inline void LowPassGate::setOutputNum(LpgFilter::OutputNumber outnum) {
  filter.setOutputNum(outnum);
};

float LowPassGate::process(float sample, float vactrol_in) {
  const float vactrol_output = vactrol.process(vactrol_in);

  const float control_output = controlcircuit.process(vactrol_output * 10.0);
  filter.setVactrolResistance(control_output);
  /* std::cout << "resistance: " << control_output << std::endl; */
  float filter_output = filter.process(sample);

  return filter_output;
}

// -------------------------------------------
// Low pass gate filter emulation
// -------------------------------------------

void LpgFilter::setup(float samplerate, float initVactrolResistance,
                      float initVCAAmount, float initResonance,
                      Linearity initProcessingMode) {

  m_samplerate = samplerate;
  setVactrolResistance(initVactrolResistance);
  setVCAAmount(initVCAAmount);
  setResonance(initResonance);
  setProcessingMode(initProcessingMode);
}

inline void LpgFilter::setVactrolResistance(float rf) {
  m_vactrol_resistance = static_cast<double>(rf);
}

inline void LpgFilter::setVCAAmount(double amount) {
  // Not sure why this is necessary to be honest:
  amount = (1.0 - amount) * 5e6 + 10e3;
  m_vcaness = amount;
}

inline void LpgFilter::setProcessingMode(Linearity mode) {
  m_processing_mode = mode;
}

inline void LpgFilter::setLowpassMode(LowpassProcessing mode) {
  m_lowpass_mode = mode;
}

inline void LpgFilter::setResonance(double resonance) {
  m_resonance = resonance;
};

inline void LpgFilter::setOutputNum(OutputNumber outnum) {
  m_output_num = outnum;
};

// Filter
float LpgFilter::process(float sample) {
  double c1 = (1e-09);
  double c2 = (2.2e-10);
  double c3 = (4.7e-09);

  if (m_lowpass_mode == LpgFilter::LowpassProcessing::Off) {
    c3 = 0.;
  }

  // r1 = 1e3;
  double rf = m_vactrol_resistance;

  // rf = 30e3;
  double r3 = m_vcaness;

  double max_res = 1.0 * (2.0 * c1 * r3 + (c2 + c3) * (r3 + rf)) / (c3 * r3);
  double a = mkutils::constrain(m_resonance, 0.0, max_res);
  /* double a = m_resonance; */
  /* double a = m_resonance * max_res; */
  double x = static_cast<double>(sample);

  double freq = 0.5 / m_samplerate;
  // double freq = 2*pi * (m_vactrol_resistance+1e-3)*0.5/m_samplerate;

  double a1 = 1.0 / (c1 * rf);
  double a2 = -((1.0 / rf) + (1.0 / r3)) / c1;
  /* std::cout << std::to_string(max_res) << std::endl; */

  double b1 = 1.0 / (rf * c2);
  double b2 = -2.0 / (rf * c2);
  double b3 = 1.0 / (rf * c2);
  double b4 = c3 / c2;
  double d1 = a;
  double d2 = -1.;
  m_tanh_xo = tanh(m_state_xo);

  double Dx = 1.0 / (1.0 - b2 * freq);
  double Do = 1.0 / (1.0 - a2 * freq);

  double yx, yo, yd;

  if (m_processing_mode == NonLinear) {

    double Dmas =
        1.0 /
        (1.0 - Dx * (freq * freq * b3 * Do * a1 +
                     b4 * freq * d1 * (1.0 - pow(m_tanh_xo, 2.0)) * Do * a1 +
                     b4 * d2));

    yx = (m_state_x + freq * b1 * x + freq * b3 * Do * m_state_o +
          freq * b4 *
              (m_state_d +
               (1.0 / freq) * d1 *
                   (m_tanh_xo - m_state_xo * (1.0 - pow(m_tanh_xo, 2.0)))) +
          b4 * d1 * (1.0 - pow(m_tanh_xo, 2.0)) * Do * m_state_o) *
         Dx * Dmas;

    yo = (m_state_o + freq * a1 * yx) * Do;

    yd = (m_state_d +
          (1.0 / freq) * d1 *
              (m_tanh_xo - m_state_xo * (1.0 - pow(m_tanh_xo, 2.0)))) +
         (1.0 / freq) * (d1 * ((1.0 - pow(m_tanh_xo, 2.0))) * yo + d2 * yx);

    m_state_x = mkutils::flushed(
        m_state_x + 2.0 * freq * (b1 * x + b2 * yx + b3 * yo + b4 * yd));

    m_state_o = mkutils::flushed(m_state_o + 2.0 * freq * (a1 * yx + a2 * yo));

    m_state_d = mkutils::flushed(
        -(m_state_d +
          (2.0 / freq) * d1 *
              (m_tanh_xo - m_state_xo * (1.0 - pow(m_tanh_xo, 2.0)))) -
        (2.0 / freq) * (d1 * (1.0 - pow(m_tanh_xo, 2.0)) * yo + d2 * yx));

    // Fixme: Is this correct?
    m_state_xo = yo;

  } else {
    double Dmas = 1.0 / (1.0 - Dx * (freq * freq * b3 * Do * a1 +
                                     b4 * freq * d1 * Do * a1 + b4 * d2));

    yx = (m_state_x + freq * b1 * x + freq * b3 * Do * m_state_o +
          freq * b4 * m_state_d + b4 * d1 * Do * m_state_o) *
         Dx * Dmas;

    yo = (m_state_o + freq * a1 * yx) * Do;

    yd = m_state_d + (1.0 / freq) * (d1 * yo + d2 * yx);

    m_state_x = mkutils::flushed(
        m_state_x + 2.0 * freq * (b1 * x + b2 * yx + b3 * yo + b4 * yd));

    m_state_o = mkutils::flushed(m_state_o + 2.0 * freq * (a1 * yx + a2 * yo));

    m_state_d =
        mkutils::flushed(-m_state_d - (2.0 / freq) * (d1 * yo + d2 * yx));
  }
  if (m_output_num == LpgFilter::Vout) {
    return yo;
  } else if (m_output_num == LpgFilter::VactrolOut) {
    return yx;
  } else if (m_output_num == LpgFilter::Differentiator) {
    return yd;
  }

  return 0.;
}

// -------------------------------------------
// Vactrol emulation
// -------------------------------------------

void LpgVactrol::setup(float samplerate) { m_samplerate = samplerate; }

// Vactrol emulation, process an input voltage
float LpgVactrol::process(float sample) {

  // Vactrol time-domain behaviour filter
  float a_base = 1000.0f * M_PI / (m_samplerate);
  float t_down = 3e3f; // Fall time
  t_down = 10.0f + t_down * (1.0f - 0.9f * m_state_1);
  float a_down = a_base / t_down;
  float t_up = 20.0f; // Rise time
  t_up = 1.0f + t_up * (1.0f - 0.999f * m_state_1);
  float a_up = a_base / t_up;

  float dm_state_1 = sample - m_state_1;

  float x = (sample - m_state_1) * a_up / (1.0f + a_up);
  if (dm_state_1 < 0.0f) {
    x = (sample - m_state_1) * a_down / (1.0f + a_down);
  }

  float y = x + m_state_1;
  m_state_1 = y + x;

  return y;
}

// -------------------------------------------
// Vactrol control circuit emulation
// -------------------------------------------

void LpgControlCircuit::setup(double initOffset, double initScale) {
  setOffset(initOffset);
  setScale(initScale);
}

void LpgControlCircuit::setOffset(double offset) { m_offset = offset; }
void LpgControlCircuit::setScale(double scale) { m_scale = scale; }

float LpgControlCircuit::process(float sample) {
  // Inputs
  double Vb = mkutils::constrain(static_cast<double>(sample), -10.0, 50.0);
  double offset = 0.9999 * m_offset + 0.0001;
  double scale = mkutils::constrain(m_scale, 0.0, 1.0);
  // FIXME:
  /* double scale = 0.48f; // This value is tuned for appropriate input range.
   */
  /* std::cout << "voltage base: " << Vb << std::endl; */

  // Constants
  const double A = 3.4645912;
  const double B = 1136.2129956;
  const double Ifmin = 10.1e-6;
  const double Ifmax = 40e-3;
  const double R2max = 10e3;
  const double R6max = 20e3;
  const double R7 = 33e3;
  const double R3 = 150e3;
  const double R5 = 100e3;
  const double R4 = 470e3;
  const double R8 = 4.7e3;
  const double R9 = 470.;
  const double VB = 3.9;
  const double VF = 0.7;
  const double VT = 26e-3;
  const double n = 3.9696;
  const double kl = 6.3862;
  const double G = 2e5;
  const double Vs = 15.;
  const double k0 = 1.468e2;
  const double k1 = 4.9202e-1;
  const double k2 = 4.1667e-4;
  const double k3 = 7.3915e-9;
  const double gamma = 0.0001;

  const double R6 = scale * R6max;
  const double R1 = (1. - offset) * R2max;
  const double R2 = offset * R2max;

  double alpha = 1.0 + (R6 + R7) * (1.0 / R3 + 1.0 / R5);
  double beta = ((1.0 / alpha) - 1.0) / (R6 + R7) - 1.0 / R8;

  double bound1 = 600. * alpha * n * VT / (G * (R6 + R7 - 1. / (alpha * beta)));

  double Ia = Vb / R5 + Vs / (R3 * (1.0 + R1 / R2));

  double V3 = 0.;
  double If = 0.;
  double flag = 0.;

  if (Ia <= -bound1) {
    V3 = -Ia / (alpha * beta);
  } else if (Ia < bound1) {
    const double x =
        G * Ia * (R6 + R7 - 1. / (alpha * beta)) / (alpha * n * VT);
    const double w = k0 + k1 * x + k2 * pow(x, 2.) + k3 * pow(x, 3.);
    V3 = -(alpha / G) * n * VT * w - Ia / (alpha * beta);
  } else {
    V3 = kl * alpha / G * n * VT - Ia * (R6 + R7);
  }

  double Ifbound1 = alpha * (Ifmin - beta * V3);
  double Ifbound2 = VB / (R6 + R7);
  double Ifbound3 = (gamma * G * VB + alpha * R9 * (VB * beta + Ifmax)) /
                    (gamma * G * (R6 + R7) + R9);

  if (Ia <= Ifbound1) {
    If = Ifmin;
  } else if (Ia <= Ifbound2) {
    If = beta * V3 + Ia / alpha;
  } else if (Ia <= Ifbound3) {
    If = gamma * G * (Ia * (R6 + R7) - VB) / (alpha * R9) - beta * VB +
         Ia / alpha;
  } else {
    If = Ifmax;
  }

  double Rf = B + A / pow(If, 1.4);
  double out = Rf;

  return static_cast<float>(out);
}

} // namespace lpg
