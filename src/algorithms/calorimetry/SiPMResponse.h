// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Dima Libman

/*
 *  Analytic SiPM response model (header only, no framework dependencies).
 *
 *  N       detected photons (photo-electrons before saturation)
 *  lambda  = N / Npix                      mean number of photons per pixel
 *  f(N)    = Npix (1 - e^-lambda)          fired pixels (saturation)
 *          + (pXT e^-lambda + cAP) N       crosstalk on free pixels + afterpulses
 *          + nDark                         dark counts in the gate
 *  S(N)    = df/dN = e^-lambda (1 + pXT (1 - lambda)) + cAP
 *
 *  Var[f]  = S^2 N                          photon statistics (Poisson), propagated
 *          + Npix (e^-lambda - (1+lambda) e^-2lambda)   pixel occupancy (binomial)
 *          + (pXT e^-lambda + cAP) N        crosstalk / afterpulse counting
 *          + sigmaNoise^2 + nDark           electronics noise + dark counts
 *
 *  Inversion ("Step B"): solve f(N) = measured by Newton iterations N <- N - (f(N)-measured)/S(N)
 *  Calibration (common, multiplicative) uncertainties enter the variance of the reconstructed N as
 *          (dGain Nt / S)^2 + (dNpix Npix dN/dNpix)^2,   dN/dNpix = (1 - e^-l - l e^-l + pXT l^2 e^-l)/S
 */

#pragma once

#include <algorithm>
#include <cmath>

namespace eicrecon::sipm {

struct SiPMParameters {
  double photonsPerGeV{0}; // detected photo-electrons per GeV (light yield x PDE)
  double numPixels{0};     // effective number of pixels
  double crossTalk{0};     // prompt crosstalk probability per fired pixel
  double afterpulse{0};    // afterpulse probability x recovery factor (p_AP * p_rec)
  double darkCounts{0};    // mean dark counts in the integration gate [p.e.]
  double noiseSigma{0};    // electronics noise [p.e.]
};

// mean measured signal [pixels] for N detected photons
inline double response(double N, const SiPMParameters& p) {
  const double el = std::exp(-N / p.numPixels);
  return p.numPixels * (1. - el) + (p.crossTalk * el + p.afterpulse) * N + p.darkCounts;
}

// slope dS/dN [pixels per photon]
inline double slope(double N, const SiPMParameters& p) {
  const double l = N / p.numPixels, el = std::exp(-l);
  return el * (1. + p.crossTalk * (1. - l)) + p.afterpulse;
}

// invert the response: measured pixels -> detected photons (Newton, starts from the pure-saturation inverse)
inline double invert(double measured, const SiPMParameters& p) {
  const double M = p.numPixels;
  double N       = (measured < 0.95 * M) ? -M * std::log(1. - measured / M) : 3. * M;
  for (int i = 0; i < 60; ++i) {
    const double d = (response(N, p) - measured) / std::max(slope(N, p), 1e-6);
    N -= d;
    if (std::fabs(d) < 1e-6 * std::max(1., std::fabs(N))) {
      break;
    }
  }
  return N;
}

// variance [GeV^2] of the energy reconstructed from one readout at true energy E [GeV]
// dGain, dNpix: relative calibration uncertainties (e.g. 0.002)
inline double energyVariance(double E, const SiPMParameters& p, double dGain, double dNpix) {
  const double N  = p.photonsPerGeV * std::max(E, 1e-6);
  const double M  = p.numPixels;
  const double l  = N / M, el = std::exp(-l);
  const double S  = slope(N, p);
  const double Nt = response(N, p);
  double varN     = (S * S * N + M * (el - (1. + l) * el * el) + (p.crossTalk * el + p.afterpulse) * N +
                 p.noiseSigma * p.noiseSigma + p.darkCounts) /
                (S * S);
  varN += std::pow(dGain * Nt / S, 2);
  varN += std::pow(dNpix * M * (1. - el - l * el + p.crossTalk * l * l * el) / S, 2);
  return varN / (p.photonsPerGeV * p.photonsPerGeV);
}

} // namespace eicrecon::sipm
