// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Dima Libman

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

#include "algorithms/calorimetry/SiPMResponse.h"

using namespace eicrecon::sipm;
using Catch::Matchers::WithinRel;

// 3 x S14160-6015PS summed, as configured for the B0ECal
static const SiPMParameters p15{3 * 13500 * 0.30, 3 * 159565., 0.005, 0.003 * 0.571,
                                0.9,              std::sqrt(3.) * 5.};
// 1 x S14160-6010PS
static const SiPMParameters p10{13500 * 0.16, 359011., 0.005, 0.002 * 0.625, 0.3, 5.};

TEST_CASE("SiPM response inversion is exact", "[SiPMResponse]") {
  for (double E : {0.01, 0.1, 1., 10., 40., 100., 170.}) {
    for (const auto& p : {p15, p10}) {
      const double N = p.photonsPerGeV * E;
      REQUIRE_THAT(invert(response(N, p), p), WithinRel(N, 1e-9));
    }
  }
}

TEST_CASE("SiPM slope is the derivative of the response", "[SiPMResponse]") {
  for (double E : {0.1, 10., 100.}) {
    const double N = p15.photonsPerGeV * E, h = 1e-3 * N;
    const double numeric = (response(N + h, p15) - response(N - h, p15)) / (2 * h);
    REQUIRE_THAT(slope(N, p15), WithinRel(numeric, 1e-5));
  }
}

TEST_CASE("SiPM energy resolution matches the analytic noise budget", "[SiPMResponse]") {
  // 10 GeV, 0.2% gain and N_MAX calibration: 3x15PS 0.383%, 10PS 0.723%
  REQUIRE_THAT(std::sqrt(energyVariance(10., p15, 0.002, 0.002)) / 10., WithinRel(0.00383, 0.01));
  REQUIRE_THAT(std::sqrt(energyVariance(10., p10, 0.002, 0.002)) / 10., WithinRel(0.00723, 0.01));
  // the two readouts are equally good at the crossover E* = 40.7 GeV
  REQUIRE_THAT(energyVariance(40.74, p15, 0.002, 0.002),
               WithinRel(energyVariance(40.74, p10, 0.002, 0.002), 0.01));
}
