// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Wouter Deconinck

#pragma once

#include <string>
#include <vector>

#include <Evaluator/DD4hepUnits.h>
#include <edm4eic/unit_system.h>

namespace eicrecon {

struct CalorimeterHitDigiConfig {

  std::vector<double> eRes;
  double tRes;

  // single hit energy deposition threshold
  double threshold{1.0 * dd4hep::keV};

  // readout settings
  std::string readoutType{"simple"};
  double lightYield{0. / edm4eic::unit::GeV};
  double photonDetectionEfficiency{0.0}; // (light collection efficiency) x (quantum efficiency)
  unsigned long long numEffectiveSipmPixels{0};
  // optional SiPM noise for readoutType "sipm" (all default to 0 = previous behaviour)
  double sipmCrossTalk{0.0};       // prompt crosstalk probability per fired pixel
  double sipmAfterpulse{0.0};      // afterpulse probability x recovery factor
  double sipmDarkCounts{0.0};      // mean dark counts in the integration gate [p.e.]
  double sipmNoiseSigma{0.0};      // electronics noise [p.e.]
  bool sipmOccupancyFluctuation{false}; // binomial fluctuation of the number of fired pixels
  bool sipmPoissonPhotons{false};       // Poisson number of photons reaching the SiPM (else fixed)

  // digitization settings
  unsigned int capADC{1};
  double capTime{1000}; // dynamic range in ns
  double dyRangeADC{1};
  unsigned int pedMeanADC{0};
  double pedSigmaADC{0};
  double resolutionTDC{1};
  std::string corrMeanScale{"1.0"};

  // signal sums
  std::string readout{""};
  std::vector<std::string> fields{};
};

} // namespace eicrecon
