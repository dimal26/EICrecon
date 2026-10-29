// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Dima Libman

#pragma once

#include <Evaluator/DD4hepUnits.h>
#include <string>

namespace eicrecon {

struct CalorimeterSiPMCombinerConfig {

  // ---- readout A (high light yield, saturates early, e.g. 3 x 15 um SiPMs summed) ----
  double photonsPerGeVA{0}; // detected photo-electrons per GeV (light yield x PDE x number of sensors)
  double numPixelsA{0};     // effective pixels (x number of sensors)
  double crossTalkA{0};
  double afterpulseA{0};
  double darkCountsA{0};
  double noiseSigmaA{0};

  // ---- readout B (low light yield, large dynamic range, e.g. 1 x 10 um SiPM) ----
  double photonsPerGeVB{0};
  double numPixelsB{0};
  double crossTalkB{0};
  double afterpulseB{0};
  double darkCountsB{0};
  double noiseSigmaB{0};

  // calibration uncertainties, used only to compute the weights and the hit energy error
  double gainUncertainty{0.002};
  double numPixelsUncertainty{0.002};

  // "weighted": inverse-variance average of A and B (no free parameter)
  // "switch"  : A below switchEnergy, B above
  // "A" / "B" : one readout only
  std::string mode{"weighted"};
  double switchEnergy{40.7 * dd4hep::GeV};

  // sampling fraction applied to the combined energy
  // (the input CalorimeterHitReco collections must use sampFrac = 1)
  double sampFrac{1.0};
};

} // namespace eicrecon
