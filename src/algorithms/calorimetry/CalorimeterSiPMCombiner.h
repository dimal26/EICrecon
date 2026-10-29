// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Dima Libman

/*
 *  Combine two SiPM readouts of the same calorimeter cells into one hit collection.
 *
 *  Each input is the output of CalorimeterHitReco (sampFrac = 1) for a CalorimeterHitDigi
 *  with readoutType = "sipm". Its energy is therefore proportional to the number of fired
 *  pixels, i.e. still saturated. For every cell this algorithm
 *    1. converts the energy back to fired pixels,
 *    2. inverts the SiPM response (saturation + crosstalk + afterpulses + dark counts),
 *    3. combines readout A and B (inverse-variance weights, or a switch energy),
 *    4. applies the sampling fraction and stores the expected energy uncertainty.
 */

#pragma once

#include <algorithms/algorithm.h>
#include <edm4eic/CalorimeterHitCollection.h>
#include <string>
#include <string_view>

#include "algorithms/calorimetry/CalorimeterSiPMCombinerConfig.h"
#include "algorithms/calorimetry/SiPMResponse.h"
#include "algorithms/interfaces/WithPodConfig.h"

namespace eicrecon {

using CalorimeterSiPMCombinerAlgorithm =
    algorithms::Algorithm<algorithms::Input<edm4eic::CalorimeterHitCollection,
                                            edm4eic::CalorimeterHitCollection>,
                          algorithms::Output<edm4eic::CalorimeterHitCollection>>;

class CalorimeterSiPMCombiner : public CalorimeterSiPMCombinerAlgorithm,
                                public WithPodConfig<CalorimeterSiPMCombinerConfig> {

public:
  CalorimeterSiPMCombiner(std::string_view name)
      : CalorimeterSiPMCombinerAlgorithm{name,
                                         {"inputHitCollectionA", "inputHitCollectionB"},
                                         {"outputHitCollection"},
                                         "Invert the SiPM response of two readouts of the same "
                                         "cells and combine them into one energy."} {}

  void init() final;
  void process(const Input&, const Output&) const final;

private:
  enum class Mode { kWeighted, kSwitch, kOnlyA, kOnlyB };
  Mode m_mode{Mode::kWeighted};
  sipm::SiPMParameters m_parA, m_parB;
};

} // namespace eicrecon
