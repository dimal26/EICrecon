// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Dima Libman

#pragma once

#include "algorithms/calorimetry/CalorimeterSiPMCombiner.h"
#include "services/algorithms_init/AlgorithmsInit_service.h"
#include "extensions/jana/JOmniFactory.h"

namespace eicrecon {

class CalorimeterSiPMCombiner_factory
    : public JOmniFactory<CalorimeterSiPMCombiner_factory, CalorimeterSiPMCombinerConfig> {

public:
  using AlgoT = eicrecon::CalorimeterSiPMCombiner;

private:
  std::unique_ptr<AlgoT> m_algo;

  PodioInput<edm4eic::CalorimeterHit> m_hitsA_input{this};
  PodioInput<edm4eic::CalorimeterHit> m_hitsB_input{this};
  PodioOutput<edm4eic::CalorimeterHit> m_hits_output{this};

  ParameterRef<double> m_photonsPerGeVA{this, "photonsPerGeVA", config().photonsPerGeVA};
  ParameterRef<double> m_numPixelsA{this, "numPixelsA", config().numPixelsA};
  ParameterRef<double> m_crossTalkA{this, "crossTalkA", config().crossTalkA};
  ParameterRef<double> m_afterpulseA{this, "afterpulseA", config().afterpulseA};
  ParameterRef<double> m_darkCountsA{this, "darkCountsA", config().darkCountsA};
  ParameterRef<double> m_noiseSigmaA{this, "noiseSigmaA", config().noiseSigmaA};
  ParameterRef<double> m_photonsPerGeVB{this, "photonsPerGeVB", config().photonsPerGeVB};
  ParameterRef<double> m_numPixelsB{this, "numPixelsB", config().numPixelsB};
  ParameterRef<double> m_crossTalkB{this, "crossTalkB", config().crossTalkB};
  ParameterRef<double> m_afterpulseB{this, "afterpulseB", config().afterpulseB};
  ParameterRef<double> m_darkCountsB{this, "darkCountsB", config().darkCountsB};
  ParameterRef<double> m_noiseSigmaB{this, "noiseSigmaB", config().noiseSigmaB};
  ParameterRef<double> m_gainUncertainty{this, "gainUncertainty", config().gainUncertainty};
  ParameterRef<double> m_numPixelsUncertainty{this, "numPixelsUncertainty",
                                              config().numPixelsUncertainty};
  ParameterRef<std::string> m_mode{this, "mode", config().mode};
  ParameterRef<double> m_switchEnergy{this, "switchEnergy", config().switchEnergy};
  ParameterRef<double> m_sampFrac{this, "sampFrac", config().sampFrac};

  Service<AlgorithmsInit_service> m_algorithmsInit{this};

public:
  void Configure() {
    m_algo = std::make_unique<AlgoT>(GetPrefix());
    m_algo->level(static_cast<algorithms::LogLevel>(logger()->level()));
    m_algo->applyConfig(config());
    m_algo->init();
  }

  void Process(int32_t /* run_number */, uint64_t /* event_number */) {
    m_algo->process({m_hitsA_input(), m_hitsB_input()}, {m_hits_output().get()});
  }
};

} // namespace eicrecon
