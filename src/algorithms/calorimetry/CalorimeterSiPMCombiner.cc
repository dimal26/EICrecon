// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Dima Libman

#include "algorithms/calorimetry/CalorimeterSiPMCombiner.h"

#include <edm4eic/CalorimeterHit.h>
#include <cmath>
#include <cstdint>
#include <format>
#include <map>
#include <stdexcept>
#include <unordered_map>

namespace eicrecon {

void CalorimeterSiPMCombiner::init() {

  if (m_cfg.mode == "weighted") {
    m_mode = Mode::kWeighted;
  } else if (m_cfg.mode == "switch") {
    m_mode = Mode::kSwitch;
  } else if (m_cfg.mode == "A") {
    m_mode = Mode::kOnlyA;
  } else if (m_cfg.mode == "B") {
    m_mode = Mode::kOnlyB;
  } else {
    error("Invalid mode \"{}\" (use weighted, switch, A or B)", m_cfg.mode);
    throw std::runtime_error(std::format("Invalid mode \"{}\"", m_cfg.mode));
  }

  if (m_cfg.photonsPerGeVA <= 0 || m_cfg.numPixelsA <= 0 || m_cfg.photonsPerGeVB <= 0 ||
      m_cfg.numPixelsB <= 0) {
    error("photonsPerGeV and numPixels must be positive for both readouts");
    throw std::runtime_error("CalorimeterSiPMCombiner: invalid SiPM parameters");
  }
  if (m_cfg.sampFrac <= 0) {
    error("sampFrac must be positive");
    throw std::runtime_error("CalorimeterSiPMCombiner: invalid sampFrac");
  }

  m_parA = {m_cfg.photonsPerGeVA, m_cfg.numPixelsA,  m_cfg.crossTalkA,
            m_cfg.afterpulseA,    m_cfg.darkCountsA, m_cfg.noiseSigmaA};
  m_parB = {m_cfg.photonsPerGeVB, m_cfg.numPixelsB,  m_cfg.crossTalkB,
            m_cfg.afterpulseB,    m_cfg.darkCountsB, m_cfg.noiseSigmaB};

  debug("mode {}, switchEnergy {} GeV, sampFrac {}", m_cfg.mode, m_cfg.switchEnergy / dd4hep::GeV,
        m_cfg.sampFrac);
}

void CalorimeterSiPMCombiner::process(const CalorimeterSiPMCombiner::Input& input,
                                      const CalorimeterSiPMCombiner::Output& output) const {

  const auto [hitsA, hitsB] = input;
  auto [outHits]            = output;

  // cellID -> index, ordered map keeps the output order reproducible
  std::map<uint64_t, std::size_t> idxA;
  std::unordered_map<uint64_t, std::size_t> idxB;
  for (std::size_t i = 0; i < hitsA->size(); ++i) {
    idxA[(*hitsA)[i].getCellID()] = i;
  }
  for (std::size_t i = 0; i < hitsB->size(); ++i) {
    idxB[(*hitsB)[i].getCellID()] = i;
  }
  // cells seen only by readout B
  std::map<uint64_t, std::size_t> onlyB;
  for (const auto& [id, i] : idxB) {
    if (idxA.count(id) == 0) {
      onlyB[id] = i;
    }
  }

  // energy of the hit [GeV] is proportional to the fired pixels: invert the SiPM response
  auto reconstruct = [](const edm4eic::CalorimeterHit& h, const sipm::SiPMParameters& p) {
    const double pixels = h.getEnergy() / dd4hep::GeV * p.photonsPerGeV;
    return sipm::invert(pixels, p) / p.photonsPerGeV; // GeV
  };

  auto combine = [&](const edm4eic::CalorimeterHit* hA, const edm4eic::CalorimeterHit* hB) {
    const edm4eic::CalorimeterHit& ref = hA != nullptr ? *hA : *hB;

    const double eA = hA != nullptr ? reconstruct(*hA, m_parA) : 0.;
    const double eB = hB != nullptr ? reconstruct(*hB, m_parB) : 0.;
    const double vA =
        sipm::energyVariance(eA, m_parA, m_cfg.gainUncertainty, m_cfg.numPixelsUncertainty);
    const double vB =
        sipm::energyVariance(hB != nullptr ? eB : eA, m_parB, m_cfg.gainUncertainty,
                             m_cfg.numPixelsUncertainty);

    double energy = 0, variance = 0;
    if (hA == nullptr) { // only B available
      energy = eB, variance = vB;
    } else if (hB == nullptr) { // only A available
      energy = eA, variance = vA;
    } else {
      switch (m_mode) {
      case Mode::kOnlyA:
        energy = eA, variance = vA;
        break;
      case Mode::kOnlyB:
        energy = eB, variance = vB;
        break;
      case Mode::kSwitch:
        if (eA < m_cfg.switchEnergy / dd4hep::GeV) {
          energy = eA, variance = vA;
        } else {
          energy = eB, variance = vB;
        }
        break;
      case Mode::kWeighted: {
        // weights from the expected variances at the readout-A energy (A is unbiased after inversion)
        const double wA = 1. / sipm::energyVariance(eA, m_parA, m_cfg.gainUncertainty,
                                                    m_cfg.numPixelsUncertainty);
        const double wB = 1. / sipm::energyVariance(eA, m_parB, m_cfg.gainUncertainty,
                                                    m_cfg.numPixelsUncertainty);
        energy   = (wA * eA + wB * eB) / (wA + wB);
        variance = 1. / (wA + wB);
        break;
      }
      }
    }

    energy /= m_cfg.sampFrac;
    const double energyError = std::sqrt(variance) / m_cfg.sampFrac;
    trace("cellID {}: E_A {} GeV, E_B {} GeV -> {} +- {} GeV", ref.getCellID(), eA, eB, energy,
          energyError);

    auto hit = outHits->create(ref.getCellID(), static_cast<float>(energy * dd4hep::GeV),
                               static_cast<float>(energyError * dd4hep::GeV), ref.getTime(),
                               ref.getTimeError(), ref.getPosition(), ref.getDimension(),
                               ref.getSector(), ref.getLayer(), ref.getLocal());
    hit.setRawHit(ref.getRawHit()); // keeps the MC links of readout A for clustering / truth
  };

  for (const auto& [id, iA] : idxA) {
    const auto hA = (*hitsA)[iA];
    auto itB      = idxB.find(id);
    if (itB != idxB.end()) {
      const auto hB = (*hitsB)[itB->second];
      combine(&hA, &hB);
    } else {
      combine(&hA, nullptr);
    }
  }
  for (const auto& [id, iB] : onlyB) {
    const auto hB = (*hitsB)[iB];
    combine(nullptr, &hB);
  }
}

} // namespace eicrecon
