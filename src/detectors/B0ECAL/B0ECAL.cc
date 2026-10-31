// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 - 2025 Whitney Armstrong, Sylvester Joosten, Chao Peng, David Lawrence, Wouter Deconinck, Kolja Kauder, Nathan Brei, Dmitry Kalinkin, Derek Anderson, Michael Pitt

#include <Evaluator/DD4hepUnits.h>
#include <JANA/JApplicationFwd.h>
#include <JANA/Utils/JTypeInfo.h>
#include <cmath>
#include <string>
#include <variant>
#include <vector>

#include "extensions/jana/JOmniFactoryGeneratorT.h"
#include "factories/calorimetry/CalorimeterClusterRecoCoG_factory.h"
#include "factories/calorimetry/CalorimeterClusterShape_factory.h"
#include "factories/calorimetry/CalorimeterHitDigi_factory.h"
#include "factories/calorimetry/CalorimeterHitReco_factory.h"
#include "factories/calorimetry/CalorimeterIslandCluster_factory.h"
#include "factories/calorimetry/CalorimeterSiPMCombiner_factory.h"
#include "factories/calorimetry/CalorimeterTruthClustering_factory.h"

extern "C" {
void InitPlugin(JApplication* app) {

  using namespace eicrecon;

  InitJANAPlugin(app);

  // ---------------------------------------------------------------------------------------
  // SiPM readout: every crystal is read by 3 x S14160-6015PS (15 um) + 1 x S14160-6010PS (10 um)
  // Light reaching one 6x6 mm2 sensor: 13500 photons / GeV (PbWO4 light yield x collection).
  // The 3 x 15PS are summed: equivalent to one sensor with 3x the light and 3x the pixels
  // (same occupancy), dark counts x3 and electronics noise x sqrt(3).
  // Photon statistics come from the SiPM simulation itself -> eRes = 0 (no double counting).
  // ---------------------------------------------------------------------------------------
  constexpr double B0ECal_photonsPerSensor = 13500. / dd4hep::GeV;
  constexpr double B0ECal_PDE15            = 0.30;
  constexpr double B0ECal_PDE10            = 0.16;
  constexpr unsigned long long B0ECal_Npix15 = 159565; // S14160-6015PS
  constexpr unsigned long long B0ECal_Npix10 = 359011; // S14160-6010PS
  constexpr double B0ECal_crossTalk        = 0.005;
  constexpr double B0ECal_afterpulse15     = 0.003 * 0.571; // p_AP x p_rec
  constexpr double B0ECal_afterpulse10     = 0.002 * 0.625;
  constexpr double B0ECal_darkCounts       = 3000e3 * 100e-9; // DCR [Hz] x gate [s] = 0.3 p.e.
  constexpr double B0ECal_noiseSigma       = 5.0;             // p.e. per sensor

  // readout A: 3 x 15PS summed (keeps the standard collection names -> MC links used downstream)
  app->Add(new JOmniFactoryGeneratorT<CalorimeterHitDigi_factory>(
      "B0ECalRawHits", {"EventHeader", "B0ECalHits"},
      {"B0ECalRawHits", "B0ECalRawHitLinks", "B0ECalRawHitAssociations"},
      {
          .eRes                      = {0.0, 0.0, 0.0},
          .tRes                      = 0.0 * dd4hep::ns,
          .threshold                 = 5.0 * dd4hep::MeV,
          .readoutType               = "sipm",
          .lightYield                = 3 * B0ECal_photonsPerSensor,
          .photonDetectionEfficiency = B0ECal_PDE15,
          .numEffectiveSipmPixels    = 3 * B0ECal_Npix15,
          .sipmCrossTalk             = B0ECal_crossTalk,
          .sipmAfterpulse            = B0ECal_afterpulse15,
          .sipmDarkCounts            = 3 * B0ECal_darkCounts,
          .sipmNoiseSigma            = std::sqrt(3.) * B0ECal_noiseSigma,
          .sipmOccupancyFluctuation  = true,
          .sipmPoissonPhotons        = true,
          .capADC                    = 16384,
          .dyRangeADC                = 170 * dd4hep::GeV,
          .pedMeanADC                = 100,
          .pedSigmaADC               = 0, // electronics noise is in sipmNoiseSigma
          .resolutionTDC             = 1e-11,
          .corrMeanScale             = "1.0",
          .readout                   = "B0ECalHits",
      },
      app));
  // readout B: 1 x 10PS
  app->Add(new JOmniFactoryGeneratorT<CalorimeterHitDigi_factory>(
      "B0ECalRawHits10PS", {"EventHeader", "B0ECalHits"},
      {"B0ECalRawHits10PS", "B0ECalRawHitLinks10PS", "B0ECalRawHitAssociations10PS"},
      {
          .eRes                      = {0.0, 0.0, 0.0},
          .tRes                      = 0.0 * dd4hep::ns,
          .threshold                 = 5.0 * dd4hep::MeV,
          .readoutType               = "sipm",
          .lightYield                = B0ECal_photonsPerSensor,
          .photonDetectionEfficiency = B0ECal_PDE10,
          .numEffectiveSipmPixels    = B0ECal_Npix10,
          .sipmCrossTalk             = B0ECal_crossTalk,
          .sipmAfterpulse            = B0ECal_afterpulse10,
          .sipmDarkCounts            = B0ECal_darkCounts,
          .sipmNoiseSigma            = B0ECal_noiseSigma,
          .sipmOccupancyFluctuation  = true,
          .sipmPoissonPhotons        = true,
          .capADC                    = 16384,
          .dyRangeADC                = 170 * dd4hep::GeV,
          .pedMeanADC                = 100,
          .pedSigmaADC               = 0,
          .resolutionTDC             = 1e-11,
          .corrMeanScale             = "1.0",
          .readout                   = "B0ECalHits",
      },
      app));

  // ADC -> energy proportional to fired pixels (still saturated), sampFrac applied after the combination
  app->Add(new JOmniFactoryGeneratorT<CalorimeterHitReco_factory>(
      "B0ECalRecHits15PS", {"B0ECalRawHits"}, {"B0ECalRecHits15PS"},
      {
          .capADC          = 16384,
          .dyRangeADC      = 170. * dd4hep::GeV,
          .pedMeanADC      = 100,
          .pedSigmaADC     = 0,
          .resolutionTDC   = 1e-11,
          .thresholdFactor = 0.0,
          .thresholdValue  = 1.0, // 1 ADC = 10.4 MeV
          .sampFrac        = "1.0",
          .readout         = "B0ECalHits",
          .sectorField     = "sector",
      },
      app));
  app->Add(new JOmniFactoryGeneratorT<CalorimeterHitReco_factory>(
      "B0ECalRecHits10PS", {"B0ECalRawHits10PS"}, {"B0ECalRecHits10PS"},
      {
          .capADC          = 16384,
          .dyRangeADC      = 170. * dd4hep::GeV,
          .pedMeanADC      = 100,
          .pedSigmaADC     = 0,
          .resolutionTDC   = 1e-11,
          .thresholdFactor = 0.0,
          .thresholdValue  = 1.0,
          .sampFrac        = "1.0",
          .readout         = "B0ECalHits",
          .sectorField     = "sector",
      },
      app));

  // invert the SiPM response and combine 3x15PS (A) with 10PS (B) -> standard B0ECalRecHits
  // default: switch at E* = 40.7 GeV (resolution crossover for 0.2% gain and N_pix calibration)
  app->Add(new JOmniFactoryGeneratorT<CalorimeterSiPMCombiner_factory>(
      "B0ECalRecHits", {"B0ECalRecHits15PS", "B0ECalRecHits10PS"}, {"B0ECalRecHits"},
      {
          .photonsPerGeVA       = 3 * B0ECal_photonsPerSensor * B0ECal_PDE15 * dd4hep::GeV,
          .numPixelsA           = 3. * B0ECal_Npix15,
          .crossTalkA           = B0ECal_crossTalk,
          .afterpulseA          = B0ECal_afterpulse15,
          .darkCountsA          = 3 * B0ECal_darkCounts,
          .noiseSigmaA          = std::sqrt(3.) * B0ECal_noiseSigma,
          .photonsPerGeVB       = B0ECal_photonsPerSensor * B0ECal_PDE10 * dd4hep::GeV,
          .numPixelsB           = 1. * B0ECal_Npix10,
          .crossTalkB           = B0ECal_crossTalk,
          .afterpulseB          = B0ECal_afterpulse10,
          .darkCountsB          = B0ECal_darkCounts,
          .noiseSigmaB          = B0ECal_noiseSigma,
          .gainUncertainty      = 0.002,
          .numPixelsUncertainty = 0.002,
          .mode                 = "switch", // 3x15PS below switchEnergy, 10PS above; or "weighted", "A", "B"
          .switchEnergy         = 40.7 * dd4hep::GeV,
          .sampFrac             = 0.998,
      },
      app));
  app->Add(new JOmniFactoryGeneratorT<CalorimeterTruthClustering_factory>(
      "B0ECalTruthProtoClusters", {"B0ECalRecHits", "B0ECalRawHitLinks"},
      {"B0ECalTruthProtoClusters"}, app));
  app->Add(new JOmniFactoryGeneratorT<CalorimeterIslandCluster_factory>(
      "B0ECalIslandProtoClusters", {"B0ECalRecHits"}, {"B0ECalIslandProtoClusters"},
      {
          .adjacencyMatrix{},
          .peakNeighbourhoodMatrix{},
          .readout{},
          .sectorDist = 5.0 * dd4hep::cm,
          .localDistXY{},
          .localDistXZ{},
          .localDistYZ{},
          .globalDistRPhi{},
          .globalDistEtaPhi{},
          .dimScaledLocalDistXY          = {1.8, 1.8},
          .splitCluster                  = false,
          .minClusterHitEdep             = 1.0 * dd4hep::MeV,
          .minClusterCenterEdep          = 30.0 * dd4hep::MeV,
          .transverseEnergyProfileMetric = "globalDistEtaPhi",
          .transverseEnergyProfileScale  = 1.,
          .transverseEnergyProfileScaleUnits{},
      },
      app));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
      "B0ECalClustersWithoutShapes",
      {
          "B0ECalIslandProtoClusters", // edm4eic::ProtoClusterCollection
          "B0ECalRawHitLinks",         // edm4eic::MCRecoCalorimeterHitLink
          "B0ECalRawHitAssociations"   // edm4eic::MCRecoCalorimeterHitAssociationCollection
      },
      {"B0ECalClustersWithoutShapes", // edm4eic::Cluster
       "B0ECalClusterLinksWithoutShapes",
       "B0ECalClusterAssociationsWithoutShapes"}, // edm4eic::MCRecoClusterParticleAssociation
      {.energyWeight = "log", .sampFrac = 1.0, .logWeightBase = 3.6, .enableEtaBounds = false},
      app));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "B0ECalClusters", {"B0ECalClustersWithoutShapes", "B0ECalClusterLinksWithoutShapes"},
      {"B0ECalClusters", "B0ECalClusterLinks", "B0ECalClusterAssociations"},
      {.energyWeight = "log", .logWeightBase = 3.6}, app));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterRecoCoG_factory>(
      "B0ECalTruthClustersWithoutShapes",
      {
          "B0ECalTruthProtoClusters", // edm4eic::ProtoClusterCollection
          "B0ECalRawHitLinks",        // edm4eic::MCRecoCalorimeterHitLink
          "B0ECalRawHitAssociations"  // edm4eic::MCRecoCalorimeterHitAssociationCollection
      },
      {"B0ECalTruthClustersWithoutShapes", // edm4eic::Cluster
       "B0ECalTruthClusterLinksWithoutShapes",
       "B0ECalTruthClusterAssociationsWithoutShapes"}, // edm4eic::MCRecoClusterParticleAssociation
      {.energyWeight = "log", .sampFrac = 1.0, .logWeightBase = 6.2, .enableEtaBounds = false},
      app));

  app->Add(new JOmniFactoryGeneratorT<CalorimeterClusterShape_factory>(
      "B0ECalTruthClusters",
      {"B0ECalTruthClustersWithoutShapes", "B0ECalTruthClusterLinksWithoutShapes"},
      {"B0ECalTruthClusters", "B0ECalTruthClusterLinks", "B0ECalTruthClusterAssociations"},
      {.energyWeight = "log", .logWeightBase = 6.2}, app));
}
}
