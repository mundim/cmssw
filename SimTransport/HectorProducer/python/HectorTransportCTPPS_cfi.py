import FWCore.ParameterSet.Config as cms

from SimG4Core.Application.hectorParameter_cfi import *
LHCTransport = cms.EDProducer("CTPPSHectorProducer",
    HepMCProductLabel = cms.string('generatorSmeared'),  ## HepMC source to be processed
    CTPPSTransport = cms.bool(True), 
    Verbosity = cms.bool(False),
    CTPPSHector = cms.PSet(
        HectorEtaCut,
        Beam1 = cms.string('SimTransport/HectorProducer/data/LHCB1_Beta0.60_6.5TeV_CR142.5_v6.503.tfs'),
        Beam2 = cms.string('SimTransport/HectorProducer/data/LHCB2_Beta0.60_6.5TeV_CR142.5_v6.503.tfs'),
        CrossingAngleBeam1  = cms.double(142.5), #in mrad
        CrossingAngleBeam2  = cms.double(142.5), #in mrad
        BeamLineLengthCTPPS = cms.double(250.0),
        CTPPSf = cms.double(203.827),    ##in meters
        CTPPSb = cms.double(203.827),    ##in meters
        smearEnergy = cms.bool(True),       ## if False: no Energy smearing(i.e. sigmaEnergy =0.0)
        sigmaEnergy = cms.double(1.11e-4),     ## beam energy dispersion (GeV); if =0.0 the default(=0.79) is used
        smearAng = cms.bool(True),       ## if False: no Angle smearing(i.e. sigmaSTX(Y) =0.0)
        sigmaSTX = cms.double(30.03),     ## x angle dispersion at IP (urad); if =0.0 the default(=30.23) is used
        sigmaSTY = cms.double(30.03),      ## y angle dispersion at IP (urad); if =0.0 the default(=30.23) is used
        CrossAngleCorr = cms.bool(True),
        BeamEnergy = cms.double(6500.0),
        VtxMeanX       = cms.double(0.),
        VtxMeanY       = cms.double(-0.15),
        VtxMeanZ       = cms.double(0.),
        MomentumMin = cms.double(3.000),	
        BeamPositions = cms.vstring("212.551:0.,3.442",
                                    "215.077:0.,3.394",
                                    "219.551:0.,3.336",
                                    "-212.551: 0.,3.311",
                                    "-215.077: 0.,3.261",
                                    "-219.551: 0.,3.201")
#        BeamPositions = cms.vstring("203.826, 3.514,-0.00642",
#                                    "212.551, 3.232,-0.163",
#                                    "215.077, 3.151,-0.192",
#                                    "-203.826, 2.930, 1.309",
#                                    "-212.551, 2.496, 1.415",
#                                    "-215.077, 2.370, 1.446")
    )
)
