import FWCore.ParameterSet.Config as cms
def customise(process):

	if hasattr(process.VtxSmeared,"X0"):
           VertexX = process.VtxSmeared.X0
           VertexY = process.VtxSmeared.Y0
           VertexZ = process.VtxSmeared.Z0

        if hasattr(process.VtxSmeared,"MeanX"):
           VertexX = process.VtxSmeared.MeanX
           VertexY = process.VtxSmeared.MeanY
           VertexZ = process.VtxSmeared.MeanZ

	process.load("SimG4Core.Application.g4SimHits_cfi")
	process.g4SimHits.Generator.HepMCProductLabel   = 'LHCTransport'
	process.g4SimHits.Generator.MinEtaCut        = -13.0
	process.g4SimHits.Generator.MaxEtaCut        = 13.0

	process.g4SimHits.SteppingAction.MaxTrackTime = cms.double(2000.0)
 	process.g4SimHits.StackingAction.MaxTrackTime = cms.double(2000.0)

	process.load('SimTransport.PPSProtonTransport.HectorTransport_cfi')
	process.LHCTransport.PPSHector.VtxMeanX  = VertexX
        process.LHCTransport.PPSHector.VtxMeanY  = VertexY
	process.LHCTransport.PPSHector.VtxMeanZ  = VertexZ
	
	process.transport_step = cms.Path(process.LHCTransport)
	#process.psim = cms.Sequence(cms.SequencePlaceholder("randomEngineStateProducer")*process.LHCTransport*process.g4SimHits)
	process.pgen = cms.Sequence(cms.SequencePlaceholder("randomEngineStateProducer")*process.pgen*process.LHCTransport)

        for path in process.paths:
            print path

	return(process)
