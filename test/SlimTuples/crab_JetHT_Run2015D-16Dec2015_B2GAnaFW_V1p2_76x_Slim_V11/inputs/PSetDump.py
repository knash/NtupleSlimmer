import FWCore.ParameterSet.Config as cms

process = cms.Process("slimntuple")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/user/lcorcodi/BstarToTW_M-1400_RH_TuneCUETP8M1_13TeV-madgraph-pythia8/crab_BstarToTW_M-1400_RH_TuneCUETP8M1_13TeV-madgraph-pythia8/160318_162851/0000/B2GEDMNtuple_2.root')
)
process.jetsAK8 = cms.EDProducer("SlimUserData_jetsAK8",
    jer = cms.string('nominal'),
    jes = cms.string('nominal')
)


process.jetsAK8jerdown = cms.EDProducer("SlimUserData_jetsAK8",
    jer = cms.string('down'),
    jes = cms.string('nominal')
)


process.jetsAK8jerup = cms.EDProducer("SlimUserData_jetsAK8",
    jer = cms.string('up'),
    jes = cms.string('nominal')
)


process.jetsAK8jesdown = cms.EDProducer("SlimUserData_jetsAK8",
    jer = cms.string('nominal'),
    jes = cms.string('down')
)


process.jetsAK8jesup = cms.EDProducer("SlimUserData_jetsAK8",
    jer = cms.string('nominal'),
    jes = cms.string('up')
)


process.Filter = cms.EDFilter("SlimUserData_Filter")


process.edmNtuplesOut = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('f1')
    ),
    dropMetaData = cms.untracked.string('ALL'),
    fileName = cms.untracked.string('B2GEDMNtuple_slim.root'),
    outputCommands = cms.untracked.vstring('drop *', 
        'keep *_*_*_slimntuple', 
        'drop edmTriggerResults_*_*_slimntuple', 
        'keep *_subjetsCmsTopTag_subjetCmsTopTagCSV_*', 
        'keep *_subjetsAK8_subjetAK8CSV_*', 
        'keep *_eventUserData_*_*', 
        'drop *_eventUserData_v*_*', 
        'keep *_eventInfo_*_*', 
        'keep *_met_*_*', 
        'drop *_met_Px_*', 
        'drop *_met_Py_*', 
        'keep *_pdfWeights*_*_*')
)


process.p = cms.Path(process.Filter+process.jetsAK8)


process.f1 = cms.Path(process.Filter)


process.endPath = cms.EndPath(process.edmNtuplesOut)


process.MessageLogger = cms.Service("MessageLogger",
    FrameworkJobReport = cms.untracked.PSet(
        FwkJob = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        optionalPSet = cms.untracked.bool(True)
    ),
    categories = cms.untracked.vstring('FwkJob', 
        'FwkReport', 
        'FwkSummary', 
        'Root_NoDictionary'),
    cerr = cms.untracked.PSet(
        FwkJob = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        FwkReport = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True),
            reportEvery = cms.untracked.int32(1000)
        ),
        FwkSummary = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True),
            reportEvery = cms.untracked.int32(1)
        ),
        INFO = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        Root_NoDictionary = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000)
        ),
        noTimeStamps = cms.untracked.bool(False),
        optionalPSet = cms.untracked.bool(True),
        threshold = cms.untracked.string('INFO')
    ),
    cerr_stats = cms.untracked.PSet(
        optionalPSet = cms.untracked.bool(True),
        output = cms.untracked.string('cerr'),
        threshold = cms.untracked.string('WARNING')
    ),
    cout = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    debugModules = cms.untracked.vstring(),
    debugs = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    default = cms.untracked.PSet(

    ),
    destinations = cms.untracked.vstring('warnings', 
        'errors', 
        'infos', 
        'debugs', 
        'cout', 
        'cerr'),
    errors = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    fwkJobReports = cms.untracked.vstring('FrameworkJobReport'),
    infos = cms.untracked.PSet(
        Root_NoDictionary = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        optionalPSet = cms.untracked.bool(True),
        placeholder = cms.untracked.bool(True)
    ),
    statistics = cms.untracked.vstring('cerr_stats'),
    suppressDebug = cms.untracked.vstring(),
    suppressInfo = cms.untracked.vstring(),
    suppressWarning = cms.untracked.vstring(),
    warnings = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    )
)


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

