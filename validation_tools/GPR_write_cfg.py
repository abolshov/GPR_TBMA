import FWCore.ParameterSet.Config as cms

process = cms.Process("write")
process.load("CondCore.CondDB.CondDB_cfi")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)
process.source = cms.Source("EmptySource")
from CondCore.DBCommon.CondDBSetup_cfi import CondDBSetup
process.PoolDBOutputService = cms.Service("PoolDBOutputService",
    CondDBSetup,
    # Writing to oracle needs the following shell variable setting (in zsh):
    # export CORAL_AUTH_PATH=/afs/cern.ch/cms/DB/conddb
    # string connect = "oracle://cms_orcoff_int2r/CMS_COND_ALIGNMENT"
    timetype = cms.untracked.string('runnumber'),
    # connect = cms.string('sqlite_file:5mradMuonGamma.db'),
    connect = cms.string('sqlite_file:geom_wGPR.db'),
    # untracked uint32 authenticationMethod = 1
    toPut = cms.VPSet(cms.PSet(
        record = cms.string('GlobalPositionRcd'),
        tag = cms.string('GlobalPositionRcd')
    ))
)
process.GlobalPositionRcdWrite = cms.EDAnalyzer("GlobalPositionRcdWrite",

    # setting 'useEulerAngles' to 'True' lets the module interpret the angles
    # 'alpha', 'beta' and 'gamma' as Euler angles. This is the original behavior
    # of the module. If 'useEulerAngles' is 'False' the angles 'alpha', 'beta'
    # and 'gamma' are interpreted as rotations around the X, Y and Z axes,
    # respectively.
    useEulerAngles = cms.bool(False),

    tracker = cms.PSet(
        y = cms.double(0.0),
        x = cms.double(0.0),
        z = cms.double(0.0),
        alpha = cms.double(0.0),
        beta = cms.double(0.0),
        gamma = cms.double(0.0)
    ),
    muon = cms.PSet(
        y = cms.double(8.29054e-02),
        x = cms.double(6.73283e-02),
        z = cms.double(-9.85197e-03),
        alpha = cms.double(4.24823e-05),
        beta = cms.double(2.72998e-05),
        gamma = cms.double(-1.32678e-04)
    ),
    ecal = cms.PSet(
        y = cms.double(0.0),
        x = cms.double(0.0),
        z = cms.double(0.0),
        alpha = cms.double(0.0),
        beta = cms.double(0.0),
        gamma = cms.double(0.0)
    ),
    hcal = cms.PSet(
        y = cms.double(0.0),
        x = cms.double(0.0),
        z = cms.double(0.0),
        alpha = cms.double(0.0),
        beta = cms.double(0.0),
        gamma = cms.double(0.0)
    ),
    calo = cms.PSet(
        y = cms.double(0.0),
        x = cms.double(0.0),
        z = cms.double(0.0),
        alpha = cms.double(0.0),
        beta = cms.double(0.0),
        gamma = cms.double(0.0)
    )
)

process.p = cms.Path(process.GlobalPositionRcdWrite)
#process.CondDBSetup.DBParameters.messageLevel = 2
