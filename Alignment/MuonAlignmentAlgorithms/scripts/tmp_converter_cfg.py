# sqlite2xml conversion
from Alignment.MuonAlignment.convertSQLitetoXML_cfg import *

process.source = cms.Source("EmptySource",
    numberEventsInRun = cms.untracked.uint32(1),
    firstRun = cms.untracked.uint32(1)
)

process.inertGlobalPositionRcd.connect = "sqlite_file:inertGlobalPositionRcd.StdTags.746p3.DBv2.db"
process.inertGlobalPositionRcd.toGet =  cms.VPSet(cms.PSet(record = cms.string('GlobalPositionRcd'), tag = cms.string('inertGlobalPositionRcd')))

process.PoolDBESSource.connect = "sqlite_file:pp_test_v3_CSC_1_01.db"
process.MuonGeometryDBConverter.outputXML.fileName = "pp_test_v3_CSC_1_01_inertGPR_relToNone.xml"

process.MuonGeometryDBConverter.outputXML.relativeto = "none"

process.MuonGeometryDBConverter.outputXML.suppressDTBarrel = True
process.MuonGeometryDBConverter.outputXML.suppressDTWheels = True
process.MuonGeometryDBConverter.outputXML.suppressDTStations = True
process.MuonGeometryDBConverter.outputXML.suppressDTChambers = False
process.MuonGeometryDBConverter.outputXML.suppressDTSuperLayers = False
process.MuonGeometryDBConverter.outputXML.suppressDTLayers = False

process.MuonGeometryDBConverter.outputXML.suppressCSCEndcaps = True
process.MuonGeometryDBConverter.outputXML.suppressCSCStations = True
process.MuonGeometryDBConverter.outputXML.suppressCSCRings = True
process.MuonGeometryDBConverter.outputXML.suppressCSCChambers = False
process.MuonGeometryDBConverter.outputXML.suppressCSCLayers = False

