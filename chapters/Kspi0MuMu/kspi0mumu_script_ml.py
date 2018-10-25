#! /usr/bin/env python
from ROOT import *
#from Bender.MainMC import *
import BenderTools.TisTos

from BenderAlgo.extrafunctions import *# trackHits, triggerBlock, muIDDetailedInfo, more_muon_things, triggerBlockDaughters

from BenderAlgo.BenderK0SPi0MuMu import * 
from BenderAlgo.BenderV0 import *
from Configurables import DaVinci

from Gaudi.Configuration import NTupleSvc, EventSelector
from   Bender.MainMC   import * # SUCCESS, AlgoMC

import LinkerInstances.eventassoc as MCLinker
import os,sys
MCParticle = cpp.LHCb.MCParticle
Track = cpp.LHCb.Track

#SAMPLE = "MCUpgrade_MU"
#SAMPLE = "MC12_Down"
#SAMPLE = "MC12_K3pi"
SAMPLE = "Data16_Down"

## Overwrite SAMPLE if the script argument is specified (VC-style)
if len(sys.argv) > 1:
    if "_" in sys.argv[1]: SAMPLE = str(sys.argv[1])#"MC12_Up"
print SAMPLE

ISGRID = 0

TRIGGER = 0
SIMULATION = int("MC" in SAMPLE)
RERUN_STP_new = 0 * SIMULATION
RERUN_STP_21 = 0 * (not RERUN_STP_new)*SIMULATION 
## ML and other dummies: you have to change this if not ganga stuff (then INTERACTIVE cant be zero)
INTERACTIVE = 10000*(not ISGRID) # 0 for all
_DEBUG = 1
TUPLE_PATH = "/afs/cern.ch/user/m/mlucioma/cmtuser/Erasmus_v13r4/Phys/Ks2Pi0MuMuTuples/python/Ks2Pi0MuMuTuples/"
#TUPLE_PATH = os.environ["HOME"] + "/vol5/"
job = "" + (not INTERACTIVE )*("_" + sys.argv[-1])
if ".py" in job: job = ""

_stpnew = RERUN_STP_new or ("Data15" in SAMPLE)##("Data16" in SAMPLE) ##ML:does this have to be changed?

DataDic = {}
DataDic["MC12_Strip_Up"] = "MC_2012_34112402_Beam4000GeV2012MagUpNu2.5Pythia8_Sim08f_Digi13_Trig0x409f0045_Reco14a_Stripping20r0p1Filtered_KSPI0MUMU.STRIP.DST" + job + ".py"
DataDic["MC12_Strip_Down"] = "MC_2012_34112402_Beam4000GeV2012MagDownNu2.5Pythia8_Sim08f_Digi13_Trig0x409f0045_Reco14a_Stripping20r0p1Filtered_KSPI0MUMU.STRIP.DST" + job +".py"
DataDic["Data12_Down"] = "Stp21_Leptonic_MagDown2012" + job+ ".py"
DataDic["Data12_Up"] = "Stp21_Leptonic_MagUp2012" + job + ".py"
DataDic["MC12_Down"] = "MC12_MD_Pythia8.py"#"MC_2012_34112401_Beam4000GeV2012MagDownNu2.5Pythia8_Sim08e_Digi13_Trig0x409f0045_Reco14a_Stripping20NoPrescalingFlagged_ALLSTREAMS.DST.py"
DataDic["MC12_Up"]= "MC12_MU_Pythia8.py"
DataDic["MC12_MB_Down"] = "MC_2012_30000000_Beam4000GeV2012MagDownNu2.5Pythia8_Sim08f_Digi13_Trig0x409f0045_Reco14a_Stripping20NoPrescalingFlagged_ALLSTREAMS.DST.py"
DataDic["MCUpgrade_MD"] = "Upgrade_MD.py"
DataDic["MCUpgrade_MU"] = "Upgrade_MU.py"
DataDic["Data15_Down"] = "S24_MD" + job + ".py"
DataDic["Data15_Up"] = "S24_MU"+job+".py"
DataDic["MC12_K3pi"] = "MC_2012_Beam4000GeV2012MagUpNu2.5Pythia8_Sim08i_Digi13_Trig0x409f0045_Reco14c_Stripping21NoPrescalingFlagged_34102408_ALLSTREAMS.DST.py"
DataDic["Data16_Down"] = "LHCb_Collision16_Beam6500GeVVeloClosedMagUp_RealData_Reco16_Stripping26_90000000_DIMUON.DST.py"
DataTypeDic = {}
for key in DataDic.keys():
    if "12" in key : DataTypeDic[key] = "2012"
    elif "11" in key: DataTypeDic[key] = "2011"
    elif "15" in key: DataTypeDic[key] = "2015"
    elif "16" in key: DataTypeDic[key] = "2016"
    elif "Upgrade" in key: DataTypeDic[key] = "Upgrade"

dataType = DataTypeDic[SAMPLE]

ParticlePath = {}
dataFormats = {}

for key in DataDic.keys():
    if "_Strip_" in key:
        ParticlePath[key] = "Kspi0mumu.Strip"
        dataFormats[key] = "MDST"
    elif "MCUpgrade" in key:
        ParticlePath[key] = ""
        dataFormats[key] = "LDST"
    elif "MC12" in key:
        ParticlePath[key] = ""
        dataFormats[key] = "DST"
    elif "Data15" in key :
        ParticlePath[key] = "Dimuon"
        dataFormats[key] = "DST"
    elif "Data16" in key :
        ParticlePath[key] = "Dimuon"
        dataFormats[key] = "DST" #ML: check this
    elif "Data" in key :
        ParticlePath[key] = "Leptonic"
        dataFormats[key] = "MDST"

dataFormat = dataFormats[SAMPLE]
####Examples of Stp  LINE PATHS:   ## Otherwise use dst-dump -f <the eos file> to see the content including unpacking
#/Event/Phys/Lambda02PiMuLine/Particles
#/Event/Kspi0mumu.Strip/pPhys/Particles ##PACKED!!
RootsInTES = {}
RootsInTES["MDST"] = "/Event/Leptonic/"
RootsInTES["DST"] = False
RootsInTES["LDST"] = False
RootInTES = RootsInTES[dataFormat]
#
#Raw event juggler to split Other/RawEvent into Velo/RawEvent and Tracker/RawEvent
#
#from Configurables import RawEventJuggler
#juggler = RawEventJuggler( DataOnDemand=True, Input=2.0, Output=4.0 )


if RERUN_STP_new:
    
    from CommonParticles.Utils import DefaultTrackingCuts
    DefaultTrackingCuts().Cuts  = { "Chi2Cut" : [ 0, 3 ],
                                "CloneDistCut" : [5000, 9e+99 ] }
    confname='RnS'
    from StrippingConf.Configuration import StrippingConf

    from StrippingSelections import buildersConf
    confs = buildersConf()
    from StrippingSelections.Utils import lineBuilder, buildStreamsFromBuilder
    streams = buildStreamsFromBuilder(confs,confname)
    leptonicMicroDSTname   = 'Leptonic'
    charmMicroDSTname      = 'Charm'
    pidMicroDSTname        = 'PID'
    bhadronMicroDSTname    = 'Bhadron'
    mdstStreams = [ leptonicMicroDSTname,charmMicroDSTname,pidMicroDSTname,bhadronMicroDSTname ]
    dstStreams  = [ "BhadronCompleteEvent", "CharmCompleteEvent", "CharmToBeSwum", "Dimuon",
                    "EW", "Semileptonic", "Calibration", "MiniBias", "Radiative" ]

    stripTESPrefix = 'Strip'

    from Configurables import ProcStatusCheck
    
    sc = StrippingConf( Streams = streams,
                        MaxCandidates = 2000,
                        AcceptBadEvents = False,
                        BadEventSelection = ProcStatusCheck(),
                        TESPrefix = stripTESPrefix,
                        ActiveMDSTStream = True,
                        Verbose = True,
                        DSTStreams = dstStreams,
                        MicroDSTStreams = mdstStreams )
    DaVinci().appendToMainSequence( [ sc.sequence() ] )


if RERUN_STP_21:
    ## from StrippingConf.Configuration import StrippingConf, StrippingStream
##     from StrippingSettings.Utils import strippingConfiguration
##     from StrippingArchive.Utils import buildStream, cloneLinesFromStream
##     from StrippingArchive.Utils import buildStreams, cloneLinesFromStream
    
##     from StrippingArchive import strippingArchive

##     stripping='stripping21'
##     config  = strippingConfiguration(stripping)
##     archive = strippingArchive(stripping)

##     _leptonic = buildStream(stripping=config, streamName='Leptonic', archive=archive)
##     streams = []
##     _leptonic.lines[:] = [ x for x in _leptonic.lines if 'K0s2Pi0MuMuSignalLine' in x.name() ]
##     for line in _leptonic.lines :
##         print "leptonic has a line called " + line.name()

##     streams.append( _leptonic )

##     for stream in streams: 
##         for line in stream.lines:
##             line._prescale = 1.0 

##     AllStreams = StrippingStream("KspizmmFULL.Strip")
    
##     for stream in streams:
##         AllStreams.appendLines(stream.lines)

##     sc = StrippingConf( Streams = [ AllStreams ],
##                     MaxCandidates = 2000 )

    from StrippingConf.Configuration import StrippingConf
    
    #Tighten Trk Chi2 to <3
    from CommonParticles.Utils import DefaultTrackingCuts
    DefaultTrackingCuts().Cuts  = { "Chi2Cut" : [ 0, 3 ],
                                "CloneDistCut" : [5000, 9e+99 ] }

    #Now build the stream
    from StrippingConf.StrippingStream import StrippingStream
    stream = StrippingStream("Test")

    from StrippingSelections.StrippingK0s2Pi0MuMuLines import config_default
    from StrippingSelections.StrippingK0s2Pi0MuMuLines import K0s2Pi0MuMuLinesConf
    config_default['SidebandLinePrescale']=1.0
    myBuilder = K0s2Pi0MuMuLinesConf( name='K0s2Pi0MuMu', config = config_default)
    stream.appendLines( myBuilder.lines() )

    config_pipi={
       'NoMuIDLinePrescale'    : 1,#1e-03,
       'NoMuIDLinePostscale'   : 1,
       'K0s2mmLinePrescale'  : 1,
       'K0s2mmLinePostscale'  : 1,
       'K0s2mmSBLinePrescale'  : 0.1,

        'K0s2mmSBLinePostscale'  : 1,
        'minMuPT' : 300,  #MeV
        'minKsPT' : 600,  #MeV
        'minMuPT' : 0,  #MeV
        'minKsPT' : 0,  #MeV
        }
    from StrippingSelections.StrippingK0s2MuMuLines import K0s2MuMuLinesConf
    config_default['NoMuIDLinePrescale']=1.0
    myBuilder_pipi = K0s2MuMuLinesConf( name='K0s2MuMu', config = config_pipi)
    stream.appendLines( myBuilder_pipi.lines() )



    #Standard configuration of Stripping, do NOT change them
    from Configurables import  ProcStatusCheck
    filterBadEvents =  ProcStatusCheck()
    
    sc = StrippingConf( Streams = [ stream ],
                       MaxCandidates = 2000,
                       AcceptBadEvents = False,
                       BadEventSelection = filterBadEvents,
                       TESPrefix = 'Strip'
                       )


    from Configurables import AuditorSvc, ChronoAuditor
    AuditorSvc().Auditors.append( ChronoAuditor("Chrono") )
    
    from Configurables import StrippingReport
    sr = StrippingReport(Selections = sc.selections())
    
    from Configurables import AlgorithmCorrelationsAlg
    ac = AlgorithmCorrelationsAlg(Algorithms = sc.selections())
    
    DaVinci().appendToMainSequence( [ sc.sequence() ] )
    DaVinci().appendToMainSequence( [ sr ] )
    DaVinci().appendToMainSequence( [ ac ] )


DaVinci().EvtMax = 0
if RootInTES:
    DaVinci().InputType = "DST"#dataFormat
    DaVinci().RootInTES = RootInTES

print "Before the DataType stuff"
DaVinci().DataType = dataType
print "dataType=",dataType
from Configurables import CondDB
CondDB().Upgrade    = ("Upgrade" in SAMPLE)

# is this data?
DaVinci().Simulation = SIMULATION
DaVinci().Lumi = True
# if MC, add the SeqKsResolved
DaVinci().UserAlgorithms = []


TUPLE_FILE = TUPLE_PATH + "kspi0mumu_ntupleFULL" + SAMPLE + job + ".root"
NTupleSvc( Output =["T DATAFILE='"+ TUPLE_FILE + "' TYP='ROOT' OPT='NEW'"] )

## RUN THE TRIGGER
#L0SelRepSeq = GaudiSequencer("L0SelRepSeq")



#L0SelRepSeq.MeasureTime = True
#from Configurables import L0SelReportsMaker, L0DecReportsMaker
#L0SelRepSeq.Members += [ L0DecReportsMaker(), L0SelReportsMaker() ]
#DaVinci().UserAlgorithms += [ L0SelRepSeq ]

if 'KS2PI0MUMUTUPLESROOT' in os.environ.keys(): importOptions("$KS2PI0MUMUTUPLESROOT/python/datacards/" + DataDic[SAMPLE])
else:
    print 'WARNING: KS2PI0MUMUTUPLESROOT not defined'
    importOptions("../datacards/" + DataDic[SAMPLE])

if not ISGRID: Eostize(EventSelector())

gaudi = GaudiPython.AppMgr()

algos = []

#######################
## for stripping signal (for MC, we do not need this two categories)
#kspi0mumu_StrippingSignal = StrippingK0SPi0MuMu("Kspi0mumuStrippingSignal")
kspi0mumu_StrippingSignal = StrippingK0SPi0MuMu("BenderKspi0mumuSignal", Inputs = ["/Event/" +ParticlePath[SAMPLE] + "/Phys/K0s2Pi0MuMuSignalLine/Particles"] ) ### this is for MCTruth matching, for data see the version below
kspi0mumu_StrippingSignal.DST = False
kspi0mumu_StrippingSignal.decayname = "K0S --> pi0mumu"

kspi0mumu_StrippingSignal.COUNTER = {}
kspi0mumu_StrippingSignal.DEBUG = _DEBUG

kspi0mumu_StrippingSignal.COUNTER["Bender(evts) " + kspi0mumu_StrippingSignal.decayname] = 0
## this corresponds to the TES[" "] location of the pars (without the "/Particles" )
kspi0mumu_StrippingSignal.LookIn = ParticlePath[SAMPLE] + "/Phys/K0s2Pi0MuMuSignalLine"#/Particles"#"Phys/Ks2Pi0MuMuSignalLine"
kspi0mumu_StrippingSignal.PROMPT_KS_COUNTER = SIMULATION                                   
algos.append(kspi0mumu_StrippingSignal)

#######################
## for stripping SB 

kspi0mumu_StrippingSideband = StrippingK0SPi0MuMu("BenderKspi0mumuSideband", Inputs =["/Event/" + ParticlePath[SAMPLE] + "/Phys/K0s2Pi0MuMuSidebandLine/Particles"])
kspi0mumu_StrippingSideband.DST = False
kspi0mumu_StrippingSideband.decayname = "K0S --> pi0mumu"

kspi0mumu_StrippingSideband.COUNTER = {}
kspi0mumu_StrippingSideband.DEBUG = _DEBUG
kspi0mumu_StrippingSideband.PROMPT_KS_COUNTER = SIMULATION

kspi0mumu_StrippingSideband.COUNTER["Bender(evts) " + kspi0mumu_StrippingSideband.decayname] = 0
## this corresponds to the TES[" "] location of the pars (without the "/Particles" )
kspi0mumu_StrippingSideband.LookIn = ParticlePath[SAMPLE] + "/Phys/K0s2Pi0MuMuSidebandLine"
algos.append(kspi0mumu_StrippingSideband)


### For Normalization
kspipi = B2QQ("Ks2pipi", Inputs =["/Event/" +  ParticlePath[SAMPLE] + "/Phys/K0s2MuMuNoMuIDLine"*(not _stpnew) + "/Phys/Ks2PiPiForRnSLine" *( _stpnew)])
kspipi.DST = False
kspipi.decayname = "K0S --> pipi"

kspipi.COUNTER = {}
kspipi.DEBUG = _DEBUG
kspipi.PROMPT_KS_COUNTER = SIMULATION
kspipi.COUNTER["Bender(evts) " + kspipi.decayname] = 0
## this corresponds to the TES[" "] location of the pars (without the "/Particles" )
kspipi.LookIn = ParticlePath[SAMPLE] + "/Phys/K0s2MuMuNoMuIDLine"*(not _stpnew) + "/Phys/Ks2PiPiForRnSLine" *( _stpnew)
algos.append(kspipi)


### Signal as partially reconstructed dimuon (effectively a V0)
kspzmmv0 = B2QQ("Ks2pizeromm_as_V0", Inputs =["/Event/" +  ParticlePath[SAMPLE] + "/Phys/TriggerTestLine/Particles"])
kspzmmv0.DST = False
kspzmmv0.decayname = "K0S --> (pi0) mu mu"

kspzmmv0.COUNTER = {}
kspzmmv0.DEBUG = _DEBUG
kspzmmv0.PROMPT_KS_COUNTER = SIMULATION
kspzmmv0.COUNTER["Bender(evts) " + kspipi.decayname] = 0
## this corresponds to the TES[" "] location of the pars (without the "/Particles" )
kspzmmv0.LookIn = ParticlePath[SAMPLE] + "/Phys/TriggerTestLine"
algos.append(kspzmmv0)


#######################
TES = gaudi.evtsvc()
for algo in algos: gaudi.addAlgorithm(algo)

gaudi.initialize()

##############
## add your custom functions
for algo in algos:
    algo.extraFunctions = [trackHits, muIDDetailedInfo, more_muon_things]#, GiacomoStuff]#, TisTosInfo]#, triggerBlock, triggerBlockDaughters]
    # if data, set this to false
    algo.MC_INFO = SIMULATION
    ###
    if algo.LookIn[0:2] == "/P": algo.LookIn = algo.LookIn.replace("/P","P")
    algo.TRIGGER = TRIGGER
    ## custom code, just ignore
    algo.NTupleLUN = "T"
    algo.addedKeys = []
    algo.TUP = 1
    algo.dataFormat = dataFormat
    algo.COUNTER['negSq'] = 0
    algo.COUNTER["weird"] = 0
    algo.COUNTER["EVT"] = 0
    algo.COUNTER["Sel"] = 0
    algo.COUNTER["MuMuCouplesAnalyzed"] = 0

    algo.evt_of = 1e6

    ##algo.l0BankDecoder=algo.tool(cpp.IL0DUFromRawTool,'L0DUFromRawTool')
    ##algo.rawBankDecoder=algo.tool(cpp.IOTRawBankDecoder,'OTRawBankDecoder')
    algo.PVRefitter = algo.tool(cpp.IPVReFitter,"AdaptivePVReFitter")
    algo.LifeFitter = algo.tool(cpp.ILifetimeFitter,"PropertimeFitter")
    algo.Geom =  algo.tool(cpp.IDistanceCalculator,"LoKi::DistanceCalculator")
    ##algo.mid = algo.tool(cpp.IMuonIDTool, "MuonIDPlusTool") ML commented this
    
    if algo.MC_INFO : algo.matcher = gaudi.toolsvc().create("MCMatchObjP2MCRelator",interface="IP2MCP")
    if RootInTES: algo.RootInTES = RootInTES

if kspipi.MC_INFO: kspipi.extraFunctions += [BQQMCtruth]
if kspi0mumu_StrippingSignal.MC_INFO: kspi0mumu_StrippingSignal.extraFunctions += [KsPi0MuMuMCtruth]
if kspi0mumu_StrippingSideband.MC_INFO: kspi0mumu_StrippingSideband.extraFunctions += [KsPi0MuMuMCtruth]
if kspzmmv0.MC_INFO: kspzmmv0.extraFunctions += [BQQMCtruth]
gaudi.initialize()

if INTERACTIVE:
    gaudi.run(INTERACTIVE)
    TES = gaudi.evtsvc()
    DET = gaudi.detsvc()
    ## ML: you have to do stop and finalize
# 
else:
    gaudi.run(-1)
    gaudi.stop()
    gaudi.finalize()
#for i in range(INTERACTIVE):
  #  try: gaudi.run(1)
 #   except: pass
    #TES = gaudi.evtsvc()
#    DET = gaudi.detsvc()
#    sensor = DET["/dd/Structure/LHCb/BeforeMagnetRegion/VP/VPRight/Module47WithSupport/Module47/Ladder3/Sensor3"]
   # if sys.argv[-1]:
    #    gaudi.stop()
    #    gaudi.finalize()
    

#else:
#    gaudi.run(-1)
#    gaudi.stop()
#    gaudi.finalize()

## gaudi.run(900)
## for i in range(100):
##     gaudi.run(1)
##     TES = gaudi.evtsvc()
##     print i, '=============='
##     if i == 28:BREAK
##     #TES.dump()
##     #if TES["/Event/" + ParticlePath[SAMPLE] + "/Phys/K0s2Pi0MuMuSignalLine/Particles"]: BREAK


