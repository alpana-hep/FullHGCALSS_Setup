#!/bin/bash                                                                                                                                            
cd /home/kalpana/t3store3/public/FullHGcal_simulation/StandAloneSimulation/CMSSW_12_0_0/src/Samples/AnalyzerCode

source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`


echo $PWD

./analyzeHGCOctTB  /home/kalpana/t3store3/public/FullHGcal_simulation/StandAloneSimulation/CMSSW_12_0_0/src/Samples/AnalyzerCode/file_list.txt /home/kalpana/t3store3/public/FullHGcal_simulation/StandAloneSimulation/CMSSW_12_0_0/src/Samples/AnalyzerCode/out.root 20
