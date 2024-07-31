# FullHGCALSS_Setup
Analyzer script for HGCAL-SS setup
Descrition of scripts:
AnalyzeHGCOctTB.cc : Main analysis code AnalyzeHGCOctTB.h : Initialize histos here HGCNtupleVariables.h: Tree variable initialization Others are helping classes
All runlist text files are under filelist directory

## Create Skimmed files with removing noisy hits
```
git clone https://github.com/alpana-hep/FullHGCALSS_Setup.git .
make
```

Example to run code -
```
./analyzeHGCOctTB <infile.txt> <outputfile.root> <energy 20>
```
or
```
./analyzeHGCOctTB infile.txt out.root 20
```

## Analyzing and making input checks plots
