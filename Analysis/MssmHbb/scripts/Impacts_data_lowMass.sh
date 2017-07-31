#!/bin/sh
# $1 mass
# $2 signal strength range

# create workspace
rm workspace_${1}.root
combineTool.py -M T2W -i hbb_mbb$1_mssm-13TeV.txt -o workspace_${1}.root

# perform initial fit
combineTool.py -M Impacts -d workspace_${1}.root -m $1 --robustFit 1 --doInitialFit --rMin -$2 --rMax $2

# scan systematic nuisances
combineTool.py -M Impacts --named CMS_lumi_13TeV,CMS_BTagEff_13TeV,SFb -d workspace_${1}.root -m $1 --robustFit 1 --doFits --rMin -$2 --rMax $2
combineTool.py -M Impacts --named CMS_PU_13TeV,JES,JER,CMS_SFl_13TeV,PtEff -d workspace_${1}.root -m $1 --robustFit 1 --doFits --rMin -$2 --rMax $2 --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER
combineTool.py -M Impacts --named CMS_lumi_13TeV,CMS_SFl_13TeV,CMS_PU_13TeV,CMS_BTagEff_13TeV,JES,JER,PtEff,SFb -m ${1} -o impacts_${1}_sys.json -d workspace_${1}.root

# scan function nuisances
combineTool.py -M Impacts --named width,par4,bias,tail,peak -d workspace_${1}.root -m $1 --robustFit 1 --doFits --rMin -$2 --rMax $2 --minimizerAlgoForMinos Minuit2,Migrad --cminFallbackAlgo "Minuit2,Minimize,0:0.2" --cminOldRobustMinimize 0 
combineTool.py -M Impacts --named width,par4,bias,tail,peak -m ${1} -o impacts_${1}_func.json -d workspace_${1}.root

# plot impacts (separate pdf files for systematic and function nuisances)
plotImpacts.py -i impacts_${1}_sys.json -o impacts_${1}_sys
plotImpacts.py -i impacts_${1}_func.json -o impacts_${1}_func
