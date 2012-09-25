#!/bin/bash     


export model0="xww_SMHiggsWW_10fb.txt"
export model0name="SMHiggsWW"

export model1="xww_TWW_2mplus_10fb.txt"
export model1name="TWW_2mplus"

#
# generate toys according to model0 / model 1
#

lands.exe -d $model0 -M Hybrid  -m 125 --minuitSTRATEGY 0 --bMultiSigProcShareSamePDF -L $CMSSW_BASE/lib/*/libHiggsAnalysisCombinedLimit.so --bWriteToys 1 -n "$model0name" --nToysForCLsb 1000 --nToysForCLb 1 --singlePoint 1 --seed 12344  -rMin 0 -rMax 10  --freq

lands.exe -d $model1 -M Hybrid  -m 125 --minuitSTRATEGY 0 --bMultiSigProcShareSamePDF -L $CMSSW_BASE/lib/*/libHiggsAnalysisCombinedLimit.so --bWriteToys 1 -n "$model1name" --nToysForCLsb 1000 --nToysForCLb 1 --singlePoint 1 --seed 12344  -rMin 0 -rMax 10  --freq


#
# Fit the toys for different models respectively
# 

lands.exe -d $model0 -M MaxLikelihoodFit -rMin 0 -rMax 10 -m 125 --NoErrorEstimate --minuitSTRATEGY 0 --bMultiSigProcShareSamePDF -L $CMSSW_BASE/lib/*/libHiggsAnalysisCombinedLimit.so  --doExpectation 1 --loadToysFromFile ${model0name}_PseudoData_sb_seed12344.root  -n "LL_toy0_model0"

lands.exe -d $model1 -M MaxLikelihoodFit -rMin 0 -rMax 10 -m 125 --NoErrorEstimate --minuitSTRATEGY 0 --bMultiSigProcShareSamePDF -L $CMSSW_BASE/lib/*/libHiggsAnalysisCombinedLimit.so  --doExpectation 1 --loadToysFromFile ${model0name}_PseudoData_sb_seed12344.root  -n "LL_toy0_model1"

lands.exe -d $model0 -M MaxLikelihoodFit -rMin 0 -rMax 10 -m 125 --NoErrorEstimate --minuitSTRATEGY 0 --bMultiSigProcShareSamePDF -L $CMSSW_BASE/lib/*/libHiggsAnalysisCombinedLimit.so  --doExpectation 1 --loadToysFromFile ${model1name}_PseudoData_sb_seed12344.root  -n "LL_toy1_model0"

lands.exe -d $model1 -M MaxLikelihoodFit -rMin 0 -rMax 10 -m 125 --NoErrorEstimate --minuitSTRATEGY 0 --bMultiSigProcShareSamePDF -L $CMSSW_BASE/lib/*/libHiggsAnalysisCombinedLimit.so  --doExpectation 1 --loadToysFromFile ${model1name}_PseudoData_sb_seed12344.root  -n "LL_toy1_model1"


root -l hypoSeparation.C\(\"$model0name\",\"LL_toy0_model0_maxlls_tree.root\",\"LL_toy0_model1_maxlls_tree.root\",\"$model1name\",\"LL_toy1_model0_maxlls_tree.root\",\"LL_toy1_model1_maxlls_tree.root\"\)
