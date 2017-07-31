#!/usr/bin/python

'''Macro to run combine tool with Asymptotic method.
    
'''

import os,sys, shutil
from subprocess import Popen

__author__ = "Rostyslav Shevchenko"
__maintainer__ = "Rostyslav Shevchenko"
__email__ = "rostyslav.shevchenko@desy.de"

def writeHbbLimits(mass):
    '''Function to write Hbb.limits file.
    
    '''
    f = open("Hbb.limits","a")
    #file name:
    result = os.path.abspath("./") + '/higgsCombineHbb.Asymptotic.mH' + mass + '.root\n'
    #write to the file
    f.write(result)
    
def checkInput(datacard_folder):
    '''Method to check whether folder exists or not.
    
    '''
    if not os.path.exists(datacard_folder):
        raise AttributeError('Wrong input folder path: ' + datacard_folder)

def runCombineTool(combine,data_card,m,log_file = 'log_asymptotic'):
    '''Function to run combine for paticular mass.
    
    '''
    proc = Popen(combine + ' -m ' + m + ' ' + data_card + ' > ' + log_file + '_' + m + '.o',shell=True)

if __name__ == '__main__':

    #working directory with datacards and stored output:
    datacard_folder = '/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_8_0_20_patch1/src/Analysis/MssmHbb/datacards/201707/26/blinded/mssm/bias/'
    checkInput(datacard_folder)
    os.chdir(datacard_folder)
    #if combination of 7+8+13 is performed:
    combination = False
    #Clean Hbb.Limits file if exists
    if os.path.exists("Hbb.limits"): os.remove("Hbb.limits")
    #list of mass points
    mass = ['300','350','400','500','600','700','900','1100','1300']
#     mass = ['1100']
    #combine preferences to be added
    combine_add = ''
    #values of rMin and rMax to be used
    rMin = '-20'
    rMax = '20'
    for m in mass:
        data_card = 'hbb_mbb' + m + '_mssm-13TeV.txt'
        if combination:
            data_card = 'Comb_7TeV_8TeV_HIG16018_13TeV_Packed_M' + m + '_card.txt'
            rMin = '-30'
            rMax = '40'

        #adjust combine command
        if m == '300' or m == '350':
            if m == '350' and combination:
                rMin = '-20'
                rMax = '30'
            if m == '300' and 'Bernstein' in datacard_folder:
                rMin = '-30'
                rMax = '30'
            if m == '300':
                rMin = '-10'
                rMax = '10'
            if m == '350':
                rMax = '20'
                rMin = '-10'
            combine = 'combine -M Asymptotic -n Hbb --rMin=' + rMin + ' --rMax=' + rMax + ' --minimizerStrategy 0 --cminFallbackAlgo "Minuit2,Minimize,0:0.1" --cminOldRobustMinimize 0 -v5'# -t -1 --expectSignal=0' 
        else:
            combine = 'combine -M Asymptotic -n Hbb --rMin=' + rMin + ' --rMax=' + rMax + ' -v 5'# -t -1 --expectSignal=0'
        combine += ' ' + combine_add
        #run combination tool
        runCombineTool(combine,data_card,m)
        #write file with results to Hbb.Limits
        writeHbbLimits(m)
