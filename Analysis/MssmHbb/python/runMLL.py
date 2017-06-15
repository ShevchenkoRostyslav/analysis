#!/usr/bin/python

'''Macro to run combine tool with MLL fit.
    
'''

import os,sys, shutil
from subprocess import Popen

__author__ = "Rostyslav Shevchenko"
__maintainer__ = "Rostyslav Shevchenko"
__email__ = "rostyslav.shevchenko@desy.de"

def checkInput(datacard_folder):
    '''Method to check whether folder exists or not.
    
    '''
    if not os.path.exists(datacard_folder):
        raise AttributeError('Wrong input folder path: ' + datacard_folder)

def runCombineTool(combine,data_card,dir_path,m):
    '''Function to run combine for paticular mass.
    
    '''
    os.chdir(dir_path)
    proc = Popen(combine + ' -m ' + m + ' ' + data_card + ' > log_' + m + '.o',shell=True)
    os.chdir("../")
    
# def runCombineTool(combine,datacards,dirs,masses):
#     '''Macro to run combine tool.
#     
#     '''
#     for card,directory,mass in zip(datacards,dirs,masses):
#         os.chdir(directory)
#         proc = Popen(combine + ' -m ' + mass + ' ' + card + ' > log_' + mass + '.o',shell=True)
#         os.chdir('../')

def copyDatacardsToDir(datacards,dirs):
    '''Function to copy datacards to corresponded fodlers.
    
    '''
    for card,directory in zip(datacards,dirs):
        shutil.copy(card, directory)

def copyDatacardToDir(datacard,directory):
    '''Function to copy datacard to directory.
    
    '''
    shutil.copy(datacard,directory)

def createCleanDirs(dir_list):
    '''Function to create empty dirs from the list.
    
    '''
    for d in dir_list:
        print d
        createCleanDir(d)

def createCleanDir(dir_name):
    '''Function to create a clean dir.
    
    '''
    try:
        if os.path.exists(dir_name):
            shutil.rmtree(dir_name)
        os.mkdir(dir_name)
    except Exception as e:
        print(e)
        

def AdjustSignalStrength(mass, bg_only):
    '''Funciton to Adjust signla strength accroding to the mass point
    and bg_only fit or not.
    
    '''
    rmin = '-10'
    rmax = '10'
    if(bg_only):
        rmin = '-10'#'-0.001'
        rmax = '10'#'0.001'
    else:
        if mass == '300' or mass == '350':
            rmin = '-10'
            rmax = '20'
        elif mass == '400':
            rmin = '-10'
            rmax = '20'
        else:
            rmin = '-20'
            rmax = '20'
    
    r_string = '--rMin=' + rmin + ' --rMax=' + rmax
    return r_string
            

if __name__ == '__main__':

    #bg_only fit?
    bg_only = True
    #working directory with datacards and stored output:
    datacard_folder = '/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_8_0_20_patch1/src/Analysis/MssmHbb/datacards/201706/14/switch_sub_ranges/mll_forBias/'
    checkInput(datacard_folder)
    os.chdir(datacard_folder)
    #list of mass points
    mass = ['500','600','1100']
    #mass = ['300','350','400','500','600','700','900','1100','1300']
    #combine preferences to be added
    combine_add = ''#'--saveWorkspace '#'--plots --saveShapes --saveWithUncertainties -t -1'#'--freezeNuisanceGroups signal --expectSignal=0'# --freezeNuisances CMS_bkgd_qcd_13TeV'
    for m in mass:
        dir_path =  datacard_folder + '/mll_M-' + m
        data_card = 'hbb_mbb' + m + '_mssm-13TeV.txt'
        #create clean directory for particular mass:
        createCleanDir(dir_path)
        #copy data card to the corresponded dir:
        copyDatacardToDir(data_card,dir_path)
        #adjust combine command
        if m == '300' or m == '350':# or m == '500':
            combine = 'combine -M MaxLikelihoodFit --robustFit 1 --minimizerAlgoForMinos Minuit2,Migrad --cminFallbackAlgo "Minuit2,Minimize,0:0.1" --cminOldRobustMinimize 0 -v 5'
        else:
            combine = 'combine -M MaxLikelihoodFit --robustFit 1 -v 5'
        combine += ' ' + combine_add
        #special emulation for the background only fit
        if(bg_only): combine +=  ' --freezeNuisanceGroups signal'
        #add r-limits
        r_string = AdjustSignalStrength(m,bg_only)
        combine += ' ' + r_string
        #run combination tool
        print (combine)
        runCombineTool(combine,data_card,dir_path,m)
