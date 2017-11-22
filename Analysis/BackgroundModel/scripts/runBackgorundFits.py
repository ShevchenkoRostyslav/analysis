#!/usr/bin/python

'''Macro to run Bg fits for 3 sub-ranges and SR and CR.

'''

import os,sys, shutil
from subprocess import Popen

__author__ = "Rostyslav Shevchenko"
__maintainer__ = "Rostyslav Shevchenko"
__email__ = "rostyslav.shevchenko@desy.de"

if __name__ == '__main__':
	#cmssw base dir
	cmsswBase = os.environ["CMSSW_BASE"]
	#output folder
	output_f = cmsswBase + '/src/Analysis/BackgroundModel/test/finale_files_LinearErrors/'
	#regions to run
	regions = ['CR','SR']
	#preferences for regions
	data = {'CR' : 'TripleBTagReverseSelectionBtoH2016_prescale_13TeV_G4.root', 'SR' : 'TripleBTagSelectionBtoH2016_13TeV.root'}
	#preferences for sub-ranges
	subranges = ['SR1','SR2','SR3']
	preferences_sr = {'SR1_xmin' : '200', 'SR1_xmax' : '650', 'SR2_xmin' : '350', 'SR2_xmax' : '1190', 'SR3_xmin' : '500', 'SR3_xmax' : '1700', 'SR1_bins' : '45', 'SR2_bins' : '42', 'SR3_bins' : '48', 'SR1_SR_func' : 'extnovoefffixprod', 'SR2_SR_func' : 'novosibirsk', 'SR3_SR_func' : 'novosibirsk', 'SR1_CR_func' : 'extnovoeffprod', 'SR2_CR_func' : 'novosibirsk', 'SR3_CR_func' : 'novosibirsk'}
	#param_modifs   = {'CR_SR1' : '','CR_SR2' : '"peak1: max=500" "tail: max=1, min=-1"', 'CR_SR3' : '', 'SR_SR1' : '', 'SR_SR2' : '', 'SR_SR3' : '"peak1: start=462, min=350, max=550" "width1: min=45, max=65" "tail: start=-0.94"'}
	param_modifs   = {'CR_SR1' : '','CR_SR2' : '"peak1: max=500" "tail: max=1, min=-1"', 'CR_SR3' : '', 'SR_SR1' : '', 'SR_SR2' : '', 'SR_SR3' : ''}
        #param_modifs   = {'CR_SR2' : '"peak1: start=50, min=10, max=70" "width1: start=35.7" "tail: start=-0.8"', 'CR_SR3' : '"peak1: start=425.236, min=400, max=440" "tail: start=-0.89" "width1: start=60.57"', 'CR_SR1' : '', 'SR_SR1' : '', 'SR_SR2' : '', 'SR_SR3' : '"peak1: start=461.819, min=400, max=500" "width1: start=58.3072, min=50, max=65" "tail: start=-0.93615"'}
	for r in regions:
		command = 'FitBackgroundModel -t data/2016DataRereco_05_01_2017/' + data[r] + ' -o ' + output_f + r + '/'
		for sr in subranges:
			command2 = command
			command2 += preferences_sr[sr + '_' + r + '_func'] + '_' + preferences_sr[sr + '_xmin'] + 'to' + preferences_sr[sr + '_xmax']
			if r == 'CR': command2 += ' --control_region 1'
			else: command2 += ' --control_region 0'
			command2 += ' --fit_min ' + preferences_sr[sr + '_xmin'] + ' --fit_max ' + preferences_sr[sr + '_xmax'] + ' --nBinsPlot ' + preferences_sr[sr + '_bins'] + ' --nbins 5000'# + preferences_sr[sr + '_bins']
			command2 += ' -b ' + preferences_sr[sr + '_' + r +'_func']
			if(param_modifs[r + '_' + sr] != ''): command2 += ' -m ' + param_modifs[r + '_' + sr]
			print command2
			proc = Popen(command2,shell=True)
	
