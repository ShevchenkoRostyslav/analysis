#!/usr/bin/python

from ROOT import TTree,TFile,TH3F
from array import array
import math

def TranslateCosB_To_sinA(cos_bins):
    """Function to trasfer cos(beta-alpha) binning in conv. B to sin(beta-alpha) in conv. A.

    convention B: 0 < beta-alpha < Pi
    convention A: -90 < beta-alpha < 90
    """
    sin_beta_alpha_A = []
    for p in cos_bins:
        #calculate the sin
        s_beta_alpha_B = math.sqrt(1 - p**2)
        #translate convention B to A
        s_beta_alpha_A = s_beta_alpha_B
        if p < 0:
            s_beta_alpha_A *= -1
        s_beta_alpha_A = round(s_beta_alpha_A,5)
        sin_beta_alpha_A.append(s_beta_alpha_A)
    return sin_beta_alpha_A

def MakeBins(l):
#    print l
    L = array('f')
    L.append(l[0] - (l[1]-l[0])/2)
    #make bin boundaries
    for i in range(len(l)):
#        print l[i], ' len = ', len(l), 'L = ', l[i]+(l[i]-L[-1])
#        print l[i], L[-1]
        L.append(l[i]+(l[i]-L[-1]))
    #check if bin boundaries are in increasing order
    for i in range(len(L)-1):
#        print L[i], L[i+1], len(L)
        if L[i+1]-L[i] <= 0:
            raise RuntimeError("Check bin boundaries!!!")

    return L

#tanBetas
# tanBetaRange = [1.,100.]
# tanBetaStep = 1
# tanBetas = [ tanBetaRange[0]+tanBetaStep*i for i in range(int((tanBetaRange[1]-tanBetaRange[0])/tanBetaStep+1)) ]
#tanBetas
tanBetaLowRange = [0.5,2]
tanBetaLowStep  = 0.05
tanBetasLow = [ round(tanBetaLowRange[0]+tanBetaLowStep*i,3) for i in range(int((tanBetaLowRange[1]-tanBetaLowRange[0])/tanBetaLowStep)) ]
tanBetaRange = [2.,100.]
tanBetaStep = 0.5
tanBetas = [ tanBetaRange[0]+tanBetaStep*i for i in range(int((tanBetaRange[1]-tanBetaRange[0])/tanBetaStep+1)) ]
tanBetas = tanBetasLow + tanBetas

#H masses
mH = [200,250,300,350,400,500,600,700,800,900,1100,1300]
#h masses
mh = [125.]
#sin(beta - alpha)
cosB_A_Range = [-1,1]
cosB_A_step = 0.005
cosB_As = [cosB_A_Range[0] + cosB_A_step* i for i in range(int((cosB_A_Range[1] - cosB_A_Range[0])/cosB_A_step + 1))]
sinB_As = TranslateCosB_To_sinA(cosB_As)

type_boson = 'type1_mH'
production = 'production_cosB_A_-1_1_tanB_0p5-100_COMBINATION'
inputTxt = '/nfs/dust/cms/user/shevchen/SusHiScaner/output/' + production + '/txtFiles/' + type_boson + '.txt'

# names of the column
lines = [ line.split() for line in open(inputTxt, 'r')]
histNames = lines[0]
print histNames
mH_bin = MakeBins(mH)
tanBetas_bin = MakeBins(tanBetas)
# sinB_As_bin = MakeBins(sinB_As)
cosB_As_bin = MakeBins(cosB_As)

print 'TAN BETAs:'
print(tanBetas)
print 'cos(beta-alpha)s:'
print sinB_As

