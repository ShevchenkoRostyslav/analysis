'''COnfig file with preferences to run and merge Bias tests.

'''

__author__ = "Rostyslav Shevchenko"
__maintainer__ = "Rostyslav Shevchenko"
__email__ = "rostyslav.shevchenko@desy.de"

#list of mass points
mass = [500,600,700,900]#[300,350,400,500]#[1100,1300]#,[600]#[300]#[300,350,400,500,600,700,900,1100,1300]#
#number of toys
n_toys = 5000
#to speedup
n_jobs = 500
#signal strength to be uaed
signal_strengths = [0,1,2,4]
#sanity check falg. If yes - alternative = pdf will be used
sanity_check = False
#bg only?
bg_only = False
#fix turn-on or not
fix_turnon = False
#Additional out name to be used in hadd -f 
additinoal_out_name = ''

#functions to be used
pdf = { 'sr1' : 'extnovoeffprod_G4_R1', 'sr2' : 'novosibirsk_G4_R2', 'sr3' : 'novosibirsk_G4_R3' }
#pdf = { 'sr1' : 'supernovoeffprod3_G4_R1', 'sr2' : 'novosibirsk_G4_R2', 'sr3' : 'novosibirsk_G4_R3'}
# pdf = { 'sr1' : 'extnovoeffprod_G4_R1', 'sr2' : 'supernovosibirsk1_G4_R2', 'sr3' : 'supernovosibirsk1_G4_R3' }
#alternative = { 'sr1' : 'superdijeteffprod2_G4_R1', 'sr2' : 'dijetexp_G4_R2', 'sr3' : 'dijetexp_G4_R3'}
alternative = { 'sr1' : 'berneffprod8par_G4_R1', 'sr2' : 'exppolynom_G4_R2', 'sr3' : 'dijetexp_G4_R3'}
# alternative = { 'sr1' : 'berneffprod8par_G4_R1', 'sr2' : 'dijetexp_G4_R2', 'sr3' : 'dijetexp_G4_R3'}
#alternative = { 'sr1' : 'berneffprod9par_G4_R1', 'sr2' : 'dijetexp_G4_R2', 'sr3' : 'dijetexp_G4_R3'}
#alternative = { 'sr1' : 'supernovoeffprod2_G4_R1', 'sr2' : 'novosibirsk_G4_R2', 'sr3' : 'novosibirsk_G4_R3'}
# alternative = { 'sr1' : 'extnovoeffprod_G4_R1', 'sr2' : 'novosibirsk_G4_R2', 'sr3' : 'novosibirsk_G4_R3'}

#pdf = { 'sr1' : 'extnovoeffprod_R1'}
#alternative = { 'sr1' : 'berneffprod8par_R1'} #berneffprod8par_R1 #superdijeteffprod2_R1

#number of events accroding to the sub-ranges
n_observed = {'sr1' : '259399', 'sr2' : '105053', 'sr3' : '26760'}
