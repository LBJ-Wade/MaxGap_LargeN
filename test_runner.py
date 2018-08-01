from data_analy import *
from mcmc_primer import *

# Test data_analy class
#A_Tool = Data_Input('fake_data/test_10000pts.dat', save_info=True)
#A_Tool.solve_coords([0.5,.5,.5])

# Test MCMC
mcmcR = MCMC_Runner('fake_data/test_10000pts.dat')
mcmcR.run_scan()
