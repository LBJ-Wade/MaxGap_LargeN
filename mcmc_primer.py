import numpy as np
from data_analy import *
import pymultinest
import os,shutil

path = os.getcwd()

class MCMC_Runner(object):
    # Implicit assumption that datafile has been mapped into unit cube
    
    def __init__(self, file_name, output_name='1-', n_live_points=100, evid_tol=0.1,
                 sample_eff=0.3, resume=False):
        self.file_name = file_name
        self.analy = Data_Input(self.file_name)
        chains_root = 'chains/'
        
        self.pnest_pars = {}
        self.pnest_pars['n_live_points'] = n_live_points
        self.pnest_pars['evidence_tolerance'] = evid_tol
        self.pnest_pars['sampling_efficiency'] = sample_eff
        self.pnest_pars['resume'] = resume
        self.pnest_pars['outputfiles_basename'] = output_name

        self.chain_path = path + '/' + chains_root
        return

    def likelihood_val(self, cube, ndim, nparams):
        box_size = self.analy.solve_coords([cube[0], cube[1], cube[2]])
        return -box_size
    
    def prior(self, cube, ndim, nparams):
        return cube

    def run_scan(self):
        if os.path.exists(self.chain_path):
            shutil.rmtree(self.chain_path)
        os.makedirs(self.chain_path)
        os.chdir(self.chain_path)
        
        pymultinest.run(self.likelihood_val, self.prior, 3, **self.pnest_pars)
        os.chdir(path)
        return
        
    def find_sides_BF(self):
        samples = np.loadtxt(self.chains_path + output_name + 'post_equal_weights.dat')
        posterior = samples[:,-1]
        max_index = posterior.argmax()
        point = samples[max_index, :-1]
        box = self.analy.find_box(point)
        print 'Volume: ', box[0]
        print 'Sides: ', box[1]
        return

    def visualize(self):

        return
