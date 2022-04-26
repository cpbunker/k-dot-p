'''
Created on Apr 17, 2020

@author: Christian
'''

import numpy as np
import scipy as sp

from Valleytronics import hamiltonian
from Valleytronics import GetEnergyLevels
from Valleytronics import EnergyLevels


                    
################################################################################
# define condition functions for the param space search

def is_degenerate(levels, dummy1, dummy2, dummy3):   
    '''
    Condition function that returns true/false depending on whether the energy levels
    satisfy the relevant condition at this point in param space
    
    Condition here is that any 2 energy levels are degenerate
    
    :param levels: 1d np array of the energy levels at this point in param space
    :param dummy1: dummy argument so that this (test) conditon will have same args as (real) condition
            SpinInversion, so code can move back and forth between them easily
    :param dummy2:
    :param dummy3:
    returns bool
    '''
    
    # check whether any of the levels have the same value
    # iter over levels
    for i in range(len(levels)):
        
        # now check equality with all of the higher levels
        # to avoid checking against self or checking twice needlessly
        for j in range(i+1, len(levels) ):
            
            # only need one equality to meet condition
            if( abs(levels[i] - levels[j]) < 0.001 ):
                return True;
            
    # we only make it here if there is no degeneracy in the levels
    return False;

    #### end is degenerate
                 
                    
def E0_highest(levels, dummy1, dummy2, dummy3 ):  
    '''
    Condition function that returns true/false depending on whether the energy levels
    satisfy the relevant condition at this point in param space
    
    Condition here is that E0 has highest energy
    
    :param levels: 1d np array of the energy levels at this point in param space
    :param dummy1: dummy argument so that this (test) conditon will have same args as (real) condition
            SpinInversion, so code can move back and forth between them easily
    :param dummy2:
    :param dummy3:
    returns bool
    ''' 
    
    # get E0
    E0 = levels[0]
    
    # check against maximum of the levels
    if( E0 == np.amax(levels) ): # E0 is the maximum
        return True;
    
    else: # is not the maximum
        return False;    
    
    #### end E0 highest 
    
    
def SpinInversion(levels, bench_levels, phi, H):
    '''
    Gets the spin ordering in ascending order for levels at a certain point in param space. Then checks if the 
    band spin ordering has switched for each band from what it was at the point specified by bench_levels
    
    :param levels, 1d array of energy levels at the test point in param space
    :param bench_levels, 1d array of energy levels at the benchmark point in param space
    :param phi, double, angle in rads between magnetization and z hat
    
    returns, bool according to whether condition was met
    '''
    
    # get the labels associated with the (in-array) level ordering from bench_levels
    labels, dummy_colors = EnergyLevels.GetBandLabels(H, bench_levels, phi);
    
    # get the energy ordering of the levels
    # specifically we need indices that sort the levels so we can access associated labels
    indices = np.argsort(levels);
    
    # go thru indices list, which accesses labels in ascending (energy) order
    # look for up/ down, up/down spin ordering
    if( 'uparrow' in labels[indices[1]] and 'downarrow' in labels[indices[0]] ): # bottom two are inverted
        if('uparrow' in labels[indices[3]] and 'downarrow' in labels[indices[2]] ): # top two are inverted
            
            # spin ordering is inverted as we want it to be
            return True;
        
    # if we reach here, no inversion
    return False;

    #### end SpinInversion
    
    
def Switch(levels, dummy1, dummy2, dummy3):
    '''
    Same idea as spin inversion, but removes spin from the equation since I'm not sure how to denote spin bands.
    Instead looks for the bottom two bands (ie valence at the test point) having opposite order as at the test point,
    ie the lower valence band is now at a higher energy. Likewise for the conduction bands.
    
    :param levels: 1d array of energy levels at this point in param space
    :param dummy1: dummy argument so that this (test) conditon will have same args as (real) condition
            SpinInversion, so code can move back and forth between them easily
    :param dummy2:
    :param dummy3:
    '''
    
    # unpack the energy levels in terms of their ordering in the array
    # this correspons to their energy orering at the benchmark point
    E0, E1, E2, E3 = levels[0], levels[1], levels[2], levels[3];
    
    # first condition is that the (originally) lower valence band is above the (originally) higher one
    if( E0 > E1 ):
        
        # check the same condition for the conduction band
        if(E2 > E3 ):
            
            # if we get here, condition is met
            return True;
      
    # otherwise it is not  
    return False;
