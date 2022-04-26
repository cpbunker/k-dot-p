'''
Created on Mar 29, 2020

@author: Christian

energylevels.py was a module for everything related to finding and plotting or otherwise visualizing the energy levels of the 
physical system defined by the combined Hamiltonian and Exchange matrix objects defined in the hamiltonian.py program. I have
since split this module into this, EnergyLevels.py, and GetEnergyLevels.py. This module contains the more complicated
functions that are often called by the functions of GetEnergyLevels. Because of this it is another level removed from the user.

In other words, if you are looking for a problem:
    FIRST check the runner file (PbSnSe or equivalent)
    SECOND check GetEnergyLevels
    ONLY THEN look in this module

The functions in this module:
    Determine for the energy levels of the system throughout the parameter space defined by a1, a2, b1, b2
    Look for specific energy level behavior (ie spin inversion) in this parameter space. The specific conditions are
    defined in the module conditions.py
    Explore anisotropy effects by finding energy levels for different values of phi (phi = 0 --> central valley, 
    phi = 123 rads --> oblique valley)
'''

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.collections as collections

from Valleytronics import hamiltonian



################################################################################
# Functions to plot energy levels


def PlotLevel2D(H_inputs, material, phi, free1, free1_limits, free2, free2_limits, fixed1, fixed1_val, fixed2, fixed2_val):
    '''
    **** Note: This function is not currently used anywhere in the code. It works but it is not useful and has nott ben kept up to date
    
    Plots one of the energy levels of the hamiltonian determined by H_inputs along with the exchange
    matrix which is varied according to the parameters free1 and free2 while parameters fixed 1 and
    fixed2 are held constant. Note that the exchange matrix E=E(a1, a2, b1, b2) and we should always
    follow this order when assigning the free and fixed params. i.e. free1 = a1 and free2 = a2 but
    not vice versa
    
    The resulting plot is 3d and shows the energy level selected by level against the two parameters.
    Should be color coded to show when level inversion occurs
    
    :param H_inputs: tuple of the 5 params needed for the bulk hamiltonian H
                    (Eg, muBohr, Bvec, vt, vl)
    :param phi: double, angle between H field and valley
    :param free1: str, tells which param is the first free param
    :param free1_limits: tuple, gives upper and lower limits of free1 param
    :param free2: see free1
    :param free2_limits: see free1
    :param fixed1: str, tells which param is first fixed
    :param fixed1_val: double, tells value of first fixed
    :param fixed2: see fixed1
    :param fixed2_val: see fixed1
    '''
    
    # get the k.p bulk hamiltonian
    H = hamiltonian.Hamiltonian(*H_inputs);
    
    # define number of points to use in parameter space
    N_param = 100
    N_levels = H.N
    del H;
    
    # iter over the free parameters
    # use linspace to get vals in the limits of free parameter 1 and 2
    free_arr1 = np.linspace(free1_limits[0], free1_limits[1], N_param);
    free_arr2 = np.linspace(free2_limits[0], free2_limits[1], N_param);
    
    # init container variable
    E_levels = np.zeros((N_levels, N_param, N_param ));
    
    # iter over the free parameters
    for i in range(len(free_arr1)):
        
        # similarly iter over the second free parameter
        for j in range(len(free_arr2)):
            
            # assign vals to free params
            val1 = free_arr1[i];
            val2 = free_arr2[j];
            
            # get input params for the exchange matrix
            
            # check to make sure the fixed params were chosen correctly
            fixed_error = "Fixed params fixed1 = "+fixed1+" and fixed2 = "+fixed2;
            fixed_error += " were chosen incorrectly. Make sure your order obeys ";
            fixed_error += "a1, a2, b1, b2";
            if( fixed1 == "a1"):
                
                pass; # no room for error here
            
            elif( fixed1 == "a2"): # fixed2 cannot be before
                
                if( fixed2 == "a1"):                    
                    print(fixed_error);
                    return;
                    
            elif( fixed1 == "b1"): # fixed2 cannot be before
                
                if( fixed2 == "a1"): 
                    print(fixed_error);
                    return;
                    
                elif( fixed2 == "a2"):
                    print(fixed_error);
                    return;
                    
            elif( fixed1 == "b2"): # fixed2 cannot be before
                
                if( fixed2 == "a1"):
                    print(fixed_error);
                    return;
                    
                elif( fixed2 == "a2"):
                    print(fixed_error);
                    return;
                    
                elif( fixed2 == "b1"):
                    print(fixed_error);
                    return;
                
            
            # have to determine what the free params are
            # iter over all the possibilities
            if( free1 == "a1"):
                
                if( free2 == "a2"):
                    
                    # define tuple for making exchange matrix
                    # syntax is (phi, a1, a2, b1, b2)
                    E_inputs = (phi, val1, val2, fixed1_val, fixed2_val);
                    
                    # also make labels for plotting
                    xlabel, ylabel = "a1 [eV]", "a2 [eV]"
                    
                    #constuct the exchange matrix
                    E = hamiltonian.Exchange(*E_inputs);
                    
                    # add it to the hamiltonian
                    H = hamiltonian.Hamiltonian(*H_inputs);
                    H.Add(E);
                    
                    # get the eigenvectors at i = 0, j= 0 only
                    if( i == 0 and j== 0):
                        eigvecs = H.Eigvecs();
        
                    #iter thru eigenvals
                    for eig_i in range(len(eigvecs)):
            
                        # get the eigval that corresponds to this eigvec
                        eigvec = eigvecs[eig_i];
                        eigval = H.GetEigenval(eigvec);
            
                        # one energy level val for this param i
                        E_levels[eig_i, i, j] = eigval;
            
                    # delete the matrices
                    del H;
                    del E;
                                    
                elif(free2 == "b1"):
                    
                    # define tuple for making exchange matrix
                    # syntax is (phi, a1, a2, b1, b2)
                    E_inputs = (phi, val1, fixed1_val, val2, fixed2_val);
                    
                    # also make labels for plotting
                    xlabel, ylabel = "a1 [eV]", "b1 [eV]"
                    
                    #constuct the exchange matrix
                    E = hamiltonian.Exchange(*E_inputs);
                    
                    # add it to the hamiltonian
                    H = hamiltonian.Hamiltonian(*H_inputs);
                    H.Add(E);
                    
                    # get the eigenvectors at i = 0, j= 0 only
                    if( i == 0 and j== 0):
                        eigvecs = H.Eigvecs();
        
                    #iter thru eigenvals
                    for eig_i in range(len(eigvecs)):
            
                        # get the eigval that corresponds to this eigvec
                        eigvec = eigvecs[eig_i];
                        eigval = H.GetEigenval(eigvec);
            
                        # one energy level val for this param i
                        E_levels[eig_i, i, j] = eigval;
            
                    # delete the matrices
                    del H;
                    del E;
                    
                elif(free2 == "b2"):
                    
                    # define tuple for making exchange matrix
                    # syntax is (phi, a1, a2, b1, b2)
                    E_inputs = (phi, val1, fixed1_val, fixed2_val, val2);
                    
                    # also make labels for plotting
                    xlabel, ylabel = "a1 [eV]", "b2 [eV]";
                    
                    #constuct the exchange matrix
                    E = hamiltonian.Exchange(*E_inputs);
                    
                    # add it to the hamiltonian
                    H = hamiltonian.Hamiltonian(*H_inputs);
                    H.Add(E);
                    
                    # get the eigenvectors at i = 0, j= 0 only
                    if( i == 0 and j== 0):
                        eigvecs = H.Eigvecs();
        
                    #iter thru eigenvals
                    for eig_i in range(len(eigvecs)):
            
                        # get the eigval that corresponds to this eigvec
                        eigvec = eigvecs[eig_i];
                        eigval = H.GetEigenval(eigvec);
            
                        # one energy level val for this param i
                        E_levels[eig_i, i, j] = eigval;
            
                    # delete the matrices
                    del H;
                    del E;
                    
            elif(free1 == "a2"):
                
                if(free2 == "b1"):
                    
                    # define tuple for making exchange matrix
                    # syntax is (phi, a1, a2, b1, b2)
                    E_inputs = (phi, fixed1_val, val1, val2, fixed2_val);
                    
                    # also make labels for plotting
                    xlabel, ylabel = "a2 [eV]", "b1 [eV]"
                    
                    #constuct the exchange matrix
                    E = hamiltonian.Exchange(*E_inputs);
                    
                    # add it to the hamiltonian
                    H = hamiltonian.Hamiltonian(*H_inputs);
                    H.Add(E);
                    
                    # get the eigenvectors at i = 0, j= 0 only
                    if( i == 0 and j== 0):
                        eigvecs = H.Eigvecs();
        
                    #iter thru eigenvals
                    for eig_i in range(len(eigvecs)):
                        
                        # get the eigval that corresponds to this eigvec
                        eigvec = eigvecs[eig_i];
                        eigval = H.GetEigenval(eigvec);
            
                        # one energy level val for this param i
                        E_levels[eig_i, i, j] = eigval;
            
                    # delete the matrices
                    del H;
                    del E;
                    
                elif(free2 == "b2"):
                    
                    # define tuple for making exchange matrix
                    # syntax is (phi, a1, a2, b1, b2)
                    E_inputs = (phi, fixed1_val, val1, fixed2_val, val2 );
                    
                    # also make labels for plotting
                    xlabel, ylabel = "a2 [eV]", "b2 [eV]"
                    
                    #constuct the exchange matrix
                    E = hamiltonian.Exchange(*E_inputs);
                    
                    # add it to the hamiltonian
                    H = hamiltonian.Hamiltonian(*H_inputs);
                    H.Add(E);
                    
                    # get the eigenvectors at i = 0, j= 0 only
                    if( i == 0 and j== 0):
                        eigvecs = H.Eigvecs();
        
                    #iter thru eigenvals
                    for eig_i in range(len(eigvecs)):
            
                        # get the eigval that corresponds to this eigvec
                        eigvec = eigvecs[eig_i];
                        eigval = H.GetEigenval(eigvec);
            
                        # one energy level val for this param i
                        E_levels[eig_i, i, j] = eigval;
            
                    # delete the matrices
                    del H;
                    del E;
                    
            elif( free1 == "b1"):
                
                if(free2 == "b2"):
                    
                    # define tuple for making exchange matrix
                    # syntax is (phi, a1, a2, b1, b2)
                    E_inputs = (phi, fixed1_val, fixed2_val, val1, val2 );
                    
                    # also make labels for plotting
                    xlabel, ylabel = "b1 [eV]", "b2 [eV]"
                    
                    #constuct the exchange matrix
                    H = hamiltonian.Hamiltonian(*H_inputs);
                    E = hamiltonian.Exchange(*E_inputs);
                    
                    # add it to the hamiltonian
                    H.Add(E);
                    
                    # get the eigenvectors at i = 0, j= 0 only
                    if( i == 0 and j== 0):
                        eigvecs = H.Eigvecs();
        
                    #iter thru eigenvals
                    for eig_i in range(len(eigvecs)):
            
                        # get the eigval that corresponds to this eigvec
                        eigvec = eigvecs[eig_i];
                        eigval = H.GetEigenval(eigvec);
            
                        # one energy level val for this param i
                        E_levels[eig_i, i, j] = eigval;
            
                    # delete the matrices
                    del H;
                    del E;
                    
            else: # free params where chosen wrong
                print("Error: the parameter selection:")
                print("free1 = "+free1+ ", free2 = "+free2+ ", fixed1 = "+fixed1+ ", fixed2 = "+fixed2)
                print("is incorrect. Make sure your ordering obeys a1, a2, b1, b2.")
                return;
                
                    
    #### end i, j loops
    
    #### now we do the plotting
    
    # get mesh of x, y vals
    x, y = np.meshgrid(free_arr1, free_arr2);
    
    # plot the energy levels
    fig = plt.figure();
    ax = fig.gca(projection = "3d");
    
    # plot the levels we found
    for level_i in range( N_levels ):
        
        # check for degeneracy before plotting
        degenerate = False;
        flat = True;
        for level_j in range(level_i): # iter over levels to check equality/degeneracy
            if( np.array_equal(E_levels[level_j], E_levels[level_i]) ): # there is degeneracy
                degenerate = True;
                
            
        # flatness means all the row arrays are equal
        row0 = E_levels[level_i][0]
        for row in E_levels[level_i]: # iter over rows to check equality
            if( not np.array_equal(row, row0) ):
                # if we ever get here it is not flat
                flat = False;
        #also means all the columns are equal
        col0 = E_levels[ level_i, :, 0];
        for col_i in range( N_param): # iter over rows to check equality
            col = E_levels[level_i,  :, col_i ];
            if( not np.array_equal(col, col0) ):
                # if we ever get here it is not flat
                flat = False;
        COLORS = ['red', 'orangered','greenyellow','green']
        # plot as square if flat
        if(flat):       
            ax.plot_wireframe(x, y, E_levels[level_i], color = COLORS[level_i], label = "$E_{"+str(level_i)+"}$", cstride = 100, rstride = 100);
        else:
            ax.plot_wireframe(x, y, E_levels[level_i], color = COLORS[level_i], label = "$E_{"+str(level_i)+"}$", cstride = 10, rstride = 10);
            
    # finally plot the mid gap energy
    mid_gap = np.zeros((N_param, N_param)); 
    ax.plot_wireframe(x, y, mid_gap, color = 'gray', label = "Mid gap energy", rstride = 100, cstride = 100)  
        
    
    # format the plot
    plt.xlabel(xlabel);
    plt.ylabel(ylabel);
    plt.title(material + " bulk E levels", loc = 'left', fontsize = 16 );
    plt.title("$\phi$ = "+str(phi)[:4]+" [rad], "+fixed1 +" = "+str(fixed1_val)+", "+fixed2+" = "+str(fixed2_val), loc = 'right', fontsize = 10);
    ax.legend();
    
    plt.show();
    
    #### end PlotLevel2D
    
    
    
def ParamSS(H_inputs, phi, condition, material, a1_lim, a2_lim, b1_lim, b2_lim, test = "None"):
    '''
    Param Space Search: Traverses the param space define by a1_lim, a2_lim etc, usually the entire param space. Finds
    the energy levels everywhere in this space, but also specifically searches for the region where
    the conditon is met. Creates an array of bools corresponding to all param space, telling where 
    the condition is met. 
    
    Can optionally search along only one axis and hold the other cals fixed at lim[0], this is done by setting
    the keyword test to the free param.
    
    :param H_inputs: tuple of the 5 params needed for the bulk hamiltonian H
                    (Eg, muBohr, Bvec, vt, vl
    :param phi: double, angle between H field and z hat in radians
    :param condition: points to function that takes array of energy levels at one point in param space
                as an argument and returns bool according to certain condition. Condition functions are
                defined below
    :param a1_lim: tuple of limts for searching a1 space
    :param a2_lim: see above
    :param b1_lim: see above
    :param b2_lim: see above
    :param test: string keyword arg, defaults to None in which case all axes are searched. if set to name
                of one param eg a1, only that axis will be searched. For now mainly used for debugging
                
    returns tuple of all the args needed for functions that plot the results (PlotCondition1D)
    tuple of ( big 5D array of energy levels in param space, big 4D array of bool vals in param space,
                1d arrays that define param space along a1, a2, b1, b2)
    '''
    
    # define important variables for the function
    N_param = 6; # how many pts along each param space axis
    N_levels = hamiltonian.Hamiltonian(*H_inputs).N; # number of energy levels
    E_levels = np.zeros((N_levels, N_param, N_param, N_param, N_param));
    
    # create arrays to iter over
    # this is actually a little complicated because in the test cases (test != None) we 
    # only iter over one arr
    if( test == "None"):
        
        # all are real linspace arrays
        a1_arr = np.linspace(*a1_lim, N_param);
        a2_arr = np.linspace(*a2_lim, N_param);
        b1_arr = np.linspace(*b1_lim, N_param);
        b2_arr = np.linspace(*b2_lim, N_param);
        
    elif( test == "a1"):
        
        # only a1 is actually itered over
        a1_arr = np.linspace(*a1_lim, N_param);
        
        # the rest are just the same val over and over
        a2_arr = np.full(N_param, a2_lim[0]);
        b1_arr = np.full(N_param, b1_lim[0]);
        b2_arr = np.full(N_param, b2_lim[0]);
        
    elif( test == "a2"):
        
        # only a2 is actually itered over
        a2_arr = np.linspace(*a2_lim, N_param);
        
        # the rest are just the same val over and over
        a1_arr = np.full(N_param, a1_lim[0]);
        b1_arr = np.full(N_param, b1_lim[0]);
        b2_arr = np.full(N_param, b2_lim[0]);
        
    elif( test == "b1"):
        
        # only b1 is actually itered over
        b1_arr = np.linspace(*b1_lim, N_param);
        
        # the rest are just the same val over and over
        a1_arr = np.full(N_param, a1_lim[0]);
        a2_arr = np.full(N_param, a2_lim[0]);
        b2_arr = np.full(N_param, b2_lim[0]);
        
    elif( test == "b2"):
        
        # only b2 is actually itered over
        b2_arr = np.linspace(*b2_lim, N_param);
        
        # the rest are just the same val over and over
        a1_arr = np.full(N_param, a1_lim[0]);
        a2_arr = np.full(N_param, a2_lim[0]);
        b1_arr = np.full(N_param, b1_lim[0]);
    
    # get energy levels for all param vals
    # iter over the entire param space!
    
    # iter over a1
    for a1_i in range(len(a1_arr)):
        
        print(a1_i);
        
        # iter over a2
        for a2_i in range(len(a2_arr)):
            
                # iter over b1
                for b1_i in range(len(b1_arr)):
                    
                    # iter over b2
                    for b2_i in range(len(b2_arr)):
                        
                        # here is where we actually get the E levels
                        
                        #assign param vals
                        a1, a2, b1, b2 = a1_arr[a1_i], a2_arr[a2_i], b1_arr[b1_i], b2_arr[b2_i]
                        
                        # make the matrices, combine
                        H = hamiltonian.Hamiltonian(*H_inputs);
                        E_inputs = (phi, a1, a2, b1, b2);
                        E = hamiltonian.Exchange(*E_inputs);
                        H.Add(E);
                        
                        # get the eigenvectors the initial run thru
                        if( a1_i == 0 and a2_i == 0 and b1_i == 0 and b2_i == 0):
                            
                            eigvecs = H.Eigvecs();
                            
                        # get the eigenvalue for each level
                        # this is done (right) by itering thru eigvecs and getting corresponding eigval
                        for eig_i in range(len(eigvecs)):
    
                            # get the eigval that corresponds to this  (initial) eigvec
                            eigvec = eigvecs[eig_i];
                            eigval = H.GetEigenval(eigvec);
    
                            # place level in correct spot of container array
                            E_levels[eig_i, a1_i, a2_i, b1_i, b2_i ] = eigval;
                            
                        # delete the matrices
                        del H;
                        del E;
                                
    #### end over itering over param space    
    # now we check condition at each point in param space
    
    # create H at Eg = 150 meV and E at params = 0 to determine benchmarking/ naming
    H_inputs = (0.15, H_inputs[1], H_inputs[2], H_inputs[3], H_inputs[4] );
    H = hamiltonian.Hamiltonian(*H_inputs);
    E = hamiltonian.Exchange(*E_inputs); # this is subtle: just takes last E_inputs given, essentially looking at right side of plot
    H.Add(E);
    
    # get levels from this H
    Hi_E_levels = np.zeros(N_levels);
    for eig_i in range(N_levels):
        Hi_E_levels[eig_i] = H.GetEigenval(eigvecs[eig_i]);
     
    # use this as the benchmark point for labeling   
    labels, colors = GetBandLabels(H.MakeCopy(), Hi_E_levels, phi, False, material );
    
    # container variable for condition checking
    condition_arr = np.zeros((N_param, N_param, N_param, N_param));
    # so we have to iter over it again
    # iter over a1
    for a1_i in range(len(a1_arr)):
        
            # iter over a2
            for a2_i in range(len(a2_arr)):
                
                    # iter over b1
                    for b1_i in range(len(b1_arr)):
                        
                            # iter over b2
                            for b2_i in range(len(b2_arr)):
                                
                                # get all the energy levels here in param space
                                energies = E_levels[:, a1_i, a2_i, b1_i, b2_i ];
                                
                                # check the whether the levels here meet the appropriate condition
                                # this is done with the condition function, which returns a bool
                                # then the result is place in the condition_arr which holds results
                                # at all point in param space
                                condition_arr[a1_i, a2_i, b1_i, b2_i] = condition(energies, Hi_E_levels,phi, H.MakeCopy());
                                #if( condition_arr[a1_i, a2_i, b1_i, b2_i] ):
                                    #a1, a2, b1, b2 = a1_arr[a1_i], a2_arr[a2_i], b1_arr[b1_i], b2_arr[b2_i]
                                    #print(a1, a2, b1, b2);
                                    
    #### end of the condition array creation                           
                                
    return E_levels, condition_arr, a1_arr, a2_arr, b1_arr, b2_arr, phi, labels, colors;

    #### end ParamSpaceSearch 
    
    
    
def ParamSS_AB(H_inputs, phi, condition, material, a1_lim, b1_lim, A, B, test = "None"):
    '''
    Traverses the param space define by a1_lim, b1_lim, with A, B fixed. Finds
    the energy levels everywhere in this space, but also specifically searches for the region where
    the conditon is met. Creates an array of bools corresponding to all param space, telling where 
    the condition is met. 
    
    Can optionally search along only one axis and hold the other cals fixed at lim[0], this is done by setting
    the keyword test to the free param.
    
    :param H_inputs: tuple of the 5 params needed for the bulk hamiltonian H
                    (Eg, muBohr, Bvec, vt, vl
    :param phi: double, angle between H field and z hat in radians
    :param condition: points to function that takes array of energy levels at one point in param space
                as an argument and returns bool according to certain condition. Condition functions are
                defined below
    :param a1_lim: tuple of limts for searching a1 space
    :param b1_lim: see above
    :param A: fixed value of input param A = a1 - a2
    :param B: fixed value of input param B = b1 - b2
    :param test: string keyword arg, defaults to None in which case all axes are searched. if set to name
                of one param eg a1, only that axis will be searched. For now mainly used for debugging
                
    returns tuple of all the args needed for functions that plot the results (PlotCondition1D)
    tuple of ( big 5D array of energy levels in param space, big 4D array of bool vals in param space,
                1d arrays that define param space along a1, a2, b1, b2)
    '''
    
    # define important variables for the function
    N_param = 50; # how many pts along each param space axis
    N_levels = hamiltonian.Hamiltonian(*H_inputs).N; # number of energy levels
    E_levels = np.zeros((N_levels, N_param, N_param));
    
    # create arrays to iter over
    # this is actually a little complicated because in the test cases (test != None) we 
    # only iter over one arr
    if( test == "None"):
        
        # all are real linspace arrays
        a1_arr = np.linspace(*a1_lim, N_param);
        b1_arr = np.linspace(*b1_lim, N_param);
        
    elif( test == "a1"):
        
        # only a1 is actually itered over
        a1_arr = np.linspace(*a1_lim, N_param);
        
        # the rest are just the same val over and over
        b1_arr = np.full(N_param, b1_lim[0]);
        
    elif( test == "b1"):
        
        # only b1 is actually itered over
        b1_arr = np.linspace(*b1_lim, N_param);
        
        # the rest are just the same val over and over
        a1_arr = np.full(N_param, a1_lim[0]);
    
    # get energy levels for all param vals
    # iter over the entire param space!
    
    # iter over a1
    for a1_i in range(len(a1_arr)):
            
        # iter over b1
        for b1_i in range(len(b1_arr)):
                        
            # here is where we actually get the E levels
            
            #assign param vals
            a1, b1 = a1_arr[a1_i], b1_arr[b1_i]
            
            # get a2 and b2 from fixed vals
            a2 = a1 - A;
            b2 = b1 - B;
            
            # make the matrices, combine
            H = hamiltonian.Hamiltonian(*H_inputs);
            E_inputs = (phi, a1, a2, b1, b2);
            E = hamiltonian.Exchange(*E_inputs);
            H.Add(E);
            
            # get the eigenvectors the initial run thru
            if( a1_i == 0 and b1_i == 0):
                
                eigvecs = H.Eigvecs();
                
            # get the eigenvalue for each level
            # this is done (right) by itering thru eigvecs and getting corresponding eigval
            for eig_i in range(len(eigvecs)):

                # get the eigval that corresponds to this  (initial) eigvec
                eigvec = eigvecs[eig_i];
                eigval = H.GetEigenval(eigvec);

                # place level in correct spot of container array
                E_levels[eig_i, a1_i, b1_i] = eigval;
                
            # delete the matrices
            del H;
            del E;
                                
    #### end over itering over param space    
    # now we check condition at each point in param space
    
    # create H at Eg = 150 meV and E at params = 0 to determine benchmarking/ naming
    H_inputs = (0.15, H_inputs[1], H_inputs[2], H_inputs[3], H_inputs[4] );
    H = hamiltonian.Hamiltonian(*H_inputs);
    print("bench H", H);
    E = hamiltonian.Exchange(*E_inputs); # this is subtle: just takes last E_inputs given, essentially looking at right side of plot
    H.Add(E);
    
    # get levels from this H
    Hi_E_levels = np.zeros(N_levels);
    for eig_i in range(N_levels):
        Hi_E_levels[eig_i] = H.GetEigenval(eigvecs[eig_i]);
     
    # use this as the benchmark point for labeling   
    labels, colors = GetBandLabels(H.MakeCopy(), Hi_E_levels, phi, False, material);
    
    # container variable for condition checking
    condition_arr = np.zeros((N_param, N_param ));
    # so we have to iter over it again
    # iter over a1
    for a1_i in range(len(a1_arr)):
                
        # iter over b1
        for b1_i in range(len(b1_arr)):
                                
            # get all the energy levels here in param space
            energies = E_levels[:, a1_i, b1_i];
            
            # check the whether the levels here meet the appropriate condition
            # this is done with the condition function, which returns a bool
            # then the result is place in the condition_arr which holds results
            # at all point in param space
            condition_arr[a1_i, b1_i] = condition(energies, Hi_E_levels,phi, H.MakeCopy());
                                    
    #### end of the condition array creation                           
                                
    return E_levels, condition_arr, a1_arr, A, b1_arr, B, phi, labels, colors;

    #### end ParamSS_AB 
    
    
def ProcessCondition( axis, energy_arr, condition_arr, a1_arr, a2_arr, b1_arr, b2_arr, phi, labels, colors):  
    '''
    What we want is a plot along the selected axis which shows the condition being met at
    some points in parameter space. Then we can use span_where to highlight this region
    
    Takes all the returns of ParamSpaceSearch
    
    :param axis: string, tells which axis to plot along, a1, a2, b1, b2
    :param energy_arr: 5D array of energy levels everywhere in param space
    :param condition_arr: 4D array of condition bools everywhere in param space
    :param a1_arr: 1d array of domain along a1 axis
    :param a2_arr: same but for a2
    :param b1_arr: same but for b1
    :param b2_arr: same but for b2
    :param phi: double, angle between H field and z hat
    :param labels: 1d array of strings which are legend labels from GetBandLabels
    :param colors: 1d array of strings which are color names from GetBandLabels
    ''' 
    
    # some relevant values
    N_params = len(a1_arr); # number of vals along each param axis
    N_levels = len(energy_arr); # energy levels at each point in param space
    cond_thresh = int(N_params/5); # how long along the axis the condition needs to be maintained
    # should change this 
    
    # control variables based on which axis is selected
    if( axis == "a1"):
        
        x_arr = a1_arr; # will be the x axis of the plot
        x_label = "a1 [eV]"; # label of the x axis
        
        param_arr = [a2_arr, b1_arr, b2_arr]; # all of param space in one array
        param_names = ["a2", "b1", "b2"]; # names for title
        
        def GetEnergyArr( i, j, k):
            return energy_arr[:, :, i, j, k];
        
        def GetConditionArr( i, j, k):
            return condition_arr[:, i, j, k];  
        
    # control variables based on which axis is selected
    if( axis == "a2"):
        
        x_arr = a2_arr; # will be the x axis of the plot
        x_label = "a2 [eV]"; # label of the x axis
        
        param_arr = [a1_arr, b1_arr, b2_arr]; # all of param space in one array
        param_names = ["a1", "b1", "b2"]; # names for title
        
        def GetEnergyArr( i, j, k):
            return energy_arr[:, i, :, j, k];
        
        def GetConditionArr( i, j, k):
            return condition_arr[i, :, j, k];   
        
    # control variables based on which axis is selected
    if( axis == "b1"):
        
        x_arr = a2_arr; # will be the x axis of the plot
        x_label = "b1 [eV]"; # label of the x axis
        
        param_arr = [a1_arr, b1_arr, b2_arr]; # all of param space in one array
        param_names = ["a1", "a2", "b2"]; # names for title
        
        def GetEnergyArr( i, j, k):
            return energy_arr[:, i, j, :, k];
        
        def GetConditionArr( i, j, k):
            return condition_arr[i, j, :, k];  
        
    # control variables based on which axis is selected
    if( axis == "b2"):
        
        x_arr = a2_arr; # will be the x axis of the plot
        x_label = "b2 [eV]"; # label of the x axis
        
        param_arr = [a1_arr, b1_arr, b2_arr]; # all of param space in one array
        param_names = ["a1", "a2", "b1"]; # names for title
        
        def GetEnergyArr( i, j, k):
            return energy_arr[:, i, j, k, :];
        
        def GetConditionArr( i, j, k):
            return condition_arr[i, j, k, :];       
        
    # sweep over the remaining 3D param space
    # find spot where the condition is met a certain number of times consecutively
    params_found = False;
        
    # iter over the param space
    # stop sweep once we find condition is met
    i, j, k = 0, 0, 0; 
    while (i < N_params and not params_found):
        while (j < N_params and not params_found):
            while (k < N_params and not params_found):
                
                # at this point in the param space, get the energy levels and condition arrays
                # this will each be rows over the axis param
                E_levels = GetEnergyArr(i, j, k);
                conditions = GetConditionArr(i, j, k);
                
                # check that the condition is met consecutively somewhere
                cond_counter = 0;
                for my_cond in conditions: # iter over the cond array
                    
                    # if met increase the counter
                    if(my_cond):
                        cond_counter += 1;
                    else: # reset the counter
                        cond_counter = 0;
                        
                    # check if counter has crossed the threshold
                    if( cond_counter > cond_thresh):                       
                        
                        # this is a good spot in param space
                        params_found = True;
                        
                        # get the relevant fixed params
                        param0, param1, param2 = param_arr[0][i], param_arr[1][j], param_arr[2][k];
        
                k += 1;                
            j += 1;
        i += 1;
        
    #### end of the param space sweep 
    
    # check that is was successful
    if( not params_found):
        print("Error in ProcessConditionArr: the condition was not met "+str(cond_thresh)+" times consecutively anywhere in the param space.");
        return;
        
    # create static arrays for passing to PlotCondition1D function
    # this also depends on the axis chosen
    if( axis == "a1"):
        a2_arr = np.full(N_params, param0);
        b1_arr = np.full(N_params, param1);
        b2_arr = np.full(N_params, param2);
    elif( axis == "a2"):
        a1_arr = np.full(N_params, param0);
        b1_arr = np.full(N_params, param1);
        b2_arr = np.full(N_params, param2);
    elif( axis == "b1"):
        a1_arr = np.full(N_params, param0);
        a2_arr = np.full(N_params, param1);
        b2_arr = np.full(N_params, param2);
    elif( axis == "b2"):
        a1_arr = np.full(N_params, param0);
        a2_arr = np.full(N_params, param1);
        b1_arr = np.full(N_params, param2);
    
    return energy_arr, condition_arr, a1_arr, a2_arr, b1_arr, b2_arr, phi, labels, colors;  

    #### end process condition arr
    
    
def ProcessCondition_AB( axis, energy_arr, condition_arr, a1_arr, A, b1_arr, B, phi, labels, colors ):  
    '''
    What we want is a plot along the selected axis which shows the condition being met at
    some points in parameter space. Then we can use span_where to highlight this region
    
    Takes all the returns of ParamSpaceSearch
    
    :param axis: string, tells which axis to plot along, a1, a2, b1, b2
    :param energy_arr: 5D array of energy levels everywhere in param space
    :param condition_arr: 4D array of condition bools everywhere in param space
    :param a1_arr: 1d array of domain along a1 axis
    :param A: double, fixed val of A = a1 - a2
    :param b1_arr: same as a1 but for b1
    :param B: double, fixed val of b1 - b2
    :param phi: double, angle between H field and z hat
    :param labels: 1d array of strings which are legend labels from GetBandLabels
    :param colors: 1d array of strings which are color names from GetBandLabels
    ''' 
    
    # some relevant values
    N_params = len(a1_arr); # number of vals along each param axis
    N_levels = len(energy_arr); # energy levels at each point in param space
    cond_thresh = 0 # int(N_params/5); # how long along the axis the condition needs to be maintained
    
    # control variables based on which axis is selected
    if( axis == "a1"):
        
        x_arr = a1_arr; # will be the x axis of the plot
        x_label = "a1 [eV]"; # label of the x axis
        
        param_arr = b1_arr; # all of param space in one array
        param_names = ["A", "b1", "B"]; # names for title
        
        def GetEnergyArr( i ):
            return energy_arr[:, :, i];
        
        def GetConditionArr( i):
            return condition_arr[:, i ];          
        
    # sweep over the remaining 3D param space
    # find spot where the condition is met a certain number of times consecutively
    params_found = False;
        
    # iter over the param space
    # stop sweep once we find condition is met
    i = 0;
    while (i < N_params and not params_found):
                
        # at this point in the param space, get the energy levels and condition arrays
        # this will each be rows over the axis param
        E_levels = GetEnergyArr(i);
        conditions = GetConditionArr(i);
        
        # check that the condition is met consecutively somewhere
        cond_counter = 0;
        for my_cond in conditions: # iter over the cond array
            
            # if met increase the counter
            if(my_cond):
                cond_counter += 1;
            else: # reset the counter
                cond_counter = 0;
                
            # check if counter has crossed the threshold
            if( cond_counter > cond_thresh):                       
                
                # this is a good spot in param space
                params_found = True;
                
                # get the relevant fixed params
                param0 = param_arr[i];
        
        i += 1;
        
    #### end of the param space sweep 
    
    # check that is was successful
    if( not params_found):
        print("Error in ProcessConditionArr: the condition was not met "+str(cond_thresh)+" times consecutively anywhere in the param space.");
        return;
        
    # create static arrays for passing to PlotCondition1D function
    # this also depends on the axis chosen
    if( axis == "a1"):
        b1_arr = np.full(N_params, param0);
    
    return energy_arr, condition_arr, a1_arr, A, b1_arr, B, phi, labels, colors;                                                          
    
               
def PlotCondition(energy_arr, condition_arr, a1_arr, a2_arr, b1_arr, b2_arr, phi, labels, colors, material):  
    '''
    Using the results of ParamSpaceSearch ( but only when it searches only 1 axis using the test keyword)
    plots the energy levels along a given axis and highlights the region where the condition is met
    
    Mainly a benchmarking function as of now. Want to revamp or make new similar but more useful
    function soon.
    
    :param energy_arr: big 5D array of energy levels everywhere in param space
    :param condition_arr: big 4d array of bool vals everywhere in param space
    :param a1_arr: 1d array that defines a1 param space
    :param a2_arr: same as a1_arr. As of now, only one of these param arrays actually iters. The rest
                    are static, we just pull out the first values as constants
    :param b1_arr: same as above
    :param b2_arr: same as above
    :param phi: double, the angle between the H field and z hat in radians
    :param labels: 1d array of strings which are legend labels from GetBandLabels
    :param colors: 1d array of strings which are color names from GetBandLabels
    :param material: string, the name of the material, for title purposes
    '''
    
    # get relevant vars
    N_levels = len(energy_arr);
    ax = plt.subplot(); # for plotting
    
    # check which arr is the free parameter
    # this gives lets us choose the following necessary values
    param_arrs = [a1_arr, a2_arr, b1_arr, b2_arr];
    param_names = ["a1", "a2", "b1", "b2"];
    
    # find the correct free parameter
    for i in range(4): # iter over 4 params
        
        # we choose the array that is not static
        # this is the array with endpoints not equal
        if( param_arrs[i][0] != param_arrs[i][-1]): # this is the one we vary over
            
            # get the x vals and labels
            x_arr = param_arrs[i];
            xlabel = param_names[i];
            
            # get the y vals ( a little complicated since the energy level arr is so big)       
            if( i == 0): # for varying over a1
                
                # determine the fixed values
                fixed1, fixed2, fixed3 = param_names[1], param_names[2], param_names[3];
                fixed1_val, fixed2_val, fixed3_val = a2_arr[0], b1_arr[0], b2_arr[0];
                
                # get the condition values
                cond_vals = condition_arr[:, 0, 0, 0];
                
                # get the energy levels we want
                E_levels = energy_arr[:, :, 0, 0, 0];

    # now we can plot
    for eig_i in range(N_levels):           

        if(phi == 0):
            plt.plot(x_arr, E_levels[eig_i], label = labels[eig_i], color = colors[eig_i], linestyle = '--');
        else:
            plt.plot(x_arr, E_levels[eig_i], label = labels[eig_i], color = colors[eig_i]);
                    
    # format the plot
    plt.xlabel(xlabel);
    plt.ylabel("Energy [eV]");
    plt.title(material + " bulk E levels", loc = 'left', fontsize = 16 );
    plt.title(fixed1 +" = "+str(fixed1_val)+", "+fixed2+" = "+str(fixed2_val)+", "+fixed3+" = "+str(fixed3_val), loc = 'right', fontsize = 10);
     
    # plot the condition array
    # this is just so that it is accessible to the span_where function
    ax.plot(x_arr, cond_vals, scaley = False, alpha = 0);
    
    # use span_where to shade the true regions
    collect = collections.BrokenBarHCollection.span_where(x_arr, 2*plt.ylim()[0], 2*plt.ylim()[1], where = cond_vals, facecolor = 'slateblue', alpha = 0.4);
    ax.add_collection(collect);
    
    # make the legend
    chartBox = ax.get_position()
    ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.6, chartBox.height])
    ax.legend(loc='upper center', bbox_to_anchor=(1.45, 0.8), shadow=True, ncol=1)
    
    # show and return 
    plt.show()               
    return;
    #### end PlotCondition1D 
    
    
def PlotCondition_AB(energy_arr, condition_arr, a1_arr, A, b1_arr, B, phi, labels, colors, material, return_results = False):  
    '''
    Using the results of ParamSpaceSearch ( but only when it searches only 1 axis using the test keyword)
    plots the energy levels along a given axis and highlights the region where the condition is met
    
    Mainly a benchmarking function as of now. Want to revamp or make new similar but more useful
    function soon.
    
    :param energy_arr: big 5D array of energy levels everywhere in param space
    :param condition_arr: big 4d array of bool vals everywhere in param space
    :param a1_arr: 1d array that defines a1 param space
    :param A: fixed value of a1 - a2
    :param b1_arr: same as a1
    :param B: fixed value of b1 - b2
    :param phi: double, the angle between the H field and z hat in radians
    :param labels: 1d array of strings which are legend labels from GetBandLabels
    :param colors: 1d array of strings which are color names from GetBandLabels
    :param material: string, the name of the material, for title purposes
    '''
    
    # get relevant vars
    N_levels = len(energy_arr);
    ax = plt.subplot(); # for plotting
    
    # check which arr is the free parameter
    # this gives lets us choose the following necessary values
    param_arrs = [a1_arr, A, b1_arr, B];
    param_names = ["a1", "A", "b1", "B"];
    
    # find the correct free parameter
    for i in [0, 2]: # iter over a1 and b1
        
        # we choose the array that is not static
        # this is the array with endpoints not equal
        if( param_arrs[i][0] != param_arrs[i][-1]): # this is the one we vary over
            
            # get the x vals and labels
            x_arr = param_arrs[i];
            xlabel = param_names[i] + " [eV]";
            
            # get the y vals ( a little complicated since the energy level arr is so big)       
            if( i == 0): # for varying over a1
                
                # determine the fixed values
                fixed1, fixed2, fixed3 = param_names[1], param_names[2], param_names[3];
                fixed1_val, fixed2_val, fixed3_val = A, b1_arr[0], B;
                
                # get the condition values
                cond_vals = condition_arr[:, 0];
                
                # get the energy levels we want
                E_levels = energy_arr[:, :, 0];
                
    # return levels if asked
    if (return_results):
        
        return x_arr, E_levels, labels, colors;

    # now we can plot
    for eig_i in range(N_levels):           

        if(phi == 0):
            plt.plot(x_arr, E_levels[eig_i], label = labels[eig_i], color = colors[eig_i], linestyle = '--');
        else:
            plt.plot(x_arr, E_levels[eig_i], label = labels[eig_i], color = colors[eig_i]);
                    
    # format the plot
    plt.xlabel(xlabel);
    plt.ylabel("Energy [eV]");
    plt.title(material + " bulk E levels", loc = 'left', fontsize = 16 );
    plt.title(fixed1 +" = "+str(fixed1_val)[:4]+", "+fixed2+" = "+str(fixed2_val)[:4]+", "+fixed3+" = "+str(fixed3_val)[:4], loc = 'right', fontsize = 10);
     
    # plot the condition array
    # this is just so that it is accessible to the span_where function
    ax.plot(x_arr, cond_vals, scaley = False, alpha = 0);
    
    # use span_where to shade the true regions
    collect = collections.BrokenBarHCollection.span_where(x_arr, 2*plt.ylim()[0], 2*plt.ylim()[1], where = cond_vals, facecolor = 'slateblue', alpha = 0.4);
    ax.add_collection(collect);
    
    # make the legend
    chartBox = ax.get_position()
    ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.6, chartBox.height])
    leg1 = ax.legend(loc='upper center', bbox_to_anchor=(1.45, 0.8), shadow=True, ncol=1)
    
    # show and return 
    #plt.show()               
    return leg1;

    #### end PlotCondition1D  
    
    
def GetIntersections(x_arr, y_arr1, y_arr2, fineness = 5, help_debug = True):
    '''
    Fits to data sets yarr1 and yarr2 to lines with their slope at x = 0. Finds if the resulting two lines intersect.
    If so, returns point of intersection.
    
    Used for plotting intersections of photon transition energies in GetEnergyLevels.EgSweepTransition
    
    :param x_arr: 1d arr of the domain values of the data, corresponds to both y's
    :param y_arr1: 1d array of y vals for the data
    :param y_arr2: 1d array of second data set
    :param fineness: int, defaults to 5, tells how many mesh points in x to use to search for the intersection, eg 10^fineness points
    :param help_debug: bool, defaults to True, tells function to pring out helpful info for debugging
    
    returns: tuple of
            x_int, double, x coordinate of intersection, None if not found
            y_int, double, same but y coordinate
    '''
    
    # define return variables
    # if there is an intersection found these will be overwritten with numbers
    x_int = None;
    y_int = None;
    
    # write each y arr in y = mx + b form, which describes line around x[i_0]
    # sometimes when there are level intersections the transition line is piecewise of two line segments
    # so fit two lines to each, one for left side of domain one for right
    
    # get slope, intercept for y1, at left endpoint of domain
    m1_left = (y_arr1[1] - y_arr1[0])/(x_arr[1] - x_arr[0] );
    b1_left = y_arr1[0] - m1_left*x_arr[0];
    
    # get slope, intercept for y1, at right endpoint of domain
    m1_right = (y_arr1[-1] - y_arr1[-2])/(x_arr[-1] - x_arr[-2] );
    b1_right = y_arr1[-1] - m1_right*x_arr[-1];
    
    # get slope, intercept for y1, at left endpoint of domain
    m2_left = (y_arr2[1] - y_arr2[0])/(x_arr[1] - x_arr[0] );
    b2_left = y_arr2[0] - m2_left*x_arr[0];
    
    # get slope, intercept for y1, at right endpoint of domain
    m2_right = (y_arr2[-1] - y_arr2[-2])/(x_arr[-1] - x_arr[-2] );
    b2_right = y_arr2[-1] - m2_right*x_arr[-1];
    
    # print out results of this
    if( help_debug):
        print("Left side fit of line 1 is: y = "+str(m1_left)+"x + "+str(b1_left) );
        print("Right side fit of line 1 is: y = "+str(m1_right)+"x + "+str(b1_right) );
        print("Left side fit of line 2 is: y = "+str(m2_left)+"x + "+str(b2_left) );
        print("Right side fit of line 2 is: y = "+str(m2_right)+"x + "+str(b2_right) );
        
    # look (very finely) for x that satisfies both
    x_vals, delta_x = np.linspace( x_arr[0], x_arr[-1], pow(10, fineness), retstep = True); # return step size also
    delta_y =  abs(delta_x); # get max possible y val change per step as motivation for equality tolerance

    # iter over x vals looking for intersection
    # check both line segments for each line against each other
    # all lines should be either parallel or not so no false positives
    for x_val in x_vals:

        # check equality of m1_left and m2_left, also make sure potential intersection point is in the y range
        #print( m1_left*x_val + b1_left );
        if( (abs( m1_left*x_val + b1_left - (m2_left*x_val + b2_left) ) < delta_y) and (m1_left*x_val + b1_left > 0 )  ):
               
            # equality means this is the intersection
            x_int = x_val;
            y_int = m1_left*x_val + b1_left;
            
            # return coordinates of intersection
            return x_int, y_int;
        
        # now check m1 left with m2 right
        elif( (abs( m1_left*x_val + b1_left - (m2_right*x_val + b2_left) ) < delta_y) and (m1_left*x_val + b1_left > 0 ) ):
               
            # equality means this is the intersection
            x_int = x_val;
            y_int = m1_left*x_val + b1_left;
            
            # return coordinates of intersection
            return x_int, y_int;
        
        elif( (abs( m1_right*x_val + b1_right - (m2_left*x_val + b2_left) ) < delta_y) and (m1_right*x_val + b1_right > 0 ) ):
               
            # equality means this is the intersection
            x_int = x_val;
            y_int = m1_right*x_val + b1_right;
            
            # return coordinates of intersection
            return x_int, y_int;
        
        elif( (abs( m1_right*x_val + b1_right - (m2_right*x_val + b2_right) ) < delta_y) and (m1_right*x_val + b1_right > 0 ) ):
               
            # equality means this is the intersection
            x_int = x_val;
            y_int = m1_right*x_val + b1_right;
            
            # return coordinates of intersection
            return x_int, y_int;
     
    # return Nones if not found   
    return None, None;
        
    
    
def GetTransitionEnergies( levels, labels):
    '''
    Function for determining the photon transition energies of a given set of energy levels. Due to spin conservation
    on the 4 band system only 2 transitions are allowed, up --> up and down --> down. Note that all transitions are absorptions.
    
    Has to know spin of each band so will use labels generated by GetBandLabels, see below
    
    Args:
    :param levels, 1d array, contains the 4 energy level vals at this spot in param space
    :param labels, 1d array of the labels corresponding to each level (label of levels[i] = labels[i] )
    
    returns tuple of
        dE_up, energy difference between up bands
        dE_down, likewise for down bands
    ''' 
    
    # seperate levels out into up and down
    # get containers
    E_up = [];
    E_down = [];
    
    # iter over levels and classify up or down
    for eig_i in range(len(levels)):
        
        if( "uparrow" in labels[eig_i] ): # this is a spin up level
            E_up.append( levels[eig_i]);
            
        elif("downarrow" in labels[eig_i]):
            E_down.append( levels[eig_i]);
            
        else:
            raise ValueError("Incorrect energy level labels in GetTransitionEnergies.")
        
    # check that this worked
    if( len(E_down) != 2 or len(E_up) != 2 ): # there's an error
        raise ValueError("Unable to separate levels by spin");
            
    # now determine transition energies
    dE_up = abs( E_up[0] - E_up[1] );
    dE_down = abs( E_down[0] - E_down[1] );
    
    # return these results
    return dE_up, dE_down;

    #### end GetTransitionEnergies
    
    
def GetBandLabels(H, bench_levels, phi, no_phi, material, help_debug = False):
    '''
    New methodology for labelling bands that replaces what we used in GetBandLabels. 
    
    The new methodology used here is to classify bands by their slope with respect to different inputs. It is easily seen from the hamiltonian that the valence
    band always decreases with Eg and the conduction band always increases with Eg. Similarly up bands always decrease with (a2 or b2) and down bands always increase
    with (a2 or b2) where it is a2 for valence bands and b2 for conduction bands. So the method is to slightly increase the test parameters and use the sign of 
    the resulting change to classify the bands.
    
    :param H: hamiltonian object representing the state of the system before changes
    :param phi: angle between longitudinal direction and magnetization
    :param material: string, tells what material this is, because naming scheme is different for different materials
    :param no_phi: bool, tells not to add phi value to labels, as this is redundant when we are plotting valleys on separate axes
    :param help_debug: bool, defaults to false. Set to True if this function is having issues to help debug
    '''
    
    if( help_debug ):
        print("\n\n");
        print(10*"*"+" Entering NewGetBandLabels.\nOptional help_debug argument set to true so function will print run info.")
        print(H);
        print("phi = "+str(phi)+" and exchange params are ", H.Exchange.a1,H.Exchange.a2,H.Exchange.b1,H.Exchange.b2, );
        
    
    def increment(x):
        # method for increasing whatever quantity we are finding slope with respect to
        # return number increased by a small amount
        
        # get 1 percent increase
        delta = abs(x)/100; 
        
        # possible issue is delta = 0 bc x was 0
        if( delta == 0):
            delta += 0.001; # hard code this in
        
        return x + delta; # always moves in positive direction
    
    # determine naming, color scheme based on the material
    # the key difference here is that according to Eqs 9 and 10 of Bauer Zawadski SST 1992, 
    # in PbSe C+ --> C up and C- --> C down, but in PbTe C- --> C up and C+ --> C down
    # see that paper for details. Remember that the cos(theta) mixing parameter dominates
    if( "Se" in material): # follow PbSe symmetry
        labels_list = [r"$|L_{6}^{+} \uparrow \rangle$", r"$|L_{6}^{+} \downarrow \rangle$", r"$|L_{6}^{-} \uparrow \rangle$", r"$|L_{6}^{-} \downarrow \rangle$"];
        colors_list = ["orangered", "orange", "green", "greenyellow" ];
        
    elif( "Te" in material): # follow PbTe symmetry
        print('nice')
        labels_list = [r"$|L_{6}^{+} \uparrow \rangle$", r"$|L_{6}^{+} \downarrow \rangle$", r"$|L_{6}^{-} \downarrow \rangle$", r"$|L_{6}^{-} \uparrow \rangle$"];
        colors_list = ["orangered", "orange", "greenyellow", "green" ];
        
    else: # default to PbSe symmetry
        labels_list = [r"$|L_{6}^{+} \uparrow \rangle$", r"$|L_{6}^{+} \downarrow \rangle$", r"$|L_{6}^{-} \uparrow \rangle$", r"$|L_{6}^{-} \downarrow \rangle$"];
        colors_list = ["orangered", "orange", "green", "greenyellow" ];
    
    # get the original eigenvalues of H
    eigvals = H.Eigvals();
    if( help_debug):
        print("old eigvals: ", eigvals);
    
    # for each eigenvalue we want the sign of its slope with respect to Eg, a2, and b2
    # make container vars for these quantitites
    Eg_slope = np.zeros(len(eigvals));
    a2_slope = np.zeros(len(eigvals)); 
    b2_slope = np.zeros(len(eigvals));
    
    # also make container for results, which are labels and colors
    labels = np.full(len(eigvals), 'dummy_label_goes_here_but_needs_filler_text_to_work'); # has to be super long
    colors = np.full(len(eigvals), 'dummy_color_here');
    
    # now iter thru eigvals to get Eg slopes
    for eig_i in range( len(eigvals) ):

        # get the change in this eigval as we vary Eg slightly
        # this requires making a new hamiltonian object
        #print(H.Eg, increment(H.Eg) )
        H_new_inputs = increment(H.Eg), H.mu_Bohr, H.B_vec, H.vt, H.vl ; # tuple of all inputs from old H except Eg is slightly increased
        H_new = hamiltonian.Hamiltonian(*H_new_inputs); # make new hamiltonian object
           
        # must also add exchange effects to new H
        # first get old exchange matrix from old hamiltonian
        E_old = H.Exchange;
        H_new.Add(E_old);
        
        # get eigenvalues of modified H
        new_eigvals = H_new.Eigvals();
        
        # get sign of slope by subtracting like eigenvals
        Eg_slope[eig_i] = np.sign(new_eigvals[eig_i] - eigvals[eig_i]);
        
        #print helpful info if asked
        if( help_debug ):
            print("Band gap "+str(H.Eg)+" --> "+str( H_new.Eg ) );
            print(str(eig_i)+"th eigenvalue is "+str(new_eigvals[eig_i])+", was "+str( eigvals[eig_i] ) );
        
        # delete H new object to prevent any errors    
        del H_new;
    
    if( help_debug):    
        print( "Signs of slopes w/r/t Eg, (should be 1 twice, -1 twice): ", Eg_slope);    
    
    # iter again thru eigvals to get a2 slopes
    for eig_i in range(len(eigvals) ):
        
        # get the change in this eigval as we vary a2 slightly
        # this requires making a new hamiltonian object
        H_new_inputs = H.Eg, H.mu_Bohr, H.B_vec, H.vt, H.vl ; # H inputs are exactly same
        H_new = hamiltonian.Hamiltonian(*H_new_inputs); # make new hamiltonian object
        
        # now get change with a2
        # first get old exchange matrix from old hamiltonian
        E_old = H.Exchange;
        
        # modify input for a2
        E_new_inputs = E_old.phi, E_old.a1, increment(E_old.a2), E_old.b1, E_old.b2, E_old.long_valley_shift ; # same inputs as before but a1 slightly increased
        E_new = hamiltonian.Exchange(*E_new_inputs);
                    
        # now add to new H to get new eigvals
        H_new.Add(E_new);
        new_eigvals = H_new.Eigvals();

        # get sign of slope by subtraction
        a2_slope[eig_i] = np.sign(new_eigvals[eig_i] - eigvals[eig_i] );
        
        # prove modifications are correct if asked
        if( help_debug ):
            print("a2 "+str(E_old.a2)+" --> "+str(E_new.a2) );
            print(str(eig_i)+"th eigenvalue is, "+str(new_eigvals[eig_i])+", was "+str( eigvals[eig_i]) );
            
        del H_new;
    
    if( help_debug):    
        print("Signs of slopes w/r/t a2, (should be 1 once, -1 once, 0 twice): ", a2_slope);    
    
    # finally iter to get slopes w/r/t b2
    for eig_i in range(len(eigvals) ):
        
        # get the change in this eigval as we vary b2 slightly
        # this requires making a new hamiltonian object
        H_new_inputs = H.Eg, H.mu_Bohr, H.B_vec, H.vt, H.vl ; # H inputs are exactly same
        H_new = hamiltonian.Hamiltonian(*H_new_inputs); # make new hamiltonian object
        
        # now get change with b2
        # first get old exchange matrix from old hamiltonian
        E_old = H.Exchange;
        
        # modify input for b2 (slightly increase)
        E_new_inputs = E_old.phi, E_old.a1, E_old.a2, E_old.b1, increment(E_old.b2), E_old.long_valley_shift ; # same inputs as before but a1 slightly increased
        E_new = hamiltonian.Exchange(*E_new_inputs);
                    
        # now add to new H to get new eigvals
        H_new.Add(E_new);
        new_eigvals = H_new.Eigvals();

        # get sign of slope by subtraction
        b2_slope[eig_i] = np.sign(new_eigvals[eig_i] - eigvals[eig_i] );
        
        # prove modifications are correct if asked
        if( help_debug ):
            print("b2 "+str(E_old.b2)+" --> "+str(E_new.b2) );
            print(str(eig_i)+"th eigenvalue is, "+str(new_eigvals[eig_i])+", was "+str( eigvals[eig_i]) );
            
        del H_new;
    
    if( help_debug):    
        print("Signs of slopes w/r/t b2, (should be 1 once, -1 once, 0 twice): ", b2_slope);
    
    # now assign based on both slopes
    # we have used Bauer's basis of V+, V-, C+, C-, where V is L6+ and C is L6- physically
    # and the fact that the cos(theta) of the mixing angle dominates the spin mixing (see Bauer equations 8, 9, 10)
    # to determine the basis. For PbSe this turns out to be V up, V down, C up, C down, while for
    # PbTe it is V up, V down, C up, C down
    # based on the form of the complete hamiltonian the signs of slopes (Eg, a2 or b2) should correspond to eigenvalues
    # (-1, -1)  0  0  0
    # 0  (-1, 1)  0  0
    # 0  0  (1, -1)  0
    # 0  0  0  (1, 1)
    # where second number is wrt a2 for valence bands (upper two diagonals) and b2 for conduction bands (lower two diagonals)
    
    # iter over eigvals and classify accordingly
    for eig_i in range(len(eigvals) ):
        
        # V+ (V up) case
        if( Eg_slope[eig_i] == -1 and a2_slope[eig_i] == -1 and b2_slope[eig_i] == 0):
            labels[eig_i] = labels_list[0] ;
            colors[eig_i] = colors_list[0] ;  
            
        # V- (V down) case
        elif( Eg_slope[eig_i] == -1 and a2_slope[eig_i] == 1 and b2_slope[eig_i] == 0):
            labels[eig_i] = labels_list[1] ;
            colors[eig_i] = colors_list[1] ; 
            
        # C+ (C up for PbSe) case
        if( Eg_slope[eig_i] == 1 and a2_slope[eig_i] == 0 and b2_slope[eig_i] == -1):
            labels[eig_i] = labels_list[2] ;
            colors[eig_i] = colors_list[2] ;   
            
        # C- (C down for PbSe) case
        if( Eg_slope[eig_i] == 1 and a2_slope[eig_i] == 0 and b2_slope[eig_i] == 1):
            labels[eig_i] = labels_list[3] ;
            colors[eig_i] = colors_list[3] ; 
            
    # now have to put eigvals in correct order (ie asending order on far left of plot)
    # create containers for correctly ordered labels, colors
    bench_labels = np.full(len(bench_levels), 'dummy_label_goes_here_but_needs_filler_text_to_work'); # has to be super long
    bench_colors = np.full(len(bench_levels), 'dummy_color_here');
    
    # iter over bench levels (indexing is same as in bench labels and bench colors)
    for bench_i in range(len(bench_levels)):
        
        # iter thru wrong ordered levels to get matching index
        for i in range(len(eigvals)):
            
            # check for match
            if( bench_levels[bench_i] == eigvals[i]): # if so it is a match, place labels and colors accordingly
                bench_labels[bench_i] = labels[i];
                bench_colors[bench_i] = colors[i];
                
    # add phi vals to labels if asked
    if( not no_phi ): # default case, phi vals will be added
        
        # have to make a new list bc string arrays are garbage
        new_labels = [];
        for my_label in bench_labels: # iter thru old labels and add phi
            new_labels.append( my_label + ", $\phi$ = "+str(phi)[:4]+" [rad]" );   
            
        # set old name as new list 
        bench_labels = np.array(new_labels);
                        
    return bench_labels, bench_colors;
    
    #print(labels, colors);
    return eigvals, bench_labels, bench_colors;

    #### end GetBandLabels
   
    
      
                                       
                    
                    
################################################################################
# test code here
if __name__ == "__main__":
    
    # test the GetBandLabels function
    # define the H inputs for the benchmark H
    Eg = 0.15; # band gap in eV
    mu_Bohr = 5.79e-5; # bohr magneton in eV/T
    B_vec = np.zeros(3); # B field in T
    vt = 1; # transverse velocity matrix element
    vl =1; # longitudinal velocity matrix element
    H_inputs = (Eg, mu_Bohr, B_vec, vt, vl);
    H_bench = hamiltonian.Hamiltonian(*H_inputs);
    print(H_bench);
    
    # add in exchange
    phi = 0; 
    a1, a2, b1, b2 = 0.01, 0.02, 0.01, 0.02 ; # define exchange params
    E_inputs = (phi, a1, a2, b1, b2);
    E = hamiltonian.Exchange(*E_inputs );
    H_bench.Add(E); #add exhange effects to existing hamiltonian

    
    
