'''
Created on Feb 6, 2020

@author: Christian

Module to create an effective Hamiltonian using matrices.py

Then diagonalize this hamiltonian to find dispersion (kx, ky != 0) energy values)
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

from Finders import main
from Finders import functions
from Finders import matrices

def E_k(k_x, k_y, func_inputs, finder_inputs, is_tilde= False, is_verbose = False):
    '''
    Function to find the eigenvalues of the perturbation matrix H_k
    H_k is a perturbation for nonzero kx, ky so these represent the dispersion bound state energies
    
    :param k_x: double, momenta in the x
    :param k_y: double, momenta in the y
    :param func_inputs: functions.InputParams object, has all relevant experimental parameters
    :param finder_inputs: main.FinderParams object, tells how to find eigenvalues at kx, ky = 0
    :param is_tilde, bool, tells whether we are diag'zing the matrix H_tilde (see matrices.py)
    :param is_verbose: bool, tells whether to print out extra things for troubleshooting
    
    returns eig_vals, 1d np array of eigenvalues, each should be doubly degenerate
    '''
    
    if(is_verbose):
        print("Entering E_k");
        
    # determine degeneracy based on which matrix
    if( not is_tilde): # we are doing H_eff, so eigvals are doubly degenerate
        degeneracy = 2;
    else: #we are doing H_tilde, so single degeneracy
        degeneracy = 1;
        
    # get hamiltonian for dispersion
    # sort by what matrix we are diagonalizing
    if( is_tilde): # matrix is n_eigvals x n_eigvals, diagonal/antidiagonal matrix H_tilde
        H_k = matrices.H_tilde(k_x, k_y, func_inputs, finder_inputs, is_verbose);
    else: # matrix is 2n_eigvals x 2n_eigvals full perturbation matrix
        H_k = matrices.H_eff(k_x, k_y, func_inputs, finder_inputs, is_verbose);
    
    # check that it is hermitian
    H_hermitian = H_k.conj().T  # get hermitian conjugate
    # iter over array elements
    for i in range( H_k.shape[0]):
        for j in range( H_k.shape[1] ):
            if( H_k[i,j] != H_hermitian[i,j]):
                print("Error in E_k: H_eff not hermitian.");
                print( i, j, H_k[i,j],H_hermitian[i,j]);
                return;
    
    # get eigenvals (and eigenvectors) for this hamiltonian
    eig_vals, eig_vecs = np.linalg.eigh( H_k );
    
    if( is_verbose and degeneracy == 2): # print out eigenvalues by degeneracy
        
        # iter over degenerate pairs of eigenvalues
        for i in range( int( len( eig_vals )/degeneracy) ):
        
            pr_string = "Eigenvalue pair 1: "+str(eig_vals[degeneracy*i])+", "+str(eig_vals[degeneracy*i + 1]);
            print( pr_string);
    
    return eig_vals;

#  End E_k #################################################################


def Vary_k(k_start, k_stop, func_inputs, finder_inputs, is_tilde = False, is_verbose = False ):
    '''
    Function to find dispersion (k nonzero) values for bound state energies
    for different values of k   

    :param k_start: double, start of k range in nm^-1
    :param k_stop: double, end of k range in nm^-1   
    :param func_inputs: functions.InputParams object, has all relevant experimental parameters
    :param finder_inputs: main.FinderParams object, tells how to find eigenvalues at kx, ky = 0
    :param is_tilde, bool, tells whether we are diag'zing the matrix H_tilde (see matrices.py)
    :param is_verbose: bool, tells whether to print out extra things for troubleshooting 

    returns tuple of k_arr, all_E_levels
    all_E_levels is a list of 1D np arrays
    each array corresponds to a list of energies for a given level
    k_arr, 1D np array of k vals
    k_arr[i] is k val corresponding to energies located at [i] of each np array in all_E_levels
    These

    '''  
    
    # determine degeneracy based on which matrix
    if( not is_tilde): # we are doing H_eff, so eigvals are doubly degenerate
        degeneracy = 2;
    else: #we are doing H_tilde, so single degeneracy
        degeneracy = 1;
    
    n_levels = degeneracy*len( main.AllRoots(func_inputs, finder_inputs) )

    # make array of k_vals
    k_arr = np.linspace(k_start, k_stop, 10); 

    #make array for each level to store E values, store all in a list
    all_E_levels = []; 
    for n in range(n_levels):
        # each array should have len(k_arr) spots
        all_E_levels.append( np.zeros( len(k_arr) ) );

    # iter thru k vals 
    for i_k in range( len( k_arr) ):

        #determine k_x, k_y, symmetrically for now
        k_val = k_arr[i_k];
        k_x = np.sqrt( k_val*k_val/2);
        k_y = k_x;      

        # get array of energies from E_k
        energies_k = E_k(k_x, k_y,  func_inputs, finder_inputs, is_tilde, is_verbose); 
        
        # sort and place into array by level
        # remember i_k is index of k_val in array, n is index of level in all_E_levels
        for n in range(n_levels):
            all_E_levels[n][i_k] = energies_k[n] ;
             

    return k_arr, all_E_levels;  

# end Vary_k ############################################################


def PlotDispersion(k_arr, all_E_levels, exp_inputs, is_tilde = False):
    

    # determine degeneracy based on which matrix
    if( not is_tilde): # we are doing H_eff, so eigvals are doubly degenerate
        degeneracy = 2;
    else: #we are doing H_tilde, so single degeneracy
        degeneracy = 1;
        
    # determine eigenenergies at k = 0 in order to determine parity
    eigenenergies=[];
    for i in range(0, len( all_E_levels ), degeneracy ): #get all levels
    
        E_level = all_E_levels[i];
        eigenenergies.append( E_level[0]);
    
    # convert to np array
    eigenenergies = np.array( eigenenergies);

    #iter through energy arrays and plot
    for E_level in all_E_levels:
        
        #### determine labels and line formatting according to TIS status and parity
        
        # determine if TIS ( ie if -2|delta| < E < 0 at k = 0
        is_TIS = False;
        if( E_level[0] < 0 and E_level[0] > 2*exp_inputs.delta):
            is_TIS = True
            
        # determine parity of state
        is_odd = matrices.Parity_E(E_level[0], eigenenergies, exp_inputs.delta);
        
        #plot and make labels accordingly
        if(is_TIS): # make plot the line red
            
            if(is_odd): # plot odd solutions dashed
                plt.plot(k_arr, E_level, color = 'red', linestyle = 'dashed');
            else: # not dashed
                plt.plot(k_arr, E_level, color = 'red') ;
                
        else: # not a TIS, plot in black
            
            if(is_odd): # plot odd solutions dashed
                plt.plot(k_arr, E_level, color = 'black', linestyle = 'dashed');
            else: # not dashed
                plt.plot(k_arr, E_level, color = 'black');  
        
    # make titles, labels etc
    title = "Dispersion E($k_{x}, k_{y}$) symmetric in k"
    y_label = "Energy [meV]";
    x_label = "k [nm^-1]";

    # make legend
    red_line = mlines.Line2D([],[],color='red',label='Topological interface states');
    black_line = mlines.Line2D([],[], color = 'black', label = 'Propagative states')
    dashed_line = mlines.Line2D([],[], color = 'black', linestyle = 'dashed', label = 'Denotes odd states')
    plt.legend(handles=[red_line, black_line, dashed_line], loc="upper left");
    
    #make a textbox showing experimental params
    exp_params_string = ""
    
    # add info of experimental params to box as desired
    for line in exp_inputs.DeclareSelf(): #this returns a list of strings describing each exp param
        exp_params_string += line;
        
    # put the textbox on the plot
    plt.text( 0.33,-70, exp_params_string);
        
    # format the resulting plot
    plt.suptitle(title); #set title
    plt.xlabel(x_label);
    plt.ylabel(y_label);
    plt.grid(True, 'both', 'both');
    y_start, y_stop = plt.ylim(); #get y axis limits
    
    # if asked for by setting linear to True, fit the slope of each dispersion level for high k
    for i in range(0, len( all_E_levels ), degeneracy ):
        
        E_level = all_E_levels[i];
        
        # determine parity of state
        print( "is_odd = "+str(matrices.Parity_E(E_level[0], eigenenergies, exp_inputs.delta) ) );
        
        #determine the velocity from the slope
        slope = (E_level[-1] - E_level[-2])/(k_arr[-1] - k_arr[-2] );
        v_D = slope / exp_inputs.h_bar / pow(10,14);
        
        # print out 
        print("Level "+str(i)+": velocity is "+str(v_D)+" nm/s." ); 
        
        # put slope on plot as a text box
        text_string = r"$\Delta E/ \, \hbar \Delta k$ = "+str( np.ceil(v_D*100 )/100 )+"e14 nm/s"
        
        # determine location and whether to plot by looking at E relative to delta
        if( E_level[0] > exp_inputs.delta and E_level[0] < 0): # this is the upper TIS
            x_box, y_box = 0.4, 70
            plt.text( x_box, y_box, text_string);
        elif( E_level[0] < exp_inputs.delta and E_level[0] > 2*exp_inputs.delta): # this is the lower TIS
            x_box, y_box = 0.4, -130
            plt.text( x_box, y_box, text_string);
        
        
    #show all on same plot
    plt.show();

if __name__ == "__main__":
################################################################################
#   Test code here
################################################################################

     #### define experimental parameters
    delta = -25; #band gap in meV
    V = 175; # potential in meV
    d = 12; # well thickness in nm
    v_z = 4.3*pow(10,14) # matrix element in nm/s
    v_x = 2.3*pow(10,14) # matrix element in nm/s
    v_y = 2.3*pow(10,14) # matrix element in nm/s
    
    # pass to InputParams object
    exp_inputs = functions.InputParams(delta, V, d, v_z, v_x, v_y);
    
    #### define finder parameters
    my_method = 'bisect';
    my_interval = ( -V, V - 2*abs(delta) );
    
    #pass to FinderParams object
    find_inputs = main.FinderParams(method = my_method, interval = my_interval, is_verbose=True)
    
    # define range of k_vals
    k_start = 0;
    k_stop = 0.5;
    
    # determine which matrix to use
    is_tilde = true;
    
    # get dispersion energy levels as function of k
    k_vals, all_E_levels = Vary_k(k_start, k_stop, exp_inputs, find_inputs, is_tilde, True);
    
    # plot the dispersion
    PlotDispersion(k_vals, all_E_levels, exp_inputs, True);
