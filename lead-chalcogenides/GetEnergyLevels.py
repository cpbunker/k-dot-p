'''
Created on Mar 29, 2020

@author: Christian

energylevels.py was a module for everything related to finding and plotting or otherwise visualizing the energy levels of the 
physical system defined by the combined Hamiltonian and Exchange matrix objects defined in the hamiltonian.py program. I have
since split this module into this, GetEnergyLevels.py, and EnergyLevels.py. This module contains the simpler functions and often
features calls to the functions of EnergyLevels.py

In other words, if you are looking for a problem:
    FIRST check the runner file (PbSnSe or equivalent)
    SECOND check this module
    ONLY THEN look in EnergyLevels.py

The functions in this module:
    Determine and plot the energy levels of the system as a function of the bulk k.p parameters such asenergy gap Eg and composition x
    Determine and plot the energy levels of the system with only one exchange parameter free, the others fixed
    Plot the results of the conditional param space search (see docstrings of ParamSS, ProcessCondition, PlotCondition dunctions in 
    EnergyLevels and SearchProcessPlot in this module if you don't know what this means
    Do the last in the reduced dimensional case  A, B fixed
'''

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.collections as collections

from Valleytronics import hamiltonian
from Valleytronics import EnergyLevels


################################################################################
# Functions to plot energy levels

def EgSweepBothValleys(H_inputs, E_inputs, Eg_lim, material, Eg_bench, x_inputs, both ):
    '''
    Combines the process of plotting the levels in both valleys against Eg into one function, so user doesn't have to call EgSweep w/
    separate phi values anymore 
    
    :param H_inputs: tuple of the inputs for the bulk k.p hamiltonian (Eg, muBohr, Bvec, vt, vl), except
                Eg is just a dummy value since it is our dependent variable
    :param E_inputs: tuple of the inputs to construct the exchange matrix (phi, a1, a2, b1, b2)
    :param Eg_lim: tuple of the start and stop values for the Eg sweep
    :param material: string of the name of the material for the title of the plot
    :param Eg_bench, double, the energy gap to benchmark (ie determine valence/conduction) at
    :param x_inputs, tuple of params for plotting against x (composition) domain if asked
                    vs_x, bool, tells to plot against x if true
                    x_lim, tuple of start and stop point for x, very often (0,1) since those are the physical limits of x
                    temperature, double, temperature in K for determing Eg( x, T)
                    Eg_func, points to function which has dependence of Eg on x, T
    :param both: bool, tells whether to force both valley levels onto same axis (True) or keep separate (False)
    
    returns None
    '''
    
    # unpack the x inputs
    # see docstirng for what these are
    # this tuple is used if we are plotting with composition x as the independent variable
    vs_x, x_lim, temperature, Eg_func = x_inputs
     
    # determine whether plotting on seperate axes or together
    if( both): # plot together
        
        # only one axis needed
        fig, ax1 = plt.subplots();
        
        # alias the axis
        ax2 = ax1;
        
        # we will want to have phi in the legend
        no_phi = False;
        
        # also determine legend placement
        legend_coef = 0.6;
        
    else: # plot seperately
        
        # create two seperate axes for long and oblique valleys        
        fig, (ax1, ax2) = plt.subplots(ncols=2, figsize = (8,4), sharey = True );
        
        # we do not want phi's in the legend, it is redundant
        no_phi = True;
        
        # also determine legend placement
        legend_coef = 0.95;
    
    # package E inputs for long and oblique
    E_long = (0, *tuple(E_inputs[1:])) # this way of list slicing, then scattering ensures all elements of the E_inputs tuple after 1st are passed
    E_obl = (70.5*np.pi/180, *tuple(E_inputs[1:]))
    
    # determine whether to plot against Eg domain, using EgSweep, or x domain, using EgSweep_vs_x (if vs_x True)
    # find, plot levels vs E gap if that is independent variable
    if( not vs_x ): # we are plotting vs gap not composition
    
        # for the longitudinal valley
        EgSweep(H_inputs, E_long, Eg_lim, material, Eg_bench, ax1, no_phi ); 
    
        # similarly for the oblique valley    
        EgSweep(H_inputs, E_obl, Eg_lim, material, Eg_bench, ax2, no_phi );
        
    else: # composition is the independent variable, have to call EgSweep_v_x
        
        # for the longitudinal valley
        EgSweep_vs_x(H_inputs, E_long, Eg_func, x_lim, temperature, material, Eg_bench, ax1, no_phi );
        
        # for the oblique valley
        EgSweep_vs_x(H_inputs, E_obl, Eg_func, x_lim, temperature, material, Eg_bench, ax2, no_phi );
        
    # format the plot
    fig.suptitle(material + " bulk E levels \n", fontsize = 14);
    
    # add text box for exchange params
    #plt.text("a1 = "+str(E_inputs[1])[:6]+", "+"a2 = "+str(E_inputs[2])[:6]+",\n"+"b1 = "+str(E_inputs[3])[:6]+", "+"b2 = "+str(E_inputs[4])[:6]+", " , loc = "right", fontsize = 10  );
        
    # make the legend
    chartBox = ax2.get_position();
    ax2.set_position([chartBox.x0, chartBox.y0, chartBox.width*legend_coef, chartBox.height]);
    ax2.legend(loc='upper center', bbox_to_anchor=(1.45, 0.8), shadow=True, ncol=1);
    leg2 = ax2.get_legend();
    leg2.set_draggable(True); # makes the legend box draggable !
    
    # textbox below the legend for exchange param vals
    textstring = "a1 = "+str(E_inputs[1])[:6]+",\n"+"a2 = "+str(E_inputs[2])[:6]+",\n"+"b1 = "+str(E_inputs[3])[:6]+",\n"+"b2 = "+str(E_inputs[4])[:6]+", ";
    text_x, text_y, = 1.03*ax2.get_xlim()[1], ax2.get_ylim()[0] #+ ax2.get_ylim()[1])/3;
    ax2.text(text_x, text_y, textstring);
       
    # show and return 
    plt.show();                   
    return;

    #### end EgSweepBothValleys
    
                      

def EgSweep(H_inputs, E_inputs, Eg_lim, material, Eg_bench, ax, no_phi, help_debug = False):
    '''
    Finds the dependence of the energy levels on the energy gap of the material, Eg
    
    :param H_inputs: tuple, the inputs for the bulk k.p hamiltonian (Eg, muBohr, Bvec, vt, vl), except
                Eg is just a dummy value since it is our dependent variable
    :param E_inputs: tuple, the inputs to construct the exchange matrix (phi, a1, a2, b1, b2)
    :param Eg_lim: tuple, the start and stop values for the Eg sweep
    :param material: string, the name of the material for the title of the plot
    :param Eg_bench: double, the energy gap to benchmark (ie determine valence/conduction) at
    :param ax: matplotlib axis object on which to plot the levels, used in conjunction with EgSweepBothValleys
        to be able to easily plot valleys on same or different axes
    :param no_phi: bool, tells to suppress printing of phi values in legend, used for seperate axis plotting
    :param help_debug: bool, defaults to False, tells function to execute many extra print statements so that user can check for bugs
    '''
    
    # get, delete a version of the k.p bulk hamiltonian
    H = hamiltonian.Hamiltonian(*H_inputs);
    
    # print info about H
    if( help_debug):
        print("\n\n");
        print(10*"*"+" Entering EgSweep.\nOptional help_debug argument set to true so function will print run info.");
        print(H);
    
    # define number of points to use in parameter space using H
    N_param = 100;
    N_levels = H.N;
    del H;
    
    # define the return variable
    # 2d array of energy levels for different Eg
    E_levels = np.zeros((N_levels, N_param));
    
    # construct the exchange matrix
    # this is the same throughout since we are only varying Eg
    E = hamiltonian.Exchange(*E_inputs);
    phi = E_inputs[0];
    
    # vary over the energy gap
    Eg_arr = np.linspace( Eg_lim[0], Eg_lim[1], N_param );
    
    # iter over the energy gap values
    for i in range(len(Eg_arr)):
        
        # for this individual energy gap
        Eg = Eg_arr[i];
        
        # get the energy levels
        # first construct the bulk matrix
        H_inputs = (Eg, H_inputs[1], H_inputs[2], H_inputs[3], H_inputs[4] );
        H = hamiltonian.Hamiltonian(*H_inputs);
        
        # combine with the exchange matrix
        H.Add(E);
        
        # first run thru, get eigvecs
        if( i == 0):
            eigvecs = H.Eigvecs();
            
        # get eigenvalues with eigenvector method
        for eig_i in range(N_levels):
            
            # get the appropriate eigval and place in holder
            E_levels[eig_i, i] = H.GetEigenval(eigvecs[eig_i]);
            
        # delete H before the next run thru
        del H;
        
    #### end of loop over Eg vals
    
    # create H at  benchmark Eg (defaults to 150 meV) to determine band labels
    H_inputs = ( Eg_bench, H_inputs[1], H_inputs[2], H_inputs[3], H_inputs[4] );
    H = hamiltonian.Hamiltonian(*H_inputs);
    H.Add(E);
    
    # get levels from this H
    Hi_E_levels = np.zeros(N_levels);
    for eig_i in range(N_levels):
        Hi_E_levels[eig_i] = H.GetEigenval(eigvecs[eig_i]);
        
    # determine labels, colors approproately
    labels, colors = EnergyLevels.GetBandLabels(H.MakeCopy(), Hi_E_levels, phi, no_phi, material );
    
    # set title and linestyle according to phi
    if( phi == 0):
        style = '--'; # plots dashed
        ax.set_title("Longitudinal valley (--)", loc = 'left', fontsize = 10);
        ax.set_xlabel("$E_{g}$ [eV]");
        ax.set_ylabel("Energy [eV]");
        
    else: #oblique valley
        style = 'solid'; # plots solid
        ax.set_title("Oblique valley ($ -$)", loc = 'right', fontsize = 10);
        ax.set_xlabel("$E_{g}$ [eV]");
        #ax.set_ylabel("Energy [eV]");
        
    # iter over each level and plot
    if(help_debug):
        print("phi = "+str(phi)+" endpoints");
    for eig_i in range(N_levels):
        
        # plot level with appropriate label, color
        ax.plot(Eg_arr, E_levels[eig_i], color = colors[eig_i], linestyle = style, label = labels[eig_i]);
        
        # get endpoint to calculate energy differences
        if( help_debug):
            print("E_"+str(eig_i)+" = "+str(E_levels[eig_i, -1] ) );
    
    return;

    #### end EgSweep
    
    
def EgSweepTransition(H_inputs, E_inputs, Eg_lim, material, Eg_bench):
    '''
    Finds the energy levels as a function of the band gap Eg. Then finds and plots the spin-spin transition energies between levels.
    
    :param H_inputs: tuple of the inputs for the bulk k.p hamiltonian (Eg, muBohr, Bvec, vt, vl), except
                Eg is just a dummy value since it is our dependent variable
    :param E_inputs: tuple of the inputs to construct the exchange matrix (phi, a1, a2, b1, b2)
    :param Eg_lim: tuple of the start and stop values for the Eg sweep
    :param material: string of the name of the material for the title of the plot
    :param Eg_bench: double, the energy gap at point we want to label energy levels at
    '''
    
    
    #### iter over Eg, get transition energies,, and plot for both phi's
    
    # after this is done we want to get intersections
    # need containers for the lines we get
    central_lines = [];
    oblique_lines = [];
    
    # iter over both phis
    for phi in [0, 70.5*np.pi/180]:
    
        # get, delete a version of the k.p bulk hamiltonian
        H = hamiltonian.Hamiltonian(*H_inputs);
        
        # define number of points to use in parameter space using H
        N_param = 100;
        N_levels = H.N;
        del H;
        
        # define the return variable
        # 2d array of energy levels for different Eg
        E_levels = np.zeros((N_levels, N_param));
        
        # construct the exchange matrix
        # first reconstruct inputs with appropriate phi
        E_inputs = (phi, *tuple(E_inputs[1:]) );
        E = hamiltonian.Exchange(*E_inputs);
        
        # get linspace to vary over the energy gap 
        Eg_arr = np.linspace( Eg_lim[0], Eg_lim[1], N_param );
        
        # iter over the energy gap values
        for i in range(len(Eg_arr)):
            
            # for this individual energy gap
            Eg = Eg_arr[i];
            
            # get the energy levels
            # first construct the bulk matrix
            H_inputs = (Eg, H_inputs[1], H_inputs[2], H_inputs[3], H_inputs[4] );
            H = hamiltonian.Hamiltonian(*H_inputs);
            
            # combine with the exchange matrix
            H.Add(E);
            
            # first run thru, get eigvecs
            if( i == 0):
                eigvecs = H.Eigvecs();
                
            # get eigenvalues with eigenvector method
            for eig_i in range(N_levels):
                
                # get the appropriate eigval and place in holder
                E_levels[eig_i, i] = H.GetEigenval(eigvecs[eig_i]);
                
            # delete H before the next run thru
            del H;
            
        #### end of loop over Eg vals
        
        # create H at  benchmark Eg (defaults to 150 meV) to determine band labels
        H_inputs = ( Eg_bench, H_inputs[1], H_inputs[2], H_inputs[3], H_inputs[4] );
        H = hamiltonian.Hamiltonian(*H_inputs);
        H.Add(E);
        
        # get levels from this H
        Hi_E_levels = np.zeros(N_levels);
        for eig_i in range(N_levels):
            Hi_E_levels[eig_i] = H.GetEigenval(eigvecs[eig_i]);
            
        # determine labels, colors approproately
        labels, colors = EnergyLevels.GetBandLabels(H.MakeCopy(), Hi_E_levels, phi, False, material);
    
        #### now that we have levels and labels we can get transition energies as a function of Eg
        
        # container for transition energies
        E_up = np.zeros(N_param);
        E_down = np.zeros(N_param);
        
        # iter over the E levels we found before
        for level_i in range(N_param):
            
            # get levels
            levels = E_levels[:, level_i ];
            
            # get transition energies from levels
            E_up[level_i], E_down[level_i] = EnergyLevels.GetTransitionEnergies(levels, labels);   
            
        # plot transition energies as a function of the gap
        if( phi == 0 ):  # plot in black for the longitudinal valley
            plt.plot( Eg_arr, E_up, color = 'black', label = r"$|L_{6}^{+} \uparrow \rangle \rightarrow |L_{6}^{-} \uparrow \rangle, \phi ="+str(phi)[:4]+"$ [rad]", linestyle = '--' );
            plt.plot( Eg_arr, E_down, color = 'black', label = r"$|L_{6}^{+} \downarrow \rangle \rightarrow |L_{6}^{-} \downarrow \rangle, \phi ="+str(phi)[:4]+"$ [rad]");
            
            # also place resulting levels in containers for determining intersections
            central_lines.append(E_up);
            central_lines.append(E_down);
            
        else: # oblique valley should be red
            plt.plot( Eg_arr, E_up, color = 'red', label = r"$|L_{6}^{+}\uparrow \rangle \rightarrow |L_{6}^{-} \uparrow \rangle, \phi ="+str(phi)[:4]+"$ [rad]", linestyle = '--' );
            plt.plot( Eg_arr, E_down, color = 'red', label = r"$|L_{6}^{+}\downarrow \rangle \rightarrow |L_{6}^{-} \downarrow \rangle, \phi ="+str(phi)[:4]+"$ [rad]");
            
            # also place resulting levels in containers for determining intersections
            oblique_lines.append(E_up);
            oblique_lines.append(E_down);
       
            
    ####  make and show the actual plot now that we have done both valleys
    
    # get intersections of corresponding lines
    # recreate x arr for these lines
    Eg_arr = np.linspace(Eg_lim[0], Eg_lim[1], N_param);
    
    # want intersection of up, central with down, oblique
    x_int_up, y_int_up = EnergyLevels.GetIntersections( Eg_arr, central_lines[0], oblique_lines[1] );
    # and of down, central with up, oblique
    x_int_down, y_int_down = EnergyLevels.GetIntersections( Eg_arr, central_lines[1], oblique_lines[0] );
    
    # convert these points into formatted strings for the plot
    # only if intersections were found (ie returns != None)
    if( type(x_int_up) != type(None) and type(y_int_up) != type(None) ):
        label_up = "("+'{:.5f}'.format(x_int_up)+","+'{:.5f}'.format(y_int_up)+")"; # use .format()   
    else: # intersection was not found
        label_up = "("+str(x_int_up)+","+str(y_int_up)+")";
     
    # similarly for the other intersection   
    if(type(x_int_down) != type(None) and type(y_int_down) != type(None) ):
        label_down = "("+'{:.5f}'.format(x_int_down)+","+'{:.5f}'.format(y_int_down)+")"; # use .format()
    else: # intersection was not found
        label_down = "("+str(x_int_down)+","+str(y_int_down)+")";
    
    # plot the result as a scatter
    # need to figure out a better way of making these double strings
    plt.scatter(np.full(1, x_int_up), np.full(1, y_int_up), color = 'black', marker = '^', label = label_up);
    plt.scatter(np.full(1, x_int_down), np.full(1, y_int_down), color = 'black', marker = 'o', label = label_down );
    
    # format the plot
    plt.title(material+" $\gamma $ absorption spectrum");
    plt.xlabel("$E_{g}$ [eV]");
    plt.ylabel("$\Delta$E [eV]");
    
    # make the legend
    ax = plt.subplot();
    chartBox = ax.get_position();
    ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.6, chartBox.height]);
    ax.legend(loc='upper center', bbox_to_anchor=(1.45, 0.8), shadow=True, ncol=1);
    leg2 = ax.get_legend();
    leg2.set_draggable(True); # makes the legend box draggable !
    
    # textbox below the legend for exchange param vals
    textstring = "a1 = "+str(E_inputs[1])[:6]+",\n"+"a2 = "+str(E_inputs[2])[:6]+",\n"+"b1 = "+str(E_inputs[3])[:6]+",\n"+"b2 = "+str(E_inputs[4])[:6]+", ";
    text_x, text_y, = 1.4*Eg_lim[1], ax.get_ylim()[0] #+ ax2.get_ylim()[1])/3;
    ax.text(text_x, text_y, textstring, fontsize = 12);
            
    plt.show();

    #### end EgSweepTransition
    
    
def EgSweep_vs_x(H_inputs, E_inputs, Eg_func, x_lim, temp, material, Eg_bench, ax, no_phi ):
    '''
    Finds the dependence of the energy levels on the relative composition of the material, x
    Finds energy gap according to x, temperature, then basically runs the same code as the Eg sweep function
    
    :param H_inputs: tuple of the inputs for the bulk k.p hamiltonian (Eg, muBohr, Bvec, vt, vl), except
                Eg is just a dummy value since it is our dependent variable
    :param E_inputs: tuple of the inputs to construct the exchange matrix (phi, a1, a2, b1, b2)
    :param Eg_func: points to function that determine Eg at a given x, T
    :param x_lim: tuple of the start and stop values for the x sweep
    :param temp: double, the temperature in K
    :param material: string, the name of the material for the title of the plot
    :param Eg_bench: double, value of the band gap to benchmark ( label valence/conduction/spin of bands) at
    :param ax: matplotlib axes object, to plot these energy levels on
    :param no_phi: bool, tells whether to include value pf phi in legend labels (it is redundant if they are on separate axes)
    :param Eg_bench, double, the energy gap to benchmark (ie determine valence/conduction) at
    
    returns None, has already made all the plots
    '''
    
    # get, delete a version of the k.p bulk hamiltonian
    H = hamiltonian.Hamiltonian(*H_inputs);
    
    # define number of points to use in parameter space using H
    N_param = 100;
    N_levels = H.N;
    del H;
    
    # define the return variable
    # 2d array of energy levels for different Eg
    E_levels = np.zeros((N_levels, N_param));
    
    # construct the exchange matrix
    # this is the same throughout since we are only varying Eg
    E = hamiltonian.Exchange(*E_inputs);
    phi = E_inputs[0];
    
    # vary over composition x
    x_arr = np.linspace( x_lim[0], x_lim[1], N_param );
    
    # iter over the x values
    for i in range(len(x_arr)):
        
        # for this individual x, get energy gap
        x = x_arr[i];
        Eg = Eg_func(x, temp);
        
        # get the energy levels
        # first construct the bulk matrix
        H_inputs = (Eg, H_inputs[1], H_inputs[2], H_inputs[3], H_inputs[4] );
        H = hamiltonian.Hamiltonian(*H_inputs);
        
        # combine with the exchange matrix
        H.Add(E);
        
        # first run thru, get eigvecs
        if( i == 0):
            eigvecs = H.Eigvecs();
            
        # get eigenvalues with eigenvector method
        for eig_i in range(N_levels):
            
            # get the appropriate eigval and place in holder
            E_levels[eig_i, i] = H.GetEigenval(eigvecs[eig_i]);
            
        # delete H before the next run thru
        del H;
        
    #### end of loop over Eg vals
    
        
    # create H at  benchmark Eg (defaults to 150 meV) to determine band labels
    H_inputs = ( Eg_bench, H_inputs[1], H_inputs[2], H_inputs[3], H_inputs[4] );
    H = hamiltonian.Hamiltonian(*H_inputs);
    H.Add(E);
    
    # get levels from this H
    Hi_E_levels = np.zeros(N_levels);
    for eig_i in range(N_levels):
        Hi_E_levels[eig_i] = H.GetEigenval(eigvecs[eig_i]);
        
    # determine labels, colors approproately
    labels, colors = EnergyLevels.GetBandLabels(H.MakeCopy(), Hi_E_levels, phi, no_phi, material);
    
    # set title and linestyle according to phi
    if( phi == 0):
        style = '--'; # plots dashed
        ax.set_title("Longitudinal valley (--)", loc = 'left', fontsize = 10);
        ax.set_xlabel("x");
        ax.set_ylabel("Energy [eV]");
        
    else: #oblique valley
        style = 'solid'; # plots solid
        ax.set_title("Oblique valley ($ -$)", loc = 'right', fontsize = 10);
        ax.set_xlabel("x");
        #ax.set_ylabel("Energy [eV]");
        
    # iter over each level and plot
    for eig_i in range(N_levels):
        
        # plot level with appropriate label, color
        ax.plot(x_arr, E_levels[eig_i], color = colors[eig_i], linestyle = style, label = labels[eig_i]);
    
    return;

    #### end Eg sweep vs x
    
################################################################################
# done making plots against Eg, x
# moving to plots vs exchange parameters

def PlotLevelBothValleys( H_inputs, material, H_bench, a1, a2, b1, b2, limits, both ):
    '''
    Similar to EgSweepBothValleys, combines all the function calls of plotting energy levels vs a given param for both valleys into one. Plots/formats together
    and gives user option of plotting valleys on same or different axes. 
    
    :param H_inputs: tuple of the inputs required to construct the bulk hamiltonian.
                        (Eg, muBohr, Bvec, vt, vl)
    :param material: string, tells what material we are plotting energies of, for the title
    :param H_bench: hamiltonian object, hamiltonian at param space point which we want to use as benchmark point, should already have exchange added ! 
    :param a1: exchange parameter, either set to double of its fixed value, or to None if its the free param
    :param a2: see a1
    :param b1: see a1
    :param b2: see a1
    :param limits: tuple of the starting and stopping limits for the free exchange parameter
    :param both: bool, tells whether to plot two valleys on same or different axes
                        
    returns: None
    '''
    
         
    # determine whether plotting on seperate axes or together
    if( both): # plot together
        
        # only one axis needed
        fig, ax1 = plt.subplots();
        
        # alias the axis
        ax2 = ax1;
        
        # we will want to have phi in the legend
        no_phi = False;
        
        # also determine legend placement
        legend_coef = 0.6;
        
    else: # plot seperately
        
        # create two seperate axes for long and oblique valleys        
        fig, (ax1, ax2) = plt.subplots(ncols=2, figsize = (8,4), sharey = True );
        
        # we do not want phi's in the legend, it is redundant
        no_phi = True;
        
        # also determine legend placement
        legend_coef = 0.95;
    

        
    # plot levels against selected exchange param
    
    # for the longitudinal valley
    textstring = PlotLevel(H_inputs, material, H_bench, ax1, 0, a1, a2, b1, b2, limits, no_phi = no_phi);
    textstring = PlotLevel(H_inputs, material, H_bench, ax2, 70.5*np.pi/180, a1, a2, b1, b2, limits, no_phi = no_phi);
        
    # format the plot
    fig.suptitle(material + " bulk E levels \n", fontsize = 14);
    
    # add text box for exchange params
    #plt.text("a1 = "+str(E_inputs[1])[:6]+", "+"a2 = "+str(E_inputs[2])[:6]+",\n"+"b1 = "+str(E_inputs[3])[:6]+", "+"b2 = "+str(E_inputs[4])[:6]+", " , loc = "right", fontsize = 10  );
        
    # make the legend
    chartBox = ax2.get_position();
    ax2.set_position([chartBox.x0, chartBox.y0, chartBox.width*legend_coef, chartBox.height]);
    ax2.legend(loc='upper center', bbox_to_anchor=(1.45, 0.8), shadow=True, ncol=1);
    leg2 = ax2.get_legend();
    leg2.set_draggable(True); # makes the legend box draggable !
    
    # textbox below the legend for exchange param vals
    text_x, text_y, = 1.1*ax2.get_xlim()[1], ax2.get_ylim()[0] #+ ax2.get_ylim()[1])/3;
    ax2.text(text_x, text_y, textstring);
       
    # show and return 
    plt.show();                   
    return;

    #### end PlotLevelBothValleys
    
    
    
def PlotLevel(H_inputs, material, H_bench, ax, phi, a1, a2, b1, b2, limits,  return_results = False, fixedAB = False, no_phi = False ):
    '''
    Plots the energy levels as a function of one of the exchange params (a1, a2, b1, or b2). The other three
    params are kept at a fixed value. We select which param to vary over by setting it to None
    
    :param H_inputs: tuple of the inputs required to construct the bulk hamiltonian.
                        (Eg, muBohr, Bvec, vt, vl)
    :param material: string, tells what material we are plotting energies of, for the title
    :param H_bench: hamiltonian object, hamiltonian at param space point which we want to use as benchmark point, should already have exchange added ! 
    :param ax: matplotlib axis object for plotting data on, recieved from PlotLevelBothValleys (see above)
    :param phi: double, tells angle between H field and z hat in radians
    :param a1: exchange parameter, either set to double of its fixed value, or to None if its the free param
    :param a2: see a1
    :param b1: see a1
    :param b2: see a1
    :param limits: tuple of the starting and stopping limits for the free exchange parameter
    :param show: bool, whether to show the plot of the levels found here. Defaults to true. Use this option to hold off on plotting
                        central valley levels so that both can be shown on the same plot
    :param return_results: bool, tells to return level results before plotting or even formatting. Defaults to false. This option is used
                        by SearchProcessPlot to achieve plotting of both valleys together. User should never have to call this option
    :param fixedAB: bool, tells whether plotting in the case of A, B fixed (2D param space). Again, used by SearchProcessPlot but should 
                        not be changed by user. defaults to False ie regular 4d param space
    :param no_phi: bool, tells not to include phi value in legend band labels, defaults to False
                        
    returns: None in most cases, since plotting has already been done.
            If return_results == True, tuple of val_arr, E_levels, labels, colors
    
    '''
    
    # if A, B are fixed, get their vals which are disguised as a2, b2
    if( fixedAB):
        A, B = a2 , b2;
    
    # get, print, delete the k.p bulk hamiltonian
    H = hamiltonian.Hamiltonian(*H_inputs);
    print(H);
    
    # define number of points to use in parameter space
    N_param = 100
    N_levels = H.N
    del H;
    
    # get the array over the param space
    val_arr = np.linspace(limits[0], limits[1], N_param);
    
    # get container for the energy levels
    E_levels = np.zeros((N_levels, N_param) );
    
    # define E inputs here so that it persists
    E_inputs = ();
    
    # iter over this param space
    for i in range(len(val_arr) ):
        
        val = val_arr[i];
        
        # determine which one is the free param
        # it is the one set to None
        if( a1 == None):
            
            xlabel = "a1 [eV]";
            E_inputs = (phi, val, a2, b1, b2 );
            
            if( fixedAB): # have to modify to determine a2, b2
                
                # in this case a2 is actually A and b2 is B
                a2 = val - A;
                b2 = b1 - B;
                
                E_inputs = (phi, val, a2, b1, b2 );
            
        elif( a2 == None):
            
            xlabel = "a2 [eV]";
            E_inputs = (phi, a1, val, b1, b2 );
            
        elif( b1 == None):
            
            xlabel = "b1 [eV]";
            E_inputs = (phi, a1, a2, val, b2 );
            
        elif( b2 == None):
            
            xlabel = "b2 [eV]";
            E_inputs = (phi, a1, a2, b1, val );
            
        # every time thru we need to recreate H
        H = hamiltonian.Hamiltonian(*H_inputs);
            
        # construct the exchange matrix
        E = hamiltonian.Exchange(*E_inputs);
        
        # combine with H and get eigenvalues
        H.Add(E);
        
        # get the eigenvectors at i = 0 only
        # ordering will always follow how we ordered it here !
        if( i == 0):
            eigvecs = H.Eigvecs();
        
        #iter thru eigenvals
        for eig_i in range(len(eigvecs)):
            
            # get the eigval that corresponds to this eigvec
            eigvec = eigvecs[eig_i];
            eigval = H.GetEigenval(eigvec);
            
            # one energy level val for this param i
            E_levels[eig_i, i] = eigval;
            
        # delete the matrices
        del H;
        del E;
        
    #### end loop over i (free param)
    
    # create H at Eg = 150 meV and E at params = 0 to determine naming
    H_bench.Add(hamiltonian.Exchange(*E_inputs)); # this is bad code, just a stopgap measure
    # get levels from the banchmark hamiltonian
    Hi_E_levels = np.zeros(N_levels);
    for eig_i in range(N_levels):
        Hi_E_levels[eig_i] = H_bench.GetEigenval(eigvecs[eig_i]);
        
    # use these to determine labels, colors
    # this step should make labels independent of the ordering
    labels, colors = EnergyLevels.GetBandLabels(H_bench.MakeCopy(), Hi_E_levels, phi, no_phi, material);
    
    # show only for oblique phi
    if( return_results):
        return val_arr, E_levels, labels, colors;

    # plot
    for eig_i in range(N_levels):           

        if(phi == 0):
            ax.plot(val_arr, E_levels[eig_i], label = labels[eig_i], color = colors[eig_i], linestyle = '--');
        else:
            ax.plot(val_arr, E_levels[eig_i], label = labels[eig_i], color = colors[eig_i]);
        
    # format the plot
    ax.set_xlabel(xlabel);
    ax.set_ylabel("Energy [eV]");
    
    # make text box showing fixed parameters
    textstring = ''
    
    # go thru a, b params
    for param_i in range(4):
        
        param = [a1, a2, b1, b2][param_i];
        pstring = ["a1","a2","b1","b2"][param_i];
        
        # check that this is not the free param
        if (param != None):
            
            #place info in text box
            textstring += pstring+" = "+str(param)+",\n";
            
    return textstring;

    #### end plot level
    
    
def SearchProcessPlot(H_inputs, condition, axis, lim, material, H_bench, ABfixed = False, ABtuple = (0, 0) ):               
    '''
    Do the param space search, process the results, then plot the results with condition highlighted,
    all at once for a given phi. Then add comparison plot of same param space region for the other phi
    
    :param H_inputs: tuple of the 5 params needed for the bulk hamiltonian H
                    (Eg, muBohr, Bvec, vt, vl
    :param condition_arr: big 4D (2D) array of bool vals everywhere in param space
    :param axis: string, a1, a2, b1, b2, tells which param space axis to plot along
                    then we search along the others
    :param lim: tuple of where to start and stop the plotted param
    :param material: string telling which mateiral this is for, for title purposes
    :param H_bench, hamiltonian object, represents bulk H at inputs chosen to define band labels at
    :param ABfixed: bool, whether we are fixing A and B ( eliminating a2 and b2 and making the param space only 2D)
            default, set to False meaning param space is 4D
    :param ABtuple: tuple of the fixed vals for A and B to be set to, defaults to 0, 0
    '''
    
    # define phi values for each valley
    phi_central = 0;
    phi_oblique = 70.5*np.pi/180; # in radians
    
    # param space limits can all be defined the same
    a1_lim, a2_lim, b1_lim, b2_lim = lim, lim, lim, lim;
    
    if( not ABfixed): # sweeps over 4D param space
    
        # run the param space sweep for the central valley
        print("Running param space search.");
        search_results = EnergyLevels.ParamSS(H_inputs, phi_central, condition, a1_lim, a2_lim, b1_lim, b2_lim);
        
        # pass results to condition process function
        print("Processing param space search results to create plot.");
        process_results = EnergyLevels.ProcessCondition(axis, *search_results);
        
        # get vals of fixed params from process function
        param0 = process_results[2][0];
        param1 = process_results[3][0];
        param2 = process_results[4][0];
        
        # reassign based on axis
        if( axis == "a1"):
            a2, b1, b2 = param0, param1, param2;
            
            # now plot other valley energy levels at same point in param space
            # note that PlotLevel needs to take an axis arg after H_bench
            # for our purposes don't need to pass this but need a placeholder hence
            dummy_ax = 0;
            obl_x, obl_levels, obl_labels, obl_colors = PlotLevel(H_inputs, material, H_bench, dummy_ax, phi_oblique, None, a2, b1, b2, lim, return_results = True);
            for eig_i in range(len(obl_levels)):
                
                plt.plot(obl_x, obl_levels[eig_i], label = obl_labels[eig_i], color = obl_colors[eig_i]);
                
        # likewise for the other axes
        elif( axis == "a2"):
            a1, b1, b2 = param0, param1, param2;
            
            # now plot other valley energy levels at same point in param space
            # note that PlotLevel needs to take an axis arg after H_bench
            # for our purposes don't need to pass this but need a placeholder hence
            dummy_ax = 0;
            obl_x, obl_levels, obl_labels, obl_colors = PlotLevel(H_inputs, material, H_bench, dummy_ax, phi_oblique, a1, None, b1, b2, lim, return_results = True);
            for eig_i in range(len(obl_levels)):
                
                plt.plot(obl_x, obl_levels[eig_i], label = obl_labels[eig_i], color = obl_colors[eig_i]);
                
        # pass these results to the plot condition function
        #plt.show();
        print("Creating plot.");
        EnergyLevels.PlotCondition(*process_results, material);
        
    else: #sweep over only a1 and b1 with A, B fixed
        
        # run the param space sweep for the central valley
        # for A, B fixed the central valley energy levels are always static !!
                      
        # likewise for the oblique valley
        search_results = EnergyLevels.ParamSS_AB(H_inputs, phi_oblique, condition, a1_lim, b1_lim, *ABtuple );
        
        # pass results to process function
        process_results = EnergyLevels.ProcessCondition_AB(axis, *search_results);
        
        # reassign based on axis
        if( axis == "a1"):
            
            # get vals of fixed params from process function
            A, b1, B = process_results[3], process_results[4][0], process_results[5]
            
            # now plot central valley energy levels at same point in param space
            cen_x, cen_levels, cen_labels, cen_colors = PlotLevel(H_inputs, material, phi_central, None, A, b1, B, lim, H_bench, return_results = True, fixedAB=True);
            for eig_i in range(len(cen_levels)):
                
                plt.plot(cen_x, cen_levels[eig_i], label = cen_labels[eig_i], color = cen_colors[eig_i], linestyle = '--');
        
        # pass process results to the plot condition function
        leg1 = EnergyLevels.PlotCondition_AB(*process_results, material);   
        
        leg1.set_title("$E_{g} = $ "+str(H_inputs[0])[:5]+" [eV]" ); 
        plt.show();  

#### end Search Process Plot
                    
################################################################################
# test code here
if __name__ == "__main__":
    
    # define the H input params
    Eg = 0.05; # band gap in eV
    mu_Bohr = 5.79e-5; # bohr magneton in eV/T
    B_vec = np.zeros(3); # B field in T
    vt = 1; # transverse velocity matrix element
    vl =1; # longitudinal velocity matrix element
    H_inputs = (Eg, mu_Bohr, B_vec, vt, vl);
    
    # define the E input params
    phi_central = 0*np.pi/180; # angle between H field and z hat for the central valley IN RADIANS
    phi_oblique = 70.5*np.pi/180; # angle between H field and z hat for the oblique valleys IN RADIANS
    fixed1, free1, fixed2, free2 = "a1", "a2", "b1", "b2"; # which params to vary and fix
    free1_limits = (0.01, 0.1); # how to vary free params
    free2_limits = (0.01, 0.02);
    fixed1_val, fixed2_val = 0.01, 0.02; # reassign fixed vals as their numerical values
    
