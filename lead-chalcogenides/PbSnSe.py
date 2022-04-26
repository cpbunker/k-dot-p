'''
Created on Mar 30, 2020

@author: Christian

Runner file to find and plot energy levels of PbSnSe using the 
hamiltonian.py, EnergyLevels.py, GetEnergyLevels.py modules.
'''

import numpy as np
import scipy as sp

from Valleytronics import hamiltonian
from Valleytronics import GetEnergyLevels
from Valleytronics import conditions


################################################################################
# define the characteristics of PbSnSe that define the bulk hamiltonian
# this is a place where the user will want to make changes
# ie inputting the values of composition x, temperature T, B field, and matrix elements vl and vt

MATERIAL = "PbSnSe"; # global var to pass to plotting functions for titleing plots

# determine the energy gap according to x and T
def Eg_func( x, T):
    '''
    Band gap in eV is function of T, temperature in K, and x, Sn composition
    This function is obtained from  Kryzman et al SPIE 2019
    '''
    # return band gap in eV
    return (72 - 710*x + 2440*pow(x,2) - 4750*pow(x,3) - 22.5 + np.sqrt(506 + 0.1*pow(T - 4.2, 2) ) )/1000;

# input the temperature and Sn composition
# and thus get the gap
x = 0.1; # unitless
T = 10; # in K
Eg = Eg_func(x, T);

            
# other input params of the bulk hamiltonian
B_vec = np.zeros(3); # B field in T
vt = 4*pow(10,5); # transverse velocity matrix element
vl = 4*pow(10,5); # longitudinal velocity matrix element

# define bohr magneton, fundamental constant, only change if you're changing unit systems
mu_Bohr = 5.79e-5; # bohr magneton in eV/T

# package up the bulk H inputs as a tuple
H_inputs = (Eg, mu_Bohr, B_vec, vt, vl);
H = hamiltonian.Hamiltonian(Eg, mu_Bohr, B_vec, vt, vl);


################################################################################
# define the inouts to, and general limits of the inputs to, the exchange matrix
# changing these parameters will change value of fixed params in plots or domains of plots

# define the phi values for differnt points in the Brillouin zone
# phi is the angle between applied field and chosen z hat direction 
# don't change these, they are built in to the geometry of the problem
phi_central = 0*np.pi/180; # gamma bar point aka central valley, in radians
phi_oblique = 70.5*np.pi/180; # m bar point aka oblique valley, in radians

# define baseline values for each exchange parameter
# we can and will overwrite these by redefining each variable later in the code
# but if you don't care about a variable enough to redefine it, this is what it will fall back to
a1 = 0.01; # in eV
a2 = 0.01; # in eV
b1 = 0.01; # in eV
b2 = 0.01; # in eV
A = 0.01; # in eV
B = 0.01; # in eV

# likewise define limits for each of these quantities for if we plot over them
# again, these are baseline vals, we will often redefine
a1_lim = 0.01, 0.1;
a2_lim = 0.01, 0.1;
b1_lim = 0.01, 0.1;
b2_lim = 0.01, 0.1;



################################################################################
# define the benchmark hamiltonian
# this amounts to determing at what input values we want to call lower bands valence
# bands and higher bands conduction bands
# definitely something the user may want to change

# define the benchmark energy level
Eg_bench = 0.15; # 150 meV

# define resst of bulk hamiltonian
H_bench_inputs = (Eg_bench, mu_Bohr, B_vec, vt, vl); #uses baseline vals, but user is free to redefine
H_bench = hamiltonian.Hamiltonian(*H_bench_inputs)


################################################################################
# make plots of the energy levels vs band gap
# if you want to make the plots in this section, set the plot_flag_1 variable below to true
plot_flag_1 = False;

if( plot_flag_1):
    
    # redefine exchange param vals if you want used in this code block by uncommenting next line
    # a1, a2, b1, b2 = 0.01, 0.015, 0.035, 0.06;
    
    # from this construct the inputs to the exchange matrix   
    E_inputs = (0, a1, a2, b1, b2);# 0 is a placeholder
        
    #### plot the levels against the energy gap
    
    # determine limits of this plot
    Eg_lim = - 0.1, 0.1; # in eV
    
    # tell code we are not plotting against x with the fololowing tuple
    x_inputs = (False, False, False, False );

    # make energy level plots for both valleys
    # these will appear in the same figure, we can force them onto the same axes by setting both = True
    GetEnergyLevels.EgSweepBothValleys(H_inputs, E_inputs, Eg_lim, MATERIAL, Eg_bench, x_inputs, both = True);
    
    #### plot the levels against the Sn composition
    
    # define input params for plotting against x instead of Eg
    vs_x = True; #tells to plot against x 
    x_lim = (0,1); # sets limits of x
    temperature = T # temperature in K inherits from what we used above, or we can redefine
    Eg_of_x = Eg_func; # function to convert from x to gap
    x_inputs = (vs_x, x_lim, temperature, Eg_of_x );
    
    # redefine E inputs back to the central valley and plot for this valley
    # can use the same E_inputs with placeholder phi
    # as above, we can force both valleys' levels onto same axis using both = True, or separate them using both = False
    GetEnergyLevels.EgSweepBothValleys(H_inputs, E_inputs, Eg_lim, MATERIAL, Eg_bench, x_inputs, both = True );
    
    
################################################################################
# make plots of the energy levels vs a single exchange param with the others fixed
# if you want to make the plots in this section, set the plot_flag variable below to true
plot_flag_2 = False;

if( plot_flag_2):
    
    # if you want to plot aganst a different param, change a1 below to your chosen param  
    b2 = None; # plots against the param set to None
    lim = 0.01, 0.1; # set plotting limits
    
    # redefine the other fixed exchange params if you want
    # don't redefine the one you set to None !
    # a2, b1, b2 = 0.05, 0.05, 0.05 ;
    
    # now make plots for both valleys against the same exchange parameter
    # can plot together or separetely using the show argument as above
    GetEnergyLevels.PlotLevelBothValleys(H_inputs, MATERIAL, H_bench, a1, a2, b1, b2, lim, both = False)


################################################################################
# plot the results of a 4D param space search
# to make the plots is this section, change the plot flag below to true
plot_flag_3 = True;

if( plot_flag_3):
    
    # define the param to plot against and its limits
    # unlike in PlotLevel, the other parameters are not fixed
    axis = "a2"; # can change to a2, b1, b2 to plot along those
    lim = 0.01, 0.1; # start and stop of domain for plot
    
    # define the condition to check for in the param space search
    # conditions are all defined and documented in conditions.py
    condition = conditions.E0_highest;
    
    # do the search and plot the results
    GetEnergyLevels.SearchProcessPlot(H_inputs, condition, axis, lim, MATERIAL, H_bench);
    
    

################################################################################
# do a param space search as above but with A, B fixed
# the result is that the param space is only 2D
# to make the plots is this section, change the plot flag below to true
plot_flag_4 = False;

if( plot_flag_4):
    
    # define the param to plot against and its limits
    # unlike in PlotLevel, the other parameters are not fixed
    axis = "a1"; # can change to b1
    lim = 0.01, 0.1; # start and stop of domain for plot
    
    # define the condition to check for in the param space search
    # conditions are all defined and documented in conditions.py
    condition = conditions.E0_highest;
    
    # redefine fixed A, B if you want
    # A, B = 0.01, 0.01
    
    # do the search and plot the results
    try:
        GetEnergyLevels.SearchProcessPlot(H_inputs, condition, axis, lim, MATERIAL, ABfixed = True, ABtuple = (A, B) )
    except:
        pass;
    

##########################################################################################
# again do a param space search with A, B fixed
# want to get an idea of good values of A, B
# start at low A, B and gradually increase
plot_flag_5 = False;

# define limits for A, B, and number of mesh points
A_lim, B_lim = (0.01, 0.1), ( 0.01, 0.1 );
N_A, N_B = 5, 5;

if( plot_flag_5):
    
    for A in np.linspace(*A_lim, N_A):
        
        for B in np.linspace(*B_lim, N_B):
        
            # tell user what A, B vals we are at
            print("A = "+str(A)+", B = "+str( B));
    
            # plot results of param space search
            # sometimes we will not find anywhere where the condition is met to this function will have an error
            # want to proceed over the A, B loop anyway
            # side step this by using try, finally structure
            try: # so that we continue if this code gives an error
                GetEnergyLevels.SearchProcessPlot(H_inputs, conditions.Switch , axis, lim, MATERIAL, True, (A, B ) );
                
            except: # do nothing, go to the next A, B vals
                pass;
        
        





