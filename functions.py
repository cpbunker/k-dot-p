'''
Created on Jan 17, 2020

@author: Christian
'''

import numpy as np
import scipy.constants
import utils

################################################################################
# 
# test math functions that the finders can call and get the roots of
#
################################################################################


def f1 (x):
    # roots are 1 and -2
    return (x-1)*(x+2);

def f2(x):
    #inf number of roots
    return np.cos(3.14159*x);

def f3(x):
    #no roots, no crossings
    return np.exp(x);

def f4(x):
    #no roots, asymptotic crossings
    
    #catch division by zero
    if(x==0):
        x= pow(10,-9)
    return 1/x ;

def f5(x):
    #bunch of roots, complicated behavior
    
    return np.log(x )*(x-10)*(x-3)*(x-3)*(x-25);

################################################################################
# Actual functions from quantum well problem for rootfinding
################################################################################

class InputParams:
    
    def __init__(self, delta, V, d, v_z, v_x = None, v_y = None): 
        '''
        Object for passing experimental inputs as params to the various transcendental functions
        
        :param delta: double, band gap in meV
        :param V: double, well potential in meV
        :param d: double, well thickness in nm
        :param v_z: double, fermi velocity in z, units nm/s (hence large size)
                input into k*p matrix element hbar*v_z which is in in nm*meV
        :param v_x, v_y, doubles, velocities in x and y directions in anisotropic case, units nm/s
                optional, set to None in which case they will be defaulted to v_z (isotropic)
        '''
        
        # set main attributes
        self.delta = delta
        self.V = V
        self.d = d
        self.v_z = v_z
        
        #set secondary attributes
        self.h_bar = scipy.constants.hbar*(1/ (1.6*pow(10,-19)) )*pow(10,3); # hbar, put in units meV*s
        # set other velocities
        if( v_x == None and v_y == None): # isotropic case
            self.v_x = v_z;
            self.v_y = v_z;
            self.v_perp = v_z;
        elif( v_x != None and v_y != None and v_x == v_y): #anisotropic case, isotropic in x and y
            self.v_x = v_x;
            self.v_y = v_y;
            self.v_perp = v_x;
        elif( v_x != None and v_y != None): # completely anisotropic
            self.v_x = v_x;
            self.v_y = v_y;
            self.v_perp = None; # there is no isotropy in perpendicular direction    
        else: # should never get here  
            print("Error in InputParams: v_x, v_y not set correctly");      

    def DeclareSelf(self):
        '''
        Function for printing out all the info about this object
        for user to see to make things clear and for troubleshooting
        '''
        
        print(self);
        print("delta = "+str(self.delta)+" meV");
        print("V = "+str(self.V)+"meV");
        print("d = "+str(self.d)+" nm");
        print("v_z = "+str(self.v_z)+" nm/s"); 
                
        # control string representation using matplotlib mathtext format
        string_list = []
        string_list.append(r"$\Delta$ = "+str(self.delta)+" meV"+"\n");
        string_list.append(r"V = "+str(self.V)+"meV"+"\n");
        string_list.append(r"d = "+str(self.d)+" nm"+"\n");
        string_list.append(r"$v_{z}$ = "+str(self.v_z/pow(10,14) )+"e14 nm/s"+"\n");
        
        return string_list;
        
# end InputParams class #####################################################################


def pos_E_even( E, inputs  ):
    '''
    Transcendl function we have to solve to get energies of the quantum well problem
    valid in the regime E < 0 and E < -2*abs(delta), even solutions only, see pos_E_odd    
    Unit scales are meV for energies, nm for distances
    
    :param E: array like, eigenenergy, the input of the equation
    :param inputs: user defined InputParams object which contains the values to set the experimental parameters equal to
                    defines V in meV, delta in meV, d in nm, and v_z in nm/s
                    optional, initialized to default values of InputParams() for convenience
    
    returns array like, value of function found at the (various) E values
    '''
    
    #define input parameters, constants
    delta = inputs.delta; # band gap in meV
    V = inputs.V; # well potential in meV
    d = inputs.d; # well thickness in nm
    h_bar = scipy.constants.hbar*(1/ (1.6*pow(10,-19)) )*pow(10,3); # hbar, put in units meV*s
    hv = inputs.v_z*h_bar ; # matrix element in nm*meV
    
    # define intermediate terms
    m_A = abs(delta)/pow((hv/h_bar),2); # mass like term
    k_z = np.sqrt( ( ( 2*m_A*E/pow(h_bar,2)) )*(1+E/ (2*abs(delta)) ) ); # wave number
    rho = np.sqrt( ( m_A/(abs(delta)*pow(h_bar,2)) )*(E+V)*(-E-2*abs(delta) + V) );
    
    #put it all together
    return rho*(E+2*abs(delta)) - np.tan(k_z*d/2)*(E+2*abs(delta)-V)*k_z;
    
#  End pos_E_even #########################################################################


def pos_E_odd( E, inputs  ):

    '''
    Transcendental function that gives eigenenergy solutions in the regime E > 0, E < -2*delta
    Gives odd solutions only! See pos_E_even
    Unit scales are meV for energies, nm for distances

    :param E: array like, energy of electron, meV
    :param inputs: user defined InputParams object which contains the values to set the experimental parameters equal to
                    defines V in meV, delta in meV, d in nm, and v_z in nm/s
                    optional, initialized to default values of InputParams() for convenience
                        
    returns array like, evaluation of function at E
    '''    
    
    #define input parameters, constants
    delta = inputs.delta; # band gap in meV
    V = inputs.V; # well potential in meV
    d = inputs.d; # well thickness in nm
    h_bar = scipy.constants.hbar*(1/ (1.6*pow(10,-19)) )*pow(10,3); # hbar, put in units meV*s
    hv = inputs.v_z*h_bar ; # matrix element in nm*meV

    # define intermediate terms
    m_A = abs(delta)/pow((hv/h_bar),2); # mass like term
    k_z = np.sqrt( ( ( 2*m_A*E/pow(h_bar,2)) )*(1+E/ (2*abs(delta)) ) ); # wave number
    rho = np.sqrt( ( m_A/(abs(delta)*pow(h_bar,2)) )*(E+V)*(-E-2*abs(delta) + V) );
    
    #if( E == 0):
        #print ("E = 0 in pos_e_odd");
        #print(-rho*(E+2*abs(delta) ) - k_z*(E+2*abs(delta) - V) );
    
    # combine everything
    return -rho*(E+2*abs(delta) ) - k_z*(E+2*abs(delta) - V)*(1/np.tan( k_z*d/2) );

# End pos_E_odd ############################################################


def neg_E_even( E, inputs  ):
    '''

    Transcendental function that gives eigenenergy solutions in the regime
    -2*abs(delta) < E < 0. Gives even solutions only, see neg_E_odd for odd ones
    Unit scales are meV for energies, nm for distances    

    :param E:    
    :param inputs: user defined InputParams object which contains the values to set the experimental parameters equal to
                    defines V in meV, delta in meV, d in nm, and v_z in nm/s
                    optional, initialized to default values of InputParams() for convenience
                    
    returns array like, evaluation of function at E
    '''  

    #define input parameters, constants
    delta = inputs.delta; # band gap in meV
    V = inputs.V; # well potential in meV
    d = inputs.d; # well thickness in nm
    h_bar = scipy.constants.hbar*(1/ (1.6*pow(10,-19)) )*pow(10,3); # hbar, put in units meV*s
    hv = inputs.v_z*h_bar ; # matrix element in nm*meV

    
    # define intermediate terms
    m_A = abs(delta)/pow((hv/h_bar),2); # mass like term
    rho = np.sqrt( ( m_A/(abs(delta)*pow(h_bar,2)) )*(E+V)*(-E-2*abs(delta) + V) );
    kappa = (1/hv)*np.sqrt( -E*(2*abs(delta) + E)) # complex wave number

    return -rho*(E+2*abs(delta) ) - kappa*(E+2*abs(delta) - V)*(np.tanh( kappa*d/2) );

# end neg_E_even ###################################################################


def neg_E_odd( E, inputs ):
    '''

    Transcendental function that gives eigenenergy solutions in the regime
    -2*abs(delta) < E < 0. Gives odd solutions only, see neg_E_even for even ones
    Unit scales are meV for energies, nm for distances    

    :param E:    
    :param inputs: user defined InputParams object which contains the values to set the experimental parameters equal to
                    defines V in meV, delta in meV, d in nm, and v_z in nm/s
                    optional, initialized to default values of InputParams() for convenience
                    
    returns array like, evaluation of function at E
    '''  

    #define input parameters, constants
    delta = inputs.delta; # band gap in meV
    V = inputs.V; # well potential in meV
    d = inputs.d; # well thickness in nm
    h_bar = scipy.constants.hbar*(1/ (1.6*pow(10,-19)) )*pow(10,3); # hbar, put in units meV*s
    hv = inputs.v_z*h_bar ; # matrix element in nm*meV
    
    # define intermediate terms
    m_A = abs(delta)/pow((hv/h_bar),2); # mass like term
    rho = np.sqrt( ( m_A/(abs(delta)*pow(h_bar,2)) )*(E+V)*(-E-2*abs(delta) + V) );
    kappa = (1/hv)*np.sqrt( -E*(2*abs(delta) + E)) # complex wave number
    
    return -rho*(E+2*abs(delta) ) - kappa*(E+2*abs(delta) - V)*(1/np.tanh( kappa*d/2) );

# end neg_E_even ###################################################################


################################################################################
# Workspace for graphing functions, excluded from module part
################################################################################
    
if (__name__ == "__main__"):
    
    # select the function to plot
    my_func = pos_E_even;
    
    # control the plotting parameters
    inputs = InputParams(d=12);
    start = 0;
    stop = 125;
    x_label = " Energy [meV]"
    y_label = " f(E)"
    title = "Quantum well function, E>0, even solutions"
    
    # plot the function
    utils.PlotF( my_func, start, stop, x_label, y_label, y_label, title, (inputs,) );
