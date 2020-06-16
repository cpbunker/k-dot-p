'''
Created on Feb 3, 2020

@author: Christian

Functions to create and manipulate matrix representations
'''

import numpy as np
import time
import scipy.constants
import scipy.integrate
from Finders import main
from Finders import functions
from Finders import utils

################################################################################
# 
# Functions for the basis |Psi_n, psi_n> = <F1n,F2n,F1n,-F1n>
#
################################################################################


def ShowMatrix(old_matrix, precision = 2):    
    '''
    Function for printing a matrix nicely
    
    :param old_matrix: square np array, matrix we want to print
    :param precision: int, optional, how many decimal places to show, default 2
    '''

    # get shape of old matrix
    n = np.shape( old_matrix)[0]
    new_matrix = np.full( (n,n), np.complex(0,0) );   

    # iter thru and replace to correct precision
    for i in range(n): # iter over rows

        for j in range(n): #iter over columns
            
            # get old val and split to real and imag
            old_val = old_matrix[i,j];
            old_real = np.real(old_val);
            old_imag = np.imag(old_val);
            
            # round each part accordingly
            new_real = int( old_real*pow(10,precision) ) / pow(10, precision);            
            new_imag = int( old_imag*pow(10,precision) ) / pow(10, precision);   
            
            # combine result
            new_val = np.complex(new_real, new_imag);         

            new_matrix[i,j] = new_val 
    
    #float_formatter = "{:.4f}".format ;
    #np.set_printoptions(formatter={'float_kind':float_formatter});         

    print(new_matrix);
    
#  End ShowMatrix ####################################################################
    
       
def DeltaH(k_x,k_y):
    '''
    creates and returns the k_x, k_y dependent hamiltonian DeltaH
    for finding the dispersion E(k_x, k_y) perturbatively
    
    :param k_x, k_y: doubles, wavenumbers (momentas) of particle in x, y directions respectively, in nm^(-1)
    '''
    
    #### define constants here
    h_bar = 1; # scipy.constants.hbar*(1/ (1.6*pow(10,-19)) )*pow(10,3); # hbar, put in units meV*s; # reduced Planck's const
    m_0 = 1; # mass term, meV/ c^2, c in nm/s
    P = 1; # Kane k.p matrix element, untiless number
    
    # make hamiltonian
    H = np.full((4,4), (0+0j));
    
    # fill hamiltonian
    H[0][3] = (h_bar*P/m_0)*(k_x - (0+1j)*k_y);
    H[1][2] = (h_bar*P/m_0)*(k_x + (0+1j)*k_y);
    H[2][1] = (h_bar*P/m_0)*(k_x - (0+1j)*k_y);
    H[3][0] = (h_bar*P/m_0)*(k_x + (0+1j)*k_y);
    
    return H;

#  End DeltaH #######################################################
  
    
def Parity_E( E, eigenenergies, delta):
    '''
    Determine parity of a eigenenergy based on the positioning of the bound states around delta
    Rule is that first state above delta is odd, first below is even
    
    :param E: energy to check the parity of
    :param eigenenergies: all the enrgy solutions for this configuration, 
            we must find our energy's relative position in this list to determine its parity
    :param delta: experimental param 1/2 band gap in well, starting point for parity determination
    '''

    # see if E is above or below delta 

    # regime 1: E > delta
    if(E > delta): # in this case first state above delta is odd                

        # get position of this E in the regime
        n_between = 0;
        for energy in eigenenergies:
            
            #check if it is between E and delta
            if( energy < E and energy > delta ):
                n_between += 1;

        # get parity based on n_between
        # if 0, 2, 4... (even #) states are between, E is odd
        if( pow(-1, n_between) > 0): # (-1)^n positive, so n_between is even, so E is odd!
            is_odd = True;   
            
        else: # (-1)^n negative, so n_between is odd, so E is even!
            is_odd = False;        

    # regime 2: E < delta
    elif(E <= delta): # in this case first state below delta is even
            
        # get position of this E in the regime
        n_between = 0;
        for energy in eigenenergies:           

            #check if it is between E and delta
            if( energy > E and energy <= delta ):
                n_between += 1;

        # get parity based on n_between
        # if 0, 2, 4... (even #) states are between, E is even
        if( pow(-1, n_between) > 0): # (-1)^n positive, so n_between is even, so E even
            is_odd = False;
            
        else: # (-1)^n negative, so n_between is odd, so E is odd
            is_odd = True;

    # print a bnch of stuff for troubleshooting
    if False:
        print("E = "+str(E) );
        print("n_between = "+str(n_between) );
        print("(-1)^n = "+str( pow(-1, n_between) ) )
        
    return is_odd;

# end ParityE ######################################################


def NormalizeF(F1, F2, k_z, kappa, exp_inputs, is_verbose = True):

    # normalization condition is from integrating over all z    
    def integrand_well(z): # this is only for -d/2 < z < d/2
        
        # check domain
        if( abs(z) <= exp_inputs.d/2 ): # we are inside the well       
            return np.conj( F1(k_z*z) )*F1(k_z*z) + np.conj( F2(k_z*z) )*F2(k_z*z); 
        
    def integrand_out(z): # integrate over wavefunction outside the well
        
        # check we are in the right domain
        if( z > exp_inputs.d/2 ): # using symmetry we only look at positive z
            return np.exp( -2*kappa*(z - exp_inputs.d/2) );
        
    # get values of the two necessary integrals
    I_well = utils.CheckIntegral( scipy.integrate.quad( integrand_well, -exp_inputs.d/2, exp_inputs.d/2), "NormalizeF: I_well");
    I_out = utils.CheckIntegral( scipy.integrate.quad( integrand_out, exp_inputs.d/2, pow(10,6)), "NormalizeF: I_out");
    
    # get value of A
    result = 2*( pow(F1(k_z*exp_inputs.d/2),2) + pow(F2(k_z*exp_inputs.d/2), 2))*I_out;
    result += I_well
    A = np.sqrt( 1/result);
    
    return A;

#  End Normalize F ####################################################################


def InnerPsi (k_x, k_y, m, n, eigenenergies, exp_inputs, is_verbose = False):
    '''
    Find the matrix element of matrix by the inner product <Psi_m, | matrix | phi_n >
    :param k_x, double, x momentum of particle
    :param k_y, double, y momentum of particle
    :param m: int, index of brah vector
    :param n: int, index of ket vector
    :param eigenenergies, 1D np array, has kx, ky = 0 energy values which are needed for k_z
    :param exp_inputs, functions.InputParams object, tells experimental values necessary for computing k_z
    :param is_verbose: bool, tells whether to print out extra info for troubleshooting
        
    returns complex number corresponding to the inner product
    '''  
    
    PLUSMINUS = 1
    
    if(is_verbose):
        print("Entering InnerPsi");   

    # get relevant energies
    E_m, E_n = eigenenergies[m-1], eigenenergies[n-1] # again recall m, n start at 1
    print( eigenenergies);
    #return

    # determine parity of each energy
    # rule is first state above -abs(delta) is odd, first below is even
    is_m_odd = Parity_E(E_m, eigenenergies, exp_inputs.delta);
    is_n_odd = Parity_E(E_n, eigenenergies, exp_inputs.delta);  

    # define wavevector, kappa and envelope functions for Psi_m, check if imag or not
    if( E_m < 0 and E_m > 2*exp_inputs.delta): # it is imag

        k_z_m = (1/exp_inputs.h_bar*exp_inputs.v_z)*np.sqrt(-E_m*(2*abs(exp_inputs.delta) + E_m) ); #define wavevector
        kappa_m = (1/exp_inputs.h_bar*exp_inputs.v_z)*np.sqrt( (E_m+exp_inputs.V)*(-E_m- 2*abs(exp_inputs.delta) + exp_inputs.V) );
        
        #check parity to assign envelope functions
        if( is_m_odd): #then F1_m is odd
            F1_m = np.sinh; #just get function pointer
            F2_m = np.cosh;  
            
        else: #F1_m is even
            F1_m = np.cosh;
            F2_m = np.sinh;      

    else: # it is not imag

        k_z_m = (1/exp_inputs.h_bar*exp_inputs.v_z)*np.sqrt(E_m*(2*abs(exp_inputs.delta) + E_m) );#define wavevector
        kappa_m = (1/exp_inputs.h_bar*exp_inputs.v_z)*np.sqrt( (E_m+exp_inputs.V)*(-E_m- 2*abs(exp_inputs.delta) + exp_inputs.V) );

        #check parity to assign envelope functions
        if( is_m_odd): #then F1_m is odd
            F1_m = np.sin; #just get function pointer
            F2_m = np.cos;
            
        else: #then F1_m is even
            F1_m = np.cos;
            F2_m = np.sin;

    # similarly for phi_n

    if(E_n < 0 and E_n > 2*exp_inputs.delta): # it is imag

        k_z_n = (1/exp_inputs.h_bar*exp_inputs.v_z)*np.sqrt(-E_n*(2*abs(exp_inputs.delta) + E_n) ); #define wavevector
        kappa_n = (1/exp_inputs.h_bar*exp_inputs.v_z)*np.sqrt( (E_n+exp_inputs.V)*(-E_n- 2*abs(exp_inputs.delta) + exp_inputs.V) );
        
        #check parity to assign envelope functions
        if( is_n_odd): #then F1_n is odd
            F1_n = np.sinh; #just get function pointer
            F2_n = np.cosh;
            
        else: # then F1_n is even
            F1_n = np.cosh;
            F2_n = np.sinh;

    else: # it is not imag

        k_z_n = (1/exp_inputs.h_bar*exp_inputs.v_z)*np.sqrt(E_n*( 2*abs(exp_inputs.delta) + E_n) );#define wavevector
        kappa_n = (1/exp_inputs.h_bar*exp_inputs.v_z)*np.sqrt( (E_n+exp_inputs.V)*(-E_n- 2*abs(exp_inputs.delta) + exp_inputs.V) );

        #check parity to assign envelope functions
        if( is_n_odd): #then F1_n is odd
            F1_n = np.sin; #just get function pointer
            F2_n = np.cos;
            
        else: # then F1_n is even
            F1_n = np.cos;
            F2_n = np.sin;
            
    # get normalization constants
    A_m = NormalizeF(F1_m, F2_m, k_z_m, kappa_m, exp_inputs, is_verbose);
    A_n = NormalizeF(F1_n, F2_n, k_z_n, kappa_n, exp_inputs, is_verbose);

    #### we are finally ready to integrate
    #### break up the integral into real and imaginary parts
    
    # real part
    def integrand_real(z): #remember this has to be piecewise in z
     
        if (abs(z) <= exp_inputs.d/2): # we are in the well

            # get envelope functions
            return A_m*F2_m(k_z_m*z)*A_n*F1_n(k_z_n*z) - PLUSMINUS* A_m*F1_m(k_z_m*z)*A_n*F2_n(k_z_n*z);

        elif( abs(z) > exp_inputs.d/2): #outside well
            
            # get sign of z
            sign_z = np.sign(z);

            #redefine F1m, F1n, F2m, F2n as decaying exponentials
            def f1_m(z):
                return np.exp(-kappa_m*( abs(z) - exp_inputs.d/2) );
            def f1_n(z):
                return np.exp(-kappa_n*( abs(z) - exp_inputs.d/2) );
            def f2_m(z):
                return np.exp(-kappa_m*( abs(z) - exp_inputs.d/2) );
            def f2_n(z):
                return np.exp(-kappa_n*( abs(z) - exp_inputs.d/2) );
            
            # define normalization constants B_m B_n according to parity
            # for B_m
            B_m = A_m*F1_m( sign_z*k_z_m*exp_inputs.d/2);

            # for B_n
            B_n = A_n*F1_n( sign_z*k_z_n*exp_inputs.d/2);
                
            # now ready to return
            #B_n, B_m = 1, 1
            return B_m*f2_m(z)*B_n*f1_n(z) - PLUSMINUS* B_m*f1_m(z)*B_n*f2_n(z);
                
    
    # imag part               
    def integrand_imag(z): #remember this has to be piecewise in z

        if ( abs(z) <= exp_inputs.d/2): # we are in the well 

            # get envelope functions
            return A_m*F2_m(k_z_m*z)*A_n*F1_n(k_z_n*z) + PLUSMINUS* A_m*F1_m(k_z_m*z)*A_n*F2_n(k_z_n*z);

        elif( abs(z) > exp_inputs.d/2): #outside well
            
            # get sign of z
            sign_z = np.sign(z);

            #redefine F1m, F1n, F2m, F2n as decaying exponentials
            def f1_m(z):
                return np.exp(-kappa_m*( abs(z) - exp_inputs.d/2) );
            def f1_n(z):
                return np.exp(-kappa_n*( abs(z) - exp_inputs.d/2) );
            def f2_m(z):
                return np.exp(-kappa_m*( abs(z) - exp_inputs.d/2) );
            def f2_n(z):
                return np.exp(-kappa_n*( abs(z) - exp_inputs.d/2) );
            
            # define normalization constants B_m B_n according to parity
            # for B_m
            B_m = A_m*F1_m( sign_z*k_z_m*exp_inputs.d/2);

            # for B_n
            B_n = A_n*F1_n( sign_z*k_z_n*exp_inputs.d/2);
            
            if( z< -6.5 and z > -6.55): 
                print(20*"*")
                print(A_n, B_n, A_n*B_n,z)
                print(A_m, B_m, A_m*B_m,z)
                print(20*"*");
                
            # now ready to return
            #B_n, B_m = 1, 1
            return B_m*f2_m(z)*B_n*f1_n(z) + PLUSMINUS* B_m*f1_m(z)*B_n*f2_n(z);               
    
    #### get value of real and imag integrals
    INF = 10000# exp_inputs.d/2; # for integration
    I_real = utils.CheckIntegral( scipy.integrate.quad(integrand_real, -INF, INF), "InnerPsi (real integral)" ) 
    I_imag = utils.CheckIntegral( scipy.integrate.quad(integrand_imag, -INF, 0), "InnerPsi (imag integral)" ) 

    # combine real and imag parts into complex inner product, include k terms and prefactor
    prefactor = exp_inputs.h_bar*exp_inputs.v_z;
    matrix_element = prefactor*np.complex(k_x*I_real, k_y*I_imag) #*prefactor;  
    
    # check that element is nonzero only for opposite parity
    if( is_m_odd == is_n_odd ): # same parity so should be zero
        if( abs((matrix_element)) > 0.001 ): # ie not nonzero
            print("Error in InnerPsi: |"+str(m)+"> odd "+str(is_m_odd)+" and |"+str(n)+"> odd "+str(is_n_odd)+
                  " but matrix element = "+str(matrix_element) );
            print("Psi_m: ",F1_m, F2_m);
            print("phi_n: ",F1_n, F2_n);
            #utils.PlotF( integrand_imag,  -80, 80, my_title = str(m)+", "+str(n))
    
    def tempf(z):
        return F2_m(k_z_n*z);      
            
    if (m==1 and n==2):
        print("m == 1 and n == 2",20*"*");
        print("F1_n = "+str(F1_n));
        print("F2_n = "+str(F2_n));
        print("F1_m = "+str(F1_m));
        print("F2_m = "+str(F2_m));
        print(k_z_n)
        utils.PlotF( tempf,  -10, 10, my_title = str(m)+", "+str(n))

    return matrix_element;

#  End InnerPsi #######################################################
    
    
def H_ij( i, j, k_x, k_y, eigenenergies, exp_inputs, check_time = False, is_verbose = False):
    '''
    Finds the m, nth matrix element of H_eff, which is a function of k_x, k_y
    This matrix element comes from integrating over the inner product of the perturbation matrix DeltaH
    in the basis | Psi_n, psi_n > (envelope functions)
    
    This function is used to build an entire n x n, k dependent matrix
    
    This function will be called many, many times when finding dispersion!!! 
    as a result it is important to time optimize it!!!
    
    :param i: int, index of brah vector
    :param j: int, index of ket vector
    :param k_x, k_y: doubles, wavenumbers (momentas) of particle in x, y directions respectively, in nm^(-1)
    :param eigenergies, np array of the k_x = k_y = 0 energies found by GetRoots
    :param exp_inputs, functions.InputParams object, tells experimental values necessary for computing k_z
    :param check_time: bool, whether to time how long this process takes
    :param is_verbose: bool, tells whether to print out extra info for troubleshooting
    
    returns double or complex number representing matrix element, or None if an error occured
    '''
    
    #begin clock if asked
    if( check_time):
        t_start = time.time();
        
    # get size of matrix
    size = 2*len(eigenenergies); # how big the overall H matrix is ( size x size )
    
    # get the delta H matrix
    H = DeltaH(k_x, k_y)
    
    #### check what to do based on i and j ####
    
    #### first figure out which quadrant / sub matrix
    
    #quadrants clockwise from top left are 1,2,3,4
    quadrant = 0; #init at 0

    # run thru four possibilities, i is row, j is column
    if( i <= size/2 and j<= size/2): # quad 1
        quadrant = 1;
    elif( i <= size/2 and j > size/2): # quad 2
        quadrant = 2;
    elif( i > size/2 and j > size/2 ): # quad 3
        quadrant = 3;
    elif( i > size/2 and j <= size/2 ): # quad 4
        quadrant = 4;
        
    #### now return val based on quadrant
    if (quadrant == 1 ): #diagonal quadrant, either eigenvalue or zero
        
        if( i == j): # need the eigenvalue            
            return_val = eigenenergies[i-1];  # matrix indices start at 1 but eigenenrgies don't      
        else:
            return_val = 0;
        
    elif (quadrant == 3 ): #diagonal quadrant, either eigenvalue or zero
        
        if( i == j): # need the eigenvalue       
            return_val = eigenenergies[i - 1 - int(size/2) ]; # matrix indices start at 1 but eigenenrgies don't    
        else:
            return_val = 0;
        
    elif( quadrant == 2): # off diagonal quadrant, need inner product
        return_val = InnerPsi(k_x, k_y, i, int(j - size/2), eigenenergies, exp_inputs, is_verbose); # j needs to be brought down, what should z be?
    
    elif( quadrant == 4): # off diagonal quadrant, need inner product
        return_val = np.conj (InnerPsi(k_x, k_y, j, int(i - size/2), eigenenergies, exp_inputs, is_verbose) ); # this quad is conjugate of quad 2
    
    else: # messed up selecting quadrant
        print( "Error in H_ij, i = "+str(i)+", j = "+str(j)+", incorrect quadrant. \n")
        return_val = None;
        
    if(check_time):
        t_stop = time.time();
        print ("H_ij run time is "+str(t_stop - t_start)+" seconds. \n");
        
    return return_val;
        
#  End H_ij ################################################################


def H_tilde(k_x, k_y, func_inputs, finder_inputs, is_verbose = False):
    '''
    Constructs a 2D np array, of dimension ( no of eigvals), which represents the bulk hamiltonian for nonzero kx, ky
    Skips perturbation step and assumes antidiagonals dominate
    
    :param k_x: double, value of momenta in x
    :param k_y: double, value of momenta in y
    :param func_inputs: functions.InputParams() object, contains experimental inputs
    :param finder_inputs: main.FinderParams() object, tells how to find roots
    :param is_verbose: bool, tells whether to print out extra info for troubleshooting
    '''
    
    if(is_verbose):
        print("Entering H_tilde");
        func_inputs.DeclareSelf();
    
    # get eigenenergies of hamiltonian from AllRoots, in ascending order
    eigenenergies = np.sort( main.AllRoots(func_inputs, finder_inputs) );
    
    # get size of effective Hamiltonian
    size = len(eigenenergies);
    
    if(is_verbose):
        print("Eigenergies are: "+str(eigenenergies) );
        print("Matrix size is "+str(size) );
    
    # make empty hamiltonian, fill with complex numbers
    H = np.full( (size, size), (0+0j) );  
    
    # determine anisotropy in x, y plane
    is_isotropic = False;
    if( func_inputs.v_x == func_inputs.v_y): # we do have isotropy
        matrix_element = func_inputs.h_bar * func_inputs.v_perp;
        is_isotropic = True;
    
    #get each matrix element
    for i in range(1,size+1): # iter over rows, start at 1
        for j in range(1,size+1): #iter over columns, start at 1
            
            # reassign with correct matrix element
            # remember we fill diagonal with eigenenergies and antidiagonal with kx, ky terms
            if( i == j): # diagonal
                # remember to correct in array indices for starting at 1
                H[i-1,j-1] = eigenenergies[i-1];
            
            elif( i == size - j + 1): # antidiagonal
                
                # row 1,3,5 etc (odd have kx - ky)
                if( pow(-1,i) < 0): # odd row, kx-ky
                    # rememebr to correct for indices starting at 1
                    H[i-1,j-1] = matrix_element*np.complex(k_x, -k_y);
                else: # even row, kx+ky
                    # rememebr to correct for indices starting at 1
                    H[i-1,j-1] = matrix_element*np.complex(k_x, k_y); 
                    
            else: # not diagonal or antidiagonal
                H[i-1,j-1] = 0;        
            
            if(is_verbose):
                print("i = "+str(i)+", j = "+str(j)+", H matrix element is " + str(H[i-1,j-1]) );
            
    return H;

#  end H_tilde ###################################################################### 


def H_landau( B, func_inputs, finder_inputs, is_verbose = False):
    '''
    Constructs a 2D np array, of dimension ( no of eigvals), which represents the bulk hamiltonian in a B field
    Skips perturbation step and assumes antidiagonals dominate
    
    :param k_x: double, value of momenta in x
    :param k_y: double, value of momenta in y
    :param func_inputs: functions.InputParams() object, contains experimental inputs
    :param finder_inputs: main.FinderParams() object, tells how to find roots
    :param is_verbose: bool, tells whether to print out extra info for troubleshooting
    '''
    
    if(is_verbose):
        print("Entering H_landau");
        func_inputs.DeclareSelf();
    
    # get eigenenergies of hamiltonian from AllRoots, in ascending order
    eigenenergies = np.sort( main.AllRoots(func_inputs, finder_inputs) );
    
    # get size of effective Hamiltonian
    size = len(eigenenergies);
    
    if(is_verbose):
        print("Eigenergies are: "+str(eigenenergies) );
        print("Matrix size is "+str(size) );
    
    # make empty hamiltonian, fill with doubles
    H = np.full( (size, size), 0.0 );   
    
    #get each matrix element
    for i in range(1,size+1): # iter over rows, start at 1
        for j in range(1,size+1): #iter over columns, start at 1
            
            # reassign with correct matrix element
            # remember we fill diagonal with eigenenergies and antidiagonal w/ B dependent terms
            if( i == j): # diagonal
                # remember to correct in array indices for starting at 1
                H[i-1,j-1] = eigenenergies[i-1];
            
            elif( i == size - j + 1): # antidiagonal
                
                H[i-1,j-1] = np.sqrt( 2*exp_inputs.h_bar*B) 
                    
            else: # not diagonal or antidiagonal
                H[i-1,j-1] = 0;        
            
            if(is_verbose):
                print("i = "+str(i)+", j = "+str(j)+", H matrix element is " + str(H[i-1,j-1]) );
            
    return H;

#  end H_tilde ###################################################################### 
    

def H_eff(k_x, k_y, func_inputs, finder_inputs, is_verbose = False):
    '''
    Constructs a 2D np array, of dimension 2*( no of eigvals), which represents effective hamiltonian for nonzero kx, ky
    Created by taking k-dependent DeltaH matrix as a perturbation
    
    :param k_x: double, value of momenta in x
    :param k_y: double, value of momenta in y
    :param func_inputs: functions.InputParams() object, contains experimental inputs
    :param finder_inputs: main.FinderParams() object, tells how to find roots
    :param is_verbose: bool, tells whether to print out extra info for troubleshooting
    '''
    
    if(is_verbose):
        print("Entering H_eff. \n)");
        func_inputs.DeclareSelf();
    
    # get eigenenergies of hamiltonian from AllRoots, in ascending order
    eigenenergies = np.sort( main.AllRoots(func_inputs, finder_inputs) );
    
    # get size of effective Hamiltonian
    size = 2*len(eigenenergies);
    
    if(is_verbose):
        print("Eigenergies are: "+str(eigenenergies) );
        print("Matrix size is "+str(size) );
    
    # make empty hamiltonian, fill with complex numbers
    H = np.full( (size, size), (0+0j) );   
    
    # call H_ij to get each matrix element
    for i in range(1,size+1): # iter over rows, start at 1
        for j in range(1,size+1): #iter over columns, start at 1
            
            # reassign with correct matrix element
            # remember to correct in array indices for starting at 1
            H[i-1,j-1] = H_ij(i, j, k_x, k_y, eigenenergies, func_inputs, is_verbose = is_verbose) ;
            
            if(is_verbose):
                print("i = "+str(i)+", j = "+str(j)+", H matrix element is " + str(H[i-1,j-1]) );
            
    return H;

#  end H_eff ###################################################################### 
    
    
if __name__ == "__main__":
    ############################################################################
    #  Test code here
    ############################################################################
    
    #### define experimental parameters
    delta = -25; #band gap in meV
    V = 175; # potential in meV
    d = 12; # well thickness in nm
    v_z = 4.5*pow(10,14) # matrix element in nm/s
    
    # pass to InputParams object
    exp_inputs = functions.InputParams(delta, V, d, v_z);
    
    #### define finder parameters
    my_method = 'bisect';
    my_interval = ( -V, V - 2*abs(delta) );
    
    #pass to FinderParams pbject
    find_inputs = main.FinderParams(method = my_method, interval = my_interval)
    
    # define k vals
    k_perp = 0.4
    k_x = np.sqrt(k_perp/2);
    k_y = k_x;
    
    # get effective hamiltonian as matrix in envelope wavefunction (Psi) basis
    H_small = H_tilde(k_x, k_y, exp_inputs, find_inputs, True);
    H_large = H_eff(k_x, k_y, exp_inputs, find_inputs, True)
    
    with np.printoptions( precision = 3, suppress = True):
        print( H_small);
        evals_small, evecs_small = np.linalg.eigh( H_small);
        print(evals_small);
        print(20*"*");
        print( H_large);
        evals_large, evecs_large = np.linalg.eigh( H_large);
        print(evals_large);
    
    
    
