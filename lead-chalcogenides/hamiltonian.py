'''
Created on Mar 27, 2020

@author: Christian

Defines the two band, two spin k.p hamiltonian for bulk IV VI materials
Follows theory in Bauer, Zawadski 1992 paper

This module does not run any code itself. Just used to define the objects that energylevels.py
uses to compute and plot the energy levels of a system with given inputs. 

Everything here should work and the user is directed to only change things in the runner files
(PbSnSe.py etc). However if you get an error with hermicity, finding eigenvectors, eigenvalues 
etc it will direct you to come here
'''

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


################################################################################
# define generic matrix class
# no physics here, just use to define helpful math methods that we want the objects
# Exchange and Hamiltonian (which carry the physics) to inherit from

class Matrix:
    '''
    Generic square matrix class from which the different variants of hamiltonian matrices
    present in the problem can inherit useful methods for plotting, getting eigenvals etc
    Attributes assumed:
        self.mat: 2d np array, represents the matrix
        self.N: int, size of the matrix (assumed to be square)
        
    Hamiltonian class inherits from this
    Exchange class inherits from this
        
    '''
    
    ################################################################################
    # overloaded methods
    
    def __plus__(self, other): 
        '''
        Overloads the + operator, ie when you write H + E where H, E are matrix objects, calls this function
        :param other: matrix object, that we add to self
        
        returns, new 2d array represinting elementwise addition of the 2 matrices
        Does not!! change matrix objects themselves
        '''
        
        # check size compatibility
        if( np.shape(self) != np.shape(other) ):
            print("Matrix addition error: matrices have different shapes.\n")
            
        # go forward with actual addition
        return self.mat + other.mat;   


    ################################################################################
    # methods for doing calculations with the matrix
    
    def CheckHermicity(self):
        '''
        Makes sure we haven't messed up by making the matrix non-hermitian
        Just used for debugging
        '''
        
        # init bool
        is_H = True;
        
        #get the hermitian transpose
        HermT = np.conj( self.mat.T);
        
        # check equality by itering thru elements
        for i in range(self.N):
            for j in range(self.N):
                
                # set tolerance
                tolerance = pow(10, -6);
                
                # check real parts
                if( abs(np.real(self.mat[i,j]) - np.real(HermT[i,j]) ) > tolerance):
                    # matrix elements are not equal!
                    is_H = False;
                    
                # same for imag parts
                if( abs(np.imag(self.mat[i,j]) - np.imag(HermT[i,j]) ) > tolerance):
                    # matrix elements are not equal!
                    is_H = False;
                
                
        if( is_H): # it is hermitian
            return True;
        
        else: # not hermitian
            
            print("The "+str(type(self))+" matrix object is not hermitian. ");
            print("There is most likely a problem in the definition of the .mat attributes of either the Hamiltonian or Exchange objects. ");
            print("Check the __init__ functions of Hamiltonian and Exchange for any issues. ");
            raise ValueError("Matrix is not hermitian. ");
            return False;
        
    #### end Check Hermicity
    
    
    def Diagonalize(self):
        '''
        Diagonalizes the square matrix (2D array) stored in the self.mat attribute
        Diagonalization is used to determine which eigenvalues are considered spin up and spin down depending on the basis
        This function assumes the matrix to be hermitian
        
        Args: none
        
        Returns a new 2D array representing the diagonalized matrix. Does not change self
        '''
        
        # find the transformation matrix S
        # this is just a matrix with eigvecs as its columns
        S = np.linalg.eigh( self.mat)[1]; #2nd return of this function is matrix of eigvecs
        
        # find inverse of S
        S_inv = np.linalg.inv(S);
        
        # transform self.mat to obtain diagonalized matrix D
        D = np.matmul( np.matmul( S_inv, self.mat), S);
        
        return D;
        
    
    def Eigvals(self, i= -1):
        '''
        Get the eigenvalues of self.mat, assuming it is hermitian. Raises ValueError if not hermitian.
        Non hermicity usually means a problem in the definitions of self.mat in either the
        Hamiltonian or Exchange object __init__ functions.
        
        :param i, int, tells which eigenvalue to grab
                defaults to -1 in which case it returns all
        '''
        
        if(self.CheckHermicity()): # make sure matrix is hermitian before proceeding
        
            # get eigenvals using np
            eigvals, eigvecs = np.linalg.eigh(self.mat);
            #eigvals.sort()
            
            # return all or one eigenvalues as asked
            if( i == -1): # want to return them all
                return eigvals;
            
            else:
                return eigvals[i];
            
    #### end Eigvals
        
    
    def Eigvecs(self, i = -1):
        '''
        Get the normalized eigenvectors of self.mat, assuming it is hermitian. Raises ValueError if not hermitian.
        Non hermicity usually means a problem in the definitions of self.mat in either the
        Hamiltonian or Exchange object __init__ functions. 
        
        :param i, int, tells which eigenvalue to grab
                defaults to -1 in which case it returns all
        '''
        
        #check that the matrix is still hermitian
        if( self.CheckHermicity()):
            
            #get eigvecs using np
            eigvals, eigvecs = np.linalg.eigh(self.mat);
            
            # get eigenvectors in same order as eigenvalues
            good_eigvecs = []; # return variable
            for j in range( len(eigvals) ):
                
                # eigenvector is the column
                good_eigvecs.append( eigvecs[:,j] );
                
            # return all or one eigenvalues as asked
            if( i == -1): # want to return them all
                return np.array(good_eigvecs);
            
            else:
                return good_eigvecs[i];
            
        #### end Eigvecs
        
    def GetEigenval(self, eigvec):
        '''
        This is one of the most important functions in this module, but is also prone to a lot of problems. basically, using np.linalg.eigh() to
        get the eigenvalues of an array is not always the best, because it simply returns the eigvals in ascending order. So if you want to see
        interesting physics ie level crossing as a parameter of the matrix is tuned, you will get bad plots that show levels "bouncing" off each
        other as the eigenvalues switch places in the array that np.linalg.eigh returns.
        
        The solution is to organize the eigenvalues not in ascending order but according to the order of their corresponding eigenvectors. The
        eigenvectors are given some fixed (ascending order of corresponding eigenvalues) order at some point in parameter space but this is not
        changed as the matrix input parameters are varied.
        
        The only problem is that the eigenvectors do not (in general) keep their exact values as we change the matrix either. The method I use
        here is to classify eigenvectors with the same signs elementwise as equal (assuming the test eigenvector wasn't actually equal to anything).
        This method hasn't failed me but it is not on 100% solid footing.
        
        :param eigvec: 1d array, representing the static eigenvector from the benchmark point in param space where we determined the initial ordering
                of the levels. Will be compared to the actual eigenvectors of this matrix.
                
        Returns double of the corresponding eigenval is successful. Raises ValueError if not.
        '''
        
        # get the index of the eigvec
        vecs = self.Eigvecs();
        
        # iter over the true eigenvectors to find if there is one with true equality to the test one
        index = None
        for i in range(len(vecs)):
            
            # check equality
            if( np.allclose(eigvec, vecs[i]) ):
                index =i;
                
        # try again if no true equality was found
        # classify them as the same if they have the same signs elementwise
        if( index == None):
            for i in range(len(vecs)):
                
                # check equality by matching signs only
                same_signs = True
                for j in range(len(eigvec)): # iter thru vector elements
                    if( np.sign( eigvec[j]) != np.sign(vecs[i,j] ) ):
                        same_signs = False;
                        
                # and if they do have the same signs, pick this index
                if( same_signs):
                    index = i;  
                    
        # try again
        # classify them as the same if they have always opposite signs elementwise
        if( index == None):
            for i in range(len(vecs)):
                
                # check equality by matching signs only
                same_signs = True
                for j in range(len(eigvec)): # iter thru vector elements
                    if( np.sign( -eigvec[j]) != np.sign(vecs[i,j] ) ):
                        same_signs = False;
                        
                # and if they do have the same signs, pick this index
                if( same_signs):
                    index = i;                           
                
        # check that that worked
        # if it didn't print some helpful info
        if( index == None): # the above process failed
            print("Error in GetEigenval(): the test eigenvector "+str(eigvec) );
            print("Did not match any of the true eigenvectors of the matrix.");
            print("Eigenvectors are: \n", vecs);
            print("You may have to resort to using the method Eigval() instead of GetEigenval() although this will result in no level crossing in the plots.");
            raise ValueError("No corresponding eigenvalue");
                
        # if it worked, get the appropriate eigenvalue
        else:
            return self.Eigvals(index);
    
    
    ################################################################################
    # methods for plotting/showing/etc
    
    def Show(self, precision = 3, imag = False):
        '''
        Method shows matrix as a matplotlib plot to avoid dealing with dumb print out
        restrictions of most IDE's
        
        :param precision, int, tells how many decimal places to include
        :param imag, bool, tells whether to print imag component of elemnts as well
        '''
        
        # determine whether to print imag parts of the elements or not
        if(imag): #this means user wants us to print no matter what so we will
            pass;
        
        else: #user has selected/defaulted no imag parts, but we double check
            
            for i in range(self.N): # over rows
            
                for j in range(self.N): # over columns
                    
                    if (np.imag( self.mat[i,j]) != 0):
                        
                        # the matrix elements have nonzero imaginary components
                        imag = True;
            
        
        # create grid in xy space to print matrix elements to
        x = np.linspace(0, 1, self.N);
        y = np.linspace(1, 0, self.N);
        
        # now iter over elements of the matrix
        for i in range(self.N): # over rows
            
            for j in range(self.N): # over columns
                
                # get the matrix element
                elem = self.mat[i,j] ;
                
                # get specified precision and make into a string
                elem_real = str( int( np.real(elem)*pow(10, precision) )/pow(10, precision) ); # real part
                elem_imag = str( int( np.imag(elem)*pow(10, precision) )/pow(10, precision) ); # imag part
                
                # determine whether to print imag part
                if imag: # matrix has imag part and we should print it
                    elem_string = elem_real +"+ "+elem_imag+"i";
                else:
                    elem_string = elem_real;
                
                # plot as a string
                plt.text(x[j], y[i], elem_string, fontsize = 9);
                
                
        
        # show with axis removed
        plt.axis('off');
        plt.show();
        
        return None;
    
    #### end Show
    
    
    

################################################################################
# define a class that represents the bulk k.p hamiltonian for this problem
# taken straight from Bauer Zawadski SST 1992 Eq 13

class Hamiltonian(Matrix):
    '''
    This object represents the k.p hamiltonian for this problem
    See Bauer Zawadski SST 1992 section 2.1 for the underlying physics
    
    Inherits methods from Matrix class which allow for showing, plotting, eigenvalue
    calculations etc
    
    Here we handle how set the inputs of the problem/ how to create the matrix
    '''
    
    ################################################################################
    # overloaded methods
    
    def __init__(self, Eg, mu_Bohr, B_vec, vt, vl):
        '''
        Create the object and its matrix representation (self.mat) from th input physical parameters
        
        :param Eg: double, band gap in eV
        :param mu_Bohr: double, bohr magneton in eV/T
        :param B_vec: 1d np array of [Bx, By, Bz] in T, typically only Bz is nonzero
        :param vt: velocity matrix element in transverse (xy) plane in m/s
        :param vl: velocity matrix element in longitudinal (z) direction in m/s
        '''
    
        
        # define input param attributes from init args
        self.Eg = Eg; # in eV
        self.mu_Bohr = mu_Bohr; # in eV/T
        self.B_vec = B_vec; # in T
        self.vt = vt; # in m/s
        self.vl = vl; # in m/s
        
        # break up the B vector for simplicity
        Bx, By, Bz = B_vec[0], B_vec[1], B_vec[2];
        
        # set m and g terms to free electron values for now
        m0 = 0.51e6 /pow(3e8, 2); # free electron mass in ev/c^2
        mt_plus = m0; # far band effective mass contributions
        mt_minus = m0;
        ml_plus = m0;
        ml_minus = m0;
        
        g0 = 2; # free electron spin gvalue, dimensionless
        gt_plus = g0; # far band spin g-value contributions
        gt_minus = g0;
        gl_plus = g0;
        gl_minus = g0;
        
        # define the P vector
        # this is coupling of momentum p to vector potential A
        Px = 0;
        Py = 0;
        Pz = 0;
        # and the +- combinations of it
        P_plus = np.complex(Px, Py);
        P_minus = np.complex(Px, -Py);
        
        # from all this we can define the band edge operators
        V_plus = -Eg/2 - (pow(Px, 2) + pow(Py, 2))/(2*mt_plus) - pow(Pz,2)/(2*ml_plus) + (1/2)*gl_plus*mu_Bohr*Bz;
        V_minus = -Eg/2 - (pow(Px, 2) + pow(Py, 2))/(2*mt_plus) - pow(Pz,2)/(2*ml_plus) - (1/2)*gl_plus*mu_Bohr*Bz;
        C_plus = Eg/2 + (pow(Px, 2) + pow(Py, 2))/(2*mt_minus) + pow(Pz,2)/(2*ml_plus) + (1/2)*gl_minus*mu_Bohr*Bz;
        C_minus = Eg/2 + (pow(Px, 2) + pow(Py, 2))/(2*mt_minus) + pow(Pz,2)/(2*ml_minus) - (1/2)*gl_minus*mu_Bohr*Bz;
        
        # empty holder for the exchange matrix when we add it
        self.Exchange = None;
        
        # now we should be able to define the actual matrix
        self.N = 4; # size of the matrix
        self.mat = np.array([
            [V_plus, (1/2)*gt_plus*mu_Bohr*Bx, vl*Pz, np.sqrt(2)*vt*P_minus ],
            [(1/2)*gt_plus*mu_Bohr*Bx, V_minus, np.sqrt(2)*vt*P_plus, -vl*Pz],
            [vl*Pz, np.sqrt(2)*vt*P_minus, C_plus, (1/2)*gt_minus*mu_Bohr*Bx],
            [np.sqrt(2)*vt*P_plus, -vl*Pz, (1/2)*gt_minus*mu_Bohr*Bx, C_minus],
            ]);
        
        
        
    #### end init
    
        
    def __str__(self):
        '''
        Overload the string representation of the object
        Basically just have it print all its input params
        '''
        
        # making it all into 1 big string due to the fact that the 
        # __str__ method is required to return a string
        printstring = '';
        
        # header
        printstring += 'Inputs of the bulk IV VI material k.p hamiltonian:\n';
        
        # go through the experimential inputs
        printstring +=  "Band gap Eg = "+str(self.Eg)+" eV\n";
        printstring += "Bohr magneton mu_Bohr = "+str(self.mu_Bohr)+" eV/T\n";
        printstring += "B field [Bx, By, Bz] = "+str(self.B_vec)+" T\n";
        printstring += "Transverse velocity vt = "+str(self.vt)+" m/s\n";
        printstring += "Longitudinal velocity vl = "+str(self.vl)+ "m/s\n";
        
        return printstring;
    
    #### end str
    
    
    def MakeCopy(self):
        '''
        Use this function when you want to need an H object to pass to a function but don't want to manipulate the actual instance of
        the object that you're working with
        '''
        
        #def H inputs and create copy
        H_inputs = (self.Eg, self.mu_Bohr, self.B_vec, self.vt, self.vl);
        copyH = Hamiltonian(*H_inputs);
        
        # the orig instance matrix might have been modified post-init 
        # ie by the addition of an Exchange matrix, so copy this separately
        copyH.mat = self.mat;
        
        # and have to pass on old exchange separately
        copyH.Exchange = self.Exchange;
        
        return copyH;
    
    
    ####################################################################
    # matrix manipulation methods
    
    def Add(self, exchange):
        '''
        Adds the elements of the exchange matrix to the elements of the hamiltonian
        Uses overloaded + but is different as it also modifies/replaces self.mat
        
        :param exchange: the exchange matrix object we are combining the hamiltonian with
        
        returns None, modifies self
        '''
        
        #check type of exchange matrix
        dummy_exchange = Exchange(0,0,0,0,0);
        if( type(exchange) == type(dummy_exchange) ): #we can proceed
            
            # this function adds the result to the .mat attribute of the hamiltonian
            # rather than just returning the result as the overloaded + does
            self.mat = self.mat + exchange.mat;
            
            # also store the exchange matrix in the .Exchange attribute so we can access its ind'l attributes
            self.Exchange = exchange;
            
            return;
            
        else:
            print("Error in Hamiltonian.Add(): second matrix must be an exchange matrix\n");
            raise ValueError("Incorrect usage of Add() method in hamiltonian.py");
 
 #### end Hamiltonian class   
    
    
################################################################################
# define a class that represents the exchange matrix for this problem
# again the theory is expounded in Bauer Zawadski SST 1992 with the
# exchange matrix defined in Eq 35

class Exchange(Matrix):
    '''
    Represents the exchange matrix for the IV VI valley problem. The physics 
    of this is expounded in Bauer Zawadski SST 1992 section 2.2
    This term encapsulates some anisotropy effects in that the central and
    oblique valleys have different orientations and therefore different phi
    values. 
    
    Remember that we are dealing with the internal auxillary field H rather than
    an external magnetic field B
    
    This is also where we bring in the spin matrix elements a1, a2, b1, b2
    It is these parameters that we want to explore
    '''
    
    ################################################################################
    # overloaded methods
    
    def __init__(self, phi, a1, a2, b1, b2, long_valley_shift = 0):
        '''
        Create the exchange matrix and its input parameters
        
        :param phi: double, angle between H and z hat (111) direction
        :param a1: spin matrix element, depends on <R|J|R>, we just set different values and check results
        :param a2: spin matrix element, depends on <S+-|J|S+->, "
        :param b1: spin matrix element, depends on <Z|J|Z>, "
        :param b2: spin matrix element, depends on <X+-|J|X+->, "
        '''
        
        # define basic attributes from args
        self.phi = phi
        self.a1 = a1
        self.a2 = a2
        self.b1 = b1
        self.b2 = b2
        self.long_valley_shift = long_valley_shift;
        
        # get linear combos of matrix elements
        A = a1 - a2;
        B = b1 - b2;
        
        # define the actual matrix
        self.N = 4; # size of the matrix
        self.mat = np.array( [
            [A*np.cos(phi), a1*np.sin(phi), 0, 0],
            [a1*np.sin(phi), -A*np.cos(phi), 0, 0 ],
            [0, 0, B*np.cos(phi), -b1*np.sin(phi)],
            [0, 0, -b1*np.sin(phi), -B*np.cos(phi)]
            ])
        
        # this block is totally garbage test code just to try something
        # Mandal Springholz et al paper reports Bi doping opens up gap at longitudinal valley but not oblique
        # for us this is equivalent to adding a linear shift along the diagonal when phi = 0
        if( phi == 0 and long_valley_shift != 0): # only bother doing this if the shift is nonzerp
            
            # construct shift matrix with +/- shift along diagonals
            # valence should be shifted down, conduction up
            shift_mat = np.array( [
                [-long_valley_shift/2, 0, 0, 0],
                [0, -long_valley_shift/2, 0, 0],
                [0, 0, long_valley_shift/2, 0],
                [0, 0, 0, long_valley_shift/2]
                ])
            
            self.mat += shift_mat;
        
    def __str__(self):
        '''
        Controls how the object is printed etc
        This will just print out all the input params for debugging purposes
        '''
        
        # has to return a string
        ret_string = '';
        
        # go thru and add info about all the inputs
        ret_string += "phi = "+str(self.phi)+"\n";
        ret_string += "a1 = "+str(self.a1)+"\n";
        ret_string += "a2 = "+str(self.a1)+"\n";
        ret_string += "b1 = "+str(self.a1)+"\n";
        ret_string += "b2 = "+str(self.a1)+"\n";
        
        return ret_string;
        
    
    
#### end Exchange class
    

    
    

if __name__ == "__main__":
    
    ############################################################################
    # test code here
    
    def increment(x):
        print("************** in increment")
        # method for increasing whatever quantity we are finding slope with respect to
        delta = abs(x)/100; # get 1 percent change
        print(x,"-->", x+delta)
        return x + delta; # always moves in positive direction
    
    # define the input params
    Eg = 10.0; # band gap in eV
    mu_Bohr = 5.79e-5; # bohr magneton in eV/T
    B_vec = np.zeros(3); # B field in T
    vt = 1; # transverse velocity matrix element
    vl =1; # longitudinal velocity matrix element
    
    # create first hamiltonian object
    H = Hamiltonian(Eg, mu_Bohr, B_vec, vt, vl);
    #H.Show();
    
    # add exchange
    phi = 0;
    a1, a2 = 2, 1; # so A = 1
    b1, b2 = 3, 1; # so B = 2
    E = Exchange(phi, a1, a2, b1, b2);
    H.Add(E);
    print("Old H w/ exchange, should have eigenvalues -4, -6, 7, 3");
    print("Code eigenvalues: ", H.Eigvals() );
    H.Show();
    
    # make copy without exchange
    H_new_inputs = H.Eg, H.mu_Bohr, H.B_vec, H.vt, H.vl ; # H inputs are exactly same
    H_new = Hamiltonian(*H_new_inputs); # make new hamiltonian object
    #H_new.Show();
    
    # get in old exchange effects and modify
    E_old = H.Exchange;
    E_new_inputs = E_old.phi, E_old.a1, increment(E_old.a2), E_old.b1, E_old.b2 ; # same inputs as before but a1 slightly increased
    E_new = Exchange(*E_new_inputs);
    H_new.Add(E_new);
    print("New H w/ modified exchange, should have eigenvalues ");
    print("Code eigenvalues: ", H_new.Eigvals() );
    H_new.Show();
    
    

    
    
   
 