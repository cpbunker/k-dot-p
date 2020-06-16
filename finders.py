'''
Created on Jan 17, 2020

@author: Christian
'''

################################################################################
# 
# Functions to find roots of math expressions using various methods
#
################################################################################

import scipy.misc

def Bisect( f, start, stop, tolerance, is_verbose, arguments=() ): 
    '''
    function to rootfind over an interval by bisecting it
    :param f: function, math expression that we find roots (zeros) of
    :param start: double, beginning of interval
    :param stop: double, end of interval
    :param tolerance: double, how close to the root to be
    :param is_verbose: bool, tells to print out everything for debugging
    :param arguments: tuple, can be passed to f if f requires add'l arguments
    
    returns double root
    '''
    
    if (is_verbose):
        print("Entering Bisect function.\n");
    
    #### check endpoint conditions #### 
    
    # otherwise the endpoints of the interval must have opposite sign
    # usually this will be guaranteed by IntervalControl
    # however when doing an AsymptoteCheck it usually will fail here
    '''
    if (f(start, *arguments)*f(stop, *arguments) > 0): #ie have the same sign
        print("Error in bisect function: \n");
        print("No sign change over interval \n");
        return None;
        '''
    
    # iniial guess and step marker
    guess = (stop+start)/2; #avg of endpoints
    n_steps=0;
    #if(is_verbose): # print out initial guess
        #print( "guess " + str(n_steps)+ ": root = "+str(guess)+"\n");
    
    
    # keep guessing to find the root
    #iterate as long as f(guess) is not close to zero
    # but dont iterate forever
    while (abs( f(guess, *arguments) ) > tolerance):
        
        # update step counter
        n_steps += 1;
        
        # hold on to the last guess before reassignment
        old_guess = guess ;

        if(f(start, *arguments)*f(guess, *arguments) < 0): # guess is above root, lower it by bisection
            #reassignments
            stop=guess; #adjust interval
            guess = (stop+start)/2;
            
            #if(is_verbose and n_steps < 20):
                #print( "guess " + str(n_steps)+ ": root = "+str(guess)+"\n");
            
        elif(f(guess, *arguments)*f(stop, *arguments) < 0): # guess is below root, raise it by bisection
            #reassignments
            start=guess;
            guess = (stop+start)/2;
            
            #if(is_verbose and n_steps < 20):
                #print( "guess " + str(n_steps)+ ": root = "+str(guess)+"\n");
                
        # prevent infinite loop (eg asymptotic behavior)
        # check that ther have been way too many steps
        # and that the guess isn't moving anymore
        if( n_steps > 5000 and (guess - old_guess < tolerance)):
            if( is_verbose):
                print("Bisect sees asymptotic behavior, root rejected.");
            return None;
        
    #### end while

    return guess;

################################################################ end bisect


def Newton( f, start, stop, tolerance, is_verbose= False ): 

    '''
    Uses Newton's method to find the root
    Pros:
    Cons:

    :param f: function, what we want to find root of
    :param start: double, where interval begins
    :param stop: double, where interval ends
    :param tolerance: double, how close to get to the real root
    :param is_verbose: bool, tells whether to print out for debugging
    '''

    if(is_verbose):
        print("Entering Newton function");

    # initial guess
    guess = (stop - start)/2;    

    # loop while we get closer to root
    n_step = 0; #counter for steps

    # continue until error is small enough
    while( f(guess) > tolerance):

        if(is_verbose):
            print("Guess "+n_step+": x = "+guess+"\n");      

        #update guess
        guess -= f(guess)/scipy.misc.derivative(f,guess) 
        #step is bigger if it is far away from f=0 (f(guess) large)
        #step is smaller if slope is large (f'(guess) large)                                      

        #update step counter
        n_step += 1;                                       

        # End while       

    return guess;

# End Newton function ####################################################################


def Secant( f, start, stop, tolerance, is_verbose = False ):
    '''
    Use secant method to find the root     

    :param f: function, what we are finding roots for
    :param start: double, interval begins
    :param stop: double, interval ends
    :param tolerance: double, how close to f=0 to get
    :param is_verbose: bool, tells whether to print extra statements for debugging
    '''

    if(is_verbose):
        print("Entering Secant function");

    # sign change necessary on interval
    if (f(start)*f(stop) > 0): 
        print("Unable to find root by secant: no sign change in interval.\n");

        return None; # end the search   

    # reassign to temp variables
    new_start = start;
    new_stop = stop;   

    # guess and step counter
    guess = new_start - f(new_start)*(new_stop - new_start)/( f(new_stop) - f(new_start) );
    n_step = 0;

    # iter until error is small enough
    while( f(guess) > tolerance ):
        
        # verbose output
        if( is_verbose):
            print("Guess "+n_step+": x = "+guess+"\n");

        # narrow interval according to conditions
        if( f(guess)*f(new_start) < 0):
            new_start += 0; # newstart unchanged
            new_stop = guess; # interval steps backward           

        elif( f(guess)*f(new_stop) < 0):
            new_start = guess; # interval steps forward
            new_stop += 0; #newstop unchanged            

        elif( f(guess) == 0):
            return guess; #found exact root       

        # update guess and step counter
        guess = new_start - f(new_start)*(new_stop - new_start)/( f(new_stop) - f(new_start) );   
        n_step += 1;

        # End while loop      

    return guess;

# End Secant function ###############################################################

 

