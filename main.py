'''
Created on Jan 17, 2020

@author: Christian
'''

import scipy.misc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

import utils
import finders
import functions
import matrices



################################################################################
# 
# Main program of the rootfinding project
# Controls user input and makes calls to rootfinding functions
#
################################################################################



class FinderParams: 
    '''
    Object that contains the inputs that control the rootfinding process   
    Sets default inputs but allows for user changesm
    '''
    
    def __init__(self, method='bisect', f = functions.f1 , is_verbose=False , interval=(-100,100), 
                 length_scale=10, n_roots=10, ):
        '''
        
        :param method: string, tells which finders.py function to call
        :param f, function, tells which functions.py math function to find roots of
        :param is_verbose: bool, tells to print everything out for debugging
        :param interval: tuple, start and stop points of rootfinding interval
        :param length_scale: double, sense of how far apart crossings/roots tend to be
        :param n_roots: int, how many roots to look for

        '''
        
        if( is_verbose):
            print("Init FinderParams object \n");

        self.method = method;
        self.f = f;
        self.is_verbose = is_verbose; 
        self.interval = interval;
        self.length_scale = length_scale;
        self.n_roots = n_roots;

        
    def Run(self, arguments):
        '''
        Pass input attributes as parameters to GetRoots function
        Basically all the params of the finder are encapsulated in the self attributes
        All the inputs to the transcendental function (experimental values like delta, etc) are in arguments
        
        :param arguments: a InputParams() object (see functions.py)
                            defines delta, V, d, v_z
        
        returns list, results
        '''
        
        if( self.is_verbose):
            print( "Entering FinderParams.Run()");
            
        #call GetRoots function
        # this returns a list
        my_results = GetRoots(self.f, self.method, self.interval[0], self.interval[1], self.length_scale, 
                              self.n_roots, self.is_verbose, arguments);
        
        return my_results;
        
    def DeclareSelf(self): 
        '''
        Displays all the current inputs to the user
        '''
        
        # concat the first part of the declare
        my_string = "The program is set to use the "+str( self.method)+ " method with " ;
            
        if(self.is_verbose):
            my_string += "verbose output ";
        else:
            my_string += "simple output ";
            
        my_string += " to find the roots of "+ str(self.f) + " defined in functions.py. \n" ;
        
        #actually print
        print( my_string);    
        print( "The input settings are:\n");
        print( "Interval: ["+str(self.interval[0])+","+str(self.interval[1])+"]\n");
        print( "length scale: "+str(self.length_scale) +"\n");
        print( "Number of roots: "+ str( self.n_roots) + "\n");
        
        
    def UserInput(self):
        '''
        Allows the user to change the input attributes
        '''
        
        if (self.is_verbose):
            print( "Entering FinderParams.UserInput \n");
            
        #  change method, make sure an available method is asked for
        system_ok = False;
        while(not system_ok): # prompts changes to method until accepted method is inputted
            
            # get user input
            new_method=input( "What rootfinding method should be used? (bisect/secant/newton) \n");

            #  check that it will work with program
            if(new_method=="bisect" or new_method =="secant" or new_method == "newton"):
                self.method = new_method;
                system_ok = True;
            else:
                print( "Not an accepted method.\n");
                
        #  change function, make sure an available function is asked for
        system_ok = False;
        while(not system_ok): # prompts changes to method until accepted method is inputted

            #get key corresponding to user inputted function
            f_key= float( input( "What function should be selected? \n") );
            
            if( utils.SelectFunction(f_key)):
                self.f = utils.SelectFunction(f_key) ;
                system_ok = True;
            else:
                print( "Not an accepted function");
                     
        #  change interval, make sure the new value makes sense
        system_ok = False;
        while(not system_ok): # prompts changes to new length_scale until it is acceptable to store it in length_scale
            
            # get user input
            new_start = float(input( "Where should the interval start?") ); #remember user input is always a string
            new_stop = float(input( "Where should the interval stop?") ); 
            
            #  check that it will work with program
            if( new_stop > new_start):
                self.interval = (new_start, new_stop);
                system_ok = True;
            else:
                print( "Not an accepted interval.\n");
                
        #  change length_scale, make sure the new value makes sense
        system_ok = False;
        while(not system_ok): # prompts changes to new length_scale until it is acceptable to store it in length_scale
            
            # get user input
            new_length_scale = float(input( "What should the length_scale be set to?\n") ); #remember user input is always a string
            
            #  check that it will work with program
            if(new_length_scale < pow(10,5) and new_length_scale > 0): #ensures 0<tolerance<1
                self.length_scale = new_length_scale;
                system_ok = True;
            else:
                print( "Not an accepted length_scale.\n");
        
        #  change n_roots, make sure the new value makes sense
        system_ok = False;
        while(not system_ok): # prompts changes to new_n_roots until it is acceptable to store it in length_scale
            
            # get user input
            new_n_roots = float( input( "How many roots should be found? \n") );

            #  check that it will work with program
            if(new_n_roots > 1):
                self.n_roots = new_n_roots;
                system_ok = True;
            else:
                print( "Not a valid number of roots. \n");
        
        #  change verbose setting
        system_ok = False;
        while(not system_ok): # prompts changes
            
            # get user input
            self.is_verbose = utils.yesno("Do you want verbose output, helpful for troubleshooting?");
            
            system_ok = True;

# end of FinderParams class ####################################################


def IntervalControl( f, start, stop, length_scale, is_verbose, arguments ): 
    '''
    
    :param f: function to find crossing intervals of
    :param start: double, beginning of the overall interval 
    :param stop: double, end of the overall interval
    :param length_scale: double, rough sense of how far apart signs changes / roots are 
            for determining step size
    :param is_verbose: bool, tells to print out everything for debugging
    :param arguments, tuple, addl arguments that can be passed to f if required
    
    returns list cross_ints, list of tuples of endpoints of all the crossing intervals
    '''
    
    if( is_verbose):
        print( "Entering IntervalControl function.");
        print("Checking interval "+str( ( start, stop) ) );
        print("f(start) = "+str( f(start, *arguments) ) );
    
    # define return variable
    cross_ints=[]; #list of tuples of points that define an interval with a sign change    
         
    #subinterval endpoints that we need to find
    step = length_scale/10; #arbitrary small step size
    newstart = start;
    newstop = start+step;
      
    #check if start of interval is a root, if so add it to cross_int
    '''
    if ( abs(f(start, *arguments) ) < length_scale*pow(10,-5) ):
        if( is_verbose):
            print("f("+str(start)+") = 0, start of interval may be a root. \n")
        cross_ints.append( (start,start) );    '''

    #step through entire interval, sectioning off subintervals
    while ( newstop < stop ):

        #look for sign change
        if (f(newstart, *arguments)*f(newstop, *arguments) < 0): #this subinterval has a sign change
            
            # good tuple to add to cross ints
            cross_ints.append( (newstart,newstop) );  
        
        # advance to the next interval
        newstart += step; # start steps forward
        newstop += step; # stop steps forward
    
        ##### end while loop
        
    if( is_verbose):
        print("Found "+str(len( cross_ints)) +" crossing intervals. \n")

    return cross_ints;

######################################################### end IntervalControl


def AsymptoteCheck(f, start, stop, length_scale, is_verbose, arguments ):  #not in use, not fully functional
    '''
    Function that takes a sign crossing interval and checks for asymptotic behavior
    Note we are looking specifically for sign crossing, ie asymptote is reflected across y=0
    
    :param f: function, what we are finding root of
    :param start: double, beginning of interval
    :param stop: double, end of interval
    :param length_scale, double, rough sense of size of interval so we can mesh it finely
    :param is_verbose: bool, tells whether to print extra stuff for debugging
    :param arguments, tuple, addl arguments that can be passed to f if required
    
    returns bool saying whether interval had an symptote or not
    '''
    
    if( is_verbose):
        print("Entering AsymptoteCheck function.\n");
           
    #screen out infinite asymptotes
    # this will probably be an intensive portion of the code
    asym_flag = 0
    
    #### test 1: check height of max
    
    # get height of max by rootfinding on the derivative
    def f_prime(x): # function that gets derivative
        return scipy.misc.derivative(f, x, length_scale/100 );
    
    # rootfind for the derivative
    if( is_verbose):
        print( "Finding maximum on interval.\n");
    max_x =  finders.Bisect(f_prime, start, stop, length_scale/100, is_verbose );
    
    # check that max height exists (Bisect returns None if there's no sign change)
    if(max_x != None):
    
        #check that height of max exceeds threshhold
        max_height = pow(10, 6);
        if( f(max_x) > max_height ):
        
            if(is_verbose):
                print("Max height = "+str(max_height)+" exceeded, f(x) = "+str(f(max_x))+"\n");
            
                # update flag
                asym_flag += 1;
        
    #### return according to whether tests were passed
    if( asym_flag == 0): #good interval
        
        if(is_verbose):
            print("No asymptote found .\n");
        return False;
    
    elif( asym_flag == 1): #all tests passed, definite asymptote
        
        if(is_verbose):
            print("All tests passed, definite asymptote on interval.\n");
        return True;

#  end AsymptoteCheck function #####################################################################
    

def GetRoots( f,  method, start, stop, length_scale, n_roots, is_verbose, arguments ): 
    '''
    
    :param f: function, what we're finding roots of
    :param method: string, tells which finder.py function to call
    :param start: double, begins overall interval
    :param stop: double, ends overall interval
    :param length_scale: double, sense of size of each crossing interval so we can set tolerance
    :param n_roots: int, tells how many roots to find
    :param is_verbose: bool, tells to print everything for debugging
    :param arguments, tuple, addl arguments that can be passed to f if required
    
    returns np array, roots, all the roots found
    '''
    
    
    if( is_verbose):
        print( "Entering GetRoots function.\n");    
    
    # create list of roots, return object
    roots=np.array([]);
    
    # break up interval into sections which cross y=0
    # returns list of tuples which define intervals which should cross x=0
    cross_ints= IntervalControl(f, start, stop, length_scale, is_verbose, arguments); 
    
    # sort by method
    if( method == "bisect"):
        
        # go until either:
        # we have as many roots as asked for
        # we run out of cross ints 
        i = 0; # index of cross int tuples
        while( len(roots) < n_roots and i < len(cross_ints) ):
            
            #iter through tuples in cross ints
            #grab local start and stop from cross ints
            start=cross_ints[i][0];
            stop=cross_ints[i][1];
            
            if( is_verbose):
                print( "Rootfinding on subinterval ["+str(start)+","+str(stop)+"] \n");
                
                
            # set tolerance from length scale
            tolerance = length_scale * pow(10,-5);
        
            #get root and add to list
            my_root = finders.Bisect(f, start, stop, tolerance, is_verbose, arguments);
            
            #for bisection we have to screen out asymptotes
            # if bisect sees asymptotic behavior it returns None
            if( my_root != None): # root is good
                roots = np.append( roots, np.array(my_root) ); # done this way to avoid pointer reference issues
                
            # update index
            i += 1;
                    
        #### end while         
      
    # were enough roots found?      
    if (len(roots) < n_roots):
        print("Not as many roots found on interval as requested. \n")
    

    return roots;

######################################################### end GetRoots


def AllRoots( func_inputs, finder_inputs):
    '''
    
    :param func_inputs: a functions.InputParams() object which contains the experimental values to be passed to the transcendental functions
                        attributes: delta, band gap in meV
                                    V, well potential in meV
                                    d, well thickness in nm
                                    v_z, Fermi velocity in nm/s
    :param finder_inputs: a FinderParams() object which contains the settings that should be passed to the finder function
                        attributes: :param method: string, tells which finders.py function to call
                                    :param f, function, tells which functions.py math function to find roots of
                                    :param is_verbose: bool, tells to print everything out for debugging
                                    :param interval: tuple, start and stop points of rootfinding interval
                                    :param length_scale: double, sense of how far apart crossings/roots tend to be
                                    :param n_roots: int, how many roots to look for
                                    
    returns roots_all, array with all E eigenenergies on specified interval
    '''
    
    if( finder_inputs.is_verbose):
        print("Entering AllRoots.\n");
        
    
    # get all roots for E on E_start < E < E_stop 
    # remember we need to break up interval since different functions give eigenenergies in different regimes
    E_2delta = -2*abs(func_inputs.delta); # markes start of regime 2, -2*abs(delta) < E < 0
    E_0 = 0; # markes stop of regime 2, -2*abs(delta) < E < 0 
    
    # make sure that our inputted intervals span this special one
    if( finder_inputs.interval[0] > E_2delta or finder_inputs.interval[1] < E_0):
        print("Error in AllRoots(): \nSpecified interval "+str(finder_inputs.interval)+
              " does not span special regime "+str(E_2delta)+" < E < "+str(E_0)+"\n");
    
    # now grab functions
    f1_even = functions.pos_E_even; # even functions in regime 1, -inf < E < -2*abs(delta)
    f1_odd  = functions.pos_E_odd ; # odd  functions in regime 1, -inf < E < -2*abs(delta)    
    f2_even = functions.neg_E_even; # even functions in regime 2, -2*abs(delta) < E < 0
    f2_odd  = functions.neg_E_odd ; # odd  functions in regime 2, -2*abs(delta) < E < 0 
    f3_even = functions.pos_E_even; # even functions in regime 3, 0 < E < inf
    f3_odd  = functions.pos_E_odd ; # odd  functions in regime 3, 0 < E < inf
    
    #### treat each regime seperately. For each:
    # Init a new set of finder inputs for this specific regime
    # get even roots
    # get odd roots and combine
    
    #### regime 1:

    
    # define interval as tuple
    int1 = (finder_inputs.interval[0], E_2delta)
    
    if( finder_inputs.is_verbose):
        print("Rootfinding on regime 1, "+str(finder_inputs.interval[0])+" < E < "+str(E_2delta)+"\n")
    
    # init new finder param mostly based on old one
    # syntax is method, f , is_verbose , interval, length scale, n_roots
    # only changes are even function, new interval and make n_roots plenty big
    inputs1_even = FinderParams( finder_inputs.method, f1_even, finder_inputs.is_verbose, int1,
                             finder_inputs.length_scale, finder_inputs.n_roots);
                             
    # get even roots with the new inputs
    # this will be an np array
    # must pass func inputs as a tuple!
    roots1_even = inputs1_even.Run( (func_inputs,) ); #this method calls GetRoots()

    # repeat above for odd functions
    # only changes are odd function, new interval and make n_roots plenty big
    inputs1_odd = FinderParams( finder_inputs.method, f1_odd, finder_inputs.is_verbose, int1,
                             finder_inputs.length_scale, finder_inputs.n_roots);
                             
    # get odd roots with the new inputs
    # this will be an np array
    # must pass func inputs as a tuple!
    roots1_odd = inputs1_odd.Run( (func_inputs,) );
    
    #combine the two
    roots1 = np.append( roots1_even, roots1_odd );
    
    #### regime 2:
    
    # define interval as tuple
    int2 = ( E_2delta, E_0);
        
    if( finder_inputs.is_verbose):
        print("Rootfinding on regime 2, "+str(E_2delta)+" < E < "+str(E_0)+"\n")
        
    # init new finder param mostly based on old one
    # syntax is method, f , is_verbose , interval, length scale, n_roots
    # only changes are even function, new interval and make n_roots plenty big
    inputs2_even = FinderParams( finder_inputs.method, f2_even, finder_inputs.is_verbose, int2,
                             finder_inputs.length_scale, finder_inputs.n_roots);
                             
    # get even roots with the new inputs
    # this will be an np array
    # must pass func inputs as a tuple!
    roots2_even = inputs2_even.Run( (func_inputs,) );
    
    # repeat above for odd functions
    # only changes are odd function, new interval and make n_roots plenty big
    inputs2_odd = FinderParams( finder_inputs.method, f2_odd, finder_inputs.is_verbose, int2,
                             finder_inputs.length_scale, finder_inputs.n_roots);
                             
    # get odd roots with the new inputs
    # this will be an np array
    # must pass func inputs as a tuple!
    roots2_odd = inputs2_odd.Run( (func_inputs,) );
    
    #combine the two
    roots2 = np.append( roots2_even, roots2_odd );
    
    #### regime 3:
    
    # define interval as tuple
    int3 = ( E_0, finder_inputs.interval[1]);
        
    if( finder_inputs.is_verbose):
        print("Rootfinding on regime 3, "+str(E_0)+" < E < "+str(finder_inputs.interval[1])+"\n");
    
    # init new finder param mostly based on old one
    # syntax is method, f , is_verbose , interval, length scale, n_roots
    # only changes are even function, new interval and make n_roots plenty big
    inputs3_even = FinderParams( finder_inputs.method, f3_even, finder_inputs.is_verbose, int3,
                             finder_inputs.length_scale, finder_inputs.n_roots);
                             
    # get even roots with the new inputs
    # this will be an np array
    # must pass func inputs as a tuple!
    roots3_even = inputs3_even.Run( (func_inputs,) );
    
    # repeat above for odd functions
    # only changes are odd function, new interval and make n_roots plenty big
    inputs3_odd = FinderParams( finder_inputs.method, f3_odd, finder_inputs.is_verbose, int3,
                             finder_inputs.length_scale, finder_inputs.n_roots);
                             
    # get odd roots with the new inputs
    # this will be an np array
    # must pass func inputs as a tuple!
    roots3_odd = inputs3_odd.Run( (func_inputs,) );
    
    #combine the two
    roots3 = np.append( roots3_even, roots3_odd )
    
    #### We have found all the roots !
    
    # now combine them into 1 big thing
    roots_all = np.append( roots1, np.append(roots2, roots3, ));
    
    return roots_all; 
    
# end AllRoots function #####################################################################


def FormatRoots( roots, E_start, E_critical, E_stop, n_roots, is_verbose ):
    '''
    Function that takes root array found by all_roots and standardizes it by:
    Having the same number of roots above and below E_critical (using dummy roots)
    sorting the roots in ascending order
    
    :param roots: array, contains all eigenenergies on E_start < E < E_stop
    :param E_start: double, begins E domain
    :param E_critical: double, horizontal asymptote to which eigenenergies converge 
                    from above and below, and are symmetrically distributed around
    :param E_stop: double, ends E domain
    :param n_roots: int, how many roots to put above and below E-critical
    :param is_verbose: bool, tells whether to make extra print outs
    
    returns good_roots: np array of roots, with right number of them, and in ascending order
    '''
    
    if(is_verbose):
        print("Entering FormatRoots function. \n");
        
    # create arrays of standardized size to hold roots
    # if not enough roots were found, we want the rest to be filled by interval endpoint values
    # these will be either E_start or E_stop
    E_lower = np.full(n_roots, E_start);
    E_upper = np.full(n_roots, E_stop);
    
    #index vars to help us (partially) fill these arrays with found roots
    i_upper = 0;
    i_lower = 0;
        
    # separate out roots above and below the horizontal asymptote E_critical)
    for root0 in roots:
        
        #check if below E_critical
        if( root0 < E_critical ): 
            
            #then place in array
            E_lower[i_lower] = root0;
            
            #and update counter
            i_lower += 1;
        
        #otherwise above E critical
        elif( root0 > E_critical ): 
            
            #then place in array
            E_upper[i_upper] = root0;
            
            #and update counter
            i_upper += 1;
            
        else: #Should not get here
            print("Error in FormatRoots function: \nroot neither above nor below E_critical");
            
            #### end for loop
            
    #### roots have been successfuly divided and the number above/below E_critical standardized!
    
    if( is_verbose):
        print("Roots above "+str(E_critical)+": "+str(E_upper)+"\n");
        print("Roots below "+str(E_critical)+": "+str(E_lower)+"\n");
    
    # now combine back together into a new array
    good_roots = np.append(E_upper, E_lower);

    #sort the result
    good_roots.sort();
    
    return good_roots;

#  End FormatRoots() function ############################################################


def VaryRoots( func_inputs, finder_inputs, x_param, x_start, x_stop):
    '''
    
    :param func_inputs: functions.InputParams object, sets experimental values for quantum well parameters
    :param finder_inputs: main.FinderParams object, tells how to go about rootfinding
    :param x_param: string, tells which experimental input (delta, V, d, v_z) to vary as x
    :param x_start: double, where to start x domain
    :param x_stop: double, where to stop x domain
    
    returns tuple, 
    root_arrays_list which holds all the root values, as varied over x, w/ each energy level a seperate array
    an array of values for the x variable
    '''
    
    if(finder_inputs.is_verbose):
        print("Entering VaryRoots function. \n")
        
    # to get roots from AllRoots we need a finder_inputs (UserInputs object) and func_inputs (InputParams object)
    # finder inputs is same for all and is passed as an argument
    # different func inputs are created below dependein on the x_param we are varying
        
    #### first sort by which function param we are varying as x
    if( x_param == 'd'): #varying well thickness
        
        if(finder_inputs.is_verbose):
            print("Varying well thickness "+str(x_start)+" < d < "+str(x_stop)+" nm.\n")
        
        # here create different func_input objects
        # create an array of d vals
        n_d_vals = int( (x_stop - x_start)*2 );
        d_arr = np.linspace(x_start, x_stop,n_d_vals);
        
        # iter through d values
        for i in range( n_d_vals ):
            
            d_val =d_arr[i];
            
            print( d_val);
            
            # skip d = 10
            if (abs( d_val - 10) < 0.2):
                print ("\n \n \n \n \n" ); 
                print ("************************************************* d_val reassigned");
                d_val = 10.5 ;
            
            # first time through, create arrays for all the roots
            if( i == 0):
                
                #store them in a list
                root_arrays_list = [];
                
                #make 2*finder_inputs.n_roots total arrays, because there are n_roots above E_critical and n_roots below
                # there is 1 array/line on graph for every root
                for j in range( 2*finder_inputs.n_roots):
                    # there is one val in each array, for each d val
                    root_arrays_list.append( np.zeros(n_d_vals) );
                
            
            # create a func_input
            # syntax is delta, V, d, v_z
            my_func_inputs = functions.InputParams( func_inputs.delta, func_inputs.V, d_val, func_inputs.v_z  );
            
            # Get and Format roots
            roots_for_d = AllRoots(my_func_inputs, finder_inputs);
            formatted_roots = FormatRoots(roots_for_d, finder_inputs.interval[0], func_inputs.delta, finder_inputs.interval[1], finder_inputs.n_roots, finder_inputs.is_verbose)
            
            # place formatted roots in arrays in preparation for plotting
            # j is index of root in formatted_roots, and index of root array in roots_array_list
            # i is index of d_val and index these roots should be placed in w/in each root_array
            for j in range( 2*finder_inputs.n_roots):
                root_arrays_list[j][i] = formatted_roots[j];
                
        #### end for loop over i
        
        if (finder_inputs.is_verbose):
            utils.PrintArrayList( root_arrays_list, "Roots");
            
        # return root arrays and d values array
        return root_arrays_list, d_arr;
    
    #### end if statement for x_param == d
    
    # vary fermi velocity
    if( x_param == 'vz'):
        
        if(finder_inputs.is_verbose):
            print("Varying fermi velocity "+str(x_start)+" < v_z < "+str(x_stop)+" nm/s\n")
        
        # here create different func_input objects
        # create an array of v_z vals
        n_vz_vals = int( 100 );
        vz_arr = np.linspace(x_start, x_stop,n_vz_vals);
        
        # iter through vz values
        for i in range( n_vz_vals ):
            
            vz_val = vz_arr[i];
            
            print( vz_val);
            
            # skip vz = 5.5e14
            if (abs( vz_val - 5.4e14) < 2e13):
                #print ("\n \n \n \n \n" ); 
                print ("************************************************* vz_val reassigned");
                vz_val = 5.6*pow(10,14) ;
                vz_arr[i] = 5.6*pow(10,14) ;
            
            # first time through, create arrays for all the roots
            if( i == 0):
                
                #store them in a list
                root_arrays_list = [];
                
                #make 2*finder_inputs.n_roots total arrays, because there are n_roots above E_critical and n_roots below
                # there is 1 array/line on graph for every root
                for j in range( 2*finder_inputs.n_roots):
                    # there is one val in each array, for each d val
                    root_arrays_list.append( np.zeros(n_vz_vals) );
                
            
            # create a func_input
            # syntax is delta, V, d, v_z
            my_func_inputs = functions.InputParams( func_inputs.delta, func_inputs.V, func_inputs.d, vz_val  );
            
            # Get and Format roots
            roots_for_vz = AllRoots(my_func_inputs, finder_inputs);
            formatted_roots = FormatRoots(roots_for_vz, finder_inputs.interval[0], func_inputs.delta, finder_inputs.interval[1], finder_inputs.n_roots, finder_inputs.is_verbose)
            
            # place formatted roots in arrays in preparation for plotting
            # j is index of root in formatted_roots, and index of root array in roots_array_list
            # i is index of d_val and index these roots should be placed in w/in each root_array
            for j in range( 2*finder_inputs.n_roots):
                root_arrays_list[j][i] = formatted_roots[j];
                
        #### end for loop over i
        
        if (finder_inputs.is_verbose):
            utils.PrintArrayList( root_arrays_list, "Roots");
            
        # return root arrays and d values array
        return root_arrays_list, vz_arr;
    
    #### end if statement for x_param == vz
    
    # vary well depth
    if( x_param == 'V'):
        
        if(finder_inputs.is_verbose):
            print("Varying well depth "+str(x_start)+" < V < "+str(x_stop)+" meV")
        
        # here create different func_input objects
        # create an array of V vals
        n_V_vals = int( 10 );
        V_arr = np.linspace(x_start, x_stop,n_V_vals);
        
        # iter through vz values
        for i in range( n_V_vals ):
            
            V_val = V_arr[i];
            
            print( V_val);
            
            # skip vz = 5.5e14
            #if (abs( V_val - 5.4e14) < 2e13):
                #print ("\n \n \n \n \n" ); 
                #print ("************************************************* vz_val reassigned");
                #vz_val = 5.6*pow(10,14) ;
                #vz_arr[i] = 5.6*pow(10,14) ;
            
            # first time through, create arrays for all the roots
            if( i == 0):
                
                #store them in a list
                root_arrays_list = [];
                
                #make 2*finder_inputs.n_roots total arrays, because there are n_roots above E_critical and n_roots below
                # there is 1 array/line on graph for every root
                for j in range( 2*finder_inputs.n_roots):
                    # there is one val in each array, for each d val
                    root_arrays_list.append( np.zeros(n_V_vals) );
                
            
            # create a func_input
            # syntax is delta, V, d, v_z
            my_func_inputs = functions.InputParams( func_inputs.delta, V_val, func_inputs.d, func_inputs.v_z  );
            
            # redefine energy regime
            E_start = -V_val;
            E_critical = exp_inputs.delta;
            E_stop = V_val - 2*abs( exp_inputs.delta);
            
            # Get and Format roots
            roots_for_V = AllRoots(my_func_inputs, finder_inputs);
            formatted_roots = FormatRoots(roots_for_V, E_start, E_critical, E_stop, finder_inputs.n_roots, finder_inputs.is_verbose)
            
            # place formatted roots in arrays in preparation for plotting
            # j is index of root in formatted_roots, and index of root array in roots_array_list
            # i is index of d_val and index these roots should be placed in w/in each root_array
            for j in range( 2*finder_inputs.n_roots):
                root_arrays_list[j][i] = formatted_roots[j];
                
        #### end for loop over i
        
        if (finder_inputs.is_verbose):
            utils.PrintArrayList( root_arrays_list, "Roots");
            
        # return root arrays and d values array
        return root_arrays_list, V_arr;
    
    #### end if statement for x_param == V
        
    
    # vary band gap
    if( x_param == 'delta'):
        
        if(finder_inputs.is_verbose):
            print("Varying band gap "+str(x_start)+" < $ \Delta $ < "+str(x_stop)+" meV")
        
        # here create different func_input objects
        # create an array of V vals
        n_delta_vals = int( 100 );
        delta_arr = np.linspace(x_start, x_stop,n_delta_vals);
        
        # iter through vz values
        for i in range( n_delta_vals ):
            
            delta_val = delta_arr[i];
            
            print( delta_val);
            
            # skip delta = -21.5
            if (abs( delta_val + 22) < 1.5):
                #print ("\n \n \n \n \n" ); 
                print ("************************************************* delta_val reassigned");
                delta_val = -20 ;
            
            # first time through, create arrays for all the roots
            if( i == 0):
                
                #store them in a list
                root_arrays_list = [];
                
                #make 2*finder_inputs.n_roots total arrays, because there are n_roots above E_critical and n_roots below
                # there is 1 array/line on graph for every root
                for j in range( 2*finder_inputs.n_roots):
                    # there is one val in each array, for each d val
                    root_arrays_list.append( np.zeros(n_delta_vals) );
                
            
            # create a func_input
            # syntax is delta, V, d, v_z
            my_func_inputs = functions.InputParams( delta_val, func_inputs.V, func_inputs.d, func_inputs.v_z  );
            
            
            # redefine energy regime
            E_start = -exp_inputs.V;
            E_critical = delta_val;
            E_stop = exp_inputs.V - 2*abs( delta_val);
            
            # Get and Format roots
            roots_for_V = AllRoots(my_func_inputs, finder_inputs);
            formatted_roots = FormatRoots(roots_for_V, E_start, E_critical, E_stop, finder_inputs.n_roots, finder_inputs.is_verbose)
            
            # place formatted roots in arrays in preparation for plotting
            # j is index of root in formatted_roots, and index of root array in roots_array_list
            # i is index of d_val and index these roots should be placed in w/in each root_array
            for j in range( 2*finder_inputs.n_roots):
                root_arrays_list[j][i] = formatted_roots[j];
                
        #### end for loop over i
        
        if (finder_inputs.is_verbose):
            utils.PrintArrayList( root_arrays_list, "Roots");
            
        # return root arrays and d values array
        return root_arrays_list, delta_arr;
    
    #### end if statement for x_param == V
        
    else:
        print("Error in VaryRoots function: ");
        print("Cannot vary over parameter "+x_param);

# End VaryRoots function ################################################################


def PlotRoots(root_arrays_list, x_arr, x_param, x_start, x_stop, exp_inputs):
    '''
    Plots energy levels of quantum well as function of whatever x param we varied
    over in VaryRoots
    
    :param root_arrays_list: list of 1D np arrays which holds all the root 
                values, as varied over x, w/ each energy level a seperate array
    :param x_arr: 1d np array, an array of values for the x variable, whether d, V, delta etc 
    :param x_param: string, tells which experimental input (delta, V, d, v_z) to vary as x
    :param x_start: double, where to start x domain
    :param x_stop: double, where to stop x domain
    :param exp_inputs: functions.InputParams object, sets experimental values for quantum well parameters
    
    returns None, just shows plots instead
    
    '''
    
    # determine eigenenergies at x_stop to check parity
    # we have to adjust how to do this based on the x param too 
    if( x_param == 'd'):
        # work with the last x val always because this is where levels are most spaced out!
        eigen_i = -1;
    elif( x_param == 'vz'):
        # work with first E val because this is where levels are spaced
        eigen_i = 2;
    else:
        eigen_i = -1;
        
    # create eigenenrgies
    eigenenergies = [];
    for E_level in root_arrays_list:
        eigenenergies.append( E_level[eigen_i]);
        
    # convert to np array
    eigenenergies = np.array( eigenenergies);
    
    #iter through root arrays and plot
    for E_level in root_arrays_list:
        
        #### determine labels and line formatting according to parity
        
        # determine parity of state
        is_odd = matrices.Parity_E(E_level[eigen_i], eigenenergies, exp_inputs.delta);
        
        #plot accordingly
        if( is_odd): # odd solution so make the plot dashed
            plt.plot(x_arr, E_level, color = 'green', linestyle = 'dashed');
        else: #even solution plot solid
            plt.plot(x_arr, E_level, color = 'green');            
        
    # make titles, labels etc according to x_param
    title = "Quantum well eigenergies varied with "+x_param;
    y_label = "Energy [meV]";
    if(x_param == 'd'):
        x_label = "Well thickness [nm]";
    elif(x_param == 'vz'):
        title = "Quantum well eigenergies varied with v_z";
        x_label = "Fermi velocity [nm/s]";
    elif( x_param == 'V'):
        x_label = "Well depth [meV]";
    elif( x_param == 'delta'):
        x_label = "Band gap [meV]";
        title = "Quantum well energies varied with $\Delta$"
    
    # extra formatting depending on x param
    if x_param == 'd': # case that we are varying thickness
        
        # draw line to demarcate TIS region
        TIS_y_arr = np.full( len(x_arr), exp_inputs.delta)
        plt.plot(x_arr, TIS_y_arr, color='red',linewidth=41, alpha=0.2);
            
        # make legend
        even_line = mlines.Line2D([],[], color = 'green', label = 'Even states')
        odd_line = mlines.Line2D([],[], color = 'green', linestyle = 'dashed', label = 'Odd states')
        TIS_line = mlines.Line2D( [],[], color='red', linewidth=8, alpha=0.2, label="TIS region");
        plt.legend(handles=[ even_line, odd_line, TIS_line], loc = (0.05,0.7));
        
        #make a textbox showing experimental params
        exp_params_string = ""
        
        # add info of experimental params to box as desired
        for line in exp_inputs.DeclareSelf(): #this returns a list of strings describing each exp param
            # skip the line which tells what d is as we are changing d
            if "d = " not in line:
                exp_params_string += line;
            
        # put the textbox on the plot
        plt.text( 1,-150, exp_params_string);
            
        # format the resulting plot
        plt.suptitle(title); #set title
        plt.xlabel(x_label);
        plt.ylabel(y_label);
        plt.grid(True, 'both', 'both');
        x_start, x_stop = plt.xlim(); # get x axis limits
        plt.xlim(x_start -5, x_stop); # extend x axis to make room for legend
        y_start, y_stop = plt.ylim(); #get y axis limits
        plt.ylim( y_start + 3, y_stop - 3); # slightly shorten y axis
    
    # extra formatting depending on x param
    if x_param == 'vz': # case that we are varying fermi velocity
        
        # draw line to demarcate TIS region
        TIS_y_arr = np.full( len(x_arr), exp_inputs.delta)
        plt.plot(x_arr, TIS_y_arr, color='red',linewidth=41, alpha=0.2);
            
        # make legend
        even_line = mlines.Line2D([],[], color = 'green', label = 'Even states')
        odd_line = mlines.Line2D([],[], color = 'green', linestyle = 'dashed', label = 'Odd states')
        TIS_line = mlines.Line2D( [],[], color='red', linewidth=8, alpha=0.2, label="TIS region");
        plt.legend(handles=[ even_line, odd_line, TIS_line], loc = (0.7,0.7) );
        
        #make a textbox showing experimental params
        exp_params_string = ""
        
        # add info of experimental params to box as desired
        for line in exp_inputs.DeclareSelf(): #this returns a list of strings describing each exp param
            # skip the line which tells what d is as we are changing d
            if "$v_{z}$ =" not in line:
                exp_params_string += line;
            
        # put the textbox on the plot
        plt.text( 7*pow(10,14),-150, exp_params_string);
            
        # format the resulting plot
        plt.suptitle(title); #set title
        plt.xlabel(x_label);
        plt.ylabel(y_label);
        plt.grid(True, 'both', 'both');
        x_start, x_stop = plt.xlim(); # get x axis limits
        plt.xlim(x_start -5, x_stop); # extend x axis to make room for legend
        y_start, y_stop = plt.ylim(); #get y axis limits
        plt.ylim( y_start + 3, y_stop - 3); # slightly shorten y axis
    
    # extra formatting depending on x param
    if x_param == 'V': # case that we are varying well depth
        
        # draw line to demarcate TIS region
        TIS_y_arr = np.full( len(x_arr), exp_inputs.delta)
        plt.plot(x_arr, TIS_y_arr, color='red', linewidth=21, alpha=0.2);
            
        # make legend
        even_line = mlines.Line2D([],[], color = 'green', label = 'Even states')
        odd_line = mlines.Line2D([],[], color = 'green', linestyle = 'dashed', label = 'Odd states')
        TIS_line = mlines.Line2D( [],[], color='red', linewidth=8, alpha=0.2, label="TIS region");
        plt.legend(handles=[ even_line, odd_line, TIS_line], loc = "upper left" );
        
        #make a textbox showing experimental params
        exp_params_string = ""
        
        # add info of experimental params to box as desired
        for line in exp_inputs.DeclareSelf(): #this returns a list of strings describing each exp param
            # skip the line which tells what d is as we are changing d
            if "V =" not in line:
                exp_params_string += line;
            
        # put the textbox on the plot
        plt.text( 40,-300, exp_params_string);
            
        # format the resulting plot
        plt.suptitle(title); #set title
        plt.xlabel(x_label);
        plt.ylabel(y_label);
        plt.grid(True, 'both', 'both');
        x_start, x_stop = plt.xlim(); # get x axis limits
        plt.xlim(x_start -5, x_stop); # extend x axis to make room for legend
        y_start, y_stop = plt.ylim(); #get y axis limits
        plt.ylim( y_start + 3, y_stop - 3); # slightly shorten y axis
    
    # extra formatting depending on x param
    if x_param == 'delta': # case that we are varying band gap
        
        # draw line to demarcate TIS region
        # this is a function of delta so isnt just a rectangle
        # use a bar graph
                
        TIS_bar = plt.bar( x_arr, 2*x_arr , width = (x_stop-x_start+7)/len(x_arr), color='red', alpha = 0.2, linewidth = 0 );
        #plt.plot(x_arr, TIS_y_arr, color='red', linewidth=21, alpha=0.2);
            
        # make legend
        even_line = mlines.Line2D([],[], color = 'green', label = 'Even states')
        odd_line = mlines.Line2D([],[], color = 'green', linestyle = 'dashed', label = 'Odd states')
        TIS_line = mlines.Line2D( [],[], color='red', linewidth=8, alpha=0.2, label="TIS region");
        plt.legend(handles=[ even_line, odd_line, TIS_line], loc = "upper left" );
        
        #make a textbox showing experimental params
        exp_params_string = ""
        
        # add info of experimental params to box as desired
        for line in exp_inputs.DeclareSelf(): #this returns a list of strings describing each exp param
            # skip the line which tells what d is as we are changing d
            if "$\Delta$ = " not in line:
                exp_params_string += line;
            
        # put the textbox on the plot
        plt.text( -90,-180, exp_params_string);
            
        # format the resulting plot
        plt.suptitle(title); #set title
        plt.xlabel(x_label);
        plt.ylabel(y_label);
        plt.grid(True, 'both', 'both');
        x_start, x_stop = plt.xlim(); # get x axis limits
        plt.xlim(x_start -5, x_stop); # extend x axis to make room for legend
        y_start, y_stop = plt.ylim(); #get y axis limits
        plt.ylim( y_start + 3, y_stop - 3); # slightly shorten y axis
        
        
    #show all on same plot
    plt.show();
                


if __name__ == '__main__': 
################################################################################
# THIS IS THE MAIN PROGRAM. IT IS DESIGNED TO BE USER FRIENDLY, SO
# MOST OF THE INPUTS CAN BE CONTROLLED FROM THE COMMAND LINE.
# IN OTHER WORDS SIMPLY RUN IT AND RESPOND TO THE PROMPTS 
################################################################################

    print( "Welcome to the rootfinding program.\n" );
    
    #### define experimental parameters
    delta = -25; #band gap in meV
    V = 175; # potential in meV
    d = 12; # well thickness in nm
    v_z = 4.5*pow(10,14) # matrix element in nm/s
    
    # set domain of energies, in meV
    E_start  =-V; #### make sure this is a double!
    E_critical = delta
    E_stop = V - 2*abs(delta) ; #### make sure this is a double!
    
    # set n roots to find
    n_roots = 20;
    
    #  display inputs to rootfinder
    # FinderParams syntax is ( method='bisect', f = functions.f1 , is_verbose=True , interval=(-100,100), 
                # length_scale=10, n_roots=1,
    my_inputs = FinderParams('bisect',functions.f1, False, (E_start, E_stop), 10, n_roots );
    my_inputs.DeclareSelf();
        
    #  check if user accepts settings
    #user_ok = utils.yesno( "Are these settings acceptable? ");
    #if (not user_ok): #change settings if needed
        #my_inputs.UserInput();

    #  End of user changing settings #############################################
    #  Actual rootfinding below ##################################################
    
    # set experimental values using InputParams() object
    exp_inputs = functions.InputParams( delta, V, d, v_z); #default vals for now
    arguments = (exp_inputs,); # exp_inputs needs to be passed as a tuple to make use of * syntax
                                # this flexibility means we can pass as many arguments as required
    
    # get all eigienergies (even, odd, all E values) using AllRoots
    my_roots = AllRoots( exp_inputs, my_inputs );
    
    '''
    print(my_roots)
    for root in my_roots:
        if( root < 0 and root > 2*exp_inputs.delta):
            E = root;
    kappa = (1/exp_inputs.hv)*np.sqrt(- E*(E+ 2*abs(exp_inputs.delta) ) )
    rho = (1/exp_inputs.hv)*np.sqrt( (E+exp_inputs.V)*( -E - 2*abs(exp_inputs.delta) + V) )
    def chi(z): # piecewise in z
        if( abs(z) <= exp_inputs.d/2):
            return np.cosh(kappa*z);
        else:
            return np.cosh(kappa*exp_inputs.d/2)*np.exp(-rho*(abs(z) - exp_inputs.d/2) );
        
    def chi2(z):
        return chi(z)*chi(z)
    
    utils.PlotF(chi2, -1.5*exp_inputs.d/2, 1.5*exp_inputs.d/2, "z [nm]", "$\chi(z)$", my_title = 
                "Even TIS probability density for d = "+str(exp_inputs.d))
    
    '''
    
    #  Found all the roots #####################################################
    #  Now organize them for graphing ##########################################  
    formatted_roots = FormatRoots(my_roots, E_start, E_critical, E_stop, n_roots, True); 

    # print formatted roots nicely
    for r in formatted_roots:
        if( r != E_start and r != E_stop ):
            print("%.2f"%r);
    
    # plot the roots against a varied InputParams attribute (delta, V, d, v_z)
    # in this case d

    # lay out the range of d
    x_param = 'delta';
    x_start = -V/2
    x_stop = -10
    
    # get roots over range of d
    roots_over_d, d_arr = VaryRoots(exp_inputs, my_inputs, x_param, x_start, x_stop);
    
    # plot roots over range of d
    PlotRoots(roots_over_d, d_arr, x_param, x_start, x_stop, exp_inputs); 
    
