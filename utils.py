'''
Created on Jan 17, 2020

@author: Christian
'''

import functions 
from matplotlib import pyplot
import numpy as np

################################################################################
# 
# Misc. functions that help with the Rootfinding project
#
################################################################################


def yesno( question ): 
    '''
    converts user input y/n into a bool
    
    :param question: string, question that is outputted to the user
    
    returns bool
    '''
    
    
    yn1 = input( question + " (y/n) \n"); #  stores the input

    if( yn1=='y' or yn1=='Y' ): #  convert to a bool
        return True;
    
    else:
        return False;
    
    
    
def SelectFunction ( my_key ): 
    '''
    Acts as a switch case block 
    Uses dictionary to retrieve function that you want to rootfind for from functions.py
    
    :param key: int, selects the function
    '''
    
    # use a dictionary to link functions to int keys
    switcher = {
        1: functions.f1,
        2: functions.f2,
        3: functions.f3,
        4: functions.f4,
        5: functions.f5,
        6: functions.pos_E_even,
        };
    
    return switcher.get(my_key, False); #false is default (no function included in dictionary)


def PrintArrayList( arrays_list, array_title):
    '''
    Nicely formats and prints the arrays contained in the list
    
    :param arrays_list: list of np arrays
    array_title: string, tells what to call each array
    '''
    
    # run through all arrays in the list
    for i in range(len (arrays_list)):
        printstring = array_title + ' at i = '+str(i)+': ' ;
        printstring += str(arrays_list[i]) + '\n' ;
        
        print( printstring);
        
        
def CheckIntegral( scipy_tup, func_name):
    
    # get integral value and error from tuple
    I_val, error = scipy_tup[0], scipy_tup[1];
    
    #check that the error is reasonable on real integral
    # case integral is nonzero
    if( abs(I_val) > 1/100 and abs(error)*100 > abs(I_val) ):
        print("Error in "+func_name+" integral error too large.")
        return None
        
    #case that integral is near zero
    elif( abs(I_val) < 1/100 and abs(error) > 1/100 ):
        print("Error in "+func_name+" integral error too large, integral near zero.")
        return None;
    
    else: # integral is ok
        return I_val;

def PlotF(f,start,stop,x_label='x axis',y_label='y axis',lab='f(x)',\
           my_title='Graph',arguments=()):
    """
    This function plots f(x) on a specified domain
    Arguments:
        f(function): the function to be plotted
        start, stop (doubles): the range to be plotted on
        x_label(string) the x axis label
        y_label(string) the y axis label
        label(optional): curve label
        arguments(tuple): any additional arguments to be passed to f(x)
    returns:
    pyplot figure object
    """
    
    #initiate figure object, plot object
    my_fig, my_plot = pyplot.subplots();

    #if no additional arguments, args is an empty tuple, so just plot

    if arguments==():
        
        # generate x array using linspace
        x_arr = np.linspace(start, stop, 1000)
        
        # generate y array by passing it the x array
        # this requires that the function be fully vectorized
        try:
            y_arr = f(x_arr)
            
        except: # the function is not vectorized
            # iter thru instead
            y_arr = np.zeros( len(x_arr) );
            for i in range( len(x_arr) ):
                y_arr[i] = f( x_arr[i])

        my_plot.plot(x_arr,y_arr,label=lab);

    #handle plotting if there are additional arguments       
    else:

        # generate x array using linspace
        x_arr = np.linspace(start, stop, 1000)
        
        # generate y array by passing it the x array
        # this requires that the function be fully vectorized
        y_arr = f(x_arr, *arguments)

        my_plot.plot(x_arr,y_arr,label=lab);

    #format the graph
    my_plot.set( xlabel = x_label, ylabel = y_label, title = my_title );
    my_plot.grid()
    
    #show the graph
    pyplot.show();

    return my_fig;

#  End of PlotF ####################################################################


        
    
    
    
