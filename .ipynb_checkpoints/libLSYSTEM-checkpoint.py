"""
pyLSYSTEM
library for L-system modelling
2024-04-21
Georg Kaufmann
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import libLSYSTEM
d2r = np.pi / 180.


#================================#
def lsystemCreate(word,rules):
    """
    L-system
    Apply rules to word, using key->value dictionary
    Input:
    word   - input string
    rules  - rules, with key original character, value the new sequence
    Output:
    newword - new string
    """
    newword = ''
    for letter in word:
        new = rules[letter]
        newword += new
    return newword


#================================#
def lsystemIterate(word,rules,iter=3):
    """
    L-system
    perform iterations
    Input:
    word   - input string
    rules  - rules, with key original character, value the new sequence
    iter   - iteration counter
    Output:
    word   - new replaced string for final iteration
    """
    for i in range(iter):
        word = libLSYSTEM.lsystemCreate(word, rules)
    return word


#================================#
def lsystemPlot(word,length=0.2,angle=90.,scale=1.,x=0,y=0,d=90,showDot=False,linewidth=1):
    """
    L-system
    Plot final word following rules defined for moving the cursor
    Input:
    word   - input string
    length - length of forward step
    angle  - angle to turn
    x,y,d  - initial position and view angle
    Output:
    (plot)
    """
    ax = plt.axes()
    #ax.set_xlim([-10,10])
    #ax.set_ylim([-10,10])
    ax.grid()
    ax.xaxis.set_major_locator(ticker.NullLocator())
    ax.yaxis.set_major_locator(ticker.NullLocator())
    ax.set_aspect('equal', 'box')
    # define empty stack
    stack = []
    # loop over letters in word
    for letter in word:
        if letter in ["F", "G", "R", "L"]:
            xfrom=x;yfrom=y;
            x = xfrom + length*np.cos(d*d2r)
            y = yfrom + length*np.sin(d*d2r)
            ax.plot([xfrom,x],[yfrom,y],lw=linewidth,color='brown')
            if (showDot): ax.plot(x,y,lw=0,marker='.',markersize=10,color='black')
        elif letter == "f":
            x = xfrom + length*np.cos(d*d2r)
            y = yfrom + length*np.sin(d*d2r)
        elif letter == "+":
            d += angle
        elif letter == "-":
            d -= angle
        elif letter == "[":
            stack.append([x,y,d])
        elif letter == "]":
            [x,y,d] = stack.pop()
        elif letter == "<":
            length /= scale
        elif letter == ">":
            length *= scale
    return


#================================#
#================================#
def cellularAutomatonSetup(nx=101):
    """
    setup of initial step for 1D cellular automaton
    input:
        nx - grid dimension
    returns:
        grid - 1D array of initial state, holding zeros and at one position 1 as initial seed
    """
    # set all cells to state 0
    grid = np.zeros(nx,dtype=int).reshape(1,nx)
    # set one cell to state 1
    grid[0,int(nx/2)] = 1
    #grid[0,0] = 1
    return grid


#================================#
def cellularAutomatonPlot(grid):
    """
    plot of 1D cellular automaton
    input:
        grid - 2D array of states 0 and 1
        nx - grid dimension
    returns:
       (none)
    """
    fig,ax=plt.subplots(figsize=(10,10)) 
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)
    ax.tick_params(axis='both', which='both', length=0)
    plt.imshow(grid,vmin=0,vmax=1,interpolation='none',cmap='viridis')
    return


#================================#
def cellularAutomatonUpdate(grid):
    """
    Update a grid of points with states 0 and 1 by a specific update rule
    """
    nx = grid.shape[1]
    # shift with rolling over edge
    shift_left = np.roll(grid[grid.shape[0]-1,:],-1)
    shift_right = np.roll(grid[grid.shape[0]-1,:],+1)
    # shift without rolling over edge
    shift_left[nx-1] = 0
    shift_right[0] = 0
    newrow = np.where((shift_left + shift_right)==1,1,0)
    grid = np.vstack((grid,newrow))
    return grid
