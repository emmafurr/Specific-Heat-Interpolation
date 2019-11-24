# This is an example script for the final.
# Example 2

#Imports
import numpy as np
import mpmath as mp
import math
import pandas as pd
import string
import matplotlib.pyplot as plt
import scipy
from mpl_toolkits.mplot3d import Axes3D
from sympy import Symbol, poly, factor

# Natural Cubic Spline Algorithm


def main():


    # Goal:
    # Our data set:
    #----------------------------------------










    #-----------------------------------------

    print('We will use now Cubic Spline Interoplation with our data set:')
    # splineConstructor(list_x, list_fx, factor_it = 'yes')

    
    #--------------------------------
     print('We will use now Lagrange Interoplation with our data set:')
    # lagrange_polynomial(nodes, function, x, degree_specification)

    
    #--------------------------------
     print('We will use now Nevilles method with a data set:')
    # nevillesMethod(x, list_x, list_fx, Q_table = None, individual = 'no',notable = 'no')


    #--------------------------------

    # Find relative error
    #relative_error(real, approx)

    # Compare results

def relative_error(real, approx):
    y = abs(real - approx) / abs(real)
    return y

def absolute_error(real, approx):
    y = abs(real - approx)
    return y

def vec(m): z = [0]*m ; return(z)
def splineConstructor(list_x, list_fx, factor_it = 'yes'):
    n =  len(list_x)
    h = vec(n-1) ; alpha = vec(n-1) ; l = vec(n+1)
    u = vec(n) ; z = vec(n+1) ; b = vec(n) ; c = vec(n+1); d = vec(n)        
    for i in range(0, n - 1):
        h[i] = list_x[i+1] - list_x[i]  
    for i in range(1, n - 1):
        alpha[i] = (3./h[i])*(list_fx[i+1] - list_fx[i])-(3./h[i - 1])*(list_fx[i] - list_fx[i-1])
        l[0] = 1 ; u[0] = 0 ; z[0] = 0
    for i in range(1, n - 1):
        l[i] = 2*(list_x[i+1] - list_x[i-1]) - h[i - 1]*u[i - 1]
        u[i] = h[i]/l[i]
        z[i] = (alpha[i] - h[i - 1]*z[i - 1])/l[i]
        l[len(l)-1] = 1 ; z[len(z)-1] = 0 ; c[len(c)-1] = 0
    for j in reversed(range(n - 1)):      
        c[j] = z[j] - u[j]*c[j+1]
        b[j] = (list_fx[j + 1] - list_fx[j])/h[j]- h[j]*(c[j + 1] + 2*c[j])/3.
        d[j] = (c[j + 1] - c[j])/(3*h[j])   
    x = Symbol('x')
    for j in range(0,  n - 1):
        if factor_it == 'no':
            s_jx = round(list_fx[j],6) +  (round(b[j],6)*(x - list_x[j]))             + round(c[j],6)*(x - list_x[j])**2 + round(d[j],6)*(x - list_x[j])**3
            print('S _',j,'(x) =', s_jx, ' --> for ['                  , list_x[j],',', list_x[j+1],']')
            print('\nCoefficients (ai, bi, ci, di):\n\n ',list_fx[j],' , ', b[j],' , ', c[j], ' , ', d[j],'\n''------------------------------------------'                  '--------------------------\n')
        if factor_it == 'yes':
            s_jx = round(list_fx[j],6) +  factor((round(b[j],6)*(x - list_x[j])))             + round(c[j],6)*(x - list_x[j])**2 + round(d[j],6)*(x - list_x[j])**3
            print('S _',j,'(x) =', s_jx, ' --> for ['                  , list_x[j],',', list_x[j+1],']')
            print('\nCoefficients (ai, bi, ci, di):\n\n ',list_fx[j],' , ', b[j],' , ', c[j], ' , ', d[j],'\n''------------------------------------------'                  '--------------------------\n')
    return

# Neville's Method

def nevillesMethod(x, list_x, list_fx, Q_table = None, individual = 'no',notable = 'no'):
    n = np.size(list_x) - 1; 
    if (Q_table == None):
        Q_table = np.zeros((n + 1, n + 1));
        
    for i in range(0,n + 1):
        Q_table[i][0] = list_fx[i];
   
    for i in range(1, n + 1):
        for j in np.arange(1, i + 1):
            Q_table[i][j] = 0.0
            Q_table[i][j] += (((x - list_x[i - j])*Q_table[i][j - 1]- (x - list_x[i])*(Q_table[i - 1][j - 1]))/(list_x[i] - list_x[i - j]))
            if individual == 'yes' and i == j:
                print('Q_(',i,',',j,') =',Q_table[i][j])
                
    if notable == 'yes':
        return
    return Q_table[n][n], Q_table; 


# nodes = (x_0, x_1, ... x_n)
# function = function that we are approximating
# x value we use to approximate and check 
# degree specification of approximating polynomial.
def lagrange_polynomial(nodes, function, x, degree_specification):
    # nth_degree = len(nodes) - 1
    whole_polynomial = 0
    nodes = nodes[:degree_specification + 1]
    for coefficients in range(0, len(nodes)):
        current_node = nodes[coefficients]
        rest_of_nodes = nodes[:coefficients] + nodes[coefficients + 1:]
        coeffcient_polynomial_numerator = 1
        coeffcient_polynomial_denominator = 1
        for remaining_nodes in rest_of_nodes:
            coeffcient_polynomial_numerator = coeffcient_polynomial_numerator * (x - remaining_nodes)
        for remaining_nodes2 in rest_of_nodes:
            coeffcient_polynomial_denominator = coeffcient_polynomial_denominator * (current_node - remaining_nodes2)
        coeffient_polynomial = function(current_node) * (coeffcient_polynomial_numerator / coeffcient_polynomial_denominator)
        whole_polynomial += coeffient_polynomial
    return whole_polynomial



