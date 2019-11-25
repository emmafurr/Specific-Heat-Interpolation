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
from sympy import Symbol, poly, factor, expand


def main():

# We use x = 24, x = 56, x = 83
    
# Our data set:
# -----------------------------------------------------------------------
    list_x = [22,42,52,82,100]
    list_fx = [4181,4179,4186,4199,4217]

    list_x2 = [10,30,50,70,90]
    list_fx2 = [4192, 4178,4182, 4191, 4208]
    x = 25
    fx = 4180

    z = 85
    fz = 4203

# -----------------------------------------------------------------------

# ======================== x = 25 degrees Celcius ======================== 
    print('\nWe will use now Cubic Spline Interoplation with our data set:\n')
    splineConstructor(list_x2, list_fx2, factor_it = 'yes')
    print('We evaluate the necessary spline function to find approximation.\n')
    def c1(x):
        #return -8.680364*(0.0454545454545454*x - 1.0) + 0.000736*(x - 22)**3 + 4181
        return -9.30357*(0.1*x - 1.0) + 0.000576*(x - 10)**3 + 4192
    
    spline_approx = c1(25)
    print('FINAL APPROXIMATION:', spline_approx)
    print('\nWe find the relative error:\n')
    s1_error = relative_error(4180, spline_approx)
    print('Relative error for Spline: \n', s1_error)
    
# -----------------------------------------------------------------------
    print('\n\n\n\nWe will use now Lagrange Interoplation with our data set:\n')
    abc = lagrange_polynomial(list_x2, list_fx2)
    print(abc)
    print('We evaluate the highest degree Lagrange to find approximation.\n')
    def g1(x):
#        return 4181*(x - 100)*(x - 82)*(x - 52)*(x - 42)/2808000 - 4179*(x - 100)*(x - 82)*(x - 52)*(x - 22)/464000 + 2093*(x - 100)*(x - 82)*(x - 42)\
#               *(x - 22)/216000 - 4199*(x - 100)*(x - 52)*(x - 42)*(x - 22)/1296000 \
#               + 4217*(x - 82)*(x - 52)*(x - 42)*(x - 22)/3908736
        return 131*(x - 90)*(x - 70)*(x - 50)*(x - 30)/120000 - 2089*(x - 90)*(x - 70)*(x - 50)\
            *(x - 10)/480000 + 2091*(x - 90)*(x - 70)*(x - 30)*(x - 10)/320000 - \
                1397*(x - 90)*(x - 50)*(x - 30)*(x - 10)/320000 + 263*(x - 70)*(x - 50)*\
                    (x - 30)*(x - 10)/240000
    lag_approx = g1(25)
    print('FINAL APPROXIMATION:', lag_approx)
    print('\nWe find the relative error:\n')
    l1_error = relative_error(4181, lag_approx)
    print('Relative error for Lagrange: ', l1_error)
    
# -----------------------------------------------------------------------
    print('\n\n\n\nWe will use now Nevilles method with a data set:')
    nev_approx,n_table = nevillesMethod(x, list_x2, list_fx2, Q_table = None\
                   , individual = 'yes',notable = 'no')
    print(n_table)
    print('\n\n FINAL APPROXIMATION: ', nev_approx)
    print('We consider the Q_4,4 th table entry - The highest degree approximation.\n') 
    print('\nWe find the relative error:\n')
    n1_error = relative_error(4180, nev_approx)
    print('Relative error for Neville: \n', n1_error)
# =========================================================================

# ======================== x = 75 degrees Celcius ========================
    print('\nWe will use now Cubic Spline Interoplation with our data set:\n')
    splineConstructor(list_x2, list_fx2, factor_it='yes')
    print('We evaluate the necessary spline function to find approximation.\n')

    def c(x):
        #return -8.680364*(0.0454545454545454*x - 1.0) + 0.000736*(x - 22)**3 + 4181
        return 44.75002*(0.0142857142857143*x - 1.0) -\
             0.000263*(x - 70)**3 + 0.015804*(x - 70)**2 + 4191

    spline_approx = c(85)
    print('FINAL APPROXIMATION:', spline_approx)
    print('\nWe find the relative error:\n')
    s1_error = relative_error(4203, spline_approx)
    print('Relative error for Spline: \n', s1_error)

# -----------------------------------------------------------------------
    print('\n\n\n\nWe will use now Lagrange Interoplation with our data set:\n')
    abc = lagrange_polynomial(list_x2, list_fx2)
    print(abc)
    print('We evaluate the highest degree Lagrange to find approximation.\n')

    def g(x):
        #        return 4181*(x - 100)*(x - 82)*(x - 52)*(x - 42)/2808000 - 4179*(x - 100)*(x - 82)*(x - 52)*(x - 22)/464000 + 2093*(x - 100)*(x - 82)*(x - 42)\
        #               *(x - 22)/216000 - 4199*(x - 100)*(x - 52)*(x - 42)*(x - 22)/1296000 \
        #               + 4217*(x - 82)*(x - 52)*(x - 42)*(x - 22)/3908736
        return 131*(x - 90)*(x - 70)*(x - 50)*(x - 30)/120000 - 2089*(x - 90)*(x - 70)*(x - 50)\
            * (x - 10)/480000 + 2091*(x - 90)*(x - 70)*(x - 30)*(x - 10)/320000 - \
            1397*(x - 90)*(x - 50)*(x - 30)*(x - 10)/320000 + 263*(x - 70)*(x - 50) *\
            (x - 30)*(x - 10)/240000
    lag_approx = g(85)
    print('FINAL APPROXIMATION:', lag_approx)
    print('\nWe find the relative error:\n')
    l1_error = relative_error(4203, lag_approx)
    print('Relative error for Lagrange: ', l1_error)

# -----------------------------------------------------------------------
    print('\n\n\n\nWe will use now Nevilles method with a data set:')
    nev_approx, n_table = nevillesMethod(
        z, list_x2, list_fx2, Q_table=None, individual='yes', notable='no')
    print(n_table)
    print('\n\n FINAL APPROXIMATION: ', nev_approx)
    print('We consider the Q_4,4 th table entry - The highest degree approximation.\n')
    print('\nWe find the relative error:\n')
    n1_error = relative_error(4203, nev_approx)
    print('Relative error for Neville: \n', n1_error)
# =========================================================================
#  Relative error function

def relative_error(actual, approx):
    rel_error = abs(actual - approx) / abs(actual)
    return rel_error

# Spline Interpolation Method

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
            print('S _',j,'(x) =', s_jx, ' for ['                  , list_x[j],',', list_x[j+1],']')
        if factor_it == 'yes':
            s_jx = round(list_fx[j],6) +  factor((round(b[j],6)*(x - list_x[j])))             + round(c[j],6)*(x - list_x[j])**2 + round(d[j],6)*(x - list_x[j])**3
            print('S _',j,'(x) =', s_jx, ' for ['                  , list_x[j],',', list_x[j+1],']')
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


# Lagrange Interpolation Method

def lagrange_polynomial(nodes, function_nodes):
    degree_specification = len(nodes) - 1
    x = Symbol("x")
    whole_polynomial = 0
    nodes = nodes[:degree_specification + 1]
    for coefficients in range(0, len(nodes)):
        current_node = nodes[coefficients]
        rest_of_nodes = nodes[:coefficients] + nodes[coefficients + 1:]
        coeffcient_polynomial_numerator = 1
        coeffcient_polynomial_denominator = 1
        for remaining_nodes in rest_of_nodes:
            coeffcient_polynomial_numerator = coeffcient_polynomial_numerator \
                                              * (x - remaining_nodes)
        for remaining_nodes2 in rest_of_nodes:
            coeffcient_polynomial_denominator = coeffcient_polynomial_denominator \
                                                * (current_node - remaining_nodes2)
        coeffient_polynomial = function_nodes[coefficients] \
                               * (coeffcient_polynomial_numerator \
                                  / coeffcient_polynomial_denominator)
        whole_polynomial += coeffient_polynomial
    return whole_polynomial

# Call the main function.
main()
