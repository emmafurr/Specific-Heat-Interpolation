import numpy as np
import matplotlib.pyplot as plt
from final_project_script import nevillesMethod


def nevillesMethod(x, list_x, list_fx, Q_table = None, individual = 'no',notable = 'no'):
    n = np.size(list_x) - 1; 
    if (Q_table == None):
        Q_table = np.zeros((n + 1, n + 1));

    for i in range(0,n + 1):
        Q_table[i][0] = list_fx[i]

    for i in range(1, n + 1):
        for j in np.arange(1, i + 1):
            Q_table[i][j] = 0.0
            Q_table[i][j] += (((x - list_x[i - j])*Q_table[i][j - 1]- (x - list_x[i])*(Q_table[i - 1][j - 1]))/(list_x[i] - list_x[i - j]))
    if individual == 'yes' and i == j:
        print('Q_(',i,',',j,') =',Q_table[i][j])

    if notable == 'yes':
        return
    return Q_table[n][n]

# Nodes 1
list_x = [20,40,50,80,90]
list_fx = [4182,4179,4182,4198,4208]

# Nodes 2
list_x2 = [10,30,50,70,90]
list_fx2 = [4192, 4178,4182, 4191, 4208]


# temperature_x = np.arange(0, 101)
# actual_specific_heat = 1000 * np.array([4.217, 4.213, 4.210, 4.207, 4.205, 4.202, 4.200, 4.198, 4.196,
#     4.194, 4.192, 4.191, 4.189, 4.188, 4.187, 4.186, 4.185, 4.184, 4.183, 4.182, 4.182, 4.181, 4.181, 4.180, 4.180, 4.180, 4.179, 4.179, 4.179,
#       4.179, 4.178, 4.178, 4.178, 4.178, 4.178, 4.178, 4.178, 4.178, 4.178, 4.179, 4.179, 4.179, 4.179, 4.179,
#       4.179, 4.180, 4.180, 4.180, 4.180, 4.181, 4.181, 4.181, 4.182, 4.182, 4.182, 4.183, 4.183, 4.183, 4.184,
#       4.184, 4.185, 4.185, 4.186, 4.186, 4.187, 4.187, 4.188, 4.188, 4.189, 4.189, 4.190, 4.190, 4.191,
#       4.192, 4.192, 4.193, 4.194, 4.194, 4.195, 4.196, 4.196, 4.197, 4.198, 4.199, 4.200, 4.200, 4.201, 4.202,
#       4.203, 4.204, 4.205, 4.206, 4.207, 4.208, 4.209, 4.210, 4.211, 4.212, 4.213, 4.214, 4.216])

temperature_x = np.arange(0, 91, 5)
actual_specific_heat = 1000 * np.array([4.217, 4.202, 4.192, 4.1855, 4.182, 4.180, 4.178, 4.178, 4.179,
    4.181, 4.182, 4.183, 4.185, 4.188, 4.191, 4.194, 4.198, 4.203, 4.208])


# Nodes 1
# Cubic Spline
# For [20, 40]
def cubic_spline_1_0(x):
    return -6.01656*(0.05*x - 1.0) + 0.000377*(x - 20)**3 + 4182


# For [40, 50]
def cubic_spline_1_1(x):
    return 6.06624*(0.025*x - 1.0) - 0.000779*(x - 40)**3 + 0.022624*(x - 40)**2 + 4179


# For [50, 80]
def cubic_spline_1_2(x):
    return 18.52225*(0.02*x - 1.0) + 0.000206*(x - 50)**3 - 0.000745*(x - 50)**2 + 4182


# For [80, 90]
def cubic_spline_1_3(x):
    return 70.5176*(0.0125*x - 1.0) - 0.000593*(x - 80)**3 + 0.01778*(x - 80)**2 + 4198


cubic_spline_y_1_list = []
for i in range(20, 40):
    cubic_spline_y_1_list.append(cubic_spline_1_0(i))
for i in range(40, 50):
    cubic_spline_y_1_list.append(cubic_spline_1_1(i))
for i in range(50, 80):
    cubic_spline_y_1_list.append(cubic_spline_1_2(i))
for i in range(80, 90):
    cubic_spline_y_1_list.append(cubic_spline_1_3(i))

nevillesMethod_y_1_list = []
for i in range(0, 91):
    nevillesMethod_y_1_list.append(nevillesMethod(i, list_x, list_fx))


cubic_spline_x_1_list = np.arange(20, 90)


# Lagrange Polynomial
def lagrange_polynomial_1(x):
    return 697*(x - 90)*(x - 80)*(x - 50)*(x - 40)/420000 - 4179*(x - 90)*(x - 80)*(x - 50)*(x - 20)/400000 + 697*(x - 90)*(x - 80)*(x - 40)*(x - 20)/60000 - 2099*(x - 90)*(x - 50)*(x - 40)*(x - 20)/360000 + 263*(x - 80)*(x - 50)*(x - 40)*(x - 20)/87500


lagrange_polynomial_1_values = lagrange_polynomial_1(np.arange(0, 91))

# Plotting Functions
plt.ylim(4170, 4230)
plt.plot(temperature_x, actual_specific_heat, label="Actual")
plt.plot(np.arange(0, 91), lagrange_polynomial_1_values, label="Lagrange Polynomial")
plt.plot(cubic_spline_x_1_list, cubic_spline_y_1_list, label="Cubic Spline")
plt.plot(np.arange(0, 91), nevillesMethod_y_1_list, label="Neville's Method")
plt.title("Using Data From Table 1")
plt.xlabel("Temperature(°C)")
plt.ylabel("Specific Heat (J / kg-°C)")
<<<<<<< HEAD:plots.py
plt.legend(), plt.grid(True), plt.show()
=======
plt.legend()
plt.draw()
plt.show()
>>>>>>> 9a7ad9b2f2c8e9b3b505f93edac2aad6f72373b4:Final_Project_Plots_Yativ_Quackenbush_Furr_Kalpakian_Opp_Soltani.py



# Nodes 2
# Cubic Spline
# For [10, 30]
def cubic_spline_2_0(x):
    return -9.30357*(0.1*x - 1.0) + 0.000576*(x - 10)**3 + 4192


# For [30, 50]
def cubic_spline_2_1(x):
    return -7.17858*(0.0333333333333333*x - 1.0) - 0.000629*(x - 30)**3 + 0.034554*(x - 30)**2 + 4178


# For [50, 70]
def cubic_spline_2_2(x):
    return 19.375*(0.02*x - 1.0) + 0.000317*(x - 50)**3 - 0.003214*(x - 50)**2 + 4182


# For [70, 90]
def cubic_spline_2_3(x):
    return 44.75002*(0.0142857142857143*x - 1.0) - 0.000263*(x - 70)**3 + 0.015804*(x - 70)**2 + 4191


cubic_spline_y_2_list = []
for i in range(10, 30):
    cubic_spline_y_2_list.append(cubic_spline_2_0(i))
for i in range(30, 50):
    cubic_spline_y_2_list.append(cubic_spline_2_1(i))
for i in range(50, 70):
    cubic_spline_y_2_list.append(cubic_spline_2_2(i))
for i in range(70, 90):
    cubic_spline_y_2_list.append(cubic_spline_2_3(i))

cubic_spline_x_2_list = np.arange(10, 90)

nevillesMethod_y_2_list = []
for i in range(0, 91):
    nevillesMethod_y_2_list.append(nevillesMethod(i, list_x2, list_fx2))


# Lagrange Polynomial
def lagrange_polynomial_2(x):
    return 131*(x - 90)*(x - 70)*(x - 50)*(x - 30)/120000 - 2089*(x - 90)*(x - 70)*(x - 50)*(x - 10)/480000 + 2091*(x - 90)*(x - 70)*(x - 30)*(x - 10)/320000 - 1397*(x - 90)*(x - 50)*(x - 30)*(x - 10)/320000 + 263*(x - 70)*(x - 50)*(x - 30)*(x - 10)/240000


lagrange_polynomial_2_values = lagrange_polynomial_2(np.arange(0, 91))

# Plotting Functions
plt.ylim(4170, 4230)
plt.plot(temperature_x, actual_specific_heat, label="Actual")
plt.plot(np.arange(0, 91), lagrange_polynomial_2_values, label="Lagrange Polynomial")
plt.plot(cubic_spline_x_2_list, cubic_spline_y_2_list, label="Cubic Spline")
plt.plot(np.arange(0, 91), nevillesMethod_y_2_list, label="Neville's Method")
plt.title("Using Data From Table 2")
plt.xlabel("Temperature(°C)")
plt.ylabel("Specific Heat (J / kg-°C)")
