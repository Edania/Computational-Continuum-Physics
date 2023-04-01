import matplotlib.pyplot as plt
import numpy as np
import sys
#sys.path.insert(0, 'C:/Users/elsa2/Documents/Master/Computational Continuum Physics/test/build/Release')
sys.path.insert(0,'C:/Users/Elsa3/Documents/Master/Computational-Continuum-Physics/test/build/Debug')
import testing

c = 0.1
omega_0 = np.sqrt(2)
omega_1 = np.sqrt(3)
h = 0.005
#tau = 0.0005
#tau = c*h**2/2
tau = 0.00005

time = 10
length = 1
result = testing.FTCS_heat(int(10/tau),int(1/h), tau, h, c, omega_0, omega_1)

#np.savetxt('result.txt', result)

nan = np.isnan(result).any()
infty = np.isinf(result).any()
neg_infty = np.isneginf(result).any() 

x = np.linspace(0,1,int(1/h))
t = np.linspace(0,10,int(10/tau))
if(nan == False and infty == False and neg_infty == False):
    #residue at given point
    x_p = 0.3
    t_p = 10
    idx_x = int(x_p/h)
    idx_t = -1#int(t_p/tau)
    residue = -(result[idx_t,idx_x] - result[idx_t-1,idx_x])/tau + c*(result[idx_t,idx_x+1] - 2*result[idx_t,idx_x] + result[idx_t,idx_x-1])/h**2
    print(f'The residue at point ({x_p}, {t_p}) is: {residue}. The value of the point is {result[idx_t, idx_x]}')

    fig, ax = plt.subplots()
    cont = ax.contour(x,t,result)
    ax.clabel(cont)
    ax.set_title('numerical')
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    fig.savefig('lab1.pdf')

    fig, ax = plt.subplots()
    ax.contourf(x,t,result)
    ax.set_title('numerical')
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    fig.savefig('lab1_f.pdf')
    ################################
    #Testing ic and bc

    ic = np.sin(np.pi*x)**2
    bc_0 = -np.cos(omega_0*t)/omega_0+1/omega_0
    bc_1 = -np.cos(omega_1*t)/omega_1+1/omega_1
    #bc_0 = np.sin(omega_0*t)
    #bc_1 = np.sin(omega_1*t)
    #print(np.gradient(result[:,0]))
    fig, axs = plt.subplots(3)
    axs[0].plot(x, ic, label = 'inital conditions')
    axs[0].plot(x, result[0,:],'--', label = 'numerical')
    axs[0].legend()
    axs[0].set_title('Initial conditions')
    axs[1].plot(t, bc_0,label = 'bc')
    axs[1].plot(t, result[:,0], '--',label = 'numerical')
    axs[1].legend()
    axs[1].set_title('Boundary condition one')
    axs[2].plot(t, bc_1,label = 'bc')
    axs[2].plot(t, result[:,-1], '--',label = 'numerical')
    axs[2].legend()
    axs[2].set_title('Boundary condition two')
    fig.tight_layout()
    fig.savefig('lab1_conditions.pdf')
else:
    print(f'Solution unstable for tau = {tau} and h = {h}')