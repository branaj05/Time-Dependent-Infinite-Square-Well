"""
AUSTIN BRANDENBERGER
Python 3.7.7
14 November 2021
This program looks at the development of a particle in the infinite square well over time
442 Intro to Quantum Mechanics. 
"""
import numpy as np
import setuptools
import matplotlib.pyplot as plt
from matplotlib import animation
from scipy.integrate import simps

#scipy.integrate.simpson(y, x=None, dx=1.0, axis=- 1, even='avg')
alpha = 10#still needs calculated--needs to be specific. 
k0 = 100
aa = 1
a = -aa/2
b = aa/2
N = 1000
xx = np.linspace(a, b, N)
#x = np.linspace(a + aa/2, b + aa/2, N)
#x = np.linspace(-5*aa, 5*aa, N)



def psi0 (x):
    alpha = 10#still needs calculated--needs to be specific.
    shift = 0
    k0 = 1
    psi0 = (2*alpha/np.pi)**.25*np.exp(1j*k0*(x+shift)-alpha*(x+shift)**2)
    return np.real(psi0*np.conj(psi0))


def psi (x, t, n): #n is the upper bound of the sum
    #finding Psi_n
    def psi_n (x, t, n):
        #finding the coefs (Cn)
        shift = +aa/2 #running in values from -a/2 to a/2 so I have to shift
        #the values going through my psi(x) b/c that is the solution to the
        #infinite square well from 0 to a. 
        def Cn(x, n):

            psi0 = (2*alpha/np.pi)**.25*np.exp(1j*k0*(x)-alpha*(x)**2) #initial wave function
            psi_x = (2/aa)**.5 * np.sin(n*np.pi*(x+shift)/aa)#time independent solution 
            y_func = psi0*psi_x
            C_n = simps(y_func, x)#performs the integral to find the coeffs
            return C_n
    
        m = 9.1093835*10**-31 #mass of an electron
        hbar = 1.0545718*10**-34
        E = (n*np.pi*hbar)**2/(2*m*aa**2)
        psi_n = Cn(x, n)*(2/aa)**.5*np.sin(n*np.pi/aa*(x+shift))*np.exp(-1j*E*t/hbar) #psi_n, one of the infinite sum
        return psi_n  
    #finding psi_xt
    i = 1
    psi_xt = np.zeros_like(x) 
    while i < n:
        psi_xt = psi_xt + psi_n(x, t, i)#adding up psi_n's as I go through values from 1-n. 
        i += 1
    return np.real(psi_xt*np.conj(psi_xt)) #psi star psi


t = 0
n = 300
"""
while t <= 2:
    #n_energy = 1
    yy = psi(xx, t, n)
    #psp = yy*np.conj(yy)
    #plt.plot(xx, yy, 'k-', label = 'psi(x, t= 0)')
    plt.plot(xx, psi0(xx), 'r-.', label = 'psi0')
    #plt.plot(xx, psp)
    #plt.legend(loc = 'best')
    plt.title('Probability distribution at t = 0')
    plt.xlabel('Position')
    plt.ylabel('Probability')
    #plt.grid()
    plt.show()
    t += 10
"""
# animation function.  This is called sequentially
# initialization function: plot the background of each frame

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

# First set up the figure, the axis, and the plot element we want to animate

fig = plt.figure()
ax = plt.axes(xlim=(a, b), ylim=(0, 10))
line, = ax.plot([], [], lw=2)
# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,

# animation function.  This is called sequentially
def animate(i):
    x = np.linspace(a, b, N)
    y = psi(x, i, n)
    line.set_data(x, y)
    return line,

t0 = 0
tf = 1500
N = 200
dt = int((tf-t0)/N)
# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=np.arange(t0, tf, tf/N), interval=.5, blit=True)
plt.xlabel('position (x)')
plt.ylabel("Probability (psi*psi)")
"""
f = r"C://Users/austi/Desktop/InfiniteSquareWell.gif"
writergif = animation.PillowWriter(fps=30) 
anim.save(f, writer=writergif)
"""

plt.show()


