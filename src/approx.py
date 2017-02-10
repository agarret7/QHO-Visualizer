import numpy as np
from scipy.special import hermite
from scipy.misc import derivative
from scipy.misc import factorial
from scipy.linalg import expm
import itertools
from matplotlib import animation
from matplotlib import pyplot as plt
import scipy.integrate as integrate

QHOBasis = lambda n, x: 1/(np.sqrt(factorial(n) * 2**n) * np.pi**(1/4)) * np.exp(-x**2 / 2) * hermite(n)(x)
dim = 3

V = lambda x: 100 if abs(x) > 2 else 0

H = np.empty([dim, dim])
for i, j in itertools.product(range(dim), range(dim)):
    func = lambda x: np.conj(QHOBasis(i, x)) * (-1/2 * derivative(lambda xPrime: QHOBasis(j, xPrime), x, dx=1e-6, n=2) + V(x) * QHOBasis(j, x))
    H[i][j] = integrate.quad(func, -np.inf, np.inf)[0]

print(integrate.quad(lambda x: np.conj(QHOBasis(0, x)) * (-1/2 * derivative(lambda xPrime: QHOBasis(2, xPrime), x, dx=1e-6, n=2) + V(x) * QHOBasis(2, x)), -np.inf, np.inf))
print(H)
U = lambda t: expm(-1j * H * t)

fig = plt.figure()
ax1 = plt.axes(xlim=(-5, 5), ylim=(-2,2))
line, = ax1.plot([], [], lw=2)
plt.xlabel('$x$')
plt.ylabel('$\psi(x)$')
plt.xticks([],[])
plt.yticks([],[])

plotlays, plotcols = [2], ["blue","red","green"]
lines = []

initWaveState = np.zeros((dim,1))
initWaveState[0][0] = 1

newWaveState = lambda t: np.dot(U(t), initWaveState)
newWaveFunction = lambda x, t: sum(newWaveState(t)[n][0] * QHOBasis(n, x) for n in range(dim))

for index in range(3):
    lobj = ax1.plot([],[],lw=2,color=plotcols[index])[0]
    lines.append(lobj)

def init():
    for line in lines:
        line.set_data([], [])
    return lines

def animate(i):
    posValues = np.linspace(-5,5,50)
    amplitudeValues = [newWaveFunction(x, i/5) for x in posValues]

    reValues = np.real(amplitudeValues)
    imValues = np.imag(amplitudeValues)
    probValues = np.absolute(amplitudeValues)**2

    ylist = [reValues, imValues, probValues]

    for lnum, line in enumerate(lines):
        line.set_data(posValues, ylist[lnum])
        
    return lines

anim = animation.FuncAnimation(fig, animate, init_func = init, frames = 200, interval = 30, blit = True)

plt.show()
