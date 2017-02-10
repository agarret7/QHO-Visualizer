import numpy as np
import scipy.integrate as integrate
from scipy.special import hermite
from scipy.misc import factorial
from scipy.misc import derivative
from matplotlib import pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
import random as rand
import itertools

# Quantum mechanics simulator for (initially) finite dimensional
# approximations of Hilbert spaces provided a basis.

'''
Code plan and structure:

One-dimensional:
1. Given an arbitrary Hamiltonian H, find the eigenvalues E_n.
2. Estimate the matrix form of the Hamiltonian in a finite-dimensional approximation in the QHO basis.
  a. Determine the numerical solutions to the differential equations H.psi = E_n.psi. (Call it numSolution[n])
  b. Perform functional optimization over the parameter space by finding the set of coefficients,
     [c_0, c_1, ..., c_n] such that the value
     
     QHOBasisApprox = numSolution[n](x) - sum([c[n] * hilbert.eigenbasis(n, x) for n in range(hilbert.dim)])

     assumes its minimum value.
  c. Take the basis approximation and form an approximation of the Hamiltonian,

     H = np.empty(2*[len(numSolution])
     for i, j in itertools.product(range(hilbert.dim), range(hilbert.dim)):
         H[i][j] = lambda x: np.conj(QHOBasisApprox[i](x)) * QHOBasisApprox[j](x)
  d. Set the unitary time operator to be

     U = lambda x, t: np.exp(-i/hbar * H(x) * t)
  e. Evolve the wave function to be
     psi = U(x, t) * initWaveFunc
3. 


Three-dimensional:
'''

class FunctionSampler:
    # Secret occult tricks to make computing our basis functions cheap.

    def __init__(self, f, minX, maxX, numSamples):
        self.minX = minX
        self.maxX = maxX
        self.range = maxX - minX
        self.numSamples = numSamples
        
        self.domain = np.linspace(minX, maxX, numSamples)
        self.image = [f(x) for x in self.domain]

    def __call__(self, x):
        x = min(max(x, self.minX), self.maxX - 1)
        x -= self.minX
        x = round(x * (self.numSamples / self.range))

        return self.image[int(x)]

class HilbertSpace:

    def __init__(self, dim, basis = 'QHO', hamiltonianPotential = None, approxMethod = 'QHO'):
        self.dim = dim

        if hamiltonianPotential != None:
            self.V = hamiltonianPotential

        if basis == 'QHO' or approxMethod == 'QHO':
            self.QHOBasis = lambda n, x: (1/np.sqrt(2**n * factorial(n)) * np.pi**(1/4)) * np.exp(-x**2/2) * hermite(n)(x)
            QHOBasisApprox = []
            for n in range(self.dim):
                QHOBasisApprox.append(FunctionSampler(lambda x: self.QHOBasis(n,x),-15,15,2000))

            self.QHOBasis = lambda n, x: QHOBasisApprox[n](x)

        elif basis == 'Square Well' or approxMethod == 'Square Well':
            self.SWBasis = lambda n, x: np.sqrt(2) * np.sin((n+1)*np.pi*(x + 1/2)) if abs(x) < 1/2 else 0

        if approxMethod == 'QHO':
            self.H = np.empty([self.dim, self.dim])
            for i, j in itertools.product(range(self.dim), range(self.dim)):
                HFunc = lambda x: np.conj(self.QHOBasis(i, x)) * (-1/2 * derivative(lambda xPrime: self.QHOBasis(2, xPrime), x, dx=1e-6, n=2) + self.V(x) * self.QHOBasis(2, x))
                self.H[i][j] = integrate.quad(HFunc, -np.inf, np.inf)[0]
                
        if basis == 'QHO':
            self.eigenbasis = self.QHOBasis
            self.eigenvalues = lambda n: 1/2 + n
            
        elif basis == 'Square Well':
            self.eigenbasis = self.SWBasis
            self.eigenvalues = lambda n: (n+1)**2 * np.pi**2
            
class WaveFunction:

    def __init__(self, hilbertSpace, initWaveFunc = None, coeff = None):
        self.hilbert = hilbertSpace

        if coeff != None:
            self.coeff = coeff
                
        else:
            self.evaluate = lambda x: initWaveFunc(x)
            self.coeff = self.orthogonalBasisProjection(initWaveFunc)

        evaluate = lambda x, t: sum([self.coeff[n] * self.hilbert.eigenbasis(n, x) * self.phaseFactor(n, t) for n in range(self.hilbert.dim)])
        normF = self.normalize(lambda x: evaluate(x, 0))
        
        self.evaluate = lambda x, t: evaluate(x, t) / normF
        
    def phaseFactor(self, n, t):
        return np.exp(-1j * self.hilbert.eigenvalues(n) * t)

    def normalize(self, waveFunc):
        return np.sqrt(integrate.quad(lambda x: np.absolute(waveFunc(x))**2, -np.inf, np.inf)[0])

    def orthogonalBasisProjection(self, waveFunc):
        coeff = []
        for n in range(hs.dim):
            coeff.append(integrate.quad(lambda x: np.conj(self.hilbert.eigenbasis(n, x)) * waveFunc(x), -np.inf, np.inf)[0])

        return coeff

class Plotter:

    def __init__(self):
        pass

    def plotWaveFunction2d(self, waveFunction, samples = 50, frames = 250, timeFactor = 1, saveName = None):

        fig = plt.figure()
        ax1 = plt.axes(xlim=(-5, 5), ylim=(-1,1))
        line, = ax1.plot([], [], lw=2)
        plt.xlabel('$x$')
        plt.ylabel('$\psi(x)$')
        plt.xticks([],[])
        plt.yticks([],[])

        plotlays, plotcols = [2], ["blue","red","green"]
        lines = []
        for index in range(3):
            lobj = ax1.plot([],[],lw=2,color=plotcols[index])[0]
            lines.append(lobj)

        def init():
            for line in lines:
                line.set_data([], [])
            return lines

        def animate(i):
            posValues = np.linspace(-5, 5,samples)
            amplitudeValues = [waveFunction.evaluate(x, timeFactor * i/50) for x in posValues]

            reValues = np.real(amplitudeValues)
            imValues = np.imag(amplitudeValues)
            probValues = np.absolute(amplitudeValues)**2

            ylist = [reValues, imValues, probValues]

            for lnum, line in enumerate(lines):
                line.set_data(posValues, ylist[lnum])
                
            return lines

        anim = animation.FuncAnimation(fig, animate, init_func = init, frames = frames, interval = 10, blit = True)

        if saveName != None:
            anim.save(saveName + '.mp4', fps = 30, extra_args = ['-vcodec', 'libx264'])

        plt.show()

    def plotWaveFunction3d(self, waveFunction, samples = 50, frames = 250, timeFactor = 1, saveName = None):

        fig = plt.figure()
        ax1 = fig.add_axes([0,0,1,1], projection = '3d')
        ax1.set_xlabel('$x$')
        ax1.set_ylabel('Im$(\psi)$')
        ax1.set_zlabel('Re$(\psi)$')
        
        ax1.set_xlim((-5, 5))
        ax1.set_ylim((-1,1))
        ax1.set_zlim((-1,1))

        ax1.set_xticks([])
        ax1.set_yticks([])
        ax1.set_zticks([])

        plotlays, plotcols = [2], ["blue","red","green"]
        
        lines = []
        pts = []
        
        for index in [0,2]:
            lobj = ax1.plot([],[], [],lw=2,color=plotcols[index])[0]
            pobj = ax1.plot([],[], [],lw=2,color=plotcols[index])[0]
            lines.append(lobj)
            pts.append(pobj)

        def init():
            for line, pt in zip(lines, pts):
                line.set_data([], [])
                line.set_3d_properties([])

            return lines + pts

        def frame_path(i):
            i+=500
            if i < 680:
                return 0, 0
            #elif 500 <= i < 680:
                #new_i = i - 500
                #return 0, -90 + new_i/2
            elif 680 <= i < 860:
                new_i = i - 680
                return new_i/6, new_i/3
            else:
                return 30, 60

        def frame_path(i):
            if i < 500:
                return 0, -90
            elif 500 <= i < 590:
                new_i = i - 500
                return new_i/3, -90 + new_i/3
            else:
                return 30, -60

        def animate(i):
            elev, azim = frame_path(i)
            ax1.view_init(elev, azim)
            
            posValues = np.linspace(-5, 5, samples)
            amplitudeValues = [waveFunction.evaluate(x, timeFactor * i/50) for x in posValues]

            reValues = np.real(amplitudeValues)
            imValues = np.imag(amplitudeValues)
            probValues = np.absolute(amplitudeValues)**2

            ylist = [(imValues, reValues), (0, probValues)]

            for n, line in enumerate(lines):
                line.set_data(posValues, ylist[n][0])
                line.set_3d_properties(ylist[n][1])

            fig.canvas.draw()
            return lines

        anim = animation.FuncAnimation(fig, animate, init_func = init, frames = frames, interval = 10, blit = True)

        if saveName != None:
            anim.save('./video/animations/' + saveName + '.mp4', fps = 30, writer = 'avconv', codec = 'libx264')

        plt.show()

p = Plotter()
hs = HilbertSpace(dim = 2, hamiltonianPotential = lambda x: 1/2*x**2, basis = 'QHO', approxMethod = None)
psi = WaveFunction(hs, coeff = [1,1])
#psi = WaveFunction(hs, lambda x: 1 if abs(x) < 5 else 0)
print('c_n = ' + str(psi.coeff))

p.plotWaveFunction3d(psi, samples = 300, frames = 1000, timeFactor = 2, saveName = None)

'''
Ideas:
- Take a state vetor and evolve it using the unitary evolution operator defined as the
  exponentiated Hamiltonian. Recall that the Hamiltonian defines the energy eigenbasis.
  If you have an eigenbasis and eigenvalues, simply define the Hamiltonian according to
  the spectral theorem. If you instead have a Hamiltonian first, we have to decompose it
  into a finite-dimensional approximation.
  One method to do this may be to calculate the first n eigenstates and eigenvalues of
  the Hamiltonian, and then to use spectral theorem to simply define the Hamiltonian.

- Long term goals may include determining the unitary evolution of a state in a time-varying potential,
  determining the projection onto the eigenstate of an observable O (to allow for measurements, which can
  be determined by sampling our probability distribution in the corresponding observable picture)

- Mapping the position wave function to the corresponding momentum picture using Fourier transforms, and
  other methods?

- Checking Vachaspati's hypothesis that the expectation value of a wave function is related to the RMS
  value of a classical wave??? Not sure what he was thinking, perhaps ask JJ.

- Coherent states.
'''
