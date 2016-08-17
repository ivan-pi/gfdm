import numpy as np
import matplotlib.pyplot as plt
import shape_functions as sf
from tree import OneDTree

class ODESolver:
    def __init__(self,ode,grid,basis,ssize):
        self.ode = ode
        self.grid = grid
        self.tree = OneDTree(self.grid.x)
        self.basis = basis
        assert ssize >= self.basis.order+1
        self.ssize = ssize

        self.A = None
        self.b = None
        self.u = None
        self.fileName = None
        
    def Solve(self):
        self.PopulateMatrix()
        self.ApplyBoundaryConditions()
        self.PopulateVector()
        
        print self.A
        print self.b
        
        self.u = np.linalg.solve(self.A,self.b)

    def SetFileName(self,fileName):
        self.fileName = fileName 

    def WriteSolutionFile(self):
        if self.fileName is not None:
            results = np.vstack([self.grid.x,self.u]).T
            np.savetxt(self.fileName,results,delimiter=' ')

    def FindSupport(self,i):
        return self.tree.Query(self.grid.x[i],self.ssize)
         
    def PopulateMatrix(self):
        n = self.grid.n
        A = np.zeros((n,n))
        
        for i in range(1,n-1):
            A[i,i] += self.ode.coeff[0]
            
            #only temporary
            sindex = self.FindSupport(i)
            print sindex
            support = self.grid.x[sindex]
            trialf = sf.ShapeFunction(self.basis,support)
            
            operator = np.zeros(3)
            for n in range(1,self.ode.order):
                operator += self.ode.coeff[n]*trialf.Derivative(n+1)
            
            for j,val in enumerate(sindex):
                A[i,val] += operator[j]
        
        self.A = A
                
    def PopulateVector(self):
        n = self.grid.n
        
        b = np.zeros(n)
        b[0] = self.ode.LeftBcVal
        b[1:n-2] = self.ode.rhs(self.grid.x[1:n-2])
        b[n-1] = self.ode.RightBcVal
        
        self.b = b
        
    def ApplyBoundaryConditions(self):
        
        if self.ode.LeftBcIsDirichlet:
            self.A[0,0] = 1.0 
        if self.ode.RightBcIsDirichlet:
            n = self.grid.n
            self.A[n-1,n-1] = 1.0
        if self.ode.LeftBcIsNeumann:
            sindex = list(range(self.ssize))
            support = self.grid.x[sindex]
            trialf = sf.ShapeFunction(self.basis,support)
            self.A[0,sindex] = trialf.Derivative(1)
        if self.ode.RightBcIsNeumann:
            n = self.grid.n
            sindex = list(range(n-1,n-1-self.ssize,-1))
            print sindex
            support = self.grid.x[sindex]
            trialf = sf.ShapeFunction(self.basis,support)
            
            self.A[n-1,sindex] = trialf.Derivative(1)
class Grid:
    def __init__(self,a=0.0,b=1.0,n=10):
        self.a = a
        self.b = b
        self.n = n
        self.dx = (self.b-self.a)/float(self.n-1)
        self.x = None
        self.FillUniform()

    def FillUniform(self):
        if self.x is None:
            self.x = np.linspace(self.a,self.b,self.n)
    
    def FillRandom(self):
        inner_nodes = np.sort(np.random.rand(self.n-2))
        x = np.hstack([self.a,inner_nodes,self.b])
        self.x = x
        
    def FillIrregular(self,irreg):
        assert irreg < self.dx
        shift = 2.*(np.random.rand(self.n-2) - 0.5)*irreg
        self.x[1:self.n-1] += shift
        
class Ode:
    def __init__(self,order,coeff,rhs):
        self.order = order
        
        coeff.reverse()
        self.coeff = coeff
        
        self.rhs = rhs
        
        self.LeftBcIsDirichlet = False
        self.RightBcIsDirichlet = False
        self.LeftBcIsNeumann = False
        self.RightBcIsNeumann = False
        
        self.LeftBcVal = None
        self.RightBcVal = None
        
    def SetLeftDirichletBC(self,val):
        self.LeftBcVal = val
        self.LeftBcIsDirichlet = True
        
    def SetRightDirichletBC(self,val):
        self.RightBcVal = val
        self.RightBcIsDirichlet = True
        
    def SetLeftNeumannBC(self,val):
        self.LeftBcVal = val
        self.LeftBcIsNeumann = True
        
    def SetRightNeumannBC(self,val):
        self.RightBcVal = val
        self.RightBcIsNeumann = True
        
def main():

    # Specify Grid
    n = 6
    grid = Grid(0.0,1.0,n)
    grid.FillIrregular(0.01)
    print 'grid = ', grid.x
    dx = 1.0/(n-1)
    
    Th = 1.0
    
    A=np.zeros((n,n))
    A[0,0] = 1
    for i in range(1,n-1):
        A[i,i-1] = 1./dx**2
        A[i,i] = -2./dx**2-Th**2
        A[i,i+1] = 1/dx**2
    A[n-1,n-3:n] = np.array([0.5,-2,1.5])/(dx)
    
    print A
    
    b = np.zeros(n)
    b[0] = 1.0
    
    c = np.linalg.solve(A,b)
    
        
    # Specify linear ODE and RHS
    def rhs(x):
        return 0.0
    
    ode = Ode(2,[1.0,-Th**2],rhs)
    
    ode.SetLeftDirichletBC(1.0)
    ode.SetRightNeumannBC(0.0)
    
    basis = sf.Monomial(2)

    system = ODESolver(ode,grid,basis,3)
    system.SetFileName('solution.dat')
    
    system.Solve()
    system.WriteSolutionFile()
    
    x = grid.x
    u = system.u
    
    plt.plot(x,u,'-o',label='numerical')
    plt.plot(np.linspace(0,1,n),c,'x-',label='uniform grid')
    plt.xlabel('x')
    plt.ylabel('u(x)')
    plt.legend()
    plt.show()
    


if __name__ == '__main__':
    main()
