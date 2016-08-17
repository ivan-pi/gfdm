import numpy as np
import matplotlib.pyplot as plt

class Monomial:
    def __init__(self,order):
        self.order = order
        self.size = order + 1
        
    def __call__(self,x):
        return [x**i for i in range(self.order+1)]
        
    def der(self,x,d):
        if d == 0:
            p = [1.,x,x**2,x**3,x**4,x**5]
        if d == 1:
            p = [0.,1.,2.*x,3.*x**2,4.*x**3,5.*x**4]
        if d == 2:
            p = [0.,0.,2.,6.*x,12.*x**2,20.*x**3]
        if d == 3:
            p = [0.,0.,0.,6.,24.*x,60.*x**2]
        if d == 4:
            p = [0.,0.,0.,0.,24.,120.*x]
        
        return p[0:self.order+1]
                    
class ShapeFunction:
    def __init__(self,basis,support,f=None):
        self.basis = basis
        self.x = support-support[0]
        self.size = len(self.x)
        self.M = self.ShapeFunctionMatrix()
        self.P = np.linalg.pinv(self.M)
        self.f = f
        if self.f is not None:
            self.a = np.linalg.lstsq(self.M,self.f)[0]
            
    def SetValues(self,f):
        self.f = f
        self.a = np.linalg.lstsq(self.M,self.f)[0]
        
    def ShapeFunctionMatrix(self):
        M = np.zeros((self.size,self.basis.size))
        for i in range(self.size):
            M[i,:] = self.basis(self.x[i])
        return M
        
    def __call__(self,x):    
        assert self.f is not None, 'function values are not set'
        xval = self.basis(x)
        return np.sum(self.a*xval)

    def Derivative(self,order,x=None):
        if x is None:
            x = self.x[0]
        xval = self.basis.der(x,order)
        return np.dot(self.P.T,xval)

def main():
    
    basis = Monomial(2)
    x = np.array([0,1,4,6,7])
    f = np.array([3,1,2,4,5])
    
    sf = ShapeFunction(basis,x)
    sf.SetValues(f)
    
    xx = np.linspace(-1,9)
    ff = [sf(i) for i in xx]
    fder = [2.*sf.a[2] for i in xx]
    
    print np.dot(sf.P,f)
    
    der = [np.dot(sf.Derivative(2,xval),f) for xval in x] 
    print der
    
    plt.plot(x,f,'o')
    plt.plot(xx,ff,'-')
    plt.plot(0,np.dot(sf.P[0,:],f),'x')
    plt.plot(x,der,'o')
    plt.plot(xx,fder)
    plt.show()
    
if __name__ == '__main__':
    main() 
