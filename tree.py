from scipy import spatial
import numpy as np

class OneDTree:
    def __init__(self,x):
        dummy = np.zeros_like(x)
        self.tree = spatial.KDTree(zip(x,dummy))

    def Query(self,pts,k):
        return self.tree.query((pts,0),k)[1]
            
def main():
    x = np.linspace(0,1,101)        
    pt = x[5:6]

    mytree = OneDTree(x)
    index = mytree.Query(pt,6)
    print pt
    print index
    print x[index]

if __name__ == '__main__':
    main()
