import numpy as np

class Lagrange():
    '''
    A class for Lagrange interpolation, integration, and derivation, using GLL (Gauss Lobatto Legendre) collocation points.
    
    gll funct.:            returning the GLL collocation points and their corresponding weigths.
    lagrange_poly func.:   interpolating using Lagrange polynomials at GLL
    integ func.:           integrating using Lagrange polynomials at GLL
    deriv func.            taking the derivative of Lagrange polynomials Legendre polynomials at GLL using
    legendre func.         returning the Legendre Polynomials at GLL
    
    
    The source of the following functions is:
            Prof. Heier Igel, 2020
            
            https://github.com/heinerigel/coursera/tree/master/Notebooks4Coursera/W9
    '''
    
    def __init__(self, N):
        self.N = N
        
        
    def gll(self):
        '''
            GLL collocation points and weights - This is for N = 2 to N = 7
        '''
        
        if self.N == 2:
            xi = [-1.0, 0.0, 1.0]
            wi = [0.33333333, 1.33333333, 0.33333333]
        elif self.N == 3:
            xi = [-1.0, -0.447213595499957, 0.447213595499957, 1.0]
            wi = [0.1666666667, 0.833333333, 0.833333333, 0.1666666666]
        elif self.N == 4:
            xi = [-1.0, -0.6546536707079772, 0.0, 0.6546536707079772, 1.0]
            wi = [0.1, 0.544444444, 0.711111111, 0.544444444, 0.1]
        elif self.N == 5:
            xi = [-1.0, -0.7650553239294647, -0.285231516480645, 0.285231516480645, 0.7650553239294647, 1.0]
            wi = [0.0666666666666667,  0.3784749562978470, 0.5548583770354862, 0.5548583770354862, 
                  0.3784749562978470, 0.0666666666666667]
        elif self.N == 6:
            xi = [-1.0, -0.8302238962785670, -0.4688487934707142, 0.0, 0.4688487934707142, 0.8302238962785670, 1.0]
            wi = [0.0476190476190476, 0.2768260473615659, 0.4317453812098627, 0.4876190476190476, 0.4317453812098627, 
                  0.2768260473615659, 0.0476190476190476]
        elif self.N == 7:
            xi = [-1.0, -0.8717401485096066, -0.5917001814331423, -0.2092992179024789, 0.2092992179024789, 
                  0.5917001814331423, 0.8717401485096066, 1.0]
            wi = [0.0357142857142857, 0.2107042271435061, 0.3411226924835044, 0.4124587946587038, 0.4124587946587038,
                  0.3411226924835044, 0.2107042271435061, 0.0357142857142857]
        else:
            raise NotImplementedError
                
        return np.array(xi), np.array(wi)

    
    def lagrange_poly(self, p, x):
        ''' 
            A function for calculating Lagrange polynomial p at x[-1, 1]
            
        '''
        [xi, wi] = self.gll()
        fac = 1
        for j in range(-1, self.N):
            if j != p:
                fac = fac * ((x - xi[j + 1]) / (xi[p + 1] - xi[j + 1]))
        return fac
    
                
    def integ(self):
        # A function for calculating the integration using Gauss quadrature at the collocation poitns
        [xi, wi] = self.gll()
        
        return sum(xi * wi)
    

    def deriv(self):
        ''' 
            A function for calculating the derivative of Lagrange polynomials 
            using Legendre polynomials at the collocation poitns
            
        '''
        
        [xi, wi] = self.gll()
        
        out = np.zeros([self.N + 1, self.N + 1])
        d = np.zeros([self.N + 1, self.N + 1])

        for i in range(-1, self.N):
            for j in range(-1, self.N):
                if i != j:
                    d[i + 1, j + 1] = self.legendre(xi[i + 1]) / \
                        self.legendre(xi[j + 1]) * 1.0 / (xi[i + 1] - xi[j + 1])
                if i == -1:
                    if j == -1:
                        d[i + 1, j + 1] = -1.0 / 4.0 * self.N * (self.N + 1)
                if i == self.N-1:
                    if j == self.N-1:
                        d[i + 1, j + 1] = 1.0 / 4.0 * self.N * (self.N + 1)

        # Calculate matrix with 1st derivatives of Lagrange polynomials
        for n in range(-1, self.N):
            for i in range(-1, self.N):
                sum = 0
                for j in range(-1, self.N):
                    sum = sum + d[i + 1, j + 1] * self.lagrange_poly(n, xi[j + 1])

                out[n + 1, i + 1] = sum
        return out
    
    
    def legendre(self, x):
        '''
        Returning Legendre Polynomials P_N(x) at position x[-1, 1].
    
        '''
        
        P = np.zeros(2 * self.N)

        if self.N == 0:
            P[0] = 1
        elif self.N == 1:
            P[1] = x
        else:
            P[0] = 1
            P[1] = x
        for i in range(2, self.N + 1):
            P[i] = (1.0 / float(i)) * ((2 * i - 1) * x * P[i - 1] - (i - 1) * P[i - 2])

        return(P[self.N])
