# Time Integration module 
import numpy as np

class Time_Int():
    
    '''
    A class for time integration in SEM formulation
    Two schemes used are centered finite differecen (cfd) and Newmark.
    Also, four boundary conditions are applied, including rigid, free, absorbing, 
        and periodic boundary conditions
    '''
    
    def __init__(self, K, m, c, rho, J, nep, dt, isrc, mode, bd, it=0, nw_acc_type='average'):
        
        '''
            K          :   stifness matrix
            m          :   mass vector
            c          :   velocity vector
            rho        :   density vector
            J          :   Jacobian
            nep        :   number of points in the domain
            dt         :   time step
            isrc       :   receiver position
            mode       :   integration scheme, either "cfd" (centered finite difference) or "newmark"
            bd         :   boundary conditions, in four cases: 'rigid', 'free', 'absorbing', 'periodic'
            it         :   time-step index
            nw_acc_type:   type of acceleration change in newmark method
        '''
        
        self.K           = K
        self.m           = m
        self.c           = c
        self.rho         = rho
        self.J           = J
        self.nep         = nep
        self.dt          = dt
        self.it          = it
        self.isrc        = isrc
        self.mode        = mode
        self.bd          = bd
        self.nw_acc_type = nw_acc_type
        # new variables
        self.minv = 1/m
        self.force_weight = (m[isrc]/ (rho[isrc] * J) ) * J # force weigth where force is applied
            
    
    def source(self, nt, freq, t0, d_amp, src_type, plot=False):
        '''
        Source: delta, source time function, or extended source in time and space!
        Here only a source time function is injected!
        
        nt    :  number of time steps
        f0    :  central frequency 
        t0    :  time shift of source time function
        d_amp :  amplitude of dirac source  
        
        '''

        self.src_type = src_type
        src = np.zeros(self.nep)
        
        if src_type == 'dirac':
            self.f[self.isrc] = d_amp * self.force_weight
        
        if src_type == 'gauss':
            time = np.arange(0, nt * self.dt, self.dt)
            self.src = - 8 * freq * (time - t0) * np.exp(-(4*freq)**2*(time - t0)**2)
            
        if (plot == True) and (src_type == 'gauss'):
            plt.figure()
            plt.plot(time, self.src)
            plt.xlabel('time (s)')
            plt.ylabel('Amplitude (m)')
            plt.title('Source time function')
            #plt.savefig('src.png')
         
        return time, self.src
            
        
    def scheme_init(self):
        # initialization of vectors for both schemes
        
        if self.mode == 'cfd':
            self.u = np.zeros(self.nep)
            self.uold = np.zeros(self.nep)
            self.unew = np.zeros(self.nep)            
            # force
            self.f = np.zeros(self.nep)
            
        if self.mode == 'newmark':
            self.u    = np.zeros(self.nep)
            self.unew = np.zeros(self.nep)
            self.v    = np.zeros(self.nep)
            self.vnew = np.zeros(self.nep)
            self.a    = np.zeros(self.nep)
            self.anew = np.zeros(self.nep)
            # force + boundary
            self.f = np.zeros(self.nep)
            self.C = np.zeros(self.nep) 
            # type of acceleration
            if self.nw_acc_type == 'average':
                self.beta  = 1/4
                self.gamma = 1/2 
            if self.nw_acc_type == 'linear':
                self.beta  = 1/6
                self.gamma = 1/2 
            
    
    def boundary(self, F):
        # Boundary conditions
        
        if self.bd == 'free': 
            pass
 
        if self.bd == 'rigid':
            self.minv[0]  = 0
            self.minv[-1] = 0 
            
            if self.mode == 'cfd':
                self.u[0]     = 0
                self.u[-1]    = 0
                self.uold[0]  = 0
                self.uold[-1] = 0
            if self.mode == 'newmark':
                self.u[0]  = 0
                self.u[-1] = 0
                self.v[0]  = 0
                self.v[-1] = 0
                self.a[0]  = 0
                self.a[-1] = 0
            
        if self.bd == 'periodic': 
            tmp   =  0.
            tmp   =  F[0]
            F[0]  =  F[0] + F[-1]
            F[-1] =  F[-1] + tmp
        
        if self.bd == 'absorbing': 
            
            if self.mode == 'cfd':
                b1 = self.rho[0]  * self.c[0]  * (self.u[0]  - self.uold[0])  / self.dt    
                b2 = self.rho[-1] * self.c[-1] * (self.u[-1] - self.uold[-1]) / self.dt
                self.f[0] = -b1
                self.f[-1]= -b2
                
            if self.mode == 'newmark':
                b1 = self.rho[0]  * self.c[0]    
                b2 = self.rho[-1] * self.c[-1] 
                self.C[0]  = -b1 
                self.C[-1] = -b2
                
            
    def run(self):
        # here, time integration begins
        # source injection
        if self.src_type == 'gauss':
            self.f[self.isrc] = self.src[self.it] * self.force_weight
        
        if self.mode == 'cfd':
            
            # Force or K*u
            F = (self.f - self.K@self.u)
            self.unew = self.dt**2 * self.minv * F + 2 * self.u - self.uold
            
            self.boundary(F)
            
            self.uold = self.u
            self.u = self.unew 
            
            return self.u
                
            
        if self.mode == 'newmark':
                
            # predictor
            self.unew = self.u + self.dt*self.v + (self.dt**2) * (0.05 - self.beta) * self.a
            self.vnew = self.v + self.dt * (1 - self.gamma) * self.a
            #self.a_new = 0
        
            # wave solution
            F = self.f + self.C*self.vnew - self.K@self.unew
            self.anew = self.minv * F
                
            # corrector
            self.vnew = self.vnew + self.gamma * self.dt * self.anew
            self.unew = self.unew + self.beta * (self.dt**2) * self.anew
            
            # boundary
            self.boundary(F)
            
            self.u = self.unew
            self.v = self.vnew
            self.a = self.anew
                        
            return self.u
        
