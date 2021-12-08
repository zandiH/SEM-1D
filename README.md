# SEM-1D
A simple program for simulating 1-D wave equation using spectral element method.<br> The time integration schemes used are center finite difference and Newmark schemes, for 4 types of boundary conditions, including free surface or zero stress, rigid or zero displacement, periodic, and absorbing boundary conditions.<br>
The input parameters are:<br>
<pre> 
        L           :  domain size (m)
        nt          :   number of time points
        c0          :  velocity, for homogenous case (m/s)
        c           :  velocity vector, for heterogenous case
        rho0        :  density,  for homogenous case (kg/m^3)
        rho         :  density vector
        freq        :  central frequency of the source
        N           :  degree of lagrange polynomials
        xsrc        :  source position (m)
        xrec        :  receiver position (m)       
        dt          :  time step
        mode        :  time integration scheme, including 'cfd' (centered finite difference) 
                       and 'newmark'
        bd          :  boundary conditions, including 'rigid', 'free', 'absorbing', 'periodic'
        nw_acc_type :  type of acceleration change in newmark method, including 'average','linear'
                       defualt is 'average'<br>
</pre>

              
