# SEM-1D
A simple programming for simulating 1-D wave equation using spectral element method.<br> The time integration schemes used are center finite difference and Newmark methods, for 4 types of boundary conditions, including free surface or zero stress, rigid or zero displacement, periodic, and absorbing boundary conditions.
The input parameters include, 
Input parameters are:<br>
<pre> 
            &emsp;&emsp;c0          :   velocity vector, for homogenous case
            &emsp;&emsp;c           :   velocity vector, for heterogenous case
            &emsp;&emsp;rho         :   density vector
            &emsp;&emsp;L           :   domain size
            &emsp;&emsp;N           :   degree of lagrange polynomials
            &emsp;&emsp;ne          :   number of elements
            &emsp;&emsp;dt          :   time step
            &emsp;&emsp;isrc        :   receiver position
            &emsp;&emsp;mode        :   time integration scheme, including 'cfd' (centered finite difference) or 'newmark'
            &emsp;&emsp;bd          :   boundary conditions, including 'rigid', 'free', 'absorbing', 'periodic'
            &emsp;&emsp;nw_acc_type :   type of acceleration change in newmark method, including 'average','linear'; defualt is 'average'<br>
</pre>
            The Lagrange class - A class for Lagrange polynomials's interpolation, integration, and derivation, 
            using GLL (Gauss Lobatto Legendre) collocation points - "Source code" is https://github.com/heinerigel/coursera/tree/master/Notebooks4Coursera/W9
              
