import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec 

def myplot(X, u, time, src, seis, xsrc, xrec, bd):
    #
    
    plt.ion()
    fig1  = plt.figure(figsize=(9, 6))
    gs = gridspec.GridSpec(2, 2, width_ratios=[1,1], hspace=0.3, wspace=0.3)

    ax1 = plt.subplot(gs[0:2])
    us, = ax1.plot(xsrc, 0, 'r*', markersize=10)
    ur, = ax1.plot(xrec, 0, 'r^', markersize=10)
    up, = ax1.plot(X, u, label = '%s boundary'%bd)
    ax1.set_xlim([X[0], X[-1]])
    #ax1.set_ylim()
    ax1.set_xlabel('distance (m)')
    ax1.set_ylabel('Amplitude')
    ax1.legend(loc='lower center')

    ax2 = plt.subplot(gs[2])
    ur2, = ax2.plot(xrec, 0, 'r^', markersize=10)
    up2_, = ax2.plot(time, seis, 'b', label='Seismogram recorded at the receiver')
    ax2.set_xlim([0, time[-1]])
    ax2.set_xlabel('time (s)')
    ax2.set_ylabel('Amplitude')
    ax2.legend(loc='lower center', ncol=2, fontsize=7)
    
    ax3 = plt.subplot(gs[3])
    up3, = ax3.plot(time, src, 'b', label='Source Time function')
    ax3.set_xlim([0, time[-1]])
    ax3.set_xlabel('time (s)')
    ax3.set_ylabel('Amplitude')
    ax3.legend(fontsize=7)

    #plt.show()
    
    return ax1, up, ax2, up2_
