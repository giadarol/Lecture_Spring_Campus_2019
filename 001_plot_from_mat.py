import numpy as np

import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1.inset_locator import mark_inset


import mystyle as ms
import myfilemanager as mfm

large_plot_srange = [-5000, 5000]
large_plot_vrange = [-12, 12]
zoom_vrange = [-1, 1]
zoom_srange = [-60, 60]

gamma = 7000.
nemitt = 2.5e-6
epsx = nemitt/gamma
epsy = nemitt/gamma


ob = mfm.myloadmat_to_obj('optics_collection.mat')

betas_plot = np.logspace(np.log10(0.100001), np.log10(10.999), 200) 

for ii, betastar in enumerate(betas_plot[::-1]):

    ilow = np.max(np.where(ob.betastar_vect<betastar)[0])
    ihigh = ilow+1

    weight_low = (ob.betastar_vect[ihigh]-betastar)/(ob.betastar_vect[ihigh]-ob.betastar_vect[ilow])
    weight_high = 1. - weight_low

    # find_ip 
    s_ip = ob.s_ip
    s = ob.s
    s_quad = ob.s_quad 
    k1l_quad = weight_low * ob.k1l_mat[ilow] + weight_high * ob.k1l_mat[ihigh]
    sigmax = weight_low * ob.sigmax_mat[ilow] + weight_high * ob.sigmax_mat[ihigh]
    
    plt.close('all')
    ms.mystyle_arial(fontsz=16, dist_tick_lab=5)
    fig1 = plt.figure(1, figsize=(8*1.5,6*1.1))
    fig1.set_facecolor('white')
    ax1 = plt.subplot2grid(shape=(4,1), loc=(1,0), rowspan=3, fig=fig1)
    ax0 = plt.subplot2grid(shape=(4,1), loc=(0,0), rowspan=1, fig=fig1, sharex=ax1)
    
    axins = ax1.inset_axes([0.7, 0.02, 0.28, 0.28]) 
    
    for ax in [ax1, axins]:
        ax.fill_between(s-s_ip, -3*sigmax*1000, 3*sigmax*1000, alpha=0.5, color='b')
    #    ax.fill_between(s-s_ip, -3*sigmay*1000, 3*sigmay*1000, alpha=0.5, color='r')
    
    axins.set_xlim(zoom_srange)
    axins.set_xticklabels('')
    axins.set_yticklabels('')
    
    ax1.set_xlim(large_plot_srange)
    
    ax1.set_ylim(large_plot_vrange)
    axins.set_ylim(zoom_vrange)
    axins.axvspan(-0.5, 0.5, color='k', alpha=.4)

    ax1.indicate_inset_zoom(axins, alpha=.8, edgecolor='k')
    
    L_bar = 30.
    
    ax0.bar(x=s_quad-s_ip, height=k1l_quad, width=L_bar, linewidth=0, alpha=.5)
    ax0.set_ylim(-0.07, 0.07)
    ax0.set_yticks([-0.05, 0, 0.05])
    
    ax1.grid(True)
    ax0.grid(True)
    # This can help
    # http://akuederle.com/matplotlib-zoomed-up-inset

    fig1.suptitle('Beam size at interaction point = %.3f mm \n(beta* = %.2f m)'%(3*np.sqrt(betastar*epsx)*1e3, betastar))

    fig1.savefig('frame_%04d.png'%ii, dpi=200)


import os 
os.system(' '.join([
    'ffmpeg',
    '-i frame_%04d.png',
    '-c:v libx264 -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2,setpts=5*PTS"',
    '-profile:v high -level:v 4.0 -pix_fmt yuv420p -crf 22 -codec:a aac optics.mp4']))
