import numpy as np

import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1.inset_locator import mark_inset

from cpymad.madx import Madx

import mystyle as ms
import interp_opt_functions_with_orbit as iof

optics_repo = '/afs/cern.ch/eng/lhc/optics/runII/2018'


large_plot_srange = [-5000, 5000]
large_plot_vrange = [-12, 12]
zoom_vrange = [-1, 1]
zoom_srange = [-60, 60]

mad = Madx()

mad.options.echo = False
mad.options.warn = False
mad.options.info = False

gamma = 7000.
nemitt = 2.5e-6
epsx = nemitt/gamma
epsy = nemitt/gamma

mad.input(f"""Beam,particle=proton,sequence=lhcb1,energy=6500.0,NPART=1.2E11,
     sige=1.1e-4,sigt=0.075,ex={epsx},ey={epsy};""")
mad.call(optics_repo + "/lhc_as-built.seq")

betastar_vect = []
sigmax_mat = []
k1l_mat = []

for i_opt in range(1, 31):
    if i_opt==29:
        continue
    mad.call(optics_repo + '/PROTON/opticsfile.%d'%i_opt)
    
    mad.use('lhcb1')
    mad.twiss()
    
    # find_ip 
    i_ip = list(mad.table.twiss.name).index('ip5:1')
    s_ip = mad.table.twiss.s[i_ip]
    
    s, keyword, name, betx, bety, dx, dy, dpy, alfay, alfax, x, y, px, py = \
        iof.interp_opt_functions_with_orbit(mad.table.twiss)
    
    sigmax = np.sqrt(betx*epsx)
    sigmay = np.sqrt(bety*epsy)
    
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
    
    mask_quad = np.abs(mad.table.twiss.k1l)>0
    k1l_quad = mad.table.twiss.k1l[mask_quad]
    s_quad = mad.table.twiss.s[mask_quad]
    L_bar = 30.
    
    ax0.bar(x=s_quad-s_ip, height=k1l_quad, width=L_bar, linewidth=0, alpha=.5)
    ax0.set_ylim(-0.07, 0.07)
    ax0.set_yticks([-0.05, 0, 0.05])
    
    ax1.grid(True)
    ax0.grid(True)
    # This can help
    # http://akuederle.com/matplotlib-zoomed-up-inset

    betastar = mad.table.twiss.betx[i_ip]

    betastar_vect.append(betastar)
    sigmax_mat.append(sigmax)
    k1l_mat.append(k1l_quad)

    fig1.suptitle('Beam size at interaction point = %.f micrometers \n(beta* = %.0f cm)'%(3*np.sqrt(betastar*epsx)*1e6, 100*betastar))
    fig1.savefig('opt%03d.png'%i_opt, dpi=200)

