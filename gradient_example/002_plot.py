import numpy as np
from scipy.interpolate import interp2d

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D #Needed to allow projection='3d'
from matplotlib import cm

import myfilemanager as mfm

from cpymad.madx import Madx

ob = mfm.myloadmat_to_obj('twiss_res.mat')
ob1 = mfm.myloadmat_to_obj('gradient_discent.mat')

mad = Madx()

with open('sequence.seq', 'r') as fid:
    mad.input(fid.read())
mad.use('toyring')

target_betax_max = ob1.target_betax_max
target_betax_min = ob1.target_betax_min
cost_func = ob1.cost_func

point_list = ob1.point_list
val_list = ob1.val_list
grad_list = ob1.grad_list

grad_list[0, :] = grad_list[1, :] # pathing a small bug :-)

XX, YY = np.meshgrid(ob.str_factor, ob.str_factor)

n_iter_range = [0, 140]
plot_grad = True
n_iter_range = [140, 141]
plot_grad = False


make_avi = True

scale_grad_xy = 0.2
scale_grad_all = 20.
i_cut = 14
L_bar = 10

vx = -ob1.grad_list[0, 0]*scale_grad_xy 
vy = -ob1.grad_list[0, 1]*scale_grad_xy
mod_grad = np.sqrt(vx**2 + vy**2)

mod_0 = np.sqrt(vx**2 + vy**2 + mod_grad**2) 

for n_iter_plot in range(n_iter_range[0], n_iter_range[1]):

    ff1 = point_list[n_iter_plot, 0]
    ff2 = point_list[n_iter_plot, 1]
    
    mad.input('kqf := %e'%(ff1*ob.k1l_quad))
    mad.input('kqd := %e'%(-ff2*ob.k1l_quad))
    tw = mad.twiss()
    
    plt.close('all')
    
    import mystyle as ms
    ms.mystyle(fontsz=14, dist_tick_lab=None, mpl_v1_style=False)
 
    fig1 = plt.figure(1, figsize=(1.7*8,1.*6))
    ax = plt.subplot2grid(shape=(5,5), loc=(0,0), 
            rowspan=5, colspan=3, fig=fig1,
            projection='3d')
    ax30 = plt.subplot2grid(shape=(5,5), loc=(1,3), 
            rowspan=1, colspan=2, fig=fig1)
    ax3 = plt.subplot2grid(shape=(5,5), loc=(2, 3),
            rowspan=3, colspan=2, fig=fig1,
            sharex=ax30)

    for signq in [1., -1]:
        mask_quad = signq*mad.table.twiss.k1l>0
        k1l_quad = mad.table.twiss.k1l[mask_quad]
        s_quad = mad.table.twiss.s[mask_quad]
        ax30.bar(x=s_quad, height=k1l_quad/ob.k1l_quad,
            width=L_bar, linewidth=0, alpha=.5)
    ax30.grid(True)
    ax30.set_ylim(-1.9, 1.9)
    ax30.set_ylabel('Quad.\nstrength')
    ax30.set_xticklabels('')
    
    ax3.plot(tw.s, np.sqrt(tw.betx), color='b', linewidth=2)
    ax3.axhline(y=np.sqrt(target_betax_max), 
            linestyle='-', linewidth=2, color='r')
    ax3.axhline(y=np.sqrt(target_betax_min), 
            linestyle='-', linewidth=2, color='r')
    ax3.set_xlim(0, 500)
    ax3.set_ylim(0, 18)
    ax3.grid(True)
    ax3.set_xlabel('s [m]')
    ax3.set_ylabel('Beam size [mm]')
    
    fig2 = plt.figure(2)
    ax2 = fig2.add_subplot(1,1,1)
    plt.pcolormesh(ob.str_factor, ob.str_factor, cost_func.T, vmax=1.)
    plt.colorbar()
    ax2.plot(point_list[:, 0], point_list[:,1], color='k', linewidth=2)
    
    ax2.axis('equal')
    
   
    
    cost_func[cost_func>1]=1.
    surf = ax.plot_surface(XX[i_cut:, i_cut:], YY[i_cut:, i_cut:], 
            cost_func[i_cut:, i_cut:].T, alpha=0.4, #cmap=cm.coolwarm,
                           linewidth=0, antialiased=True)
    ax.plot(point_list[:n_iter_plot+1, 0], point_list[:n_iter_plot+1,1], 
            val_list[:n_iter_plot+1], '.-',
            color='k', linewidth=2, markersize=10)
    
    if plot_grad:
        vx = -ob1.grad_list[n_iter_plot, 0]*scale_grad_xy 
        vy = -ob1.grad_list[n_iter_plot, 1]*scale_grad_xy
        mod_grad = np.sqrt(vx**2 + vy**2)

        mod = np.sqrt(vx**2 + vy**2 + mod_grad**2)

        vect_to_plot = []
        ax.quiver(point_list[n_iter_plot, 0], point_list[n_iter_plot,1], 
            val_list[n_iter_plot]*1.01, 
            scale_grad_all*vx/mod*mod_0, 
            scale_grad_all*vy/mod*mod_0,
            -scale_grad_all*mod_grad/mod*mod_0, color='red') 
    
    
    ax.set_xlabel('Q1 strength', labelpad=10)
    ax.set_ylabel('Q2 strength', labelpad=10)
    ax.set_zlabel('Cost function')
    
    fig1.subplots_adjust(top=1., left=.0, wspace=.6, bottom=.12)
    ax.azim=59
    ax.elev=16
    
    ax.azim=104
    ax.elev=12
    
    ax.azim=121
    ax.elev=12.5
    
    
    fig1.suptitle('\nIteration %d'%n_iter_plot)
    
    fig1.savefig('grad_3d_iter%03d.png'%n_iter_plot, dpi=200)

plt.show()


if make_avi:
    
    import os 
    os.system(' '.join([
        'ffmpeg',
        '-i grad_3d_iter%03d.png',
        '-c:v libx264 -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2,setpts=7*PTS"',
        '-profile:v high -level:v 4.0 -pix_fmt yuv420p -crf 22 -codec:a aac gradient.mp4']))
