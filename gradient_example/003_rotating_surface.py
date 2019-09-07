import numpy as np
from scipy.interpolate import interp2d

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D #Needed to allow projection='3d'
from matplotlib import cm

import myfilemanager as mfm

from cpymad.madx import Madx

ob = mfm.myloadmat_to_obj('twiss_res.mat')
ob1 = mfm.myloadmat_to_obj('gradient_discent.mat')

cost_func = ob1.cost_func


XX, YY = np.meshgrid(ob.str_factor, ob.str_factor)

make_avi = True

i_cut = 14

   
plt.close('all')

import mystyle as ms
ms.mystyle(fontsz=14, dist_tick_lab=None, mpl_v1_style=False)

fig1 = plt.figure(1, figsize=(1.*8,1.*6))
ax = plt.subplot2grid(shape=(5,1), loc=(0,0), 
        rowspan=5, colspan=3, fig=fig1,
        projection='3d')


cost_func[cost_func>1]=1.
surf = ax.plot_surface(XX[i_cut:, i_cut:], YY[i_cut:, i_cut:], 
        cost_func[i_cut:, i_cut:].T, alpha=0.4, #cmap=cm.coolwarm,
                       linewidth=0, antialiased=True)


ax.set_xlabel('Q1 strength', labelpad=10)
ax.set_ylabel('Q2 strength', labelpad=10)
ax.set_zlabel('Cost function')

fig1.subplots_adjust(top=1., left=.03, wspace=.6, bottom=.03, right=.93)
ax.azim=59
ax.elev=16

ax.azim=104
ax.elev=12

ax.azim=121
ax.elev=12.5

for ii, dd in enumerate(range(0, 361, 2)):
    ax.azim = 121 + dd
    fig1.savefig('rot_iter%03d.png'%ii, dpi=200)

plt.show()


if make_avi:
    import os 
    os.system(' '.join([
        'ffmpeg',
        '-i rot_iter%03d.png',
        '-c:v libx264 -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2,setpts=3*PTS"',
        '-profile:v high -level:v 4.0 -pix_fmt yuv420p -crf 22 -codec:a aac rot.mp4']))
