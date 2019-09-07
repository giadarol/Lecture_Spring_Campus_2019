import numpy as np
from scipy.interpolate import interp2d

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D #Needed to allow projection='3d'
from matplotlib import cm

import myfilemanager as mfm

ob = mfm.myloadmat_to_obj('twiss_res.mat')
ob1 = mfm.myloadmat_to_obj('gradient_discent.mat')


target_betax_max = ob1.target_betax_max
target_betax_min = ob1.target_betax_min
cost_func = ob1.cost_func

point_list = ob1.point_list
val_list = ob1.val_list

XX, YY = np.meshgrid(ob.str_factor, ob.str_factor)

plt.close('all')

fig2 = plt.figure(2)
ax2 = fig2.add_subplot(1,1,1)
plt.pcolormesh(ob.str_factor, ob.str_factor, cost_func.T, vmax=1.)
plt.colorbar()
ax2.plot(point_list[:, 0], point_list[:,1], color='k', linewidth=2)

ax2.axis('equal')

fig1 = plt.figure(1, figsize=(1.4*8,1.4*6))
ax = fig1.gca(projection='3d')

i_cut = 14

cost_func[cost_func>1]=1.
surf = ax.plot_surface(XX[i_cut:, i_cut:], YY[i_cut:, i_cut:], 
        cost_func[i_cut:, i_cut:].T, alpha=0.5, #cmap=cm.coolwarm,
                       linewidth=0, antialiased=True)
ax.plot(point_list[:, 0], point_list[:,1], val_list, '.-',
        color='k', linewidth=2, markersize=10)
ax.set_xlabel('Q1 strength')
ax.set_ylabel('Q2 strength')
ax.set_zlabel('Cost function')

fig1.subplots_adjust(bottom=.02, top=.98, left=.02, right=.98)
ax.azim=59
ax.elev=16


plt.show()
