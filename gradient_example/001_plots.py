import numpy as np
from scipy.interpolate import interp2d

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D #Needed to allow projection='3d'
from matplotlib import cm

import myfilemanager as mfm

ob = mfm.myloadmat_to_obj('twiss_res.mat')

cost_func = np.sqrt((ob.betax_max-200)**2/200**2 + (ob.betax_min-45)**2/45**2)

i_max, j_max = np.unravel_index(np.argmin(cost_func), cost_func.shape)

start_point = (1.5, 1.2)
start_point = (1.8, 1.6)
#start_point = (1.1, 1.7)

Gx_map, Gy_map = np.gradient(cost_func)

Gx_map_norm = Gx_map / np.sqrt(Gx_map**2 + Gy_map**2) 
Gy_map_norm = Gy_map / np.sqrt(Gx_map**2 + Gy_map**2) 

XX, YY = np.meshgrid(ob.str_factor, ob.str_factor)
cost_interp = interp2d(XX, YY, cost_func.T, kind='linear')
Gx_interp = interp2d(XX, YY, Gx_map.T, kind='linear')
Gy_interp = interp2d(XX, YY, Gy_map.T, kind='linear')

gamma = 0.5

N_steps = 200
point_list = []
point = np.array(start_point)
grad_list = []
val_list = []
for ii in range(N_steps):
    #Gx = Gx_interp(point[0], point[1])[0]
    #Gy = Gy_interp(point[0], point[1])[0]

    ix = np.argmin(np.abs(ob.str_factor-point[0]))
    iy = np.argmin(np.abs(ob.str_factor-point[1]))

    Gx = Gx_map[ix, iy]
    Gy = Gy_map[ix, iy]

    gx = Gx #/ np.sqrt(Gx**2 +Gy**2)
    gy = Gy #/ np.sqrt(Gx**2 +Gy**2)

    point -= gamma*np.array([gx, gy])
    point_list.append(point.copy())

    grad_list.append(np.array([gx, gy]))

    val_list.append(cost_interp(point[0], point[1])[0])
    #val_list.append(cost_func[ix, iy])

point_list = np.array(point_list)
grad_list = np.array(grad_list)
val_list = np.array(val_list)

plt.close('all')

fig2 = plt.figure(2)
ax2 = fig2.add_subplot(1,1,1)
plt.pcolormesh(ob.str_factor, ob.str_factor, cost_func.T, vmax=1.)
plt.colorbar()
ax2.scatter(ob.str_factor[i_max], ob.str_factor[j_max], color='r')
ax2.scatter(start_point[0], start_point[1], color='y')
ax2.plot(point_list[:, 0], point_list[:,1], color='k', linewidth=2)
ax2.quiver(XX.flatten(), YY.flatten(), -100*Gx_map_norm.T.flatten(), -100*Gy_map_norm.T.flatten())

i_point = -1
ax2.quiver(point_list[i_point, 0], point_list[i_point, 1],
    -10*grad_list[i_point, 0], -10*grad_list[i_point, 1], color='r')

ax2.axis('equal')

fig1 = plt.figure(1, figsize=(1.4*8,1.4*6))
ax = fig1.gca(projection='3d')

i_cut = 14

cost_func[cost_func>1]=1.
surf = ax.plot_surface(XX[i_cut:, i_cut:], YY[i_cut:, i_cut:], 
        cost_func[i_cut:, i_cut:].T, alpha=0.5, #cmap=cm.coolwarm,
                       linewidth=0, antialiased=True)
ax.plot(point_list[:, 0], point_list[:,1], val_list, color='k', linewidth=2)
ax.scatter(ob.str_factor[i_max], ob.str_factor[j_max], cost_func[i_max, j_max], color='r')
ax.set_xlabel('Q1 strength')
ax.set_ylabel('Q2 strength')
ax.set_zlabel('Cost function')

fig1.subplots_adjust(bottom=.02, top=.98, left=.02, right=.98)
ax.azim=59
ax.elev=16


plt.show()
