import numpy as np
import matplotlib.pyplot as plt

from cpymad.madx import Madx


L_halfcell = 50.
phase_adv_cell = np.pi/3
n_cells_arc = 11 #Must be odd
n_arcs = 4
n_dip_half_cell = 3


mad = Madx()

circum = n_arcs*(n_cells_arc)*L_halfcell*2.
n_dip_total = n_arcs*(n_cells_arc)*n_dip_half_cell*2


# Bending angle
b_ang = 2.*np.pi/n_dip_total

# Strength quadrupole
focal_length = L_halfcell/np.sin(phase_adv_cell*1.111111) # I add 10% to avoid integer resonance
k1l_quad = 1./focal_length

# Start building sequence
sequence = ''

sequence+='''
CAV: RFCAVITY, L := 0.;
!Quad strengths (to be matched)
kqf:=%e;'''%k1l_quad
sequence+='''
kqd:=%e;'''%(-k1l_quad)

sequence+='''
!Sext strengths (to be matched)
ksf:=0.01;''' # match needs to start from a non zero value
sequence+='''
ksd:=0.01;''' # match needs to start from a non zero value

sequence+='''
mb: multipole,knl=%e;'''%b_ang
sequence+='''
qf: multipole,knl:={0,kqf};
qd: multipole,knl:={0,kqd};
sf: multipole,knl:={0,0,ksf};
sd: multipole,knl:={0,0,ksd};'''


sequence+='''
circum = %e;
toyring: sequence, refer=centre, l=circum; 
CAV1:CAV, at=0.;
'''%circum

s_start_cell = 0.
for i_arc in range(n_arcs):

	#Half Arc left
	for i_cell in range(n_cells_arc):
		
		sequence += 'qf_arc%dl_cell%d: qf, at=%e;\n'%(i_arc, i_cell, s_start_cell)
		sequence += 'sf_arc%dl_cell%d: sf, at=%e;\n'%(i_arc, i_cell, s_start_cell)
		for i_bend in range(n_dip_half_cell):
			sequence += 'mb_arc%dl_cell%d_%d: mb, at=%e;\n'%(i_arc, i_cell, i_bend, 
									s_start_cell+(i_bend+1)*L_halfcell/(n_dip_half_cell+1))
		sequence += 'qd_arc%dl_cell%d: qd, at=%e;\n'%(i_arc, i_cell, s_start_cell+L_halfcell)
		sequence += 'sd_arc%dl_cell%d: sd, at=%e;\n'%(i_arc, i_cell, s_start_cell+L_halfcell)
		for i_bend in range(n_dip_half_cell):
			sequence += 'mb_arc%dl_cell%d_%d: mb, at=%e;\n'%(i_arc, i_cell, i_bend+n_dip_half_cell, 
									s_start_cell+(i_bend+1)*L_halfcell/(n_dip_half_cell+1)+L_halfcell)
		s_start_cell += L_halfcell*2


	

sequence+='''
end_machine: marker at=circum;
endsequence;
'''

mad.input(sequence)
mad.input('beam, particle = proton, sequence=toyring, energy = 6500., NPART=1.05E11, sige= 2.5e;')
mad.use('toyring')

str_factor_list = np.linspace(.7, 2, 20)

plt.close('all')
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)

betax_max = np.zeros((len(str_factor_list), len(str_factor_list)))
betax_min = np.zeros((len(str_factor_list), len(str_factor_list)))

for ii1, ff1 in enumerate(str_factor_list):
    mad.input('kqf := %e'%(ff1*k1l_quad))
    for ii2, ff2 in enumerate(str_factor_list):
        mad.input('kqd := %e'%(-ff2*k1l_quad))
        try:
            tw = mad.twiss()
            plt.plot(tw.s, tw.betx)
            betax_max[ii1, ii2] = np.max(tw.betx)
            betax_min[ii1, ii2] = np.min(tw.betx)

        except:
            print('Skipped')

plt.figure(); plt.pcolormesh(betax_max.T, vmax=500); plt.colorbar()
plt.figure(); plt.pcolormesh(betax_min.T, vmax=200); plt.colorbar()

cost_func = np.sqrt((betax_max-200)**2 + 4*(betax_min-40)**2)

i_max, j_max = np.unravel_index(np.argmin(cost_func), cost_func.shape)

plt.figure(); pcolormesh(cost_func.T, vmax = 100)
plt.plot(i_max, j_max, '.', color='r', markersize = 10)

import scipy.io as sio
sio.savemat('twiss_res.mat', {
    'k1l_quad': k1l_quad,
    'str_factor': str_factor_list,
    'betax_max': betax_max,
    'betax_min': betax_min,
    }, oned_as='row')


plt.show()



