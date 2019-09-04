import numpy as np

import matplotlib.pyplot as plt

from cpymad.madx import Madx
import interp_opt_functions_with_orbit as iof

optics_repo = '/afs/cern.ch/eng/lhc/optics/runII/2018'

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
mad.call(optics_repo + '/PROTON/opticsfile.30')

mad.use('lhcb1')
mad.twiss()

# find_ip 
i_ip = list(mad.table.twiss.name).index('ip5:1')
s_ip = mad.table.twiss.s[i_ip]

s, keyword, name, betx, bety, dx, dy, dpy, alfay, alfax, x, y, px, py = \
    iof.interp_opt_functions_with_orbit(mad.table.twiss)

plt.close('all')
fig1 = plt.figure(1)
ax1 = fig1.add_subplot(1,1,1)

ax1.fill_between(s-s_ip, -np.sqrt(betx), np.sqrt(betx), alpha=0.5)

# This can help
# http://akuederle.com/matplotlib-zoomed-up-inset

plt.show()

# This can help
# http://akuederle.com/matplotlib-zoomed-up-inset

plt.show()
