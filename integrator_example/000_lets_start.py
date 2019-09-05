import numpy as np
import matplotlib.pyplot as plt

import mystyle as ms

g1 = -.2e3
g2 = -0.e1

def f(x):
   return g1 * x + g2 * x * x

def U(x):
    return -g1 * x**2/2 - g2 * x**3 / 3

def K(v):
    return v**2 / 2

def rk_integ_step(p, f, Dt):
    def F(p):
        x = p[0]
        v = p[1]
        return np.array([v, f(x)])
    k1 = F(p)
    k2 = F(p + Dt / 2 * k1)
    k3 = F(p + Dt / 2 * k2)
    k4 = F(p + Dt * k3)
    p += Dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

def symp_integ_step(p, f, Dt):
    p[0] += p[1] * Dt/2
    p[1] += f(p[0]) * Dt
    p[0] += p[1] * Dt/2

def simulate(f, p0, T_sim, Dt, integ_step):
    t_integ = np.arange(0, T_sim, Dt)
    p = p0.copy()
    p_vect = []
    for ii, tt in enumerate(t_integ):
        p_vect.append(p.copy())
        integ_step(p, f, Dt)
    p_vect = np.array(p_vect)
    return t_integ, p_vect

p0 = np.array([1., 0.])

flag_inset = True

T_sim = 3.2
markersize = 10
linewidth = 0
alpha = 0.7

# T_sim = 6000
# markersize = 0
# linewidth = 2
# alpha = 1

compare_against_harmonic = True

Dt_rk = 0.01
t_integ_rk, p_vect_rk = simulate(f, p0, T_sim, Dt_rk, rk_integ_step)

Dt_sym = 0.02
t_integ_sym, p_vect_sym = simulate(f, p0, T_sim, Dt_sym, symp_integ_step)


plt.close('all')
ms.mystyle_arial(fontsz=16, dist_tick_lab=5)
fig1 = plt.figure(1, figsize=(8*2, 6))
fig1.set_facecolor('w')
ax1 = fig1.add_subplot(2,1,1)
ax2 = fig1.add_subplot(2,1,2, sharex=ax1)

if flag_inset:
    axins = ax1.inset_axes([0.65, 0.15, 0.24, 0.8])
else:
    axins = None

for ax in [ax1, axins]:

    if ax is None:
        continue

    if compare_against_harmonic:
        ax.plot(t_integ_rk, np.cos(np.sqrt(-g1)*t_integ_rk), 'k', linewidth=2, alpha=1.)
    ax.plot(t_integ_sym, p_vect_sym[:, 0], linestyle='-', linewidth=linewidth, markersize=markersize,
        marker='.', alpha=alpha, color='blue')
    ax.plot(t_integ_rk, p_vect_rk[:, 0], linestyle='-', linewidth=linewidth, markersize=markersize,
        marker = '.', alpha=alpha, color='green')

ax1.set_ylim(-1.2, 1.2)
if flag_inset:
    axins.set_xlim(2.8, 2.95)
    axins.set_ylim(-1.1, -.5)
    axins.set_xticklabels('')
    axins.set_yticklabels('')
    # axins.grid(True)

    box, lines = ax1.indicate_inset_zoom(axins, alpha=.8, edgecolor='k')
    for ll in lines:
        ll.set_visible(True)
        ll.set_linewidth(1.5)
    box.set_linewidth(1.5)

ax2.plot(t_integ_sym, U(p_vect_sym[:, 0]) + K(p_vect_sym[:, 1]), label='Symplectic',
    linestyle='-', linewidth=linewidth, markersize=markersize,
    marker='.', alpha=alpha, color='blue')

ax2.plot(t_integ_rk, U(p_vect_rk[:, 0]) + K(p_vect_rk[:, 1]),
    linestyle='-', linewidth=linewidth, markersize=markersize, label='Runge-Kutta',
    marker='.', alpha=alpha, color='green')

E0 = U(p_vect_rk[0, 0]) + K(p_vect_rk[0, 0])

if compare_against_harmonic:
    ax2.plot(t_integ_rk, 0*t_integ_rk + E0, color='k', linewidth=2)

ax2.set_ylim(bottom=0, top = 1.1*E0)
ax1.set_xlim(0, T_sim*.95)

ax1.grid(True)
ax2.grid(True)

ax2.set_xlabel('Turn')

ax1.set_ylabel('x')
ax2.set_ylabel('Energy')
leg = ax2.legend(loc='lower left', prop={'size':16})
plt.show()

fig1.savefig(f'nturns_{T_sim:.0f}'+{True:'with_inset', False:''}[flag_inset]+'.png', dpi=200)

