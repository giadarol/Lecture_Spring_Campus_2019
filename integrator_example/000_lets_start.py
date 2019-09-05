import numpy as np
import matplotlib.pyplot as plt

g1 = -.5e3
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
T_sim = 3000

Dt_rk = 0.01
t_integ_rk, p_vect_rk = simulate(f, p0, T_sim, Dt_rk, rk_integ_step)

Dt_sym = 0.03
t_integ_sym, p_vect_sym = simulate(f, p0, T_sim, Dt_sym, symp_integ_step)


plt.close('all')
fig1 = plt.figure(1)
ax1 = fig1.add_subplot(2,1,1)
ax2 = fig1.add_subplot(2,1,2, sharex=ax1)

ax1.plot(t_integ_sym, p_vect_sym[:, 0], linestyle='-')
ax1.plot(t_integ_rk, p_vect_rk[:, 0], linestyle='-')

ax2.plot(t_integ_sym, U(p_vect_sym[:, 0]) + K(p_vect_sym[:, 1]))
ax2.plot(t_integ_rk, U(p_vect_rk[:, 0]) + K(p_vect_rk[:, 1]))

ax2.set_ylim(bottom=0)

plt.show()
