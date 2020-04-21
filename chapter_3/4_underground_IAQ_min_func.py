# -*- coding: utf-8 -*-
# quelques couleurs
rouge_A, vert1_A, vert2_A, vert3_A, vert4_A, gris1_A = '#C60C2E', '#005157', '#627D77', '#9EB28F', '#C5E5A4', '#595A5C'
coule = [rouge_A, vert1_A, vert2_A, vert3_A, vert4_A, gris1_A]

from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize


# Crank-Nicolson scheme for air polluton equation in the underground
def fc_crank_nicolson(Cp, C, dt, alpha, beta, tau, delta, N, Ce):
    term_C = dt * (alpha * N ** 2 + (beta * N + tau) * (Ce - C) - delta * C)
    term_Cp = dt * (alpha * N ** 2 + (beta * N + tau) * (Ce - Cp) - delta * Cp)
    return -Cp + C + 0.5 * term_C + 0.5 * term_Cp


minimisation = False


def min_func(x):
    alpha, beta, delta, tau = x
    dt = 0.25  # 1/h
    time_meas = np.loadtxt("time")
    sim_time = max(time_meas)
    traffic = np.loadtxt("traffic")
    PM10 = np.loadtxt("PM10")
    PM10 = PM10 * 100
    Ce = 15
    Cinit = PM10[0]
    Ce = 15
    C = Cinit
    t, i = min(time_meas), 0
    time, concentration = [], []
    while t <= max(time_meas):
        time.append(t)
        concentration.append(C)
        N = traffic[i]
        # solve for C+, concentration at the next time step
        Cp = fsolve(fc_crank_nicolson, C, args=(C, dt, alpha, beta, tau, delta, N, Ce))
        C = Cp
        t += dt
        i += 1

    PM10 = np.asarray(PM10)
    concentration = np.asarray(concentration)
    dC = abs(concentration - PM10)
    if minimisation == True:
        return np.mean(dC)
    else:
        return time, concentration, np.mean(dC)


minimisation = True
#           alpha      beta     delta      tau
bnds = ((0.2, 0.5), (0.05, 0.2), (0.01, 10), (0.1, 1))
k = 0.5
x0 = [k * (bnds[0][0] + bnds[0][1]),
      k * (bnds[1][0] + bnds[1][1]),
      k * (bnds[2][0] + bnds[2][1]),
      k * (bnds[3][0] + bnds[3][1])]

# sol = minimize(min_func, x0, bounds=bnds, method='L-BFGS-B', tol=1e-3)
# sol = minimize(min_func, x0, bounds=bnds, method='SLSQP', tol=1e-3)
sol = minimize(min_func, x0, bounds=bnds, method='TNC', tol=1e-3)

alpha, beta, delta, tau = sol.x[0], sol.x[1], sol.x[2], sol.x[3]
minimisation = False
time, concentration, mean_error = min_func([alpha, beta, delta, tau])
print("alpha\t, beta\t, delta\t, tau")
print(sol.x)
print("mean error", round(mean_error[0], 2), " [µg/m3]")

PM10 = np.loadtxt("PM10")
PM10 = PM10 * 100
plt.clf()
plt.xlabel(r"Time [h]")
plt.ylabel(r"Underground $PM_{10}$ [µg/m$^3$]")
plt.plot(time, concentration, color=coule[-1], linestyle="--", alpha=0.75, marker='', label='model')
plt.plot(time, PM10, color=coule[1], linestyle="", alpha=0.25, marker='o', markersize='2.5',
         label='scaled measurements')
plt.legend()
plt.savefig("./PM10_underground.pdf", dpi=200, bbox_inches='tight')
