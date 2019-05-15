# from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
import math


def f_gw(x, tc, t, B = 16.6):
    f = np.zeros(len(t))
    for i in range(len(t)):
        if tc > t[i]:
            f[i] = B*x**(-5/8)*(tc-t[i])**(-3/8)
        elif tc == t[i]:
            f[i] = 1e-5
        else:
            f[i] = -B*x**(-5/8)*(t[i]-tc)**(-3/8)
        # print(i)
    return f


def h_full(func, x, args, phi, Ampl=1, phi_c=0):
    f = func(x, *args)
    h = np.zeros(len(f))
    for i in range(len(f)):
        if f[i] > 0:
            h[i] = Ampl*f[i]**(2/3)*np.cos(phi[i]+phi_c)
        elif f[i] == 0:
            h = 1e-5
        else:
            h[i] = -Ampl*(-f[i])**(2/3)*np.cos(phi[i]+phi_c)

    h[np.where(h > 0.5)[0]] = 0
    return h


def phi(x, tc, t, B=16.6):
    phi = np.zeros(len(t))
    # t = t[np.where(t < tc)]
    for i in range(len(t)):
        if tc > t[i]:
            phi[i] = -2*3/8*np.pi*B*x**(-5/8)*math.pow((tc-t[i]), (5/8))
        elif tc == t[i]:
            phi[i] = 1e-5
        else:
            phi[i] = 2*3/8*np.pi*B*x**(-5/8)*math.pow((t[i]-tc), (5 / 8))
    # for i in range(len(t)):
    #     phi[i] = 2*np.pi*B*x**(-5/8)*(-3/8)*(tc-t[i])**(5/8)
        # print(i)
    return phi


def h_simplified(x, t, tc, phi_c=0, B=16.6):
    h = (B*x**(-5./8.)*(tc-t)**(-3./8.))**(2./3.)*np.cos(-39.11*x**(-5./8.)*(tc-t)**(5./8.)+phi_c)
    return h


def filter(data, tc, t, x, std):
    SN = 0
    norm = 0
    args = [tc, t]
    phivalues = phi(x, tc, t)
    # print('phivalues:\n', phivalues)
    s = h_full(f_gw, x, args, phivalues)
    # print('s:\n', s)
    # for i in range(len(t)):
    #     SN += (data[i]*s[i])**2
    #     norm += (s[i]/std)**2
    SN = sum((data*s)**2)
    norm = sum((s/std)**2)
    return SN/norm


def plot_signal(t, s):
    plt.figure()
    plt.plot(t, s, color = 'k')
    plt.xlabel('Time [s]')
    plt.ylabel('Signal: strain')


def chirp_mass(m1,m2):
    return (m1*m2)**(3/5)/(m1+m2)**(1/5)


def template(m1, m2):
    x = chirp_mass(m1, m2)

    tc = 0.48  # the moment of merger
    B = 16.6  # constant
    t = np.linspace(0, 0.45, 10000)  # time

    gw_frequency = B*x**(-5/8)*(tc-t)**(-3/8)
    amplitude = gw_frequency**(2/3)

    t_h = np.linspace(-450, 0, 10000)
    t_merge_h = 10
    phase = 2*np.pi*B*x**(-5/8)*(-3/8)*(t_merge_h - t_h)**(5/8)
    f = B * x ** (-5 / 8) * (t_merge_h - t_h) ** (-3 / 8)
    h = f ** (2 / 3) * np.cos(phase)

    return t, gw_frequency, t_h, h, x


def plot_chirplines():
    tc = 0.5
    for x in [20, 25, 30, 35, 40]:
        t = np.linspace(0.25, tc*0.9999, 100)
        f = f_gw(x, tc, t)
        plt.plot(t, f, label='x = '+ str(x))

    plt.axvline(tc, color='k', label='Merging time', linestyle='--')
    plt.legend()
    plt.xlabel('Time [s]')
    plt.ylabel('Frequency [Hz]')
    plt.semilogy()

    plt.figure(2)
    for tc in [0.3, 0.5, 1, 2 ,5, 10]:
        x = 30
        t = np.linspace(0, tc*0.9999, 100)
        f = f_gw(x, tc, t)
        plt.plot(t, f, label='tc = ' + str(tc))
        plt.axvline(tc, color='grey', linestyle='--', alpha=0.3)

    plt.legend()
    plt.xlabel('Time [s]')
    plt.ylabel('Frequency [Hz]')
    plt.semilogy()


def plot_grav_waves(xrange, tc, phi_c=0):
    t = np.linspace(-200, tc*0.99, 10000)

    # plt.figure()
    for x in xrange:
        h = h_simplified(x, t, tc, phi_c)
        plt.plot(t, h, label='x = ' + str(x))

    plt.legend()
    plt.xlabel('time')
    plt.ylabel('strain, h(t)')


# xrange = [15, 30, 45]
# plot_grav_waves(xrange, tc)
# plot_grav_waves([30], tc)
# plot_grav_waves([30], tc, phi_c=np.pi/2)

file = np.loadtxt('AllWithNoise.dat', unpack=True)
times = file[0, :]
s1 = file[1, :]
s2 = file[2, :]
s3 = file[3, :]

# tc_min = max(times)+1
# tc_max = tc_min+len(times)
std = 0.2
m = []

# f = filter(s3, tc_min, times, x, std)

x = 90
for tc in times:
    f = filter(s1, tc, times, x, std)
    if tc % 200 == 0:
        print('Iterating merge time: ', tc)
    # plt.plot(times, f)
    m.append(f)

plot_signal(times, s1)

plt.figure()
plt.plot(times, m, color='grey')
plt.xlabel('Time [s]')
plt.ylabel('Matched filter output')