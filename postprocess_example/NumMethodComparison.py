from scipy import fftpack
import numpy as np
import h5py
import matplotlib as mpl
mpl.use('GTK3Cairo')
import matplotlib.pyplot as plt

fig1 = plt.figure(figsize=(14,2))
ax1 = fig1.add_axes((0.15, 0.32, 0.8, 0.6))
fig2 = plt.figure(figsize=(14,2))
ax2 = fig2.add_axes((0.15, 0.32, 0.8, 0.6))
fig3 = plt.figure(figsize=(14,2))
ax3 = fig3.add_axes((0.15, 0.32, 0.8, 0.6))
fig4 = plt.figure(figsize=(14,2))
ax4 = fig4.add_axes((0.15, 0.32, 0.8, 0.6))
fig5 = plt.figure(figsize=(7,2))
ax5 = fig5.add_axes((0.15, 0.32, 0.8, 0.6))
fig6 = plt.figure(figsize=(3.2,2))
ax6 = fig6.add_axes((0.32, 0.32, 0.6, 0.6))
fig7 = plt.figure(figsize=(3.2,2))
ax7 = fig7.add_axes((0.32, 0.32, 0.6, 0.6))
fig8 = plt.figure(figsize=(3.2,2))
ax8 = fig8.add_axes((0.32, 0.32, 0.6, 0.6))
fig9 = plt.figure(figsize=(7,2))
ax9 = fig9.add_axes((0.2, 0.32, 0.7, 0.6))

f1_Abq= np.loadtxt('d82_4_N30_8Hz_F0p10_dt1eM5_n30_U2.txt')
f2_Abq= np.loadtxt('d82_4_N30_76Hz_F0p50_dt1eM5_n30_U2.txt')
f3_Abq= np.loadtxt('d82_4_N30_96Hz_F0p60_dt1eM5_n30_U2.txt')
f4_Abq= np.loadtxt('d82_4_N30_35Hz_F1p40_dt1eM5_n30_U2.txt')
t1_Abq = np.transpose(f1_Abq[..., 0])
wEnd1_Abq = np.transpose(f1_Abq[..., 1])
t2_Abq = np.transpose(f2_Abq[..., 0])
wEnd2_Abq = np.transpose(f2_Abq[..., 1])
t3_Abq = np.transpose(f3_Abq[..., 0])
wEnd3_Abq = np.transpose(f3_Abq[..., 1])
t4_Abq = np.transpose(f4_Abq[..., 0])
wEnd4_Abq = np.transpose(f4_Abq[..., 1])

f1_NB = h5py.File("d82_4_NB_0.1N_8Hz.h5", 'r')
f2_NB = h5py.File("d82_4_NB_0.5N_76Hz.h5", 'r')
f3_NB = h5py.File("d82_4_NB_0.6N_96Hz.h5", 'r')
f4_NB = h5py.File("d82_4_NB_1.4N_35Hz.h5", 'r')
t1_NB = np.array([0])
t2_NB = np.array([0])
t3_NB = np.array([0])
t4_NB = np.array([0])
wEnd1_NB = np.array([0])
wEnd2_NB = np.array([0])
wEnd3_NB = np.array([0])
wEnd4_NB = np.array([0])
for n in range(6000):
	t1_NB = np.append(t1_NB, f1_NB['/u'+str(n+1)].attrs['t'][0])
	t2_NB = np.append(t2_NB, f2_NB['/u'+str(n+1)].attrs['t'][0])
	t3_NB = np.append(t3_NB, f3_NB['/u'+str(n+1)].attrs['t'][0])
	t4_NB = np.append(t4_NB, f4_NB['/u'+str(n+1)].attrs['t'][0])
	wEnd1_NB = np.append(wEnd1_NB, f1_NB['/u'+str(n+1)][-3])
	wEnd2_NB = np.append(wEnd2_NB, f2_NB['/u'+str(n+1)][-3])
	wEnd3_NB = np.append(wEnd3_NB, f3_NB['/u'+str(n+1)][-3])
	wEnd4_NB = np.append(wEnd4_NB, f4_NB['/u'+str(n+1)][-3])

f1_RK = h5py.File("d82_4_RK4_0.1N_8Hz_dt1eM5.h5", 'r')
f2_RK = h5py.File("d82_4_RK4_0.5N_76Hz_dt1eM5.h5", 'r')
f3_RK = h5py.File("d82_4_RK4_0.6N_96Hz_dt1eM5.h5", 'r')
f4_RK = h5py.File("d82_4_RK4_1.4N_35Hz_dt1eM5.h5", 'r')
t1_RK = f1_RK["/t"][...]
wEnd1_RK = f1_RK["/x"][...,-3]
t2_RK = f2_RK["/t"][...]
wEnd2_RK = f2_RK["/x"][...,-3]
t3_RK = f3_RK["/t"][...]
wEnd3_RK = f3_RK["/x"][...,-3]
t4_RK = f4_RK["/t"][...]
wEnd4_RK = f4_RK["/x"][...,-3]

inverse_arr = 1./(np.arange(0,len(wEnd1_RK))+1)
rmsErr_f1_NB = np.sqrt(inverse_arr*np.cumsum((wEnd1_Abq-wEnd1_NB)**2))
rmsErr_f2_NB = np.sqrt(inverse_arr*np.cumsum((wEnd2_Abq-wEnd2_NB)**2))
rmsErr_f3_NB = np.sqrt(inverse_arr*np.cumsum((wEnd3_Abq-wEnd3_NB)**2))
rmsErr_f4_NB = np.sqrt(inverse_arr*np.cumsum((wEnd4_Abq-wEnd4_NB)**2))
rmsErr_f1_RK = np.sqrt(inverse_arr*np.cumsum((wEnd1_Abq-wEnd1_RK)**2))
rmsErr_f2_RK = np.sqrt(inverse_arr*np.cumsum((wEnd2_Abq-wEnd2_RK)**2))
rmsErr_f3_RK = np.sqrt(inverse_arr*np.cumsum((wEnd3_Abq-wEnd3_RK)**2))
rmsErr_f4_RK = np.sqrt(inverse_arr*np.cumsum((wEnd4_Abq-wEnd4_RK)**2))
print(rmsErr_f1_NB[-1])
print(rmsErr_f2_NB[-1])
print(rmsErr_f3_NB[-1])
print(rmsErr_f4_NB[-1])
print(rmsErr_f1_RK[-1])
print(rmsErr_f2_RK[-1])
print(rmsErr_f3_RK[-1])
print(rmsErr_f4_RK[-1])

ax1.plot(t1_Abq, wEnd1_Abq, ls='-', lw=2, color='black', label='Abaqus')
ax1.plot(t1_NB[::10], wEnd1_NB[::10], ls='None', marker='x', markersize=8, color='red', label='NB')
ax1.plot(t1_RK[::10], wEnd1_RK[::10], ls='None', marker='s', markersize=5, color='blue', label="RK4")
ax2.plot(t2_Abq, wEnd2_Abq, ls='-', lw=2, color='black', label='Abaqus')
ax2.plot(t2_NB, wEnd2_NB, ls='None', marker='x', markersize=8, color='red', label='NB')
ax2.plot(t2_RK, wEnd2_RK, ls='None', marker='s', markersize=5, color='blue', label="RK4")
ax3.plot(t3_Abq, wEnd3_Abq, ls='-', lw=2, color='black', label='Abaqus')
ax3.plot(t3_NB, wEnd3_NB, ls='None', marker='x', markersize=8, color='red', label='NB')
ax3.plot(t3_RK, wEnd3_RK, ls='None', marker='s', markersize=5, color='blue', label="RK4")
ax4.plot(t4_Abq, wEnd4_Abq, ls='-', lw=2, color='black', label='Abaqus')
ax4.plot(t4_NB[::4], wEnd4_NB[::4], ls='None', marker='x', markersize=8, color='red', label='NB')
ax4.plot(t4_RK[::4], wEnd4_RK[::4], ls='None', marker='s', markersize=5, color='blue', label="RK4")
ax5.plot(t4_Abq, wEnd4_Abq, ls='-', lw=2, color='black', label='Abaqus')
ax5.plot(t4_NB[::4], wEnd4_NB[::4], ls='None', marker='x', markersize=8, color='red', label='NB')
ax5.plot(t4_RK[::4], wEnd4_RK[::4], ls='None', marker='s', markersize=5, color='blue', label="RK4")

ax1.set_xlim(0,1.0)
ax2.set_xlim(0.8,1.0)
ax3.set_xlim(0.8,1.0)
ax4.set_xlim(0,0.8)
ax5.set_xlim(3.6,4.0)
ax1.set_ylim(-0.03,0.03)
ax2.set_ylim(-0.024,0.024)
ax3.set_ylim(-0.05,0.05)
ax4.set_ylim(-1.3,1.4)
ax5.set_ylim(-1.3,1.4)
ax1.set_xlabel("Time (s)", fontsize=22)
#ax1.set_ylabel("Displacement (mm)", fontsize=22)
ax1.set_ylabel("$w_{2,30}$ (mm)", fontsize=22)
ax2.set_xlabel("Time (s)", fontsize=22)
ax2.set_ylabel("$w_2$ (mm)", fontsize=22)
ax3.set_xlabel("Time (s)", fontsize=22)
ax3.set_ylabel("$w_2$ (mm)", fontsize=22)
ax4.set_xlabel("Time (s)", fontsize=22)
ax4.set_ylabel("$w_2$ (mm)", fontsize=22)
ax5.set_xlabel("Time (s)", fontsize=22)
ax5.set_ylabel("$w_2$ (mm)", fontsize=22)
ax1.tick_params(labelsize=18)
ax2.tick_params(labelsize=18)
ax3.tick_params(labelsize=18)
ax4.tick_params(labelsize=18)
ax5.tick_params(labelsize=18)
ax1.legend(loc=1, fontsize=16)


### Frequency analysis
N = 4000
f_s = N/4.

wEndHat4_Abq = fftpack.fft(wEnd4_Abq[-N:])
freq = fftpack.fftfreq(N, 1./f_s)
mask = np.where(freq > 0)
wEndHat4_Abq = wEndHat4_Abq[mask]
freq = freq[mask]

wEndHat4_NB = fftpack.fft(wEnd4_NB[-N:])
wEndHat4_NB = wEndHat4_NB[mask]

wEndHat4_RK = fftpack.fft(wEnd4_RK[-N:])
wEndHat4_RK = wEndHat4_RK[mask]

#ax1.plot(t1_NB[::10], wEnd1_NB[::10], ls='None', marker='x', markersize=8, color='red', label='NB')
#ax1.plot(t1_RK[::10], wEnd1_RK[::10], ls='None', marker='s', markersize=5, color='blue', label="RK4")

ax6.plot(freq, abs(wEndHat4_Abq), ls='-', lw=2, color='black')
ax6.set_xlim(0,60)
ax6.set_xlabel("Frequency (Hz)", fontsize=22)
ax6.set_ylabel(r"$\hat{w}_2$ (Hz)", fontsize=22)
ax6.set_xticks([0,20,40,60])
ax6.set_yticks([0,200,400])
ax6.tick_params(labelsize=18)

ax7.plot(freq, abs(wEndHat4_NB), ls='None', marker='x', markersize=8, color='red')
ax7.set_xlim(0,60)
ax7.set_xlabel("Frequency (Hz)", fontsize=22)
ax7.set_ylabel(r"$\hat{w}_{2,30}$ (Hz)", fontsize=22)
ax7.set_xticks([0,20,40,60])
ax7.set_yticks([0,200,400])
ax7.tick_params(labelsize=18)

ax8.plot(freq, abs(wEndHat4_RK), ls='None', marker='s', markersize=5, color='blue')
ax8.set_xlim(0,60)
ax8.set_xlabel("Frequency (Hz)", fontsize=22)
ax8.set_ylabel(r"$\hat{w}_2$ (Hz)", fontsize=22)
ax8.set_xticks([0,20,40,60])
ax8.set_yticks([0,200,400])
ax8.tick_params(labelsize=18)

ax9.plot(freq, abs(wEndHat4_Abq)+100, ls='-', lw=2, color='black')
ax9.plot(freq, abs(wEndHat4_NB)+50, ls='None', marker='x', markersize=5, color='red')
ax9.plot(freq, abs(wEndHat4_RK), ls='None', marker='s', markersize=3, color='blue')
ax9.set_xlim(0,50)
ax9.set_xlabel("Frequency (Hz)", fontsize=22)
ax9.set_ylabel(r"$\hat{w}_2$ (Hz)", fontsize=22)
ax9.set_xticks([0,10,20,30,40,50])
ax9.set_yticks([0,200,400])
ax9.tick_params(labelsize=18)


plt.show()
