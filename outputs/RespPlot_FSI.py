from scipy import fftpack
from scipy import signal
import numpy as np
import h5py
import matplotlib.pyplot as plt

fig1 = plt.figure(figsize=(4,2.5))
ax1 = fig1.add_axes((0.2, 0.26, 0.76, 0.7))
fig2 = plt.figure(figsize=(4,2.5))
ax2 = fig2.add_axes((0.2, 0.26, 0.76, 0.7))
fig3 = plt.figure(figsize=(4,5))
ax3 = fig3.add_axes((0.23, 0.13, 0.73, 0.8))
fig4 = plt.figure(figsize=(6,4))
ax4 = fig4.add_axes((0.23, 0.2, 0.68, 0.7))
### solitonic resonance from structural excitation
#ax4 = fig4.add_axes((0.23, 0.2, 0.30, 0.7))
#ax4a = fig4.add_axes((0.58, 0.2, 0.30, 0.7))

f1 = h5py.File("P4_subsonic17InPlane_A2.0_f1.0_RK4.h5")
f2 = h5py.File("P4_subsonic17InPlane_F0.3_f190.0_RK4.h5")
f3 = h5py.File("P4_subsonic17InPlane_F0.06_f190.0_rand_RK4.h5")

t1 = f1["/t"][...]
wOut1 = f1["/u"][...,6*20-3]
uOut1 = f1["/u"][...,6*20-6]

Nt = round(len(t1)/4); print(Nt)
t_tot = round(f1["/t"][-1]/4,8); print("t_tot: ", t_tot)
f_s = Nt/t_tot

wOutPart1 = wOut1[-Nt:]
wOutPartHat1 = fftpack.fft(wOutPart1)
wOutPartHat1 = wOutPartHat1/max(abs(wOutPartHat1))
#wOutPartHat = fftpack.fft(wOutPart*signal.blackman(Nt))
uOutPart1 = uOut1[-Nt:]
uOutPartHat1 = fftpack.fft(uOutPart1)
freq1 = fftpack.fftfreq(Nt, 1./f_s)
mask1 = np.logical_and(freq1 > 10, freq1 < 100000)
#mask = np.where(freq > 0)


t2 = f2["/t"][...]
wOut2 = f2["/u"][...,6*20-3]
uOut2 = f2["/u"][...,6*20-6]

Nt = round(len(t2)/4); print(Nt)
t_tot = round(f2["/t"][-1]/4,8); print("t_tot: ", t_tot)
f_s = Nt/t_tot

wOutPart2 = wOut2[-Nt:]
wOutPartHat2 = fftpack.fft(wOutPart2)
wOutPartHat2 = wOutPartHat2/max(abs(wOutPartHat2))
freq2 = fftpack.fftfreq(Nt, 1./f_s)
mask2 = np.logical_and(freq2 > 10, freq2 < 100000)


t3 = f3["/t"][...]
wOut3 = f3["/u"][...,6*20-3]
uOut3 = f3["/u"][...,6*20-6]

Nt = round(len(t3)/4); print(Nt)
t_tot = round(f3["/t"][-1]/4,8); print("t_tot: ", t_tot)
f_s = Nt/t_tot

wOutPart3 = wOut3[-Nt:]
wOutPartHat3 = fftpack.fft(wOutPart3)
wOutPartHat3 = wOutPartHat3/max(abs(wOutPartHat3))
freq3 = fftpack.fftfreq(Nt, 1./f_s)
mask3 = np.logical_and(freq3 > 10, freq3 < 100000)

ax1.plot(t1, uOut1, lw=2, color='black')
ax2.plot(t1, wOut1, lw=2, color='blue')
#ax3.plot(freq[mask], abs(uOutPartHat[mask]), lw=2, color='red')
ax4.plot(freq1[mask1], abs(wOutPartHat1[mask1]), lw=2, color='magenta', label="quasistatic")
ax4.plot(freq2[mask2], abs(wOutPartHat2[mask2]) - 0.15, lw=2, color='blue', label="harmonic")
ax4.plot(freq3[mask3], abs(wOutPartHat3[mask3]) - 0.3, lw=2, color='red', label="random")

#ax1.set_xlim(0,0.01)
#ax1.set_ylim(-8,24)
ax1.set_xlabel("Time (s)", fontsize=20)
ax1.set_ylabel("$u_{1}$ (mm)", fontsize=20)
#ax1.set_xticks([5, 5.5, 6])
#ax1.set_yticks([-8, 0, 8, 16, 24])
ax1.tick_params(axis='x', labelsize=16)
ax1.tick_params(axis='y', labelsize=16)

#ax2.set_xlim(0,0.01)
#ax2.set_ylim(-1,1)
ax2.set_xlabel("Time (s)", fontsize=20)
ax2.set_ylabel("$w_{2}$ (mm)", fontsize=20)
#ax2.set_xticks([5, 5.5, 6])
ax2.tick_params(axis='x', labelsize=16)
ax2.tick_params(axis='y', labelsize=16)

#ax3.set_xlim(0,500000)
#ax3.set_ylim(0,1)
ax3.set_xlabel("Frequency (Hz)", fontsize=20)
ax3.set_ylabel("$\hat{u}_{1}$", fontsize=20)
#ax3.set_xticks([10, 26.281, 40, 60])
#ax3.set_xticklabels([10, 26.281,40, 60])
ax3.tick_params(axis='x', labelsize=16)
ax3.tick_params(axis='y', labelsize=16)

ax4.set_xlim(20,1220)
#ax4.set_ylim(0,40000)
ax4.set_xlabel("Frequency (Hz)", fontsize=20)
ax4.set_ylabel("$\hat{w}_{2,30}$", fontsize=20)
#ax4.set_xticks([0, 1500, 3000])
#ax4.set_xticklabels([10, 26.3, 40, 60])
ax4.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax4.tick_params(axis='x', labelsize=16)
ax4.tick_params(axis='y', labelsize=16)
#ax4.set_yscale('log')
ax4.legend(loc=2, fontsize=14, handlelength=2.0)

#ax4a.set_xlim(86000,90000)
#ax4a.set_ylim(0,40000)
#ax4a.set_xticks([87000, 89000])
#ax4a.tick_params(axis='x', labelsize=16)
#ax4a.tick_params(axis='y', labelsize=16)

plt.show()
