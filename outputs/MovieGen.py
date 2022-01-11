import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import matplotlib as mpl

filename = "P1_d01_F0.0_f0.0_RK4.h5"
#filename = "P2_d02_F0.0_f0.0_RK4.h5"
#filename = "P4_d04_F1.4_f35.0_RK4.h5"
#filename = "P4_d04_F0.1_f8.0_RK4.h5"
fskip = 100
f = h5py.File(filename, 'r')

Nt= len(f["/t"])
N_frame = int(Nt/fskip)
frame_range = range(0,N_frame)
#t_range = [2.0, 3.0]
#dt = f["/t"][1] - f["/t"][0]
#frame_range = range(int(t_range[0]/dt/fskip), int(t_range[1]/dt/fskip)+1)

dt_frame = 1
prob_flag = f.attrs["ProblemType"][0]

if prob_flag == 1:
    def init():
        mass.set_data([], [])
        line.set_data([], [])
        time_text.set_text('')
    #    time_scaled_text.set_text('')
        return mass, line, time_text

    def animate(i, site, fskip, L):
        global ymin, ymax

        i = i*fskip
        t = f["/t"][i]
        u = f["/u"][i,0:]
        x = site
        y = L*np.sin(u-np.pi/2)
        mass.set_data(x, y)
        line.set_data(site, u[...])

        if 'ymin' in globals():
            umin = u.min(); umax = u.max()
            if umin < ymin:
                ymin = umin
            if umax > ymax:
                ymax = umax
        else:
            ymin = u.min(); ymax = u.max()

        ax1.set_ylim(1.1*ymin,1.1*ymax)
        time_text.set_text("Time = {:8.3f}".format(t))
        return mass, line, time_text

    N = f.attrs["N"][0]
    L = f.attrs["L"][0]

    fig = plt.figure(figsize=(16,8))
    ax0 = fig.add_axes((0.1, 0.7, 0.85, 0.2))
    ax1 = fig.add_axes((0.1, 0.1, 0.85, 0.5))
    ax0.set_axis_off()
    ax0.set_xlim(0,N+1)
    ax0.set_ylim(-1.2*L,1.2*L)
    ax1.set_xlim(0,N+1)
    ax1.set_xlabel("Site", fontsize=18)
    ax1.set_ylabel("Angle", fontsize=18)
    ax1.tick_params(labelsize=16)

    time_text = ax0.text(0.87, 1.05, '', fontsize=14, transform=ax0.transAxes)
    mass, = ax0.plot([], [], 'bo', ms=4)
    ax0.plot([0,N+1], [0,0], '--', lw=1, c='black')
    site = np.arange(1,N+1)
    line, = ax1.plot([], [], 'o--', ms=5, lw=2, c='black')
    ymin = 0
    ymax = 0

    ani = ani.FuncAnimation(fig, animate, init_func=init, frames=N_frame, fargs=(site, fskip, L), interval=dt_frame, blit=True)

elif prob_flag == 2:
    def init():
        mass.set_data([], [])
        line.set_data([], [])
        time_text.set_text('')
    #    time_scaled_text.set_text('')
        return mass, line, time_text

    def animate(i, site, fskip):
        global ymin, ymax

        i = i*fskip
        t = f["/t"][i]
        u = f["/u"][i,0:]
        x = site + u
        mass.set_data(x, 0.0)
        line.set_data(site, u[...])

        if 'ymin' in globals():
            umin = u.min(); umax = u.max()
            if umin < ymin:
                ymin = umin
            if umax > ymax:
                ymax = umax
        else:
            ymin = u.min(); ymax = u.max()

        ax1.set_ylim(1.1*ymin,1.1*ymax)
        time_text.set_text("Time = {:8.3f}".format(t))
        return mass, line, time_text

    N = f.attrs["N"][0]

    fig = plt.figure(figsize=(16,8))
    ax0 = fig.add_axes((0.1, 0.7, 0.85, 0.2))
    ax1 = fig.add_axes((0.1, 0.1, 0.85, 0.5))
    ax0.set_axis_off()
    ax0.set_xlim(0,N+1)
    ax0.set_ylim(-1,1)
    ax1.set_xlim(0,N+1)
    ax1.set_xlabel("Site", fontsize=18)
    ax1.set_ylabel("Displacement", fontsize=18)
    ax1.tick_params(labelsize=16)

    time_text = ax0.text(0.87, 1.05, '', fontsize=14, transform=ax0.transAxes)
    mass, = ax0.plot([], [], 'bo', ms=4)
    site = np.arange(1,N+1)
    line, = ax1.plot([], [], 'o--', ms=5, lw=2, c='black')

    ani = ani.FuncAnimation(fig, animate, init_func=init, frames=N_frame, fargs=(site, fskip), interval=dt_frame, blit=True)


elif prob_flag == 4:
    global xOrig, yOrig

    def init():
        mass.set_data([], [])
        line1.set_data([], [])
        line2.set_data([], [])
        line3.set_data([], [])
        line4.set_data([], [])
        line5.set_data([], [])
        line6.set_data([], [])
        time_text.set_text('')
    #    time_scaled_text.set_text('')
        return mass, line1, line2, line3, line4, line5, line6, time_text

    def animate(i, site, fskip):
        global y1min, y1max, y2min, y2max, y3min, y3max, y4min, y4max, y5min, y5max, y6min, y6max

        i = i*fskip
        t = f["/t"][i]
        x = xOrig + f["/u"][i,0::2]
        y = yOrig + f["/u"][i,1::2]
        u1 = f["/u"][i,0::6]
        w1 = f["/u"][i,1::6]
        u2 = f["/u"][i,2::6]
        w2 = f["/u"][i,3::6]
        u3 = f["/u"][i,4::6]
        w3 = f["/u"][i,5::6]
        mass.set_data(x, y)
        line1.set_data(site, u1[...])
        line2.set_data(site, w1[...])
        line3.set_data(site, u2[...])
        line4.set_data(site, w2[...])
        line5.set_data(site, u3[...])
        line6.set_data(site, w3[...])

        if 'y1min' in globals():
            u1min = u1.min(); u1max = u1.max()
            if u1min < y1min:
                y1min = u1min
            if u1max > y1max:
                y1max = u1max
            ax1.set_ylim(1.1*y1min,1.1*y1max)

            w1min = w1.min(); w1max = w1.max()
            if w1min < y2min:
                y2min = w1min
            if w1max > y2max:
                y2max = w1max
            ax2.set_ylim(1.1*y2min,1.1*y2max)

            u2min = u2.min(); u2max = u2.max()
            if u2min < y3min:
                y3min = u2min
            if u2max > y3max:
                y3max = u2max
            ax3.set_ylim(1.1*y3min,1.1*y3max)

            w2min = w2.min(); w2max = w2.max()
            if w2min < y4min:
                y4min = w2min
            if w2max > y4max:
                y4max = w2max
            ax4.set_ylim(1.1*y4min,1.1*y4max)

            u3min = u3.min(); u3max = u3.max()
            if u3min < y5min:
                y5min = u3min
            if u3max > y5max:
                y5max = u3max
            ax5.set_ylim(1.1*y5min,1.1*y5max)

            w3min = w3.min(); w3max = w3.max()
            if w3min < y6min:
                y6min = w3min
            if w3max > y6max:
                y6max = w3max
            ax6.set_ylim(1.1*y6min,1.1*y6max)

        else:
            y1min = u1.min(); y1max = u1.max()
            y2min = w1.min(); y2max = w1.max()
            y3min = u2.min(); y3max = u2.max()
            y4min = w2.min(); y4max = w2.max()
            y5min = u3.min(); y5max = u3.max()
            y6min = w3.min(); y6max = w3.max()

        time_text.set_text("Time = {:8.3f}".format(t))
        return mass, line1, line2, line3, line4, line5, line6, time_text

    N = f.attrs["N"][0]
    L = f.attrs["L"]

    fig = plt.figure(figsize=(16,12))
    #fig, axes = plt.subplots(2, 1, figsize=(16,8), gridspec_kw={'height_ratios': [1,3]})
    ax0 = fig.add_axes((0.1, 0.80, 0.85, 0.15))
    ax1 = fig.add_axes((0.1, 0.675, 0.85, 0.105))
    ax2 = fig.add_axes((0.1, 0.55, 0.85, 0.105))
    ax3 = fig.add_axes((0.1, 0.425, 0.85, 0.105))
    ax4 = fig.add_axes((0.1, 0.30, 0.85, 0.105))
    ax5 = fig.add_axes((0.1, 0.175, 0.85, 0.105))
    ax6 = fig.add_axes((0.1, 0.05, 0.85, 0.105))

    ax0.set_axis_off()
    ax0.set_xlim(-0.05*N*L[0], 1.05*N*L[0])
    ax0.set_ylim((L[1]+L[2])/2-(1.1*N*L[0])/20, (L[1]+L[2])/2+(1.1*N*L[0])/20)
    ax1.set_xlim(0,N+1)
    ax2.set_xlim(0,N+1)
    ax3.set_xlim(0,N+1)
    ax4.set_xlim(0,N+1)
    ax5.set_xlim(0,N+1)
    ax6.set_xlim(0,N+1)
    ax1.set_ylabel("$u_{1}$ (mm)", fontsize=16)
    ax2.set_ylabel("$w_{1}$ (mm)", fontsize=16)
    ax3.set_ylabel("$u_{2}$ (mm)", fontsize=16)
    ax4.set_ylabel("$w_{2}$ (mm)", fontsize=16)
    ax5.set_ylabel("$u_{3}$ (mm)", fontsize=16)
    ax6.set_ylabel("$w_{3}$ (mm)", fontsize=16)
    ax1.tick_params(labelsize=12)
    ax2.tick_params(labelsize=12)
    ax3.tick_params(labelsize=12)
    ax4.tick_params(labelsize=12)
    ax5.tick_params(labelsize=12)
    ax6.tick_params(labelsize=12)
    ax1.axes.xaxis.set_ticklabels([])
    ax2.axes.xaxis.set_ticklabels([])
    ax3.axes.xaxis.set_ticklabels([])
    ax4.axes.xaxis.set_ticklabels([])
    ax5.axes.xaxis.set_ticklabels([])
    ax6.set_xlabel("Site", fontsize=16)

    time_text = ax0.text(0.88, 1.2, '', fontsize=14, transform=ax0.transAxes)
    mass, = ax0.plot([], [], 'bo', ms=5)
    site = np.arange(1,N+1)
    line1, = ax1.plot([], [], 'o--', ms=5, lw=2, c='black')
    line2, = ax2.plot([], [], 'o--', ms=5, lw=2, c='black')
    line3, = ax3.plot([], [], 'o--', ms=5, lw=2, c='black')
    line4, = ax4.plot([], [], 'o--', ms=5, lw=2, c='black')
    line5, = ax5.plot([], [], 'o--', ms=5, lw=2, c='black')
    line6, = ax6.plot([], [], 'o--', ms=5, lw=2, c='black')

    xOrig = []
    yOrig = []
    for n in range(N):
        xOrig = np.append(xOrig, n*L[0]-L[3])
        xOrig = np.append(xOrig, n*L[0])
        xOrig = np.append(xOrig, n*L[0])
        yOrig = np.append(yOrig, L[2])
        yOrig = np.append(yOrig, L[1]+L[2])
        yOrig = np.append(yOrig, 0.)


    ani = ani.FuncAnimation(fig, animate, init_func=init, frames=frame_range, fargs=(site, fskip), interval=dt_frame, blit=True)


ani.save(f"{filename[:-3]}.mp4", fps=20, extra_args=['-vcodec', 'libx264'])
print(f" >> The video saved as \"{filename[:-3]}.mp4\" in ./outputs folder")

