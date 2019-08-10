import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.mlab as mlab

# infinite system algorithm graphics
def energy_entropy(name_data, name_long):
    data = np.loadtxt(name_data, dtype=float)
    y_1 = data[:, 0]
    y_2 = data[:, 1]

    x = np.loadtxt(name_long, dtype=float)

    plt.plot(x[::2], y_1[::2]/x[::2], marker=".", color="k")
    plt.xlabel(r"$L$", fontsize=25)
    plt.ylabel(r"$E/L$", fontsize=25)
    plt.xticks(size=20)
    plt.yticks(size=20)
    plt.xlim(xmin=0, xmax=500)
    plt.tight_layout()
    plt.show()

    plt.plot(x[::2][5:], y_2[::2][5:], marker=".", color="k")
    plt.xlabel(r"$L$", fontsize=25)
    plt.ylabel(r"$S$", fontsize=25)
    plt.xticks(size=20)
    plt.yticks(size=20)
    plt.xlim(xmin=50, xmax=500)
    plt.tight_layout()
    plt.show()


def cge_energy():
    values = [50, 60, 70, 80, 90, 100]
    k_s = []
    for m in values:
        data = np.loadtxt("datam=" + str(m) + ".txt", dtype=float)
        y_1 = data[:, 0]
        x = np.loadtxt("longm=" + str(m) + ".txt", dtype=float)

        k = fitting((1/x[::2])[10:], (y_1[::2]/x[::2])[10:])
        k_s.append(k[1] - 0.25 + np.log(2))

        if m == 100:
            plt.plot((1/x[::2])[10:], (y_1[::2]/x[::2])[10:], "ko", markersize=5, label="DMRG")
            plt.plot((1/x[::2])[10:], k[0] * (1/x[::2])[10:] + k[1], color="k", label="Fitting")
            plt.legend(fontsize=18)
            plt.xlabel(r"$1/L$", fontsize=25)
            plt.ylabel(r"$E/L$", fontsize=25)
            plt.xticks(size=20)
            plt.yticks(size=20)
            plt.tight_layout()
            plt.show()
            print(k)

    plt.plot(values, (10**5) * np.array(k_s), marker="o", markersize=5, color="k")
    plt.xlabel(r"$m$", fontsize=25)
    plt.ylabel(r"$(E/L)_{dmrg} - (E/L)_{0}$ $[10^{-5}]$", fontsize=25)
    plt.xticks(size=20)
    plt.yticks(size=20)
    plt.xlim(xmin=50, xmax=100)
    plt.tight_layout()
    plt.show()


def fitting(x, y):
    line = np.polyfit(x, y, deg=1)
    return line


def error():
    values = [50, 60, 70, 80, 90, 100]
    for m in values:
        data = np.loadtxt("datam=" + str(m) + ".txt", dtype=float)
        y_1 = data[:, 2]
        x = np.loadtxt("longm=" + str(m) + ".txt", dtype=float)

        plt.plot(x[::2], (10**7) * y_1[::2], label="m = " + str(m))
        plt.xlabel(r"$L$", fontsize=25)
        plt.ylabel(r"$error$ $[10^{-7}]$", fontsize=25)
        plt.xticks(size=20)
        plt.yticks(size=20)
        plt.xlim(xmin=0, xmax=500)
        plt.legend(fontsize=18)
    plt.tight_layout()
    plt.show()


# finite system algorithm of dmrg
def energy_finite(name_data):
    data = np.loadtxt(name_data, dtype=float)
    energy = []

    for i in range(10):
        energy.append(np.array(data[i*1000:(i+1)*1000][::2])-0.25 + np.log(2))
        if i > 5:
            plt.plot((10**4) * energy[i], label="m = " + str((i+1)*10))
            plt.xlabel(r"$n$ - paso", fontsize=25)
            plt.ylabel(r"$(E/L)_{dmrg} - (E/L)_{0}$ $[10^{-4}]$", fontsize=25)
            plt.xticks(size=20)
            plt.yticks(size=20)
            plt.xlim(xmin=0, xmax=500)
            plt.ylim(ymin=3.769, ymax=3.775)
            plt.legend(fontsize=18)
    plt.tight_layout()
    plt.show()


def error_finite(name_data):
    data = np.loadtxt(name_data, dtype=float)
    energy = []

    for i in range(10):
        energy.append(np.array(data[i * 1000:(i + 1) * 1000][::2]))
        if i > 5:
            plt.plot((10 ** 4) * energy[i], label="m = " + str((i + 1) * 10))
            plt.xlabel(r"Iteration", fontsize=25)
            plt.ylabel(r"$error$", fontsize=25)
            plt.xticks(size=20)
            plt.yticks(size=20)
            plt.xlim(xmin=0, xmax=500)
            plt.legend(fontsize=18)
    plt.tight_layout()
    plt.show()


# Spin symmetry infinite system
def spin_infinite():
	
    for k in [0, 1]:
        y = np.loadtxt("spinlen_" + str(k), dtype=float)
        x = np.loadtxt("spin_" + str(k) + "_infi", dtype=float)

        plt.plot(x[::2], y[::2], ".", label=r"$S_z = $" + str(k))
        plt.xlabel(r"$L$", fontsize=25)
        plt.ylabel(r"$E/L$", fontsize=25)
        plt.legend(fontsize=18)
        plt.xticks(size=20)
        plt.yticks(size=20)
        plt.xlim(xmin=0, xmax=500)
        #plt.ylim(ymin=-0.441, ymax=-0.443)
    plt.tight_layout()
    plt.show()



def exac_energy():
    data_1 = []
    data_2 = []
    data_3 = []
    for k in range(2, 13, 2):
        y = np.loadtxt(str(k) + ".dat", dtype=float)
        data_1.append(y[:,3][0])
        data_2.append(y[:,3][1])
        data_3.append(y[:,3][2])
    plt.plot(np.arange(2, 13, 2), data_1/np.arange(2, 13, 2), "-^", label="0ST")
    plt.plot(np.arange(2, 13, 2), data_2/np.arange(2, 13, 2), "-s", label="1ST")
    plt.plot(np.arange(2, 13, 2), data_3/np.arange(2, 13, 2), "-o", label="2ST")
    plt.xlabel(r"$N$", fontsize=25)	
    plt.ylabel(r"$S_z$", fontsize=25)
    plt.xticks(size=20)
    plt.yticks(size=20)
    plt.legend(fontsize=18)
    #plt.xlim(xmin=0, xmax=500)
    plt.tight_layout()
    plt.show()


def histogram():

    date = np.loadtxt("12.dat", dtype=float)
    date_1 = date[:, 1]
    media_1 = date_1.mean()

    sigma_1 = date_1.std()


    num_bins = 100
    
    # the histogram of the data
    n, bins_1, patches = plt.hist(date_1, num_bins, normed=1, facecolor="blue", alpha=1.0, edgecolor="black")

    # add a 'best fit' line
    y_1 = mlab.normpdf(bins_1, media_1, sigma_1)
   
   
    plt.ylabel('Degeneration', fontsize=25)
    plt.xlabel(r'$E$', fontsize=25)
    plt.xticks(size=20)
    plt.yticks(size=20)
    #plt.grid(True, linestyle="--", color="0.5")
    
    # Tweak spacing to prevent clipping of ylabel
    plt.subplots_adjust(left=0.15)
    plt.tight_layout()
    plt.show()


def main():
    rc("font", **{"family": "sans-serif", "sans-serif": ["Helvetica"]})
    rc("text", usetex=True)

    plt.style.use("seaborn-darkgrid")


    #energy_entropy("datam=100.txt", "longm=100.txt")
    #cge_energy()
    #error()

    #energy_finite("Energy")
    #error_finite("error")

    #spin_infinite()
    #exac_energy()
    histogram()

if __name__ == '__main__':
    main()
