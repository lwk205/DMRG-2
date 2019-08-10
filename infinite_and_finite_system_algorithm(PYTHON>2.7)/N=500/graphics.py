import numpy as np
import matplotlib.pyplot as plt

def single_graphic_infinite():
	
	for m in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:
		data = np.loadtxt("datam="+str(m)+".txt", dtype=float)
		L = np.loadtxt("longm="+str(m)+".txt", dtype=float)

		E = data[:,0]
		e = data[:,2]
		S = data[:,1]
	
		E_site = np.array([E[i]/L[i] for i in range(len(L)) if i%2 != 0])
		in_L = np.array([1/L[i] for i in range(len(L)) if i%2 != 0])
		e_p = np.array([e[i]/L[i] for i in range(len(L)) if i%2 != 0])
		plt.plot(in_L, e_p, marker="o", markersize=5)

		#k = ajuste(in_L,E_site)
		#plt.plot(in_L, E_site, "ko", markersize=5, label="iDMRG")
		#plt.plot(in_L, k[0]*in_L + k[1], color="r", label="linear fit")
		plt.xlabel(r"$1/L$", fontsize=25)
		plt.ylabel(r"$e/L$", fontsize=25)
		plt.grid(True, linestyle="--", color="0.5")
		plt.xticks(size=20)
		plt.yticks(size=20)
		plt.legend(fontsize=25)
		#plt.xlim(xmin=0, xmax=500)
	plt.show()

def ajuste(x,y):
	base_l = np.polyfit(x, y, deg=1)
	return base_l


def main():
	single_graphic_infinite()

if __name__ == '__main__':
	main()
