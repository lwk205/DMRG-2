import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
plt.style.use("seaborn-darkgrid")


def single_graphic_infinite():
	kas=[]
	mes =[50, 60, 70, 80, 90, 100]
	for m in mes:
		data = np.loadtxt("datam="+str(m)+".txt", dtype=float)
		L = np.loadtxt("longm="+str(m)+".txt", dtype=float)

		E = data[:,0]
		e = data[:,2]
		S = data[:,1]
	
		E_site = np.array([E[i]/L[i] for i in range(len(L)) if i%2 != 0])
		in_L = np.array([1/L[i] for i in range(len(L)) if i%2 != 0])
		#e_p = np.array([e[i] for i in range(len(L)) if i%2 != 0])
		#plt.plot(1/in_L, E_site, marker="o", markersize=5, label="m="+str(m))

		k = ajuste(in_L,E_site)
		kas.append(k[1])
		
		plt.plot(in_L[5:], E_site[5:], "ko", markersize=5, label="iDMRG")
		plt.plot(in_L, k[0]*in_L + k[1], color="r", label="linear fit")
		plt.ylabel(r"$(E/L)_{0}$  ", fontsize=25)
		plt.xlabel(r"$m$", fontsize=25)
		#plt.grid(True, linestyle="--", color="0.5")
		plt.xticks([0.005,0.015,0.025,0.035],size=20)
		plt.yticks(size=20)
		#plt.legend(fontsize=25)
		#plt.xlim(xmin=0, xmax=500)
	plt.tight_layout ()
	plt.show()
	"""plt.plot(mes, kas, marker="o", markersize=7, linestyle="-",color="k")
	plt.ylabel(r"$(E/L)_{0}$", fontsize=25)
	plt.xlabel(r"$m$", fontsize=25)
	plt.grid(True, linestyle="--", color="0.5")
	plt.xticks(size=20)
	plt.yticks(size=20)
	plt.show()"""
	

def ajuste(x,y):
	base_l = np.polyfit(x, y, deg=1)
	return base_l


def main():
	single_graphic_infinite()

if __name__ == '__main__':
	main()
