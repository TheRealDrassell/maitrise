import numpy as np
import matplotlib.pyplot as plt

from tdmsRead import extract

from scipy.signal import find_peaks

def derv4(y, h = 1) :
    return (-y[:-4] + 8*y[1:-3] - 8 * y[3:-1] + y[4:])/(12*h)

def dataWrap(dataIter, sr) :
    for i in dataIter :
        data = i[:]
        nbrPoint = data.size
        fin = nbrPoint/sr
        yield data, np.linspace(0, fin, nbrPoint)

def main() :
    fichier = "../alcheData/2025-06-06_JN_Gquadexp1/DNA_10mKCl/DNA_KCl_10mM_SR1e4_500kOmhsVg0_1_long.tdms"
    sr = int(1e4)

    dataIter, nbrData = extract(fichier, sr, 3, offset=100, chunks = False)
    dataZip = list(dataIter)[0]

    #data = data[100:]
    #time = time[100:]

    fig = plt.figure()
    ax = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)
    for data, time in dataZip :

        rolling = np.lib.stride_tricks.sliding_window_view(data, 11)
        var = np.var(rolling, axis=1)

        moyVar = np.mean(var)
        stdVar = np.std(var)

        peaks = np.where(var>(moyVar+stdVar*11))

        derv = derv4(data) 
        dervTime = time[2:-2]
        varTime = time[5:-5]

        ax.plot(time, data)
        ax.vlines(dervTime[peaks], 6.5e-7, 9.5e-7, ls=":", color="k")

        ax2.plot(dervTime, derv)

        ax3.plot(varTime, var)
        ax3.plot(varTime[peaks], var[peaks], ls="", marker="o", markersize=1, color="red")
        ax3.axhline(moyVar, ls=":", color="k")
        ax3.axhline(moyVar+stdVar*10, ls=":", color="k")

        """
        moyDerv = np.mean(derv)
        stdDerv = np.std(derv)

        peaks = np.where(np.abs(derv) > ( moyDerv + stdDerv*6 ))

        ax.plot(time, data)
        ax.vlines(dervTime[peaks], 6.5e-7, 9.5e-7, ls=":", color="k")

        ax2.plot(dervTime, derv)
        ax2.plot(dervTime[peaks], derv[peaks], ls="", marker="o", markersize=1, color="red")
        """



    plt.show()

    """
    i = 0
    for data, time in dataZip[1] :
        i+=1

        data = data[100:] 
        time = time[100:]

        derv = derv4(data, 1) * 3
        moy = np.mean(data)

        #ax.plot(time[2:-2], (derv + moy))
        ax.plot(time, data)

    #fig.legend()
    plt.show()
    """

if __name__ == "__main__" :
    main()
