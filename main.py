import numpy as np
import matplotlib.pyplot as plt

from tdmsRead import extract
from amplificateur import amplificateur, derv4, rollingVar, npfft, stateChange

from scipy.signal import find_peaks

def main() :
    fichier = "../alcheData/2025-06-06_JN_Gquadexp1/DNA_10mKCl/DNA_KCl_10mM_SR1e4_500kOmhsVg0_1_long.tdms"
    sr = int(1e4)

    dataIter, nbrData = extract(fichier, sr, 3, offset=100, chunks = False)
    dataZip = list(dataIter)[0]

    fig = plt.figure()
    ax = fig.add_subplot(311)
    axAmp = fig.add_subplot(312, sharex=ax)
    axAmp2 = fig.add_subplot(313, sharex=ax)
    #axHU = fig.add_subplot(325)
    #axHD = fig.add_subplot(326)

    ales = 0
    iter = 0
    for data, time in dataZip :
        iter += 1

        windowS = 5
        recAmp, indx = amplificateur(data, sr, windowS)
        rra, indx2 = rollingVar(recAmp, windowS)
        indx2 += indx

        eps = 5e-38
        #eps = np.mean(rra) + np.std(rra) * 10
        peaks, indx3 = stateChange(recAmp, windowS, eps)
        indx3 += indx

        #print(peaks.size)
        #print(eps + np.std(rra))
        ales += peaks.size
        print(iter, ales, peaks.size)

        up = []
        down = []
        indxLast = 0
        stateTime = time[4:-4]

        for j, i in enumerate(peaks) :
            dtime = stateTime[i] - stateTime[indxLast]
            indxLast = i
            
            if j%2 :
                up.append(dtime)
            else :
                down.append(dtime)

        up = np.array(up)
        down = np.array(down)
        #recFft = npfft(data, time, cutoff=15000)

        ax.plot(time, data, color="tab:blue")
        axAmp.plot(time[indx[0]:indx[1]],recAmp, color="tab:blue")
        axAmp2.plot(time[indx2[0]:indx2[1]], rra, color="tab:blue")
        axAmp2.plot(time[indx2[0]:indx2[1]][peaks], rra[peaks],ls="", marker="x", markersize=1, color="red")

        ax.vlines(time[indx2[0]:indx2[1]][peaks], data.min(), data.max(), ls=":", color="k")
        axAmp.vlines(time[indx2[0]:indx2[1]][peaks], recAmp.min(), recAmp.max(), ls=":", color="k")
        axAmp2.axhline(eps, ls=":", color="k")

        #rVar, indx = rollingVar(derv4(data), windowS)
        #ampDD = np.cumsum(derv4(data)[2:-2] * rVar )
        #axAmp2.plot(time[4:-4], ampDD)

        """
        histUp, upEdges = np.histogram(up, bins=50, range=[up.min(), 1])
        histDown, downEdges = np.histogram(down, bins=50, range=[down.min(), 1])

        zerosU = (histUp==0)
        zerosD = (histDown==0)

        axHU.plot(upEdges[1:][~zerosU],histUp[~zerosU], marker="o", markersize=2)
        axHD.plot(downEdges[1:][~zerosD],histDown[~zerosD], marker="o", markersize=2)

        axHU.semilogy()
        axHD.semilogy()
        """


    fig.tight_layout()
    plt.show()

if __name__ == "__main__" :
    main()
