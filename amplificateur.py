import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

from tdmsRead import extract

def derv4(y, h = 1) :
    return -(-y[:-4] + 8*y[1:-3] - 8 * y[3:-1] + y[4:])/(12*h)

def rollingVar(data, windowS) : 
    rollingData = np.lib.stride_tricks.sliding_window_view(data, windowS)
    rVar = np.var(rollingData, axis=1)

    indxG = np.floor( (windowS-1)/2 ).astype(int)
    indxD = -np.ceil( (windowS-1)/2 ).astype(int)

    return rVar, np.array([indxG, indxD])

def amplificateur(data, sr, windowS) :
    # derv wise ou none

    rVar, indx = rollingVar(data, windowS)
    derv = derv4(data, 1/sr)

    if windowS == 5:
        ampData = derv * rVar
    elif windowS == 3 :
        ampData = derv * rVar[1:-1]
        indx = np.array([2,-2])
    else :
        ampData = derv[indx[0]-2:indx[1]+2] * rVar

    reconstruction = np.cumsum(ampData)
    #reconstruction = ampData

    return reconstruction, indx

def stateChange(data, window, eps = 5.5e-46) :
    rData, indx = rollingVar(data, window)
    peaks, _ = find_peaks(rData, height=eps)

    return peaks, indx

def npfft(data, time, cutoff = None) :
        sp = np.fft.rfft(data)
        freq = np.fft.rfftfreq(time.size, d = 1/time.size)
        ampl = np.abs(sp)
        if cutoff != None :
            masque = freq < cutoff     # optimisaton machine, trouver le bon cutoff, déscente de gradient cutoff le plus bas pour le meilleur signal
            sp *= masque

        req = np.fft.irfft(sp)

        return req

def test() :
    fichier = "../alcheData/2025-06-06_JN_Gquadexp1/DNA_10mKCl/DNA_KCl_10mM_SR1e4_500kOmhsVg0_1_long.tdms"
    #fichier = "../alcheData/2025-07-29_JN_Gquadexp4/10mM_KCl/10mMKCl_SR1e4_4h.tdms"
    sr = int(1e4)

    #rD, indx = rollingVar(derv4(data), 5)
    #ddd = np.cumsum(derv4(data)[2:-2] * rD)

    dataIter, nbrData = extract(fichier, sr, 3, offset=100, chunks=False)
    dataZip = list(dataIter)[0]

    windowS = 5

    fig = plt.figure()
    ax = fig.add_subplot(311)
    ax2 = fig.add_subplot(312, sharex=ax)
    ax3 = fig.add_subplot(313, sharex=ax)

    figHist = plt.figure()
    axHistU = figHist.add_subplot(211)
    axHistD = figHist.add_subplot(212)

    for data, time in dataZip :

        dataAmp, indx = amplificateur(data, windowS)
        rDataAmpVar, indx2 = rollingVar(dataAmp, windowS)
        indx2 += indx

        eps = 5.5e-46
        peaks, indx3 = stateChange(dataAmp, windowS, eps)
        indx3 += indx

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

        #axHistU.plot(np.histogram(up, bins=50)[0], marker="o")
        #axHistD.plot(np.histogram(down, bins=50)[0], marker="o")
        #axHistU.hist(up, bins=30)
        #axHistD.hist(down, bins=30)

        histUp, upEdges = np.histogram(up, bins=50, range=[up.min(), 3])
        histDown, downEdges = np.histogram(down, bins=50, range=[down.min(), 3])

        zerosU = (histUp==0)
        zerosD = (histDown==0)

        axHistU.plot(upEdges[1:][~zerosU],histUp[~zerosU], marker="o", markersize=2)
        axHistD.plot(downEdges[1:][~zerosD],histDown[~zerosD], marker="o", markersize=2)

        #axHistU.plot(upEdges[1:],histUp)
        #axHistD.plot(downEdges[1:],histDown)
        #axHistU.plot(histUp[~zerosU])
        #axHistU.plot(histUp)

        #axHistU.plot(np.sort(up)[::-1])
        #axHistD.plot(np.sort(down)[::-1])
        axHistU.semilogy()
        axHistD.semilogy()

        #axHistU.plot(up, marker="o", markersize=1)
        #axHistD.plot(down, marker="o", markersize=1)


        ax.plot(time, data)
        ax.vlines(time[indx3[0]:indx3[1]][peaks], 6e-7, 9.5e-7, ls=":", color="k")

        ax2.plot(time[indx[0]:indx[1]], dataAmp, color="tab:blue")
        ax2.vlines(time[indx3[0]:indx3[1]][peaks], np.min(dataAmp), np.max(dataAmp), ls=":", color="k")

        ax3.plot(time[indx2[0]:indx2[1]], rDataAmpVar, color="tab:blue")
        ax3.plot(time[indx3[0]:indx3[1]][peaks], rDataAmpVar[peaks], ls="", marker="x", markersize=1, color="red")
        ax3.axhline(eps, ls=":", color="k")

    plt.show()

if __name__ == "__main__" :
    test()

