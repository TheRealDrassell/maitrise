import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

from tdmsRead import extract
from main import derv4

def test() :
    fichier = "../alcheData/2025-06-06_JN_Gquadexp1/DNA_10mKCl/DNA_KCl_10mM_SR1e4_500kOmhsVg0_1_long.tdms"
    sr = int(1e4)

    dataIter, nbrData = extract(fichier, sr, 3, offset=100, chunks=True)
    dataZip = list(dataIter)[0]

    #next(dataZip)

    #comparaison par std dev, si moy selon une certaine distance de std dev, changement
    # si std dev explose, chnagement d'état

    #next(dataZip)

    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex=ax)

    j = 0
    for data, time in dataZip :
        j+=1
        print(j)

        fenetre = 5
        nbr = int(data.size/fenetre)
        #print(nbr)
        #print(data.size)

        dd = data.reshape(nbr, fenetre)
        tt = time.reshape(nbr, fenetre)

        moy = np.mean(dd, axis=1)
        std = np.std(dd, axis=1)
        var = np.var(dd, axis=1)

        varmoy = np.mean(var)
        varstd = np.std(var)

        peaks = np.where(var > ( varmoy + varstd*6 ))[0]
        indxTime = peaks * fenetre

        ax.plot(time, data, color="tab:blue")
        #for i in indxTime :
        #    ax.axvline(time[i], ls=":", color="k")
        #ax.plot(time[::fenetre], moy.ravel(), color="tab:orange")
        #ax.plot(time[::fenetre], moy.ravel(), color="k", ls="", marker="o", markersize=1)
        #ax2.plot(time[::fenetre], std.ravel(), color="tab:blue", alpha=0.5)
        ax2.plot(time[::fenetre], var.ravel(), color="tab:orange")
        #ax2.plot(np.array([time[0], time[-1]]),np.ones(2)*varmoy, color="k" )
        #ax2.plot(np.array([time[0], time[-1]]),np.ones(2)*(varmoy+varstd*6), color="red", ls=":" )
        #ax2.axhline(varmoy, color="k")
        #ax2.axhline(varstd + varmoy, ls=":", color="red")
        #ax.errorbar(time[::fenetre], moy.ravel(), yerr=std, color="k", ls="", marker="o", markersize=1)
        #ax.errorbar(time[::10], moy.ravel(), yerr=std,color="k")

        #for i in range(nbr) :
            #ax.plot(tt[i], dd[i])
            #ax.plot(tt[i], np.ones(10)*moy[i], color="k")

        #break

    plt.show()

def test2() :
    fichier = "../alcheData/2025-06-06_JN_Gquadexp1/DNA_10mKCl/DNA_KCl_10mM_SR1e4_500kOmhsVg0_1_long.tdms"
    sr = int(1e4)

    dataIter, nbrData = extract(fichier, sr, 3, offset=100, chunks=False)
    dataZip = list(dataIter)[0]

    windowS = 20

    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex=ax)
    #ax3 = fig.add_subplot(413, sharex=ax)
    #ax4 = fig.add_subplot(414, sharex=ax)

    for data, time in dataZip :

        derv = derv4(data)
        dervTime = time[2:-2]

        rollingData = np.lib.stride_tricks.sliding_window_view(data, windowS)
        rollingTime = np.lib.stride_tricks.sliding_window_view(time, windowS)

        rollingDerv = np.lib.stride_tricks.sliding_window_view(np.abs(derv)*3, windowS)
        rollingDervTime = np.lib.stride_tricks.sliding_window_view(dervTime, windowS)

        rollingVar = np.var(rollingData, axis=1)
        rollingVarDerv = np.var(rollingDerv, axis=1)

        moyVar = np.mean(rollingVar)
        stdVar = np.std(rollingVar)

        ax.plot(time, data)
        ax2.plot(time[:-(windowS-1)], rollingVar)
        #ax3.plot(dervTime, derv)
        #ax4.plot(dervTime[windowS -1:], rollingVarDerv)

        ax2.axhline(moyVar, color="k")
        ax2.axhline(moyVar+stdVar*6, ls=":", color="k")

        #eps = 0.85e-15
        eps = 1e-15
        #eps = moyVar + stdVar*6
        peaks = np.where(rollingVar > eps)
        peaks2, _ = find_peaks(rollingVar, height=eps, width=2)
        ax2.axhline(eps, ls=":", color="red")
        #print(peaks2.size)

        #upDown = np.arange(peaks2.size)%2

        #gauche = time[windowS -1:][peaks2][:-1]  # aye....
        #droite = time[windowS -1:][peaks2][1:] 

        #upda = np.vstack((upDown, upDown)).T.ravel()[:-1]
        #upTime = np.vstack( (gauche, droite) ).T.ravel()
        #upTime = np.append( upTime, droite[-1] )

        #print(upTime[::-1])
        #print(peaks2)

        ax2.plot(time[:-( windowS-1 )][peaks2], rollingVar[peaks2], ls="", marker="x", markersize=1, color="red")
        #ax2.plot(time[windowS-1:][peaks2], upDown)
        #ax2.plot(upTime, upda)
        ax.vlines(time[:-( windowS-1 )][peaks2], 6e-7, 9.5e-7, ls=":", color="k")

        #break


    fig.tight_layout()
    plt.show()

def test3() :
    fichier = "../alcheData/2025-06-06_JN_Gquadexp1/DNA_10mKCl/DNA_KCl_10mM_SR1e4_500kOmhsVg0_1_long.tdms"
    #fichier = "../alcheData/2025-07-29_JN_Gquadexp4/10mM_KCl/10mMKCl_SR1e4_4h.tdms"
    sr = int(1e4)

    dataIter, nbrData = extract(fichier, sr, [3, 1], offset=100, chunks=False)
    dataZip = list(dataIter)

    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex=ax)

    data, time = next(dataZip[0])
    ax.plot(time, data)

    data, time = next(dataZip[1])
    ax2.plot(time, data)

    plt.show()

if __name__ == "__main__" :
    test2()

