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

def amplificateur(data, windowS) :

    rVar, indx = rollingVar(data, windowS)
    #rDVar, indx = rollingVar(derv4(data), windowS)

    # windowS > 5, sinon indexer rollingVar
    ampData = derv4(data) * rVar  #[ indx[0]-2: -(indx[1] - 2 + 1) ]   # et si variance de la derv?
    #ampData = derv4(data)[2:-2] * rDVar#[2:-2]  # et si variance de la derv?
    reconstruction = np.cumsum(ampData)

    return reconstruction, indx

def stateChange(data, window, eps = 5.5e-46) :
    rData, indx = rollingVar(data, window)
    peaks, _ = find_peaks(rData, height=eps)

    return peaks, indx

def fourier(data, time, cutoff = 700) :
        sp = np.fft.rfft(data)
        freq = np.fft.rfftfreq(time.size, d = 1/time.size)
        ampl = np.abs(sp)
        masque = freq < cutoff     # optimisaton machine, trouver le bon cutoff, déscente de gradient cutoff le plus bas pour le meilleur signal
        sp_masque = sp * masque
        req = np.fft.irfft(sp_masque)

        return req

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

    dataIter, nbrData = extract(fichier, sr, 3, offset=100, chunks=False)
    dataZip = list(dataIter)[0]

    fig = plt.figure()
    ax = fig.add_subplot(311)
    ax2 = fig.add_subplot(312,sharex=ax)
    ax3 = fig.add_subplot(313,sharex=ax)
    for data, time in dataZip :

        #ax.plot(time, data)

        windowS = 5
        rollingData = np.lib.stride_tricks.sliding_window_view(data, windowS)
        rollingVar = np.var(rollingData, axis=1)

        indxG = np.floor( (windowS-1)/2 ).astype(int)
        indxD = np.ceil( (windowS-1)/2 ).astype(int)

        #ampData = data[indxG:-indxD] * rollingVar
        #derv = derv4(ampData)
        #cs = np.cumsum(derv)
        ampData = derv4(data) * rollingVar
        cs = np.cumsum(ampData)

        ax.plot(time, data)
        ax2.plot(time[indxG:-(indxD)], cs)

        rollingCS = np.lib.stride_tricks.sliding_window_view(cs, windowS)
        rollingCSVAR = np.var(rollingCS, axis=1)

        ax3.plot(time[indxG*2:-2*indxD], rollingCSVAR)

        eps = 5.5e-46
        peaks, _ = find_peaks(rollingCSVAR, height=eps)
        ax3.plot(time[indxG*2:-2*indxD][peaks], rollingCSVAR[peaks], ls="", color="red", marker="x", markersize=1)
        print(peaks.size)
        ax2.vlines(time[indxG*2:-(indxD*2)][peaks], np.min(cs), np.max(cs), ls=":", color="k")
        ax.vlines(time[indxG*2:-(indxD*2)][peaks], 6e-7, 9.5e-7, ls=":", color="k")
        
        #ax.plot(time[2:-2],derv4(data))

    #ax.set_title("normal")
    #ax2.set_title("«amplifié»")

    fig.tight_layout()
    plt.show()

def test4() :
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
    test4()

