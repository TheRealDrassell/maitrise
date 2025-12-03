import numpy as np
import matplotlib.pyplot as plt

from tdmsRead import extract

from scipy.signal import find_peaks

def derv4(y, h = 1) :
    return -(-y[:-4] + 8*y[1:-3] - 8 * y[3:-1] + y[4:])/(12*h)

def rollingVar(data, windowS) :

    rollingData = np.lib.stride_tricks.sliding_window_view(data, windowS)
    rVar = np.var(rollingData, axis=1)

    return rVar

def main() :
    fichier = "../alcheData/2025-06-06_JN_Gquadexp1/DNA_10mKCl/DNA_KCl_10mM_SR1e4_500kOmhsVg0_1_long.tdms"
    sr = int(1e4)

    dataIter, nbrData = extract(fichier, sr, 3, offset=100, chunks = True)
    dataZip = list(dataIter)[0]

    fig = plt.figure()
    ax = fig.add_subplot(211)
    axVar = fig.add_subplot(212, sharex=ax)
    for data, time in dataZip :

        ax.plot(time, data, color="tab:blue")

        windowS = 5
        rVar = rollingVar(data, windowS)

        indxG = np.floor( (windowS - 1)/2 ).astype(int)
        indxD = np.ceil( (windowS - 1)/2 ).astype(int)

        #axVar.plot(time[indxG:-indxD], rVar, color="tab:blue")

        #### FFT  #####
        sp = np.fft.rfft(data)
        freq = np.fft.rfftfreq(time.size, d = 1/time.size)
        ampl = np.abs(sp)
        masque = freq < 700     # optimisaton machine, trouver le bon cutoff, déscente de gradient cutoff le plus bas pour le meilleur signal
        sp_masque = sp * masque
        req = np.fft.irfft(sp_masque)

        #### FFT   ######

        reqR = rollingVar(req, windowS)

        eps = 1e-16
        peaks2, _ = find_peaks(reqR, height=eps, width=2)

        axVar.plot(time[indxG:-indxD], reqR, color="tab:orange")
        axVar.plot(time[indxG:-indxD][peaks2], reqR[peaks2], ls="", marker="x", markersize=1, color="red")
        axVar.axhline(eps, ls=":", color="k")

        ax.vlines(time[indxG:-indxD][peaks2], 6e-7, 9.5e-7, ls=":", color="k")

        ax.plot(time, req, color="tab:orange")
        #axVar.plot(freq, ampl)

        windowS = 30
        rrVar = np.lib.stride_tricks.sliding_window_view(rVar,windowS)
        rStdVar = np.std(rrVar, axis=1)
        rMoyVar = np.mean(rVar)

        indxG2 = np.floor( (windowS - 1)/2 ).astype(int)
        indxD2 = np.ceil( (windowS - 1)/2 ).astype(int)

        #axVar.plot( time[indxG:-indxD][indxG2:-indxD2], rMoyVar + rStdVar, ls=":", color="red" )

    plt.show()

if __name__ == "__main__" :
    main()
