import numpy as np
import matplotlib.pyplot as plt
from tdmsRead import extract

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
    ax = fig.add_subplot()
    #ax2 = fig.add_subplot(212)

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
        for i in indxTime :
            ax.axvline(time[i], ls=":", color="k")
        ax.plot(time[::fenetre], moy.ravel(), color="tab:orange")
        #ax.plot(time[::fenetre], moy.ravel(), color="k", ls="", marker="o", markersize=1)
        #ax2.plot(time[::fenetre], std.ravel(), color="tab:blue", alpha=0.5)
        #ax2.plot(time[::fenetre], var.ravel(), color="tab:orange")
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


if __name__ == "__main__" :
    test()

