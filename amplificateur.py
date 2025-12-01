import numpy as np
import matplotlib.pyplot as plt
from tdmsRead import extract

def test() :
    fichier = "../alcheData/2025-06-06_JN_Gquadexp1/DNA_10mKCl/DNA_KCl_10mM_SR1e4_500kOmhsVg0_1_long.tdms"
    sr = int(1e4)

    dataIter, nbrData = extract(fichier, sr, 3, chunks=True)
    dataZip = list(dataIter)[0]

    nbr = int(1e4/10)

    #comparaison par std dev, si moy selon une certaine distance de std dev, changement
    # si std dev explose, chnagement d'état


    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    for data, time in dataZip :
        dd = data.reshape(nbr, 10)
        tt = time.reshape(nbr, 10)

        moy = np.mean(dd, axis=1)
        std = np.std(dd, axis=1)
        var = np.std(dd, axis=1)

        ax.plot(time, data)
        #ax.plot(time[::10], moy.ravel(), color="k", ls="", marker="o", markersize=1)
        ax2.plot(time[::10], std.ravel(), color="k")
        ax.errorbar(time[::10], moy.ravel(), yerr=std, color="k", ls="", marker="o", markersize=1)
        #ax.errorbar(time[::10], moy.ravel(), yerr=std,color="k")

        #for i in range(nbr) :
            #ax.plot(tt[i], dd[i])
            #ax.plot(tt[i], np.ones(10)*moy[i], color="k")

        #break

    plt.show()


if __name__ == "__main__" :
    test()

