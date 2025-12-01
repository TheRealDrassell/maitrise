import numpy as np
import matplotlib.pyplot as plt

from tdmsRead import extract

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

    dataIter, nbrData = extract(fichier, sr, 3, chunks = True)
    dataZip = list(dataIter)[0]

    #data = data[100:]
    #time = time[100:]

    fig = plt.figure()
    ax = fig.add_subplot()
    for data, time in dataZip :

        moy = np.mean(data)
        derv = derv4(data) * 3

        rollingMean = np.cumsum(data) / np.arange(1, data.size + 1)

        #ax.plot(time[2:-2], derv + moy, color="tab:blue")
        ax.plot(time, data, color="k")
        ax.plot(np.array( [time[0], time[-1]] ), np.ones(2)*moy, color="tab:orange")
        ax.plot(time, rollingMean, color="tab:blue")


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
