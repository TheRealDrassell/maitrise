import numpy as np
import matplotlib.pyplot as plt
import nptdms

def gen(chunks, sr) :
    print(chunks)
    for j,i in enumerate(chunks) :
        nbrPoint = len(i)
        timeStamp = nbrPoint/sr

        time = np.linspace( j*timeStamp, (j+1)*timeStamp, nbrPoint )
        yield i[:], time

def test() :
    fichier = "../alcheData/2025-06-06_JN_Gquadexp1/DNA_10mKCl/DNA_KCl_10mM_SR1e4_500kOmhsVg0_1_long.tdms"
    #fichier = "../alcheData/2025-06-06_JN_Gquadexp1/DNA_10mKCl/DNA_KCl_10mM_SR1e4_500kOmhsVg0_1_court.tdms"

    sr = int(1e4)

    f1 = nptdms.TdmsFile.open(fichier)
    g1 = f1.groups()[0]
    c = g1.channels()

    cc = [i.data_chunks() for i in c]

    i = cc[3]
    print(i)
    dataZipe = gen([i], sr)

    fig = plt.figure()
    ax = fig.add_subplot()
    for data, time in dataZipe :
        ax.plot(time, data)
    
    plt.show()

if __name__ == "__main__" :
    test()
