import numpy as np
import matplotlib.pyplot as plt
import nptdms

def test() :
    fichier = "../../alcheData/2025-06-06_JN_Gquadexp1/DNA_10mKCl/DNA_KCl_10mM_SR1e4_500kOmhsVg0_1_long.tdms"

    f1 = nptdms.TdmsFile.open(fichier)
    g = f1.groups()[0]
    c = g.channels()[3]

    md = nptdms.TdmsFile.read_metadata(fichier)
    print(md.groups()[0].channels()[0])

    da = md.groups()[0].channels() 

    for name, value in da[1].properties.items() :
        print(f"{name} {value}")

    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    ax.plot(c[:])
    ax.plot(c.read_data(scaled=False))
    ax2.plot( c[:] - c.read_data(scaled=False) )
    #for j, i in enumerate(c.data_chunks()) :
    #    data = i[:]
    #    ax.plot(np.arange(j*data.size, (j+1)*data.size), data)


    plt.show()


if __name__ == "__main__" :
    test()
