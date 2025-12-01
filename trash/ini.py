import numpy as np
import matplotlib.pyplot as plt
import nptdms

def main() :
    """
    fichier = nptdms.TdmsFile.read("../dataJND/DNA_KCl_10mM_SR1e4_500kOmhsVg0_1_court.tdms")

    group = fichier["AI channels"]
    channels = group.channels()
    channel = channels[0]

    #x = channel.time_track()
    y = channel[:]

    print(channel.properties.items())
    print(channel.path)


    #print(dir(channel))

    #plt.plot(y)
    #plt.show()
    """

    with nptdm.TdmsFile.open("../dataJND/DNA_KCl_10mM_SR1e4_500kOmhsVg0_1_court.tdms") as fichier :
        groupes = fichier.groups()
        channels = groupes.channels()


if __name__ == "__main__" :
    main()
