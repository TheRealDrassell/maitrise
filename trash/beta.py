import numpy as np
import matplotlib.pyplot as plt
import nptdms

def data(channel) :
    return channel[:]

def extract(fichier, i) :
    fichierTDMS = nptdms.TdmsFile.read(fichier)
    groupe = fichierTDMS.groups()[0]
    channels = groupe.channels()

    return data(channels[i])

def main() :
    #fichier = "../dataJND/DNA_KCl_10mM_SR1e4_500kOmhsVg0_1_court.tdms"
    fichier = "../alcheData/2025-06-06_JN_Gquadexp1/DNA_10mKCl/DNA_KCl_10mM_SR1e4_500kOmhsVg0_1_court.tdms"

    figure = plt.figure(1)
    ax = figure.add_subplot()

    ales = 0 

    for i in range(16) :
        data = extract(fichier, i)
        ax.plot(data)

    plt.show()

if __name__ == "__main__" :
    main()
