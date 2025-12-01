import sys, os

import numpy as np
import matplotlib.pyplot as plt
import nptdms

# fonction «scrolling window» retourne un generator de taille N M fois

# fichier de lecture de .tdms


# liste à faire

# lire un fichier TDMS
# documentation iter

# deux manières d'ouvrire
# 1 : tout ouvrir, lourds sur la mémoire -> read
# 2 : open le fichier et juste allez chercher les channels nécessaire un à la fois -> open

channel_dictionnary = {
    'B4' : 'ai0-B4', 
    'B5' : 'ai1-B5',
    'B3' : 'ai2 -B3',
    'B6' : 'ai3 -B6',
    'B7' : 'ai4 -B7',
    'B2' : 'ai5 -B2',
    'B8' : 'ai6 -B8',
    'B1' : 'ai7 -B1',
    'A4' : 'ai8 -A4',
    'A5' : 'ai9-A5',
    'A3' : 'ai10-A3',
    'A6' : 'ai11 -A6',
    'A7' : 'ai12 -A7',
    'A2' : 'ai13 -A2',
    'A8' : 'ai14 -A8',
    'A1' : 'ai15 -A1'}

def init(fichierTDMS, selectGroupe) :
    """
    - retourne tout les channels du groupe choisis du fichier donné

    arguments :
        - string fichierTDMS : fichier.tdms
        - int selectGroupe : indice du groupe, 0 à N

    return :
        - objet tdms contenant tout les channels du groupe
    """

    fiTDMS = nptdms.TdmsFile.open(fichierTDMS)
    groupe = fiTDMS.groups()[selectGroupe] # prendre pour acquis que yaura toujours un seul groupe
    channels = groupe.channels()
    return channels

def tdmsRead(fichierTDMS, select, selectGroupe, chunks) : # toujours 16?
    """
    - retourne tout les channels voulu d'un groupe

    arguments :
        - string fichierTDMS : fichier.tdms
        - auto select : slice, list ou int des indices des channels voulues
        - int selectGroupe : groupe voulue
        - bool chunks : channels en format data_chunks ou non
            - si True, les channels sont des generator à iter N fois pour avoir tout le data

    return :
        - iter dataIter : iterator contenant les channels voulue 
        - int nbrData : nombre de channels dans dataIter
    """

    channels = init(fichierTDMS, selectGroupe)
    
    if isinstance(select, slice) :
        dataIter = channels[select]

    elif isinstance(select, int) :
        dataIter = [channels[select]]

    elif isinstance(select, list) :
        dataIter = [channels[i] for i in select]

    # si possible, pouvoir index avec les noms?

    if  chunks :
        dataIter = [ i.data_chunks() for i in dataIter ]

    return iter(dataIter), len(dataIter)


def extract(fichierTDMS, sr, select=slice(0,16), offset=0, selectGroupe = 0, chunks = False) :
    """
    - wrapper de tdmsRead, retourne le data de channels voulus selon plusieurs options 

    arguments :
        - string fichierTDMS : fichier.tdms
        - int sr : sample rate
        - auto select = slice(0,16) : slice, list ou int des indices des channels voulus
        - int offset = 0, offset de combien sur la première mesure
        - int selectGroupe = 0 : indice du groupe voulus
        - bool chunks = False : retourn les channels en generator data_chunks

    return :
        - iter dataZip : iterator contenant les generator de data des channels et du temps associé
        - int nbrData : nbr de channels
    """

    def dataWrap(dataIter, sr, offset) :    # data._length
        """
        fonction interne, ne pas utiliser
        """

        for data in dataIter :

            nbrPoint = data._length
            fin = nbrPoint/sr
            yield data[offset:], np.linspace(offset/sr, fin, nbrPoint-offset)

    def dataWrapChunks(dataIterChunks, sr, offset) :  # nécessaire????  Oui? len vs ._length ?
        """
        fonction interne, ne pas utiliser
        """
        # si offset > len(data), crash, c'est un feature
        for j, data in enumerate(dataIterChunks) :
            nbrPoint = len(data)
            timeStamp = nbrPoint/sr

            offset = offset if j == 0 else 0
            time = np.linspace( j*timeStamp, (j+1)*timeStamp, nbrPoint )

            yield data[offset:], time[offset:]

    channelIter, nbrChannel = tdmsRead(fichierTDMS, select, selectGroupe, chunks)

    if not chunks :
        #dataZip = dataWrap(channelIter, sr)
        dataZip = iter([ dataWrap([i], sr, offset) for i in channelIter ])
    else :
        dataZip = iter([dataWrapChunks(i, sr, offset) for i in channelIter])

    return dataZip, nbrChannel


def test() :

    #fichier = "../alcheData/2025-07-29_JN_Gquadexp4/50mM_KCl/50mMKCl_SR1e4_3h_Vds_1.tdms"
    fichier = "../alcheData/2025-06-06_JN_Gquadexp1/DNA_10mKCl/DNA_KCl_10mM_SR1e4_500kOmhsVg0_1_long.tdms"
    #fichier = "../alcheData/2025-06-06_JN_Gquadexp1/DNA_10mKCl/DNA_KCl_10mM_SR1e4_500kOmhsVg0_1_court.tdms"

    sr = int(1e4)

    #print(help(tdmsRead))

    dataIter, nbrData = extract(fichier, sr, [0,3], offset=0, chunks=False)   # si chunk de plusieurs, vraiment bon?
    dataIter2, nbrData2 = extract(fichier, sr, [0,3], offset=100, chunks=False)   # si chunk de plusieurs, vraiment bon?
    # regarder si possible de diviser en plusieurs generator

    #for data, time in dataZip :
    #    print(data)
    #    print(time)

    #for zip1, zip2 in zip(dataIter, dataIter2) :
    #    data1, _ = next(zip1)
    #    data2, _ = next(zip2)
    #    for i, _ in zip1 :
    #        data1 = np.concatenate((data1, i))

    #    print(np.all(data1 == data2))   # True!!!
    dataZip1 = list(dataIter)[1]
    dataZip2 = list(dataIter2)[1]

    for data, time in dataZip1 :
        plt.plot(time, data)

        a = data

    for data, time in dataZip2 :
        plt.plot(time, data)

        b = data


    plt.show()


    print(np.all( a[100:] == b ))


    

    """
    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    for j, dataZip in enumerate(dataIter) :
        #data individuel
        for data, time in dataZip :
            ax.plot(time, data, color=["tab:orange", "tab:blue"][j]) #, label=f"{j}")

    for j, dataZip in enumerate(dataIter2) :
        #data individuel
        for data, time in dataZip :
            ax2.plot(time, data, color=["tab:orange", "tab:blue"][j]) #, label=f"{j}")


    #fig.legend()
    plt.show()
    """

if __name__ == "__main__" :
    test()
