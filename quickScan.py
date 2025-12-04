import numpy as np
import matplotlib.pyplot as plt
import os, sys

from tdmsRead import extract

def scan() :

    alesDir = []

    base = "../alcheData/2025-06-06_JN_Gquadexp1/"
    for root, dirs, files in os.walk(base) :
        if dirs :
            for d in dirs :
                alesDir.append(f"{base}{d}")

    curr = 3
    alesFile = os.listdir(alesDir[curr])
    tdmsFiles = [i for i in alesFile if i[::-1][:5] == "smdt."]

    for head in tdmsFiles :
        fichier = f"{alesDir[curr]}/{head}"
        sr = 1e4  # pas toujours vrai

        dataIter, nbrData = extract(fichier, sr, offset=100, chunks=False)

        fig, axs = plt.subplots(nbrData, 1, figsize=(6,10))
        #fig = plt.figure()
        #ax = fig.add_subplot()
        for j, dataZip in enumerate(dataIter) :
            for data, time in dataZip :

                moy = np.mean(data)
                std = np.std(data)
                axs[j].hist(data, bins=100, range= [moy - std*3,moy+std*3])
                #ax.plot(time, data)

        fig.suptitle(head)
        fig.tight_layout()
        plt.show()

if __name__ == "__main__" :
    scan()
