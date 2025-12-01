import numpy as np
import matplotlib.pyplot as plt
import nptdms

class tdms :
    def __init__(self, fichierNom) :
        self.fichier = nptdms.TdmsFile.open(fichierNom)
        self.groupes = self.fichier.groups()
        self.channels = self.fichier["AI channels"].channels() # jusqu'à preuve du contraire

    def showAles(self) :

        # faire les beau plots

        for i in range( len(self.channels) ) :

            pass

def main() :
    fichier = "../dataJND/DNA_KCl_10mM_SR1e4_500kOmhsVg0_1_court.tdms"

    tocha = tdms(fichier)

if __name__ == "__main__" :
    main()
