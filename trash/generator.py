import numpy as np
import matplotlib.pyplot as plt

def genEXP(x) :
    yield np.exp(x)

def test() :
    x = np.linspace(0, 1, 10000).tolist()
    #y = np.exp(x)
    y = genEXP(x)

    print(next(y))

    #fig = plt.figure()
    #ax = fig.add_subplot()

    #ax.plot(x, y)

    #plt.show()

if __name__ == "__main__" :
    test()
