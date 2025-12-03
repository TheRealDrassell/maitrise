import numpy as np
from timeit import default_timer as time
N = int(1e6)
M = int(1e2)

def main() :
    a = np.random.random(N)

    windowS = 5
    nbrWind = a.size - windowS + 1

    start = time()
    for i in range(M) :

        b = np.lib.stride_tricks.sliding_window_view(a, windowS)
    end = time()
    print((end - start)/M)

    start = time()
    for i in range(M) :

        indx = np.arange(windowS) + np.arange(nbrWind).reshape(1, nbrWind, -1)
        c = a[indx]
    end = time()
    print((end - start)/M)


 ( np.abs(sp) >= 500 )   print(np.all( b == c ))

if __name__ == "__main__" :
    main()
