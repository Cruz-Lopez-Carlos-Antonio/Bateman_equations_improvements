#Restricted partitios.

import itertools
import numpy as np

def partitions_restricted(m,n,L):
    List = [ ]                                               
    for k in range(m+1):                                                         #Limitation introduced Eq.(42) of the repository  
        List.append(k)                                                           

    for k in itertools.product(List,repeat=n):                                  #Cartesian product
        a = np.array(k)                                                         #Transform the result of the cartesian product to an array (vector)
        l = np.sum(a)                                                           # Performs the sum of each vector
        if l==m:                                                                #Kronecker's Delta
            L.append(k)
