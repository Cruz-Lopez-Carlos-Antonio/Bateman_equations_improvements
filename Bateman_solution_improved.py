import itertools
import numpy as np
import math
from decimal import*
getcontext().prec=64



def partitions_restricted(m,n,L):
    List = [ ]                                               
    for k in range(m+1):                                                          #Creation of the set  
        List.append(k)
    List2 = [ ]
    
    for k in itertools.product(List,repeat=n):                            #Creation of the set P_n
        a = np.array(k)
        l = np.sum(a)
        if l==m:
            L.append(list(a))
        
def chi(i,j,Mu,Lambd,L):
    c = 0
    for u in L:
        Aux_L = Mu.copy()
        Aux_L.remove(Mu[i])
        Aux2 = Lambd.copy()
        Aux2.remove(Lambd[i])
        a =1
        for k in range(len(Aux_L)):
            b_f = Decimal(math.comb(Aux_L[k]+u[k],Aux_L[k]))
            
            dif = (Decimal(1)/(Decimal(Lambd[i])-Decimal(Aux2[k]))**Decimal(int(u[k])))
            a = a*b_f*dif
        c = c+a
        #print(c)
    return c


def GPS(X0,DC,t):
    term1 = Decimal(X0)/DC[-1]
    Red_S = set(DC)
    Au1 = list(Red_S)
    Mu = [ ]
    for u in Au1:
        Mu.append(DC.count(u)-1)
    term2 = 1
    for z in range(len(Au1)):
        term2 =term2*(Au1[z]**(Mu[z]+1))
    s1 = 0
    for i in range(len(Au1)):
        p1 = Decimal.exp(-Decimal(Au1[i])*Decimal(t))
        Au2 = Au1.copy()
        Au2.remove(Au1[i])
        Mu_m = Mu.copy()
        Mu_m.remove(Mu[i])
        p2 = Decimal(1)
        for j in range(len(Au2)):
            p2 =p2*(Decimal(1)/((Decimal(Au2[j])-Decimal(Au1[i])))**Decimal((Mu_m[j]+1)))       
        s2 = 0
         
        for l in range(Mu[i]+1):
            L = [ ]
            z1 = Decimal((t**l))/Decimal(math.factorial(l))
            partitions_restricted(Mu[i]-l,len(Au2),L)
            z2 = chi(i,Mu[i]-l,Mu,Au1,L)
            s2 = Decimal(s2)+z1*z2

        s1 = Decimal(s1)+Decimal(p1)*p2*s2

    solution = Decimal(term1)*Decimal(term2)*s1
    return solution

Half_lifes = [2,2,3,3,3,4]
D = [ ]
for z in Half_lifes:
    D.append(Decimal(math.log(2))/Decimal(z))
Time_vector = [0,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01]


Vector_tiempo = [0]
#Time_vector = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100]

for k in Time_vector:
    solucion = GPS(6.023E23,D,k)
    print(float(solucion))

