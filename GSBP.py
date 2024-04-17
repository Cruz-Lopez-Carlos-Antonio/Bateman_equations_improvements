#General Solution based on Partitions
#Developed by:
#Carlos-Antonio Cruz-Lopez, cacl.nucl@gmail.com
#Gilberto Espinosa-Paredes, gepe@xanum.uam.mx
#Juan Luis Francois, juan.luis.francois@gmail.com


#Libraries that are required
import itertools
import numpy as np
import math
from decimal import*

#Precision
getcontext().prec=64


# Code that solves the Diophantine Equations.
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

# Function Chi        
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


#General Solution Based on Partitions
#Inputs:
#X0 = Initial condition
#DC = Decay Constants
#t = time
def GPS(X0,DC,t):
    #Quotient X0/lambda_n
    term1 = Decimal(X0)/DC[-1]

    #Transforms the list of Decay constants in a set
    Red_S = set(DC)
    Au1 = list(Red_S)
    
    
    #Compute the factors mu_k, that represents the number
    #of times that the isotope k is repeated and store it
    # in the list Mu
    Mu = [ ]
    for u in Au1:
        Mu.append(DC.count(u)-1)

    #Product of the Step 1
    term2 = 1
    for z in range(len(Au1)):
        term2 =term2*(Au1[z]**(Mu[z]+1))

    #Main Sum
    s1 = 0
    for i in range(len(Au1)):
        #Compute the exponential
        p1 = Decimal.exp(-Decimal(Au1[i])*Decimal(t))

        #Shifted method described in Section 4.5.1
        Au2 = Au1.copy()
        Au2.remove(Au1[i])
        Mu_m = Mu.copy()
        Mu_m.remove(Mu[i])

        #Product
        p2 = Decimal(1)
        for j in range(len(Au2)):
            p2 =p2*(Decimal(1)/((Decimal(Au2[j])-Decimal(Au1[i])))**Decimal((Mu_m[j]+1)))       
        s2 = 0

        #Second sum where the Chi function is valuated 
        for l in range(Mu[i]+1):
            L = [ ]
            z1 = Decimal((t**l))/Decimal(math.factorial(l))

            #Calling the function of Diophantine Equations
            partitions_restricted(Mu[i]-l,len(Au2),L)

            #Calling the Chi function
            z2 = chi(i,Mu[i]-l,Mu,Au1,L)
            s2 = Decimal(s2)+z1*z2

        s1 = Decimal(s1)+Decimal(p1)*p2*s2

    #Final product
    solution = Decimal(term1)*Decimal(term2)*s1
    return solution


######################################################################################
######################################################################################
######################################################################################
#################################### INPUT ###########################################
#Vector Half_life


Half_lifes = [2,2,3,3,3,4]

#Generation of the Decay constant vector
#Do not modify this lines
D = [ ]
for z in Half_lifes:
    D.append(Decimal(math.log(2))/Decimal(z))

#Time vector
Time_vector = [0,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01]

#Initial condition
IC = 6.023E23

#Exection of the code
for k in Time_vector:
    solucion = GPS(IC,D,k)
    print("Time",k ,"Concentration",float(solucion))
