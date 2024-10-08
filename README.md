# General Solution of the Bateman Equations using Cuachy products and Divided Differences. 
The present repository contains a set of Python algorithms that solve and improve the General Solution of the Bateman Equations, which were developed in the following article **"General Solution of Bateman Equations using Cauchy Products and the Theory of Divided Differences"**, which was submitted to the **"Annals of Nuclear Energy"** journal.

The programs are licensed under a Creative Commons Attribution 4.0 International License: http://creativecommons.org/licenses/by/4.0/

**Authors**: Carlos-Antonio Cruz-López (cacl.nucl@gmail.com), Gilberto Espinosa-Paredes (gepe@xanum.uam.mx), Juan-Luis François-Lacouture (juan.luis.francois@gmail.com)

Mathematical and algorithmical generalities of the codes are described in the following lines with the purpose to provide some insight of the developed work. Nevertheless, a more detailed and precise discussion is provided in the submitted article.

## Financial Support.
The authors appreciate the financial support received from the Consejo Nacional de Ciencia y Tecnología, CONACYT, under the program “Estancias Posdoctorales por México, 2022”, with the project entitled: “Desarrollo de modelos fenomenológicos energéticos de orden fraccional, para la optimización y simulación en reactores nucleares de potencia”, by which the present development was possible.

## Software specifications.
The AnalyticNPKE codes were written in the Python programming language in its version 3. The codes require the following libraries:

- [x] itertools
- [x] numpy
- [x] decimal
- [x] math

## Index of the Repository
1. [Mathematical description of the problem](https://github.com/Cruz-Lopez-Carlos-Antonio/Bateman_equations_improvements/blob/main/README.md#1-mathematical-description-of-the-problem).
   - [1.1 Differential mass-balance equations](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticNPKE?tab=readme-ov-file#51-sums)
   - [1.2 Laplace transform of the system.](https://github.com/Cruz-Lopez-Carlos-Antonio/Bateman_equations_improvements#12-laplace-transform-of-the-system)
   - [1.3 Partial fraction decomposition.](https://github.com/Cruz-Lopez-Carlos-Antonio/Bateman_equations_improvements/blob/main/README.md#13-partial-fraction-decomposition)
   - [1.4 Solution of Bateman Equations.](https://github.com/Cruz-Lopez-Carlos-Antonio/Bateman_equations_improvements/blob/main/README.md#14-solution-of-the-bateman-equations)
2. [Cetnar's solution.](https://github.com/Cruz-Lopez-Carlos-Antonio/Bateman_equations_improvements/blob/main/README.md#1-mathematical-description-of-the-problem)
   - [2.1 Solution based on increments.](https://github.com/Cruz-Lopez-Carlos-Antonio/Bateman_equations_improvements/blob/main/README.md#21-solution-based-on-increments)
   - [2.2 Issues related to nested sums.](https://github.com/Cruz-Lopez-Carlos-Antonio/Bateman_equations_improvements/blob/main/README.md#22-the-issues-of-the-nested-sums)
   - [2.3 Diophantine Equations](https://github.com/Cruz-Lopez-Carlos-Antonio/Bateman_equations_improvements/blob/main/README.md#23-diophantie-equations)
3. [Improvements on the Cetnar's solution using Cauchy products](https://github.com/Cruz-Lopez-Carlos-Antonio/Bateman_equations_improvements/blob/main/README.md#3-improvements-on-the-cetnars-solution-using-cauchy-products)
   - [3.1 Reduction of $\psi_{i,j}$](https://github.com/Cruz-Lopez-Carlos-Antonio/Bateman_equations_improvements#31-reduction-of-the-product-involved-in-psi_ij)
   - [3.2 Including Diophantine Equations.](https://github.com/Cruz-Lopez-Carlos-Antonio/Bateman_equations_improvements?tab=readme-ov-file#32-including-diophantine-equations)
   - [3.3 Further Simplification.](https://github.com/Cruz-Lopez-Carlos-Antonio/Bateman_equations_improvements/tree/main#33-further-simplification)
4. [Algorithmic Implementation](https://github.com/Cruz-Lopez-Carlos-Antonio/Bateman_equations_improvements/tree/main#4-algorithmic-implementation))
   - [4.1 Diophantine Equations](https://github.com/Cruz-Lopez-Carlos-Antonio/Bateman_equations_improvements#41-diophantine-equations)
   - [4.2 Cartesian product](https://github.com/Cruz-Lopez-Carlos-Antonio/Bateman_equations_improvements/blob/main/README.md#42-cartesian-product) 

## 1. Mathematical description of the problem.
### 1.1. Differential mass-balance equations. 
Bateman equations describe the time evolution of a set of nuclides in succesive transformations due to decay and transmutation processes. These transformations can be represented by the following elementary structure:

$$ X_1\buildrel\lambda_1\over\rightarrow X_2\buildrel\lambda_2\over\rightarrow\ldots\buildrel\lambda_{n-1}\over\rightarrow X_n\buildrel\lambda_n\over\rightarrow \tag{1} $$

where $X_1, X_2,...,X_n$ are the concentrations of the nuclides and $\lambda_1, \lambda_2, ..., \lambda_n$ are their decay constants. Structure given in Eq. (1) is known as a linear chain, and its time evolution can be modeled using the following set of coupled differential equations:

$$\frac{dX_i(t)}{dt}=\lambda_{i-1}X_{i-1}\left(t\right)-\lambda_{i}X_{i}(t), \tag{2}$$

with $\lambda_{i-1} X_{i-1}=0$ and $1\leq i\leq n$. 
### 1.2 Laplace transform of the system
One of the most elementary ways to solve the system given in Eq.(2) consists of using the Laplace transform. For such task it is necessary to consider the following relationship:
$$\mathcal{L}\\{f^{(n)}(t),t\\}=s^n\mathcal{L} \\{f(t),t\\}-s^{n-1}f(0)-s^{n-2}f^\prime (0)-\ldots-f^{(n-1)}(0), \tag{3}$$
where it is assummed that the Laplace transform exists, and where $f^{(n-1)}$ denotes the $n-1$ derivative of the function $f(t)$. 
Applying the last relationship to both sides of Eq. (2), the following system is obtained:

$$s{\widetilde{x}}_ i-X_i(0)=\lambda_{i-1}{\widetilde{x}}_{i-1}-\lambda_i{\widetilde{x}}_i \tag{4}$$

where ${\widetilde{x}}_ i$ and ${\widetilde{x}}_ {i-1}$ are the Laplace transforms of $X_i(t)$ and $X _ {i-1}$, respectively, and $X _ i(0)$ are the initial conditions. The system given in (4) can be rewritten as:
$${\widetilde{x}}_ i=\frac{\lambda_{i-1}{\widetilde{x}}_ {i-1}}{s+\lambda_i}+\frac{X_i\left(0\right)}{s+\lambda_i}. \tag{5}$$
Considering the case where:
$$X_i (0)=0, \ i≥2 \tag{6},$$
the system given in Eq. (5) is reduced to: 
$${\widetilde{x}}_ 1=\frac{X_1\left(0\right)}{s+\lambda_1},\ \ {\widetilde{x}}_ i=\frac{\lambda_{i-1}{\widetilde{x}}_{i-1}}{s+\lambda_i}. \tag{7}$$
After multiple replacements, starting with $\widetilde{x}_1$ and continuing forward, it follows that: 

$$\begin{matrix}
{\widetilde{x}}_ 2=&\frac{\lambda_1X_1(0)}{(s+\lambda_2)(s+\lambda_1)},\\
{\widetilde{x}}_ 3=&\frac{\lambda_2\lambda_1X_1(0)}{(s+\lambda_3)(s+\lambda_2)(s+\lambda_1)},\\
\vdots&\vdots\\
{\widetilde{x}}_ n=&\frac{\lambda_{n-1}\lambda_{n-2}\cdots\lambda_1X_1(0)}{\left(s+\lambda_n\right)\left(s+\lambda_{n-1}\right)\cdots(s+\lambda_1)} .\\
\end{matrix} 
\tag{8}$$

The last expression, ${\widetilde{x}}_n$, is the most general relationship because it provides a general formula for each nuclide. Using the product notation:

$$\lambda_{n-1}\lambda_{n-2}\cdots\lambda_1=\prod_{k=1}^{n-1}\lambda_k ,\tag{9}$$

$$\frac{1}{\left(s+\lambda_n\right)\left(s+\lambda_{n-1}\right)\cdots(s+\lambda_1)}=\prod_{i=1}^{n}\frac{1}{(s+\lambda_i),} \tag{10}$$

the last expression in Eq. (8) can be written as:
$${\widetilde{x}}_ n=\prod_{k=1}^{n-1}\lambda_k\prod_{i=1}^{n}\frac{1}{(s+\lambda_i)}. \tag{11}$$

### 1.3 Partial fraction decomposition.
The Inverse Laplace transform of Eq. (11) can be found in a straighforward way if the Eq. (10) is expressed in terms of partial fractions:
$$\prod_{i=1}^{n}\frac{1}{(s+\lambda_i)}=\sum_{i=1}^{n}\frac{c_i}{s+\lambda_i} \tag{12}$$
where the coefficients, $c_i$ must to be computed. This last can be done **if we assume that all the decay constants, $\lambda_i$ are different**. In such case, it follows that:
$$c_i=\lim_{s\rightarrow-\lambda_i}{\prod_{j=1}^{n}\frac{(s+\lambda_i)}{(s+\lambda_j)}}=\prod_{j=1,j\neq i}^{n}\frac{1}{(\lambda_j-\lambda_i)}. \tag{13}$$
Combining the last relationships and replace them in Eq. (11), it follows that:
$${\widetilde{x}}_ n=\prod_{k=1}^{n-1}\lambda_k\sum_{i=1}^{n}\frac{c_i}{s+\lambda_i}=\prod_{k=1}^{n-1}\lambda_k\sum_{i=1}^{n}{\prod_{j=1,j\neq i}^{n}\frac{1}{(\lambda_j-\lambda_i)}\frac{1}{s+\lambda_i}} \tag{14}$$
### 1.4 Solution of the Bateman Equations. 
Using Eq. (14), the inverse Laplace transform of Eq. (11) is reduced to
$$x_n(t)=\mathcal{L}^{-1}\\{{\widetilde{x}}_ n\\}=\prod_{k=1}^{n-1}\lambda_k\sum_{i=1}^{n}{\prod_{j=1,j\neq i}^{n}\frac{1}{(\lambda_j-\lambda_i)}\mathcal{L}^{-1}\\{\frac{1}{s+\lambda_i}}\\}. \tag{15}$$
Using that $\mathcal{L}^{-1}\\{\frac{1}{s+\lambda_i}\\}=\exp(-\lambda_i t)$, it follows that:
$$x_n\left(t\right)=\prod_{k=1}^{n-1}\lambda_k\sum_{i=1}^{n}\prod_{j=1,j\neq i}^{n}\frac{\exp(-\lambda_it)}{(\lambda_j-\lambda_i)}, \tag{16}$$
which is known as the Bateman solution.
## 2. Cetnar's solution.
### 2.1 Solution based on increments. 
The Eq. (16) was developed assuming that the lambda constants, $\lambda_i, 1\leq i \leq n$, are different. Nevertheless, there cases where such condition is not fulfilled. A common practice to address this issue consists of introducing small increments $\Delta_i$ in the repeated decay constants. For example, if we considered that the decay $\lambda_i$ appears $m$ times, then such decay constants are modified as:
$$\lambda_i,\ \lambda_i+\Delta_i,\lambda_i+2\Delta_i,…,\lambda_i+(m-1)\Delta_i. \tag{17}$$
Jerzy Centar (2006) extended this idea and proposed the following approximation of the Bateman solution including the increments:
$$X_n(t)\approx \chi_n(t,\Delta_i,\Delta_j)=\frac{X_1\left(0\right)}{\lambda_n}\sum_{i=1}^{n}\sum_{m=0}^{\mu_i}\exp{\left(-\left(\lambda_i+m\Delta_i\right)t\right)}$$
$$\times\left(\prod_{l=0,l\neq m}^{\mu_i}\frac{\lambda_i+l\Delta_i}{\left(l-m\right)\Delta_i}\right)\prod_{j=1,j\neq i}^{n}\prod_{k=0}^{\mu_i}\frac{\lambda_j+k\Delta_j}{\lambda_j+k\Delta_j-\lambda_i-m\Delta_i}, \tag{18}$$
where the function $\chi_n(t,\Delta_i,\Delta_j)$ represents an approximation whose accuracy depends on the way in which the increments $\Delta_i,\Delta_j$ are chosen, being $n$ the number of different nuclides (or different decay constants) and where $\mu_i$ is the number of times that the isotope of the type $i$ is repeated in the linear chain. Then, in a very ingenious procedure, Cetnar found the exact analytical solution considering the following limit:
$$\lim_{\Delta_i,\Delta_j\rightarrow0}{\chi_n(t,\Delta_i,\Delta_j})=X_n. \tag{19}$$
After a very difficult procedure, the following equation is obtained:
$$X_n\left(t\right)=\frac{X_1\left(0\right)}{\lambda_n}\sum_{i=1}^{n}{\lambda_i\alpha_i\exp(-\lambda_it)}\cdot\sum_{l=0}^{\mu_i}{\frac{\left(\lambda_it\right)^l}{l!}\cdot\psi_{i,\mu_i-l}}, \tag{20}$$
where:
$$\alpha_i=\prod_{j=1,j\neq i}^{n}\left(\frac{\lambda_j}{\lambda_j-\lambda_i}\right)^{\mu_j+1}, \tag{21}$$
and:

$$\psi_{i,j}=\sum_{h_1=0}^{j} \sum_{h_2=0}^{j}\cdots\sum_{h_{i-1}=0}^{j}\sum_{h_{i+1}=0}^{j} \cdots\sum_{h_n=0}^{j} \prod_{k=1,k\neq i}^{n} \binom{h_k+\mu_k}{\mu_k} \left(\frac{\lambda_i}{\lambda_i-\lambda_k}\right)^{h_k}\delta_{j,p} \tag{22} $$

where:

$$p=\sum_{l=1,l\neq i}^{n}h_l,  \ \ \ \ \delta_i,p=\begin{cases}
  1  &  \text{ if $j=p$ } \\
  0  &  \text{ if $j\neq p$}
\end{cases} \tag{23}$$
### 2.2 Issues related to nested sums.
Even when the solution (as well as the procedure that was followed to obtain it) is very ingenious, the way in which it is expressed is very disadvantageous. The main issue is related to the nested sum in the term $\psi_{i,j}$ in Eq. (22), which in turn involves the Kronecker's delta given in Eq. (23). 
This implies, in computational terms, that several terms of the total sum in $\psi_{i,j}$ will be multiplied by zero, and therefore they will not be taken into account. The origin of this nested sum can be explained in terms of the following relationship:
$$\sum_{h_1=0}^{j}\sum_{h_2=0}^{j}\cdots\sum_{h_n=0}^{j}{f(h_1,h_2,\ldots,h_n)}\delta_{h_1+h_2+\ldots+h_n,j}$$
$$=\sum_{h_1+h_2+\ldots+h_n=j}{f(h_1,h_2,\ldots,h_n)}, \tag{24}$$

which in turn arises due to the presence of **Cauchy products** as it can be verified in the submitted paper (Cruz-López and Espinosa-Paredes, 2024). 
One of the main contribution of our work consists of having identified the Eq. (24), which allows simplifying the nested multiple sum to a single one. Additionally, the Delta's kronecker is removed, which implies an important simplification. Nevertheless, it is necessary to analyze the set over which the sum of the right side of Eq. (24) is carried out, which is done in detail in the following section. 

### 2.3 Diophantie Equations
The sum given in Eq. (24) is carried out over the following set:
$$h_1+h_2+\ldots+h_n=j, \tag{25}$$
where the indexes $h_i$, $1\leq i\leq n$ are non-negative integers. Therefore the Eq. (25) represents all the different permutations of non-negative integers whose sum is equal to $n$. The relationship given in Eq. (25) is a Diophantine Equation, and its solution is a set of vectors over which the sum in Eq. (24) is valuated. For example, considering the case:
$$h_1+h_2+h_3=4 \tag{25},$$
which is the set of all the permutations of non-integer numbers $h_1,h_2$ and $h_3$, whose sum is equal to 4. This set is the same that:

$$\begin{Bmatrix}
(0, 0, 4),& (1, 0, 3),&(2,1,1) \\
(0, 1, 3),&(1, 1, 2),&(2, 2, 0) \\ 
(0, 2, 2),&(1, 2, 1),&(3, 0, 1) \\ 
(0, 3, 1),& (1, 3, 0),&(3, 1, 0)\\
(0, 4, 0),&(2, 0, 2),&(4, 0, 0)
\end{Bmatrix} \tag{26}$$

where each triad represents a solution of the Diophantine Equation given in Eq. (25). Then, instead of using nested sums and the Kronecker's Delta, it is possible to rewrite such term as:
$$\sum_{h_1=0}^{4}\sum_{h_2=0}^{4}{\sum_{h_3=0}^{4}{f(h_1,h_2,h_3)}\delta_{h_1+h_2+h_3,4}}$$
$$=\sum_{h_1+h_2+h_3}{f(h_1,h_2,h_3)}=\sum_{\left(0,0,4\right),\ \ \left(1,0,3\right),\ldots,(4,0,0)}{f(h_1,h_2,h_3)}$$
$$=f\left(0,0,4\right)+f\left(1,0,3\right)+f\left(0,2,2\right)+\ldots+f(4,\ 0,\ 0) \tag{27}.$$

As it can be observed in the last expression, the nested sums have been removed as well as the Kronecker's Delta. Even more, it is possible to observe that such sum only depends on the set of Diophantine Equations.

## 3. Improvements on the Cetnar's solution using Cauchy products. 
### 3.1 Reduction of the product involved in $\psi_{i,j}$
The main contribution of the submitted article consists of simplifying the Cetnar's solution, removing the nested sums as well as the Kronecker's delta. Using the findings that were described in the last section, which are based on Cauchy products, it is possible to reduce the following product that appears inside Eq. (27) as follows:
$$\prod_{k=1,k\neq i}^{n}\left(\frac{\lambda_i}{\lambda_i-\lambda_k}\right)^{h_k}=\prod_{k=1,k\neq i}^{n}\left(\frac{1}{\lambda_i-\lambda_k}\right)^{h_k}\prod_{l=1,l\neq i}^{n}\lambda_i^{h_l}$$ 
$$\prod_{k=1,k\neq i}^{n}\left(\frac{1}{\lambda_i-\lambda_k}\right)^{h_k}\lambda_i^{h_1+h_2+\ldots+h_{i-1}+h_{i+1}+\ldots+h_l} \tag{28}$$
$$=\lambda_i^j\prod_{k=1,k\neq i}^{n}\left(\frac{1}{\lambda_i-\lambda_k}\right)^{h_k}, \tag{29}$$
where it has been used that:
$$h_1+h_2+\ldots h_{i-1}+h_i+\ldots+h_n=j \tag{30}$$
### 3.2 Including Diophantine Equations.
Using the procedure that is detailed in the submitted paper, the Cetnar's solution can be rewritten as:
$$X_n\left(t\right)=\frac{X_1\left(0\right)}{\lambda_n}\sum_{i=1}^{n}\lambda_i\prod_{j=1,j\neq i}^{n}{\frac{\lambda_j^{\mu_j+1}}{\left(\lambda_j-\lambda_i\right)^{\mu_j+1}}\exp(-\lambda_it)}\sum_{l=0}^{\mu_i}{\frac{\left(\lambda_it\right)^l}{l!}\psi_{i,\mu_i-l}}, \tag{30}$$
where:
$$\psi_{i,j}=\lambda_i^j\sum_{h_1+h_2+\ldots+h_{i-1}+h_{i+1}+\ldots+h_n=j}{\ \prod_{k=1,k\neq i}^{n}\binom{h_k+\mu_k}{\mu_k} \frac{1}{\left(\lambda_i-\lambda_k\right)^{h_k}}}. \tag{31}$$
As it can be observed, Eq. (31) is identical to Eq. (20), where the only difference is the definition of $\psi_{i,j}$, which is given in terms of the Diophantine Equations. This equivalence is fundamental with the purpose to develop an adequate algorithmic implementation.. 
### 3.3 Further simplification. 
It is possible to remove the term $\lambda_i^j$ from $\psi_{i,j}$. Similarly, the sum whose index is $l$ can be arranged, factoring the term $\lambda_i^l$, and removing some restrictions (see the submited paper for more details). After this modifications it follows that:
$$X_n\left(t\right)=\frac{X_1\left(0\right)}{\lambda_n}\ \prod_{k=1}^{n}\lambda_k^{\mu_k+1}\sum_{i=1}^{n}\prod_{j=1,j\neq i}^{n}\frac{\exp(-\lambda_it)}{\left(\lambda_j-\lambda_i\right)^{\mu_j+1}}\cdot\sum_{l=0}^{\mu_i}{\frac{t^l}{l!}\chi_{i,\mu_i-l}}, \tag{32}$$
where:
$$\chi_{i,j}=\sum_{h_1+h_2+\ldots+h_{i-1}+h_{i+1}+\ldots+h_n=j}{\ \prod_{k=1,k\neq i}^{n} \binom{h_k+\mu_k}{\mu_k}\frac{1}{\left(\lambda_i-\lambda_j\right)^{h_k}}}. \tag{33}$$
This last equation represents an improvement because several redundant operations have been reduced and simplified.
## 4. Algorithmic implementation. 
### 4.1 Diophantine Equations.
The first step of the algorithmic implementation consists of computing the Diophantine Equations, i.e., building the set of solutions of the following equation:
$$h_1+h_2+\ldots h_{i-1}+h_i+\ldots+h_n=j \tag{34}.$$
As it will be explained later, this task can be precomputed and **it must carried out a single time**, therefore, it is possible to use a standard or direct algorithm even if this one is not the most efficient. An elementary way of solving Eq. (34) consists of setting the Diophantine equation in vectorial terms, defining the following function:
$$P\left(\vec{x}\right)=P\left(\left(x_1,x_2,\ldots,x_n\right)\right)=\sum_{i=1}^{n}x_i=m, \tag{35}$$
where $\vec{x}$ is a vector of dimension $n$, whose entries are $x_1,x_2,...,x_n$, which in turn are non-negative integers with:
$$0\leq x_i\leq n. \tag{36}$$
The core of the idea is building all the possible vectors who fulfill these three condition, i.e., building the following set:
$$V_{\Sigma=1}=\set{\ \vec{x}=(x_1,x_2,\ldots,x_n)|\ \vec{x}\in \mathbb{R}^n ,0\le x_i\le n,\ x_i\in \mathbb{Z}^+, P(\vec{x})=m}.\ \ \tag{37}$$
### 4.2 Cartesian product.
As it can be observed in Eq. (37), $V_{\Sigma=1}$ is a subset of $\mathbb{R}^n$, and it contains vectors whose entries are non-negative integers lower than $n$. Therefore, it is possible to build another set, $V$ as follows:
$$V=\underbrace{\set{0,\ 1,\ 2,\ ...,n}\times \set{0,\ 1,\ 2,\ ...,n}\times\cdots\times \set{0,\ 1,\ 2,\ ...,n}}_{n}, \tag{38}$$

where it follows that:
$$V_{\Sigma=1}\subset V. \tag{39}$$
In other words, it is possible to build $V_{\Sigma=1}$ from $V$, i.e.:
$$V_{\Sigma=1}=\set{\vec{x}=(x_1,x_2,\ldots,x_n)\in V|P(\vec{x})=m}.\tag{40} $$
Therefore, in a first stage it is necessary to build $V$, which can be done in a straighforward way using the itertools Python library. For example, if we want to build the following set:
$$V_2=\set{0, 1}\times \set{0,\ 1}=\set{(0, 0),\ (0,\ 1),\ (1,\ 0),\ (1,\ 1)}, \tag{41}$$
which contains four ordered pairs, it is only necessary to carry out the product of the set $\set{0, 1}$ with repetition two, as it is shown in the following Python 3 code:

**Code 1. Example of a cartesian product**
```Python
import itertools
# Cartesian product example

List=[0,1]
z=itertools.product(List,repeat=2)
for u in z:
    print(u)
```
which produces the following output:

**Output of Code 1**
```Python
(0, 0)
(0, 1)
(1, 0)
(1, 1)
```
The last procedure can be generalized as follows:

**Code 2. Generalized Cartesian product**
```Python
import itertools

def cartesian_product(n,L):
    List = [ ]                                               
    for k in range(n+1):                                                          
        List.append(k)    #Creation of the set {0, 1, ..., n}                                                           

    for k in itertools.product(List,repeat=n):        #Cartesian product              
        L.append(k)                                   
```
where the function cartesian product admits an input list L, where the cartesian product $\set{0, 1, 2,...,n}\times \cdots \times \set{0, 1, 2,...,n}$ is storaged. 
### 4.3 Discrimination ordered pairs.
It is necessary to identify two key aspects related to the cartesian product, and with the Diophantine Equations described in Section 4.1. First, the number of products is related to the extension of the vector $\vec{x}$, which translated to our problem, is equal to the $n$ number of sum given in Eq. (34), which in turn is related to the number of different isotopes in a linear chain. Now we are interested in the sum of each of the vectors of the cartesian product, and particularly that such sum be equal to a value $m$. 
Therefore, it is possible to add a restriction to the entries of the vector as:
$$0\le x_i\le m,\ \ 1\le i\le n. \tag{42}$$
This lead to the definition of a new set that representes a refinement of the set $V$, which is given by:
$$V_{\leq m}=\set{\vec{x}=(x_1,x_2,\ldots,x_n) |0\le x_i\le m,\ P(\vec{x})=m}. \tag{43}$$
For decay and transmutation networks is common to consider that $m\leq n$, and in such case the following relationship is true:
$$V_{\Sigma=1}\subset V_{\le m}\subset V. \tag{44}$$
This represents a first discrimination of the original ordered pairs. Finally, once the set of vectors $\vec{x}$ with the last condition is built, it is necessary to add the final condition given by $P(x)=m$. These modifications can be included in the **Code 2** of the Cartesian product as follows:

**Code 3. Diophantine Equations/ Partitions restricted**
```Python
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
```
The first loop runs over the interval $(0, m)$, and in the second loop the condition of the Kronecker's Delta of Eq. (23) is introduced. 

### 4.4 Example of an application of Code 3.
Considering the following Diophantine Equation:
$$x_1+x_2+\ldots+x_7=4, \tag{45}$$
which is related to a linear chain of seven different isotopes, and where one of them appears, at least, four times. In order to compute the set of solutions of Eq. (45) the Code 3 is used as follows:
**Aplication of Code 3**
```Python
List_solutions = [ ]
partitions_restricted(4,7,List_solutions)
```
The output that is obtained consists of a list (List_solutions), where all the different solutions are storaged. Adding two lines related to the print of the content:

```Python
for z in range(len(List_solutions)):
    print(List_solutions[z])
print(len(List_solutions))
```
The corresponding output consists of a list with 210 entries, each of them with a tuple of seven numbers, whose sum is equal to 4. The reader can consult such output in the following hidden section:
**Output of the solution of Diophantine Equations**


<details><summary>CLICK HERE to expand the Output of the past example of Diophantine Equations </summary>
<p>


```Python
(0, 0, 0, 0, 0, 0, 4)
(0, 0, 0, 0, 0, 1, 3)
(0, 0, 0, 0, 0, 2, 2)
(0, 0, 0, 0, 0, 3, 1)
(0, 0, 0, 0, 0, 4, 0)
(0, 0, 0, 0, 1, 0, 3)
(0, 0, 0, 0, 1, 1, 2)
(0, 0, 0, 0, 1, 2, 1)
(0, 0, 0, 0, 1, 3, 0)
(0, 0, 0, 0, 2, 0, 2)
(0, 0, 0, 0, 2, 1, 1)
(0, 0, 0, 0, 2, 2, 0)
(0, 0, 0, 0, 3, 0, 1)
(0, 0, 0, 0, 3, 1, 0)
(0, 0, 0, 0, 4, 0, 0)
(0, 0, 0, 1, 0, 0, 3)
(0, 0, 0, 1, 0, 1, 2)
(0, 0, 0, 1, 0, 2, 1)
(0, 0, 0, 1, 0, 3, 0)
(0, 0, 0, 1, 1, 0, 2)
(0, 0, 0, 1, 1, 1, 1)
(0, 0, 0, 1, 1, 2, 0)
(0, 0, 0, 1, 2, 0, 1)
(0, 0, 0, 1, 2, 1, 0)
(0, 0, 0, 1, 3, 0, 0)
(0, 0, 0, 2, 0, 0, 2)
(0, 0, 0, 2, 0, 1, 1)
(0, 0, 0, 2, 0, 2, 0)
(0, 0, 0, 2, 1, 0, 1)
(0, 0, 0, 2, 1, 1, 0)
(0, 0, 0, 2, 2, 0, 0)
(0, 0, 0, 3, 0, 0, 1)
(0, 0, 0, 3, 0, 1, 0)
(0, 0, 0, 3, 1, 0, 0)
(0, 0, 0, 4, 0, 0, 0)
(0, 0, 1, 0, 0, 0, 3)
(0, 0, 1, 0, 0, 1, 2)
(0, 0, 1, 0, 0, 2, 1)
(0, 0, 1, 0, 0, 3, 0)
(0, 0, 1, 0, 1, 0, 2)
(0, 0, 1, 0, 1, 1, 1)
(0, 0, 1, 0, 1, 2, 0)
(0, 0, 1, 0, 2, 0, 1)
(0, 0, 1, 0, 2, 1, 0)
(0, 0, 1, 0, 3, 0, 0)
(0, 0, 1, 1, 0, 0, 2)
(0, 0, 1, 1, 0, 1, 1)
(0, 0, 1, 1, 0, 2, 0)
(0, 0, 1, 1, 1, 0, 1)
(0, 0, 1, 1, 1, 1, 0)
(0, 0, 1, 1, 2, 0, 0)
(0, 0, 1, 2, 0, 0, 1)
(0, 0, 1, 2, 0, 1, 0)
(0, 0, 1, 2, 1, 0, 0)
(0, 0, 1, 3, 0, 0, 0)
(0, 0, 2, 0, 0, 0, 2)
(0, 0, 2, 0, 0, 1, 1)
(0, 0, 2, 0, 0, 2, 0)
(0, 0, 2, 0, 1, 0, 1)
(0, 0, 2, 0, 1, 1, 0)
(0, 0, 2, 0, 2, 0, 0)
(0, 0, 2, 1, 0, 0, 1)
(0, 0, 2, 1, 0, 1, 0)
(0, 0, 2, 1, 1, 0, 0)
(0, 0, 2, 2, 0, 0, 0)
(0, 0, 3, 0, 0, 0, 1)
(0, 0, 3, 0, 0, 1, 0)
(0, 0, 3, 0, 1, 0, 0)
(0, 0, 3, 1, 0, 0, 0)
(0, 0, 4, 0, 0, 0, 0)
(0, 1, 0, 0, 0, 0, 3)
(0, 1, 0, 0, 0, 1, 2)
(0, 1, 0, 0, 0, 2, 1)
(0, 1, 0, 0, 0, 3, 0)
(0, 1, 0, 0, 1, 0, 2)
(0, 1, 0, 0, 1, 1, 1)
(0, 1, 0, 0, 1, 2, 0)
(0, 1, 0, 0, 2, 0, 1)
(0, 1, 0, 0, 2, 1, 0)
(0, 1, 0, 0, 3, 0, 0)
(0, 1, 0, 1, 0, 0, 2)
(0, 1, 0, 1, 0, 1, 1)
(0, 1, 0, 1, 0, 2, 0)
(0, 1, 0, 1, 1, 0, 1)
(0, 1, 0, 1, 1, 1, 0)
(0, 1, 0, 1, 2, 0, 0)
(0, 1, 0, 2, 0, 0, 1)
(0, 1, 0, 2, 0, 1, 0)
(0, 1, 0, 2, 1, 0, 0)
(0, 1, 0, 3, 0, 0, 0)
(0, 1, 1, 0, 0, 0, 2)
(0, 1, 1, 0, 0, 1, 1)
(0, 1, 1, 0, 0, 2, 0)
(0, 1, 1, 0, 1, 0, 1)
(0, 1, 1, 0, 1, 1, 0)
(0, 1, 1, 0, 2, 0, 0)
(0, 1, 1, 1, 0, 0, 1)
(0, 1, 1, 1, 0, 1, 0)
(0, 1, 1, 1, 1, 0, 0)
(0, 1, 1, 2, 0, 0, 0)
(0, 1, 2, 0, 0, 0, 1)
(0, 1, 2, 0, 0, 1, 0)
(0, 1, 2, 0, 1, 0, 0)
(0, 1, 2, 1, 0, 0, 0)
(0, 1, 3, 0, 0, 0, 0)
(0, 2, 0, 0, 0, 0, 2)
(0, 2, 0, 0, 0, 1, 1)
(0, 2, 0, 0, 0, 2, 0)
(0, 2, 0, 0, 1, 0, 1)
(0, 2, 0, 0, 1, 1, 0)
(0, 2, 0, 0, 2, 0, 0)
(0, 2, 0, 1, 0, 0, 1)
(0, 2, 0, 1, 0, 1, 0)
(0, 2, 0, 1, 1, 0, 0)
(0, 2, 0, 2, 0, 0, 0)
(0, 2, 1, 0, 0, 0, 1)
(0, 2, 1, 0, 0, 1, 0)
(0, 2, 1, 0, 1, 0, 0)
(0, 2, 1, 1, 0, 0, 0)
(0, 2, 2, 0, 0, 0, 0)
(0, 3, 0, 0, 0, 0, 1)
(0, 3, 0, 0, 0, 1, 0)
(0, 3, 0, 0, 1, 0, 0)
(0, 3, 0, 1, 0, 0, 0)
(0, 3, 1, 0, 0, 0, 0)
(0, 4, 0, 0, 0, 0, 0)
(1, 0, 0, 0, 0, 0, 3)
(1, 0, 0, 0, 0, 1, 2)
(1, 0, 0, 0, 0, 2, 1)
(1, 0, 0, 0, 0, 3, 0)
(1, 0, 0, 0, 1, 0, 2)
(1, 0, 0, 0, 1, 1, 1)
(1, 0, 0, 0, 1, 2, 0)
(1, 0, 0, 0, 2, 0, 1)
(1, 0, 0, 0, 2, 1, 0)
(1, 0, 0, 0, 3, 0, 0)
(1, 0, 0, 1, 0, 0, 2)
(1, 0, 0, 1, 0, 1, 1)
(1, 0, 0, 1, 0, 2, 0)
(1, 0, 0, 1, 1, 0, 1)
(1, 0, 0, 1, 1, 1, 0)
(1, 0, 0, 1, 2, 0, 0)
(1, 0, 0, 2, 0, 0, 1)
(1, 0, 0, 2, 0, 1, 0)
(1, 0, 0, 2, 1, 0, 0)
(1, 0, 0, 3, 0, 0, 0)
(1, 0, 1, 0, 0, 0, 2)
(1, 0, 1, 0, 0, 1, 1)
(1, 0, 1, 0, 0, 2, 0)
(1, 0, 1, 0, 1, 0, 1)
(1, 0, 1, 0, 1, 1, 0)
(1, 0, 1, 0, 2, 0, 0)
(1, 0, 1, 1, 0, 0, 1)
(1, 0, 1, 1, 0, 1, 0)
(1, 0, 1, 1, 1, 0, 0)
(1, 0, 1, 2, 0, 0, 0)
(1, 0, 2, 0, 0, 0, 1)
(1, 0, 2, 0, 0, 1, 0)
(1, 0, 2, 0, 1, 0, 0)
(1, 0, 2, 1, 0, 0, 0)
(1, 0, 3, 0, 0, 0, 0)
(1, 1, 0, 0, 0, 0, 2)
(1, 1, 0, 0, 0, 1, 1)
(1, 1, 0, 0, 0, 2, 0)
(1, 1, 0, 0, 1, 0, 1)
(1, 1, 0, 0, 1, 1, 0)
(1, 1, 0, 0, 2, 0, 0)
(1, 1, 0, 1, 0, 0, 1)
(1, 1, 0, 1, 0, 1, 0)
(1, 1, 0, 1, 1, 0, 0)
(1, 1, 0, 2, 0, 0, 0)
(1, 1, 1, 0, 0, 0, 1)
(1, 1, 1, 0, 0, 1, 0)
(1, 1, 1, 0, 1, 0, 0)
(1, 1, 1, 1, 0, 0, 0)
(1, 1, 2, 0, 0, 0, 0)
(1, 2, 0, 0, 0, 0, 1)
(1, 2, 0, 0, 0, 1, 0)
(1, 2, 0, 0, 1, 0, 0)
(1, 2, 0, 1, 0, 0, 0)
(1, 2, 1, 0, 0, 0, 0)
(1, 3, 0, 0, 0, 0, 0)
(2, 0, 0, 0, 0, 0, 2)
(2, 0, 0, 0, 0, 1, 1)
(2, 0, 0, 0, 0, 2, 0)
(2, 0, 0, 0, 1, 0, 1)
(2, 0, 0, 0, 1, 1, 0)
(2, 0, 0, 0, 2, 0, 0)
(2, 0, 0, 1, 0, 0, 1)
(2, 0, 0, 1, 0, 1, 0)
(2, 0, 0, 1, 1, 0, 0)
(2, 0, 0, 2, 0, 0, 0)
(2, 0, 1, 0, 0, 0, 1)
(2, 0, 1, 0, 0, 1, 0)
(2, 0, 1, 0, 1, 0, 0)
(2, 0, 1, 1, 0, 0, 0)
(2, 0, 2, 0, 0, 0, 0)
(2, 1, 0, 0, 0, 0, 1)
(2, 1, 0, 0, 0, 1, 0)
(2, 1, 0, 0, 1, 0, 0)
(2, 1, 0, 1, 0, 0, 0)
(2, 1, 1, 0, 0, 0, 0)
(2, 2, 0, 0, 0, 0, 0)
(3, 0, 0, 0, 0, 0, 1)
(3, 0, 0, 0, 0, 1, 0)
(3, 0, 0, 0, 1, 0, 0)
(3, 0, 0, 1, 0, 0, 0)
(3, 0, 1, 0, 0, 0, 0)
(3, 1, 0, 0, 0, 0, 0)
(4, 0, 0, 0, 0, 0, 0)

```
</p>
</details>

### 4.5 Implementation of Chi function, and shifted method ###
Once the solution of Diophantine Equations is computed with the **Code 3**, it is possible to develop and algorithm to implement the $\chi_{i,j}$ function given in Eq. (33).
For such task a shifted method will be used, where the limitations given for the index $i$, both in the product, $k=1, k\neq i$ and in the Diophantine Equations is removed. 
#### 4.5.1 Shifted method.
The shifted method consists of redefine an original set, ommiting the element in the position $i$. It can be expressed, in mathematical terms, defining the following set:
$$\Omega_i=\set{\lambda_j|1\le j\le n,\ j\neq i}, \tag{46}$$
after that, it is necessary to use a different index notation to enumerate the elements on $\Omega_j$, which can be understood considering the equivalence between the two sets:
$$\set{\lambda_1,\lambda_2,\ldots,\lambda_{i-1},\lambda_{i+1},\ldots,\lambda_n}\rightarrow \set{\lambda_1^\ast,\lambda_2^\ast,\ldots,\lambda_{i-1}^\ast,\lambda_i^\ast,\lambda_{i+1}^\ast,\ldots,\lambda_{n-1}^\ast}, \tag{47}$$
which can be summarized as follows:

$$\Omega_j=\begin{cases}
  \lambda_j  &  \text{ if $j \< i$ } \\
  \lambda_{j+1}  &  \text{ if $j\ge i$}
\end{cases}. \tag{48}$$

Using this redefinition, it follows that:

$$\sum_{h_1+h_2+\ldots+h_{i-1}+h_{i+1}+\ldots+h_n}{f(h_1,h_2,\ldots,h_{i-1},h_{i+1},\ldots,h_n)} $$

$$=\sum_{h_1^\ast+h_2^\ast+\ldots+h_{n-1}^\ast} f(h_1^\ast,h_2^\ast,\ldots,h_{n-1}^\ast) \tag{49}$$

and:
$$\prod_{k=1,\ k\neq i\ }^{n} \binom{h_k+\mu_k}{\mu_k}\frac{1}{\left(\lambda_i-\lambda_k\right)^{h_k}}= \prod_{k=1}^{n-1} \binom{h_k^\ast+\mu_k^\ast}{\mu_k^\ast} \frac{1}{\left(\lambda_i-\lambda_k^\ast\right)^{h_k^\ast}} \tag{50}$$


This shifted methodology can be implemented in a straighforward way in Python 3, using lists and the remove method. The following code contains the way in which this shifted method can be implemented:

**Code 4**
```Python
List_lambdas = [1,2,3,4,5,6,7,8] #Original list
Lambda_3=List_lambdas.remove(List_lambdas[3]) 
#List where the element in the position 3 (counting from zero, from the left to the right) is removed
```
#### 4.5.2 The Chi function. 
The following code contains the implementation of Eq. (33), using the shifted method that was discussed in the past section:
```Python
def chi(i,j,Mu,Lambd,L):
    c = 0
    for u in L:
        Aux_L = Mu.copy()
        Aux_L.remove(Mu[i])                 #using the shift method to remove mu_i
        Aux2 = Lambd.copy()
        Aux2.remove(Lambd[i])              #using the shift method to remove lambda_i
        a =1
        for k in range(len(Aux_L)):
            b_f = Decimal(math.comb(Aux_L[k]+u[k],Aux_L[k]))           
            dif = (Decimal(1)/(Decimal(Lambd[i])-Decimal(Aux2[k]))**Decimal(int(u[k])))
            a = a*b_f*dif
        c = c+a
    return c
```


## 5. Optimized Cetnar's solution.
Using the theory described in the last sections, it is possible develop a **first optimization** of the Cetnar's solution, which will be denoted as **Optimized Cetnar's solution (OCS)**. Essentially this first improvement consists of replacing the nested sums by the sum over the partitions over the solutions of the Diophantine Equations, i.e., using the following relationship:

$$\sum_{h_1=0}^{j}\sum_{h_2=0}^{j}\cdots\sum_{h_n=0}^{j}{f(h_1,h_2,\ldots,h_n)}\delta_{h_1+h_2+\ldots+h_n,j}$$

$$=\sum_{h_1+h_2+\ldots+h_n=j}{f(h_1,h_2,\ldots,h_n)}. \tag{51}$$

Nevertheless, this optimized version does not use the pre-calculations methodology that was described before, and therefore it is limited. On the other hand, the related code of this improved solution can be found in the file "OptimizedCetnarSolution.py" that is provided in this repository.
The methodology of this first optimization can be summarized as follows:

$$\underset{\downarrow}{\buildrel{\\mathrm{Step\ 1}}\over{\fbox{$\mathrm{Determine}\ \mu_k$}}}\ $$

$$\underset{\downarrow}{\buildrel{\\mathrm{Step\ 2}}\over{\fbox{$\mathrm{Compute}\ \prod_{k=1}^{n}\lambda_k^{\mu_k+1}$}}}\ $$

$$\underset{\downarrow}{\buildrel{\\mathrm{Step\ 3}}\over{\fbox{$\mathrm{Main\ Sum}\ $}}}\\ $$

$$\buildrel{\\mathrm{Step\ 4}}\over{\fbox{$\mathrm{Product\ of\ the\ factors\ }\ $}}\\ $$

The "Main" sum, in turns, involves a second sum where the $\chi_{i,j}$ function is involved. The corresponding code is given as follows:

$$\sum_{i=1}^{n}\prod_{j=1,j\neq i}^{n}\frac{\exp(-\lambda_it)}{\left(\lambda_j-\lambda_i\right)^{\mu_j+1}}\cdot\ \sum_{l=0}^{\mu_i}\frac{t^l}{l!}\chi_{i,\mu_i-l}$$

### 5.1 Example of an application. 
As a first application of the OCS code, the following linear chain, originally proposed by Dreher (2013) will be solved:





