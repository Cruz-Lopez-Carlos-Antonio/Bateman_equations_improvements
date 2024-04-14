# Bateman equations 
The present repository contains a set of Python algorithms that solve and improve the General Solution of the Bateman Equations, which were developed in the following article **"General Solution of Bateman Equations using Cauchy Products and the Theory of Divided Differences"**, which was submitted to the **"Annals of Nuclear Energy"** journal.

The programs are licensed under a Creative Commons Attribution 4.0 International License: http://creativecommons.org/licenses/by/4.0/

Authors: Carlos-Antonio Cruz-López (cacl.nucl@gmail.com), Gilberto Espinosa-Paredes (gepe@xanum.uam.mx)

Mathematical and algorithmical generalities of the codes are described in the following lines with the purpose to provide some insight of the developed work. Nevertheless, a more detailed and precise discussion is provided in the submitted article.

## Financial Support.
The authors appreciate the financial support received from the Consejo Nacional de Ciencia y Tecnología, CONACYT, under the program “Estancias Posdoctorales por México, 2022”, with the project entitled: “Desarrollo de modelos fenomenológicos energéticos de orden fraccional, para la optimización y simulación en reactores nucleares de potencia”, by which the present development was possible.

## Software specifications.
The AnalyticNPKE codes were written in the Python programming language in its version 3. The codes require the following libraries:

- [x] itertools
- [x] numpy
- [x] decimal
## Index of the Repository
1. [Mathematical description of the problem](https://github.com/Cruz-Lopez-Carlos-Antonio/Bateman_equations_improvements/blob/main/README.md#1-mathematical-description-of-the-problem).
   - [1.1 Differential mass-balance equations](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticNPKE?tab=readme-ov-file#51-sums)
   - [1.2 Laplace transform of the system.](https://github.com/Cruz-Lopez-Carlos-Antonio/Bateman_equations_improvements#12-laplace-transform-of-the-system)
   - [1.3 Partial fraction decomposition.](https://github.com/Cruz-Lopez-Carlos-Antonio/Bateman_equations_improvements/blob/main/README.md#13-partial-fraction-decomposition)
   - [1.4 Solution of Bateman Equations.](https://github.com/Cruz-Lopez-Carlos-Antonio/Bateman_equations_improvements/blob/main/README.md#14-solution-of-the-bateman-equations)
2. [Cetnar's solution.](https://github.com/Cruz-Lopez-Carlos-Antonio/Bateman_equations_improvements/blob/main/README.md#1-mathematical-description-of-the-problem)
   - [2.1 Solution based on increments.](https://github.com/Cruz-Lopez-Carlos-Antonio/Bateman_equations_improvements/blob/main/README.md#21-solution-based-on-increments)
   - [2.2 Issues related to nested sums.](https://github.com/Cruz-Lopez-Carlos-Antonio/Bateman_equations_improvements/blob/main/README.md#22-the-issues-of-the-nested-sums)

## 1. Mathematical description of the problem.
### 1.1. Differential mass-balance equations. 
Bateman equations describe the time evolution of a set of nuclides in succesive transformations due to decay and transmutation process. These transformations can be represented by the following elementary structure:
$$X_1\buildrel\lambda_1\over\rightarrow X_2\buildrel\lambda_2\over\rightarrow\ldots\buildrel\lambda_{n-1}\over\rightarrow X_n\buildrel\lambda_n\over\rightarrow \tag{1}$$
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
$$X_n\left(t\right)=\frac{X_1\left(0\right)}{\lambda_n}\sum_{i=1}^{n}{\lambda_i\alpha_i\exp(-\lambda_it)}\cdot\sum_{m=0}^{\mu_i}{\frac{\left(\lambda_it\right)^m}{m!}\cdot\psi_{i,\mu_i-m}}, \tag{20}$$
where:
$$\alpha_i=\prod_{j=1,j\neq i}^{n}\left(\frac{\lambda_j}{\lambda_j-\lambda_i}\right)^{m_j}, \tag{21}$$
and:

$$\psi_{i,j}=\sum_{h_1=0}^{j} \sum_{h_2=0}^{j}\cdots\sum_{h_{i-1}=0}^{j}\sum_{h_{i+1}=0}^{j} \cdots\sum_{h_n=0}^{j} \prod_{k=1,k\neq i}^{n} \binom{h_k+\mu_k}{\mu_k} \left(\frac{\lambda_i}{\lambda_i-\lambda_k}\right)^{h_k}\delta_{j,p} \tag{22} $$

where:

$$p=\sum_{l=1,l\neq i}^{n}h_l,  \ \ \ \ \delta_i,p=\begin{cases}
  1  &  \text{ if $j=p$ } \\
  0  &  \text{ if $j\neq p$}
\end{cases} \tag{23}$$
### 2.2 Issues related to nested sums.
Even when the solution (as well as the procedure that was followed to obtain it) is very ingenious, the way in which it is expressed is very disadvantageous. The main issue is related to the nested sum in the term $\psi_{i,j}$ in Eq. (22), which in turns involves the Kronecker's delta given in Eq. (23). 
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
### 3.1 Reduction of the product inside $\psi_{i,j}$
The main contribution of the submitted article consists of simplifying the Cetnar's solution, removing the nested sums as well as the Kronecker's delta. Using the findings that were described in the last section, which are based on Cauchy products, it is possible to reduce the following product that appears inside Eq. (27) as follows:
$$\prod_{k=1,k\neq i}^{n}\left(\frac{\lambda_i}{\lambda_i-\lambda_k}\right)^{h_k}=\prod_{k=1,k\neq i}^{n}\left(\frac{1}{\lambda_i-\lambda_k}\right)^{h_k}\prod_{l=1,l\neq i}^{n}\lambda_i^{h_l}$$ 
$$\prod_{k=1,k\neq i}^{n}\left(\frac{1}{\lambda_i-\lambda_k}\right)^{h_k}\lambda_i^{h_1+h_2+\ldots+h_{i-1}+h_{i+1}+\ldots+h_l} \tag{28}$$
$$=\lambda_i^j\prod_{k=1,k\neq i}^{n}\left(\frac{1}{\lambda_i-\lambda_k}\right)^{h_k}, \tag{29}$$
where it has been used that:
$$h_1+h_2+\ldots h_{i-1}+h_i+\ldots+h_n=j \tag{30}$$
### 3.2 Including Diophantine Equations.
Using the procedure that is detailed in the submitted paper, the Cetnar's solution can be rewritten as:
$$X_n\left(t\right)=\frac{X_1\left(0\right)}{\lambda_n}\sum_{i=1}^{n}\lambda_i\prod_{j=1,j\neq i}^{n}{\frac{\lambda_j^{\mu_j+1}}{\left(\lambda_j-\lambda_i\right)^{\mu_j+1}}\exp(-\lambda_it)}\sum_{l=0}^{\mu_i}{\frac{\left(\lambda_it\right)^l}{l!}\psi_{i,}\mu_i-l} \tag{30}$$





