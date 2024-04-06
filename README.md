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
1. [Mathematical description of the problem](https://github.com/Cruz-Lopez-Carlos-Antonio/AnalyticNPKE#1-mathematical-description-of-the-problem).

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
$$x_n(t)=\mathcal{L}^{-1}\\{{\widetilde{x}}_ n\\}=\prod_{k=1}^{n-1}\lambda_k\sum_{i=1}^{n}{\prod_{j=1,j\neq i}^{n}\frac{1}{(\lambda_j-\lambda_i)}\mathcal{L}^{-1}\\{\frac{1}{s+\lambda_i}}\\}, \tag{15}$$
that using that $\mathcal{L}^{-1}\\{\frac{1}{s+\lambda_i}}\\}$

