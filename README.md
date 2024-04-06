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

## Mathematical description of the problem.
Bateman equations describe the time evolution of a set of nuclides appearing in a linear chain. This last structure can be represented as:
$$X_1\ \rightarrow\above\lambda_1\ \ X_2\ \rightarrow\above\lambda_2\ \ldots\rightarrow\above\lambda_(n-1)\ \ X_n\ \rightarrow\above\lambda_n\ \tag{1}$$
