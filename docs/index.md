# Welcome to PySHbundle


[![image](https://img.shields.io/pypi/v/pyshbundle.svg)](https://pypi.python.org/pypi/pyshbundle)


**PySHbundle: A Python implementation of MATLAB codes SHbundle**


-   Free software: GNU General Public License v3
-   Documentation: <https://mn5hk.github.io/pyshbundle>
    


## How to install <br>
We recommend using Mamba to install required packages <br>
`conda create pysh` <br>
`conda activate pysh` <br>
`conda install -c conda-forge mamba` <br>
`mamba install -c conda-forge numpy pandas netCDF4 scipy xarray julian scipy geopandas matplotlib rasterio salem shapely` <br>

Convert [SHbundle matlab codes](https://www.gis.uni-stuttgart.de/en/research/downloads/shbundle/) to python<br>
[Geodesy for Earth system science (GESS) research Group at ICWaR, IISc](https://ultra-pluto-7f6d1.netlify.app/)
![](https://visitor-badge.glitch.me/badge?page_id=mn5hk.mat2py)

Currently we are upgrading (under process) the package to be implementable on binder. 
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mn5hk/pyshbundle/HEAD)


### Mathematics

In this section, we present a mathematical representation of the spherical harmonics analysis. According to potential theory, the gravitational field of a body fulfils the Laplace equation $\nabla^2\phi = 0$. Laplace's equation in spherical coordinates can be written as follows: 

\begin{equation}
    \frac{1}{r^2}\frac{\partial}{\partial r}\bigg( r^2\frac{\partial \phi}{\partial r}\bigg)  
    +
    \frac{1}{r^2\sin\vartheta}\frac{\partial}{\partial \vartheta}\bigg(\sin\vartheta\frac{\partial \phi}{\partial \vartheta}\bigg) 
    +
    \frac{1}{r^2\sin^2\vartheta}\frac{\partial^2 \phi}{\partial \lambda^2}
    = 0 ,
\end{equation}


where 
$\phi$ is the potential, 
$r$ is the radius, 
$\vartheta$ is the co-latitude and 
$\lambda$ is the longitude. 

We perform a separation of variables and insert $\phi(r, \vartheta, \lambda) =f(r)g(\vartheta)h(\lambda)$ into the Laplace equation to get three independent equations:


\begin{equation}
r^2\frac{d^2f}{dr^2}+2r\frac{df}{dr} - n(n+1)f = 0,
\end{equation}



\begin{equation}
\frac{d^2g}{d\vartheta^2}
+
\frac{dg}{d\vartheta}\cot\vartheta
+
\bigg(  n(n+1) - \frac{m^2}{\sin^2\vartheta}   \bigg) g = 0 ,
\end{equation}



\begin{equation}
\frac{d^2h}{d\lambda^2} + m^2h = 0,
\end{equation}


where $m$ and $n$ are the degree and order respectively. Solving $(2), (3)$ and $(4)$, we obtain: 


\begin{equation}
f(r) \in \{r^n, r^{-(n+1)}\},
\end{equation}



\begin{equation}
g(\vartheta) \in \{P_{n,m}(\cos \vartheta), Q_{n,m}(\cos \vartheta)\} ,
\end{equation}



\begin{equation}
h(\lambda) \in \{\cos m\lambda, \sin m\lambda\}.
\end{equation}\\


Thus, the Laplace equation's solution takes the following form: 


\begin{equation}
\phi(r, \vartheta, \lambda) = \sum_{n=0}^{\infty} \sum_{m=0}^{n} 
\alpha_{n,m}
\begin{Bmatrix}
P_{n,m}(\cos\vartheta)\\
Q_{n,m}(\cos\vartheta)\\
\end{Bmatrix}
\dot{•}
\begin{Bmatrix}
\cos m\lambda\\
\sin m\lambda\\
\end{Bmatrix}
\dot{•}
\begin{Bmatrix}
r^n\\
r^{(n+1)}\\
\end{Bmatrix}
.
\end{equation}


Solutions for $f(r)$ and $h(\lambda)$ are fairly straightforward. Eq - (3) for $g(\vartheta)$ is in the form of a Legendre differential equation and its solutions are $P_{n,m}(\cos \vartheta)$ and $Q_{n,m}(\cos \vartheta)$, the associated Legendre functions of the first and second kind. We now apply two constraints to the solution:

* $\phi \rightarrow 0$ when $r \rightarrow \infty$,
* $\phi$ is limited on the sphere,

which leads us to eliminate $Q_{n,m}(\cos \vartheta)$ and $r^n$.The $4\pi$ normalization of the Associated Legendre functions [8] is utilized in our package and is given by: 


\begin{equation}
\bar{P}_{n,m}(\cos\vartheta) = P_{n,m}(\cos\vartheta)\sqrt{(2-\delta_{m0})(2n+1)\frac{(n-m)!}{(n+m)!}},
\end{equation}

where $\delta_{m0}$ is the Kronecker delta function,

\begin{equation}
P_{n,m}(t) = (1-t^2)^{\frac{m}{2}}\frac{d^mP_n(t)}{dt^m},
\end{equation}

and 

\begin{equation}
nP_n(t)=-(n-1)P_{n-2}(t) + (2n-1)tP_{n-1}(t).
\end{equation}


Spherical harmonics are the angular portion of a set of solutions to Laplace's equation. They take into account $\vartheta$ and $\lambda$. They are functions modelled on the surface of a sphere, denoted by $Y_{n,m}(\vartheta,\lambda)$. They are of three kinds: 

* Zonal harmonics: $m=0$ - they are only latitude dependent,
* Tesseral harmonics: $0 < m < n$, and 
* Sectorial harmonics: $m=n$.

Quantities like the gravitational potential, height of water column, gravity anomaly and so on are the functionals of the gravity field which are obtained by differentiating the potential $\phi$ with respect to the spherical coordinates. 

The gravitational potential anomaly $V$ is given by:


\begin{equation}
    V(r, \vartheta, \lambda) = 
    \frac{GM}{r} \sum_{n=0} ^{N_{max}} \sum_{m=0} ^{n} 
    \left(\frac{R}{r}\right) ^{n+1}
    \bar{P}_{n,m}(\cos \vartheta) [C_{n,m}\cos m\lambda+S_{n,m}\sin m\lambda].
\end{equation}


Here, $R$ refers to the radius of the Earth,
$\bar{P}_ {n,m}$ refers to the Associated Legendre functions with $4\pi$ normalization,
$C_{lm}$ and  $S_{lm}$ refer to the spherical harmonic coefficients. Similarly, another functional, the change in surface mass density, is represented by:


\begin{equation}
    \Delta\sigma(\vartheta, \lambda) = 
    \frac{a\rho_{ave}}{3} 
    \sum_{n=0}^{N_{max}}\sum_{m=0}^{n} 
    \left(\frac{R}{r}\right)^{n+1} 
    \bar{P}_{n,m}(\cos\vartheta)
    \frac{2n+1}{1+k_l}
    [C_{n,m}\cos m\lambda + S_{n,m}\sin m\lambda],
\end{equation}


where
$\rho_{ave}$ refers to the average density of the Earth in $g/cm^3$ and
$k_n$ refers to the load Love number of degree $n$.

## Credits

This package was created with [Cookiecutter](https://github.com/cookiecutter/cookiecutter) and the [giswqs/pypackage](https://github.com/giswqs/pypackage) project template.
