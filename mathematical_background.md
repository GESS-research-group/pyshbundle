# Mathematic Backround

GRACE works on the principal of gravimetric changes. Level 2 GRACE data consists of the spherical harmonic coefficients  <i>$C_{l,m}$</i> and <i>$S_{l,m}$</i>. Gravimetric potential function <i>V ( r, θ, λ )</i> can be represented by the spherical harmonic coefficients in the frequency domain with the help of the following relation `(Vishwakarma, 2017; Kaula, 1996; Chao & Gross, 1987; Wahr et. al., 1998)`:

$
\begin{equation}
    V(r, \theta, \lambda) = 
    \frac{GM}{r} \sum_{l=0} ^ {\infty} 
    \left(\frac{a}{r}\right) ^ {l}
    \sum_{m=0} ^ {l} 
    \bar{P}_{l,m}(\cos \theta)[C_{l,m}\cos m\lambda+S_{l,m}\sin m\lambda],
\end{equation}
$

where <i>G</i> is the Gravitational constant, <i>M</i> represents the total Earth mass, <i>a</i> is the average radius of the Earth, <i>$P_{l,m}$</i>  represents the the fully normalized Legendre functions of the first kind, <i>$C_{l,m}$</i> and <i>$S_{l,m}$</i> represent the fully normalized spherical harmonic coefficients, and <i>l</i> and <i>m</i> represent the degree and order, respectively.

It should be noted that <i>Equation 1</i> does not deal with the variability of gravimetric potential function over time. However, a major application of the GRACE satellite system is to retrieve the time-variable gravity information. This is acheived by taking the variation of the spherical harmonic coefficients over time. To obtain this, a long-term mean of the monthly values of the spherical harmonic coefficients is removed from the monthly spherical harmonic coefficients obtained from GRACE L2 data. This can be denoted by <i>$ \Delta C_{l,m}$</i> and <i>$ \Delta P_{l,m}$</i>. Thus, <i>Equation 1</i> can be modified to obtain the change in gravimetric potential function over time.

Since we are interested in the change in mass in our system, we need to obtain the change in density from the change in gravity potential function. It is further assumed that the redistribution of the mass of the earth takes place within a thin layer close to the surface of the Earth. Furthermore, mass redistribution takes place with a deformation. This is accounted for by the load Love numbers <i>k<sub>l</sub></i> `(Wahr et. al., 1998)`. <i>Equation 1</i> further resolves to:

$
\begin{equation}
    \Delta \sigma (\theta, \lambda) = 
    \frac{a \rho_{avg}}{3} \sum_{l=0} ^ {\infty} 
    \sum_{m=0} ^ {l} 
    \bar{P}_{l,m}(\cos \theta)
    \frac{2 l + 1}{1 + k_{l}}
    [\Delta C_{l,m}\cos m\lambda+ \Delta S_{l,m}\sin m\lambda],
\end{equation}
$

Here, <i>$\Delta \sigma (\theta, \lambda)$</i> represents the change in surface density of the Earth, and <i>$\rho_{avg}$</i> represents the average density of the Earth (<i>5517 kg / m<sup>3</sup></i>).

As the mass redistribution on Earth over a monthly time scale is dominated by the hydroogical processes, the density change <i>$\Delta \sigma (\theta, \lambda)$</i> relates to the <i>Equivalent Water Height (EWH)</i> by: <i>$\Delta \sigma (\theta, \lambda) = EWH (\theta, \lambda) . \rho_{water}$</i>. Thus, <i>equation 2</i> can be rewritten in terms of <i>EWH</i> as:

$
\begin{equation}
    EWH (\theta, \lambda) = 
    \frac{a \rho_{avg}}{3 \rho_{water}} 
    \sum_{l, m}
    \bar{P}_{l,m}(\cos \theta)
    \frac{2 l + 1}{1 + k_{l}}
    [\Delta C_{l,m}\cos m\lambda+ \Delta S_{l,m}\sin m\lambda],
\end{equation}
$

Thus, we can obtain the hydrological parameter <i>EWH</i> from GRACE Level 2 data using <i>equation 3</i>.

The accuracy and precision of the <i>EWH</i> computed depends upon the accuracy and precision of the <i>$\Delta C_{l,m}$</i> and <i>$\Delta P_{l,m}$</i>, obtained from GRACE. However, these GRACE products are both noisy and coarse in resolution `(Wahr et. al., 1998)`. A tradeoff exists between the noise and resolution of the spherical harmonic products. To capture the spherical harmonic products at a higher spatial resolution, their values at higher degree and order needs to be used. However, noise increases with the increase in degree and order, making the computed <i>EWH</i> also noisy. Similarly, if the spherical harmonics are truncated at a lower degree and order, the noise in the computed <i>EWH</i> decreases, however, the spatial resolution of the obtained <i>EWH</i> also reduces.

To improve the signal-to-noise ratio of the obtained <i>EWH</i>, various filtering techniques have been used. An ideal filter retains all of the signal while filtering out all of the noise. A popular filter used for GRACE applications is the Gaussian filter. The weights, <i>w</i>, for the Gaussian spatial averaging is given by:

$
\begin{equation}
    \omega (\psi) = 
    \frac{\beta}{2 \pi} 
    \frac{exp [-\beta (1 - \cos \psi)]}{1 - \exp ^ {-2 \beta}},
\end{equation}
$

where, $\beta = \frac{\ln (2)}{(1 - \cos(\frac{r_{fil}}{a}))'}$. Here, $r_{fil}$ is the averaging radius of the filter. Thus, the Gaussian filter obtained in the spectral domain is written as `(Wahr et. a., 1998)`:

$$
\begin{equation}
\bar{\sigma}(\theta, \lambda) = 
\frac{2 a \rho_{avg} \pi}{3} 
    \sum_{l, m} W_l
    \bar{P}_{l,m}(\cos \theta)
    \frac{2 l + 1}{1 + k_{l}}
    [\Delta C_{l,m}\cos m\lambda+ \Delta S_{l,m}\sin m\lambda],
\end{equation}
$$

<i>Equation 5</i> is similar to <i>equation 3</i>, but for an additional multiplication factor, <i>$W_l$</i>, defined as <i>$W_l = \int_0^\pi w (\psi) P_l (\cos \psi) \sin \psi d\psi$</i> and <i>$P_l = \frac{\bar{P_l}}{\sqrt {2l + 1}}$</i>. 

<i>Equation 5</i> defines a Gaussian filter that decays with only degree. However, for our GRACE spherical harmonics, the decay occurs with the location as well as with the degrees and orders. Thus, <i>equation 5</i> is further generalized as `(Wahr et. al., 1998; Devaraju, 2015)`:
$$
\begin{equation}
\bar{\sigma}(\theta, \lambda) = 
\frac{a \rho_{avg}}{12 \pi} 
    \sum_{l, m}
    \sum_{n, k} W_{lm}^{nk}
    \bar{P}_{l,m}(\cos \theta)
    \frac{2 l + 1}{1 + k_{l}}
    [\Delta C_{l,m}\cos m\lambda+ \Delta S_{l,m}\sin m\lambda],
\end{equation}
$$

where <i>$W_{lm}^{nk}$</i> represents the spectral weight in its general form. <i>Equation 6</i> is the final result we obtain after spectral harmonic synthesis and application of Gaussian filter.
