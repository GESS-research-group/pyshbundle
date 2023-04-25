# Convert Data Formats

Spherical harmonic functions or coefficients, Legendre functions and their derivatives can be arranged in different ways. There are multiple functions in SHBundle for reordering from one format to another. Some of them have been translated to Python in PySHBundle. Couple of new ones have also been added.

## Spherical Harmonics Data Formats

### clm-format

This is a standard format to store spherical harmonic coefficients in the indexed column-vector-format (abbreviatedL clm-format)

\begin{equation}
  \left( n, m, \overline{C}_{n, m}, \overline{S}_{n, m}, \left[ \sigma_{\overline{C}_{n, m}}, \sigma_{\overline{S}_{n, m}} \right] \right)
\end{equation}

The first column represents the degree $n$, the second column represents the order $m$ (both n,m are integers), followed by the coefficients $\overline{C}_{n, m}, \overline{S}_{n, m}$ and the last two columns contain their respective standard deviations $\sigma_{\overline{C}_{n, m}}, \sigma_{\overline{S}_{n, m}}$

### klm-format

This is a variation of the clm-format for compact notation with just 3 or 4 columns. The coefficients are sorted first w.r.t. degree and then the order, particularly the sine-coefficients are arranged starting first with negative orders. The following matrix represents the klm-format:

```math
\begin{bmatrix}
  0 & 0 & \overline{C}_{0, 0} & \sigma_{\overline{C}_{0, 0}} \\
  0 & 0 & \overline{C}_{0, 0} & \sigma_{\overline{C}_{0, 0}} \\
  0 & 0 & \overline{C}_{0, 0} & \sigma_{\overline{C}_{0, 0}} \\
  0 & 0 & \overline{C}_{0, 0} & \sigma_{\overline{C}_{0, 0}} \\
  0 & 0 & \overline{C}_{0, 0} & \sigma_{\overline{C}_{0, 0}} \\
  0 & 0 & \overline{C}_{0, 0} & \sigma_{\overline{C}_{0, 0}} \\
  0 & 0 & \overline{C}_{0, 0} & \sigma_{\overline{C}_{0, 0}} \\
  0 & 0 & \overline{C}_{0, 0} & \sigma_{\overline{C}_{0, 0}} \\
  0 & 0 & \overline{C}_{0, 0} & \sigma_{\overline{C}_{0, 0}} \\
  \vdots & & & \vdots \\
  N_{max} & N_{max} & \overline{C}_{N_{max}, N_{max}} & \sigma_{\overline{C}_{N_{max}, N_{max}}} \\
\begin{bmatrix}
```

### $\left | C \backslash S \right |$-format

This is another well known arrangement of Spherical Harmonic coefficients. This is a square matrix of size $n_{max}, n_{max}$.

The lower traingular terms are made of the cosine terms


### $\left / S | C \right \backslash$-format

This is yet another popular format where the sine-coefficients are flipped from left to right, to obtain a triangular arrangement which is completed by zeros.

The following figure illustrates the  $\left | C \backslash S \right |$ and  $\left / S | C \right \backslash$ format respectively.

![Spherical Harmonic Formats](img/sh_formats.png)

::: pyshbundle.cs2sc
::: pyshbundle.sc2cs
::: pyshbundle.clm2cs
::: pyshbundle.clm2sc
::: pyshbundle.klm2sc

## Reference
  - Nico Sneeuw, Matthias Weigelt, Markus Antoni, Matthias Roth, Balaji Devaraju, et. al. (2021). SHBUNDLE 2021. http://www.gis.uni-stuttgart.de/research/projects/Bundles.