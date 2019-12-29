===================
Mathematical Proofs
===================


Proof of equivalence between the matrix representation of Curie-spin dipole--dipole cross-correlated relaxation (Paramagpy implementation) and the original description by Ghose and Prestegard
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


This proof concerns equations 24-27 of the `article <https://doi.org/10.5194/mr-2019-3>`_ on the software Paramagpy, which describe the cross-correlated relaxation between the Curie spin and dipole--dipole relaxation mechanisms. The equations 24-27 follow the matrix representation of relaxation theory (equations 20-21) by `Suturina et al. <https://doi.org/10.1039/C8CP01332B>`_ . It is shown below that these equations are equivalent to the equations established previously by `Ghose and Prestegard <https://doi.org/10.1006/jmre.1997.1227>`_ and reported in slightly modified form by `Bertini et al. <https://doi.org/10.1006/jmre.2001.2378>`_. The proof begins with equations 24-27 of the manuscript and describes their rearrangement into the form given by `Bertini et al. <https://doi.org/10.1006/jmre.2001.2378>`_

For a specific example, we consider the case of a :sup:`15`\N - :sup:`1`\H group, with the Curie-spin shielding tensor :math:`\sigma` at the site of the :sup:`1`\H spin located at :math:`\vec r=[x,y,z]` and distance :math:`r` from the paramagnetic centre.

.. math::

    \sigma = \frac{1}{4\pi r^5}
    \begin{bmatrix}
    (3x^2-r^2) & 3xy & 3xz\\
    3xy & (3y^2-r^2) & 3yz\\
    3xz & 3yz & (3z^2 - r^2)
    \end{bmatrix}.
    \begin{bmatrix}
    \chi_{xx} & \chi_{xy} & \chi_{xz}\\
    \chi_{xy} & \chi_{yy} & \chi_{yz}\\
    \chi_{xz} & \chi_{yz} & \chi_{zz}
    \end{bmatrix}

We  choose the orientation of the electron--:sup:`1`\H vector to be aligned with the :math:`z`\-axis such that :math:`z\rightarrow r`\. In the case of an isotropic magnetic susceptibility, :math:`x=y=0` and the :math:`\chi` tensor is represented by a diagonal matrix with three identical elements :math:`\chi_\text{iso}`\, yielding the following simplification

.. math::

    \sigma &= \frac{1}{4\pi r^5}
    \begin{bmatrix}
    -r^2 & 0 & 0\\
    0 & -r^2 & 0\\
    0 & 0 & 2r^2
    \end{bmatrix}.
    \begin{bmatrix}
    \chi_\text{iso} & 0 & 0\\
    0 & \chi_\text{iso} & 0\\
    0 & 0 & \chi_\text{iso}
    \end{bmatrix}\\
    &= \xi_\text{DSA}
    \begin{bmatrix}
    -r^2 & 0 & 0\\
    0 & -r^2 & 0\\
    0 & 0 & 2r^2
    \end{bmatrix}\\
    \text{where} \qquad \xi_\text{DSA}&=\frac{\chi_\text{iso}}{4\pi r^5}

The nuclear dipole shielding tensor arising from the :sup:`15`\N spin can be described in the same coordinate frame for an arbitrary orientation of the bond vector :math:`\vec r_\text{HN}=[x,y,z]` with bond length :math:`r_\text{HN}` by

.. math::

    \sigma_N &= \frac{1}{B_0}\frac{\mu_0}{4\pi}\gamma_N\hbar I \left[ 3\frac{\vec r_\text{HN}\otimes \vec r_\text{HN}^T}{r_\text{HN}^5} - \frac{\mathbb{I}_3}{r_\text{HN}^3} \right]\\
    \sigma_N &= \xi_\text{DD}
    \begin{bmatrix}
    (3x^2-r_\text{HN}^2) & 3xy & 3xz\\
    3xy & (3y^2-r_\text{HN}^2) & 3yz\\
    3xz & 3yz & (3z^2 - r_\text{HN}^2)
    \end{bmatrix}\\
    \text{where} \qquad \xi_\text{DD} &= \frac{1}{B_0} \frac{\mu_0}{4\pi} \frac{\gamma_N\hbar I}{r_\text{HN}^5}

and :math:`x`, :math:`y`, :math:`z` denote the coordinates of the :sup:`15`\N spin relative to the :sup:`1`\H spin.

The effective shielding tensor at the site of the :sup:`1`\H spin, when the :sup:`15`\N partner is in the spin-up state, is given by the sum of the two tensors

.. math::

    \sigma_\uparrow &= \sigma + \sigma_N\\
    &=
    \begin{bmatrix}
    (3x^2-r_\text{HN}^2)\xi_\text{DD} -(r^2)\xi_\text{DSA}  & (3xy)\xi_\text{DD} & (3xz)\xi_\text{DD}\\
    (3xy)\xi_\text{DD} & (3y^2-r_\text{HN}^2)\xi_\text{DD} -(r^2)\xi_\text{DSA} & (3yz)\xi_\text{DD}\\
    (3xz)\xi_\text{DD} & (3yz)\xi_\text{DD} & (3z^2 - r_\text{HN}^2)\xi_\text{DD} +(2r^2)\xi_\text{DSA}
    \end{bmatrix}

Note that this matrix is symmetric. Therefore we can ignore equation 18 of the main text and only need to substitute matrix elements into equation 19. Expanding and simplifying (via symbolic processing in the program Mathematica), this yields

.. math::

    \Delta(\sigma_\uparrow)^2={}& \sigma_{xx}^2+\sigma_{yy}^2+\sigma_{zz}^2
    -\sigma_{xx}\sigma_{yy}-\sigma_{xx}\sigma_{zz}-\sigma_{yy}\sigma_{zz} \notag\\
    &+\frac{3}{4}\left[
    ((\sigma_{xy}+\sigma_{yx})^2 + (\sigma_{xz}+\sigma_{zx})^2 + (\sigma_{yz}+\sigma_{zy})^2\right]\\
    ={}& 
    \left((3x^2-r_\text{HN}^2)\xi_\text{DD} -(r^2)\xi_\text{DSA}\right)^2 \notag\\
    &+\left((3y^2-r^2)\xi_\text{DD} -(r^2)\xi_\text{DSA}\right)^2 \notag\\
    &+\left((3z^2 - r_\text{HN}^2)\xi_\text{DD} +(2r^2)\xi_\text{DSA}\right)^2 \notag\\ 
    &-\left((3x^2-r_\text{HN}^2)\xi_\text{DD} -(r^2)\xi_\text{DSA}\right)
    \left((3y^2-r^2)\xi_\text{DD} -(r^2)\xi_\text{DSA}\right) \notag\\
    &-\left((3x^2-r_\text{HN}^2)\xi_\text{DD} -(r^2)\xi_\text{DSA}\right)
    \left((3z^2 - r_\text{HN}^2)\xi_\text{DD} +(2r^2)\xi_\text{DSA}\right) \notag\\
    &-\left((3y^2-r^2)\xi_\text{DD} -(r^2)\xi_\text{DSA}\right)
    \left((3z^2 - r_\text{HN}^2)\xi_\text{DD} +(2r^2)\xi_\text{DSA}\right) \notag\\
    &+\frac{3}{4}\left[
    (6xy\xi_\text{DD})^2 
    + (6xz\xi_\text{DD})^2 
    + (6yz\xi_\text{DD})^2
    \right]\\
    ={}& 9\xi_\text{DD}^2r_\text{HN}^4 + 9\xi_\text{DSA}^2r^4 - 
    9\xi_\text{DD}(r_\text{HN}^2-3z^2)\xi_\text{DSA}r^2

The angle :math:`\theta` between the electron-nuclear vector :math:`\vec r` and the nuclear bond vector :math:`\vec r_\text{HN}` is captured by the dot product formula


.. math::

    \vec r.\vec r_\text{HN} &= |\vec r||\vec r_\text{HN}| \cos{\theta}\\
    \implies \qquad [0,0,r].[x,y,z] &= r.r_\text{HN} \cos{\theta}\\
    \implies \qquad  z &= r_\text{HN} \cos{\theta}

Using above equation to substitute $z$ yields

.. math::

    \Delta(\sigma_\uparrow)^2= 9\xi_\text{DD}^2r_\text{HN}^4 + 9\xi_\text{DSA}^2r^4 - 9\xi_\text{DD}r_\text{HN}^2(1-3\cos^2{\theta})\xi_\text{DSA}r^2

where the first two terms account for the dipolar and Curie spin auto-relaxation terms respectively, and the last term accounts for their cross-correlation. The :math:`R_2` relaxation rate can be calculated by substitution of :math:`\Delta^2` into equation 21 of the main text.

.. math::

    R_2^\text{Curie}(\sigma_\uparrow) = \frac{1}{45}\omega^2 \left[4\mathbf{J}(0) + 3\mathbf{J}(\omega)\right]
    \left( 9\xi_\text{DD}^2r_\text{HN}^4 + 9\xi_\text{DSA}^2r^4 + 
    9\xi_\text{DD}r_\text{HN}^2(3\cos^2{\theta}-1)\xi_\text{DSA}r^2 \right)

The same derivation for :math:`\sigma_\downarrow` yields the same result except for a sign change in the cross term:

.. math::

    R_2^\text{Curie}(\sigma_\downarrow) = \frac{1}{45}B_0^2\gamma_H^2 \left[4\mathbf{J}(0) + 3\mathbf{J}(\omega)\right]
    \left( 9\xi_\text{DD}^2r_\text{HN}^4 + 9\xi_\text{DSA}^2r^4 - 
    9\xi_\text{DD}r_\text{HN}^2(3\cos^2{\theta}-1)\xi_\text{DSA}r^2 \right)

Taking the difference we obtain

.. math::

    R_2^\text{Curie$\times$DD} &= R_2^\text{Curie}(\sigma_\uparrow) - R_2^\text{Curie}(\sigma_\downarrow)\\
    &= \frac{1}{45}B_0^2\gamma_H^2 \left(18\xi_\text{DD}r_\text{HN}^2(3\cos^2{\theta}-1)\xi_\text{DSA}r^2 \right)\left[4\mathbf{J}(0) + 3\mathbf{J}(\omega)\right]\\
    &=\frac{18}{45} \frac{\mu_0}{4\pi} \frac{B_0\gamma_H^2\gamma_N\hbar I}{r_\text{HN}^3}\frac{\chi_\text{iso}}{4\pi r^3}(3\cos^2{\theta}-1)\left[4\mathbf{J}(0) + 3\mathbf{J}(\omega)\right]

Substituting :math:`\chi_\text{iso}` and the spin of :sup:`15`\N as :sup:`I=1/2` yields

.. math::

    R_2^\text{Curie$\times$DD}&=\frac{9}{45} \left(\frac{\mu_0}{4\pi}\right)^2 \frac{B_0\gamma_H^2\gamma_N\hbar}{r_\text{HN}^3} \frac{\mu_B^2g^2J(J+1)}{3 k T r^3}(3\cos^2{\theta}-1)\left[4\mathbf{J}(0) + 3\mathbf{J}(\omega)\right]\\
    &=\frac{2}{15} \left(\frac{\mu_0}{4\pi}\right)^2 \frac{B_0\gamma_H^2\gamma_N\hbar}{r_\text{HN}^3} \frac{\mu_B^2g^2J(J+1)}{k T r^3}\frac{(3\cos^2{\theta}-1)}{2}\left[4\mathbf{J}(0) + 3\mathbf{J}(\omega)\right]

The differential line width can be calculated from the relaxation rate as :math:`\Delta\nu = R_2^\text{Curie$\times$DD}/\pi` and thus this equation matches equation 7 from `Bertini et al. <https://doi.org/10.1006/jmre.2001.2378>`_.




