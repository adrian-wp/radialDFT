# Theory

General Kohn-Sham equation
$$ \left[ -\frac{1}{2} \nabla^2 + V_\mathrm{eff}(\mathbf{r}) \right] \psi_i(\mathbf{r}) = \varepsilon_i \psi_i(\mathbf{r}) $$

Spherical Symmetry
$$ \psi_{nlm}(\mathbf{r}) = R_{nl}(r) Y_{lm}(\theta, \phi) \qquad R_{nl}(r) = \frac{u_{nl}(r)}{r} $$
$$ \nabla^2 \psi_{nlm} = \left( \frac{\mathrm{d}^2}{\mathrm{d}r^2} - \frac{l(l+1)}{r^2} \right) \frac{u_{nl}(r)}{r} Y_{lm}(\theta, \phi) $$

Radial Kohn-Sham equation
$$ \left[ -\frac{1}{2} \frac{\mathrm{d}^2}{\mathrm{d}r^2} + \frac{l(l+1)}{2 r^2} + V_\mathrm{eff}(r) \right] u_{nl}(r) = \varepsilon_{nl} u_{nl}(r) $$

Radial Poisson Equation for the Hartree Potential
$$ \frac{1}{r^2} \frac{\mathrm{d}}{\mathrm{d}r} \left( r^2 \frac{\mathrm{d}}{\mathrm{d}r} V_\mathrm{H}(r) \right) = -4 \pi \rho(r) $$

Electron Density
$$ \rho(r) = \sum_{n,l} \frac{f_{nl}}{4 \pi r^2} |u_{nl}(r)|^2 $$

### Logarithmic Grid

Substitution
$$ x = \ln r \qquad r = e^x \qquad u(r) = u(e^x) = \tilde u(x) $$

Derivatives
$$ \frac{\mathrm{d}u}{\mathrm{d}r} = \frac{\mathrm{d}\tilde u}{\mathrm{d}x} \cdot \frac{\mathrm{d}x}{\mathrm{d}r} = \frac{1}{r} \cdot \frac{\mathrm{d} \tilde u}{\mathrm{d}x} $$
$$ \frac{\mathrm{d}^2u}{\mathrm{d}r^2} = \frac{\mathrm{d}}{\mathrm{d}r} \left( \frac{1}{r} \cdot \frac{\mathrm{d} \tilde u}{\mathrm{d}x} \right) = \frac{1}{r^2} \left( \frac{\mathrm{d}^2 \tilde u}{\mathrm{d}x^2} - \frac{\mathrm{d} \tilde u}{\mathrm{d}x}\right) $$

Kohn-Sham equation
$$ \left[ -\frac{1}{2} \left( \frac{\mathrm{d}^2}{\mathrm{d}x^2} - \frac{\mathrm{d}}{\mathrm{d}x} \right) + \frac{l(l+1)}{2} + e^{2x} V_\textrm{eff}(e^x) \right] \tilde u_{nl}(x) = e^{2x} \varepsilon_{nl} \tilde u_{nl}(x) $$

Poisson Equation
$$ \left( \frac{\mathrm{d}^2}{\mathrm{d}x^2} + \frac{\mathrm{d}}{\mathrm{d}x} \right) \tilde V_\mathrm{H}(x) = -4 \pi e^{2x} \rho(e^x)$$

Transformation to remove first derivative and keep stencil symmetric
$$ \chi(x) = e^{-\frac{1}{2}x} \tilde u(x) \qquad \tilde u(x) = e^{\frac{1}{2}x} \chi(x) $$
$$ \tilde u'(x) = e^{\frac{1}{2}x} \left( \chi'(x) + \frac{1}{2} \chi(x) \right) $$
$$ \tilde u''(x) = e^{\frac{1}{2}x} \left( \chi''(x) + \chi'(x) + \frac{1}{4} \chi(x) \right) $$
$$ -\frac{1}{2} \chi_{nl}''(x) + \left[ \frac{1}{8} + \frac{l(l+1)}{2} + e^{2x} V_\textrm{eff}(e^x) \right] \chi_{nl}(x) = e^{2x} \varepsilon_{nl} \chi_{nl}(x) $$

Boundary Conditions
$$ r \rightarrow \infty \qquad u(r) \rightarrow 0 $$
$$ r \rightarrow 0 \qquad u(r) = C r^{l+1} = C e^{(l+1)x} $$
$$ \tilde u'(x) = C (l+1) e^{(l+1)x} = e^{\frac{1}{2}x} \left( \chi'(x) + \frac{1}{2} \chi(x) \right) $$
$$ \chi'(x) = C (l+1) e^{(l+\frac{1}{2})x} - \frac{1}{2} \chi(x) $$
$$ \chi'(x_\textrm{min}) = (l + \frac{1}{2}) \chi(x_\textrm{min}) $$

### Numerical Solution

Finite Differences

Greens Function

### VWN functional

libxc's implementation was used for reference[^1].

No spin polarization for LDA:

$$ r_s = \left( \frac{3}{4 \pi n(r)} \right)^\frac{1}{3} \qquad \zeta = \frac{n_\uparrow - n_\downarrow}{n} = 0 \qquad f(\zeta) = 0 $$

Exchange part:

$$ \varepsilon_x(r_s, 0) = \varepsilon^\mathrm{P}_x(r_s, 0) = -\frac{3}{2}\pi \left( \frac{4}{9} \pi \right)^\frac{1}{3} r_s $$

Correlation part:

$$ \varepsilon_c(r_s, 0) = \varepsilon^\mathrm{P}_c(r_s, 0) = A \left\{ \ln \frac{x^2}{X(x)} + \frac{2b}{Q} \arctan \frac{Q}{2x+b} - \frac{bx_0}{X(x_0)} \left[ \ln\frac{(x-x_0)^2}{X(x)} + \frac{2(b + 2x_0)}{Q} \arctan \frac{Q}{2x + b} \right] \right\} $$

Helper functions:

$$ x = \sqrt{r_s} \qquad X(x) = x^2 + bx + c \qquad Q = (4c - b^2)^\frac{1}{2} $$

Constants (paramagnetic case):

$$ A^\mathrm{P} = 0.0621814 \qquad x_0 = -0.10498 \qquad b = 3.72744 \qquad c = 12.9352 $$

Energy and potential in radial coordinates:

$$ E_{xc}[n(r)] = \int_0^\infty 4 \pi r^2 n(r) \varepsilon_{xc}[n(r)] \, \mathrm{d}r $$

$$ V_{xc}(r) = \varepsilon_{xc}[n(r)] + n(r) \frac{\mathrm{d} \varepsilon_{xc}[n(r)]}{\mathrm{d}n} \qquad \frac{\mathrm{d} \varepsilon_{xc}}{\mathrm{d}n} = \frac{\mathrm{d} \varepsilon_{xc}}{\mathrm{d} r_s} \cdot \left( -\frac{r_s}{3 n} \right) $$

$$ V_{xc}(r_s) = \varepsilon_{xc} - \frac{r_s}{3} \frac{\mathrm{d} \varepsilon_{xc}}{\mathrm{d}r_s} $$

[^1]: https://github.com/ElectronicStructureLibrary/libxc/blob/4bd0e1e36347c6d0a4e378a2c8d891ae43f8c951/maple/vwn.mpl