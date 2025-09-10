# Theory

The equations in this file were mainly used as reference for the implementation and are therefore not explained in much
detail.

General Kohn-Sham equation

$$ \left[ -\frac{\hbar^2}{2m} \nabla^2 + V_\mathrm{eff}(\mathbf{r}) \right] \psi_i(\mathbf{r}) = \varepsilon_i \psi_i(\mathbf{r}) $$

Using atomic units

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

$$ w(x) = e^{-\frac{1}{2}x} \tilde u(x) \qquad \tilde u(x) = e^{\frac{1}{2}x} w(x) $$
$$ \tilde u'(x) = e^{\frac{1}{2}x} \left( w'(x) + \frac{1}{2} w(x) \right) $$
$$ \tilde u''(x) = e^{\frac{1}{2}x} \left( w''(x) + w'(x) + \frac{1}{4} w(x) \right) $$
$$ -\frac{1}{2} w_{nl}''(x) + \left[ \frac{1}{8} + \frac{l(l+1)}{2} + e^{2x} V_\textrm{eff}(e^x) \right] w_{nl}(x) = e^{2x} \varepsilon_{nl} w_{nl}(x) $$

Dirichlet boundary condition on right side

$$ r \rightarrow \infty \qquad u(r) \rightarrow 0 $$
$$ w(x_\textrm{max}) = 0 $$

Neumann boundary condition on left side

$$ r \rightarrow 0 \qquad u(r) = C r^{l+1} = C e^{(l+1)x} $$
$$ \tilde u'(x) = C (l+1) e^{(l+1)x} = e^{\frac{1}{2}x} \left( w'(x) + \frac{1}{2} w(x) \right) $$
$$ w'(x) = C (l+1) e^{(l+\frac{1}{2})x} - \frac{1}{2} w(x) $$
$$ w'(x_\textrm{min}) = (l + \frac{1}{2}) w(x_\textrm{min}) $$

### Numerical Solution

Grid (equally spaced in x, logarithmic in r)

$$ x_0 = e^{r_\textrm{min}} \qquad x_{N-1} = e^{r_\textrm{max}} \qquad x_i = x_0 + i h_x $$
$$ w_i = w(x_i)$$

Finite Differences (second derivative, second-order accuracy)

$$ w''(x_i) = \frac{w_{i-1} - 2 w_i + w_{i+1}}{h_x^2}$$ 

Transformed radial Kohn-Sham equation as general symmetric eigenvalue problem, where A is tridiagonal and B is diagonal

$$ \mathbf A \mathbf w = \varepsilon \mathbf B \mathbf w $$
$$ a_{ii} = \frac{1}{h_x^2} + \frac{1}{8} + \frac{l (l + 1)}{2} + e^{2x_i} V_\textrm{eff}(x_i) $$
$$ a_{i-1i} = a_{i+1i} = -\frac{1}{2h_x^2} $$
$$ b_{ii} = e^{2x_i} $$

Neumann Boundary Condition (divide first row by 2 to keep A symmetric)

$$ a_{00} = \frac{1}{2} \left( \frac{1}{h_x^2} + \frac{l + \frac{1}{2}}{h_x} + \frac{1}{8} + \frac{l (l + 1)}{2} + e^{2x_i} V_\textrm{eff}(x_i) \right) $$
$$ b_{00} = \frac{1}{2} e^{2x_0} $$

The tridiagonal eigenvalue solver from LAPACK does not work with a right side B

$$ \mathbf C = \mathbf B^{-\frac{1}{2}} \mathbf A \mathbf B^{-\frac{1}{2}} $$
$$ \mathbf C \mathbf y = \varepsilon \mathbf y $$
$$ \mathbf w = \mathbf B^{-\frac{1}{2}} \mathbf y $$

Green's Function solution for Hartree potential (integration with Simpson's rule)

$$ V_H(r) = \tilde V_H(x) = 4 \pi \left( e^{-x} \int_{-\infty}^x e^{3x'} \rho(e^{x'}) \, \mathrm d x' + \int_x^\infty e^{2x'} \rho(e^{x'}) \, \mathrm d x' \right)$$

### VWN functional

Based on the original VWN paper[^1]. The maple code of libxc was also used as reference[^2] during implementation. The
definition of x was taken directly from the paper and has nothing to do with the logarithmic grid.

No spin polarization for LDA

$$ r_s = \left( \frac{3}{4 \pi n(r)} \right)^\frac{1}{3} \qquad \zeta = \frac{n_\uparrow - n_\downarrow}{n} = 0 \qquad f(\zeta) = 0 $$

Exchange part

$$ \varepsilon_x(r_s, 0) = \varepsilon^\mathrm{P}_x(r_s, 0) = -\frac{3}{2}\pi \left( \frac{4}{9} \pi \right)^\frac{1}{3} r_s $$

Correlation part

$$ \varepsilon_c(r_s, 0) = \varepsilon^\mathrm{P}_c(r_s, 0) = A \left[ \ln \frac{x^2}{X(x)} + \frac{2b}{Q} \arctan \frac{Q}{2x+b} - \frac{bx_0}{X(x_0)} \left( \ln\frac{(x-x_0)^2}{X(x)} + \frac{2(b + 2x_0)}{Q} \arctan \frac{Q}{2x + b} \right) \right] $$

Helper functions

$$ x = \sqrt{r_s} \qquad X(x) = x^2 + bx + c \qquad Q = (4c - b^2)^\frac{1}{2} $$

Derivative (could be simplified further but some repetitive expressions were replaced with helper functions in the code)

$$ \varepsilon_x'(r_s) = \frac{3}{4} \left( \frac{3}{2\pi} \right)^{\frac{2}{3}} \frac{1}{r_s^2} $$
$$ \varepsilon_c'(r_s) = A \left(\frac{1}{x^2} + \left( \frac{b x_0}{X(x_0)} - 1 \right) \frac{\tilde X'(x)}{X(x)} - \frac{b x_0}{X(x_0)} \frac{1}{x^2 - x_0 x} - \left( \frac{2b}{Q} - \frac{b x_0}{X(x_0)} \frac{2 (b + 2 x_0)}{Q} \right) \frac{Q}{x (Q^2 + (2x + b)^2)} \right) $$
$$ \tilde X'(x) = \frac{\mathrm d X}{\mathrm d r_s} = 1 + \frac{1}{2} b r_s^{-\frac{1}{2}} = 1 + \frac{b}{2x} $$

Constants (paramagnetic case)

$$ A^\mathrm{P} = 0.0621814 \qquad x_0 = -0.10498 \qquad b = 3.72744 \qquad c = 12.9352 $$

Energy and potential in radial coordinates

$$ E_{xc}[n(r)] = \int_0^\infty 4 \pi r^2 n(r) \varepsilon_{xc}[n(r)] \, \mathrm{d}r $$
$$ V_{xc}(r) = \varepsilon_{xc}[n(r)] + n(r) \frac{\mathrm{d} \varepsilon_{xc}[n(r)]}{\mathrm{d}n} \qquad \frac{\mathrm{d} \varepsilon_{xc}}{\mathrm{d}n} = \frac{\mathrm{d} \varepsilon_{xc}}{\mathrm{d} r_s} \cdot \left( -\frac{r_s}{3 n} \right) $$
$$ V_{xc}(r_s) = \varepsilon_{xc} - \frac{r_s}{3} \frac{\mathrm{d} \varepsilon_{xc}}{\mathrm{d}r_s} $$

[^1]: https://doi.org/10.1139/p80-159
[^2]: https://github.com/ElectronicStructureLibrary/libxc/blob/4bd0e1e36347c6d0a4e378a2c8d891ae43f8c951/maple/vwn.mpl