cancer integrator
=================

Simulates survival curves in generalised multi-stage clonal expansion models of
Moolgavkar-Vernon-Knudsen type. This project allows us to
define a birth-death-mutation process on an arbitrary directed graph
(representing mutations), and generate survival curves using both Gillespie algorithm simulations and
fast direct integration of the characteristics.

Sources: 

 * R. Meza PNAS 2008; 
 * K.S. Crump, Risk Analysis, Vol. 25, No. 4, 2005
 * E.G. Luebeck et al 2012, doi: 10.1158/0008-5472.CAN-12-2198
 * D. Quinn, Risk Analysis, Vol. 9, Issue 3, 1989, doi: 10.1111/j.1539-6924.1989.tb01006.x

We should see an initial Armitage-Doll type curve with $p \sim t^n$, then a
mean-field regime like $p \sim t^k e^{s t}$ (source: Ivana and me 2020, also Armitage
and Doll 1958, the less famous one), and then after a sojourn time

$$T_s \sim \frac{\ln(s / \mu)}{s}$$

we should see a transition region, followed by an exponential tail (source:
Georg Luebeck 2008).

Integration methods
-------------------

This compares different numerical simulation methods for birth-death-mutation
processes: the classic Gillespie algorithm, and an experimental algorithm based on the method of
characteristics. 

The two methods:

1. Use the Gillespie algorithm for many stochastic, "exact" samples of the
master equation. e.g. for the two-hit model:

$$\Gamma = \mu_0 * N_0 + (\mu_1 + s) * N_1$$

$$x = Unif(0,\Gamma)$$

$$Pr(mut0) = (\mu_0 * N_0) / \Gamma$$

$$Pr(birth) = (s * N_1) / \Gamma$$

$$Pr(mut1) = (\mu_1 * N_1) / \Gamma $$

$$t += Exp(0, 1/\Gamma)$$

Repeat while $t < $ maximum time.

2. Use the method of characteristics to numerically integrate the generating function:
  * Generating function $\Psi =$ the Fourier transform of the probability
    distribution $P$, a function of the conjugate variables $q_j$
  * The $q_j$ are conjugate variables to the populations $N_j$
  * The Kolmogorov forward equation is of the form

$$\frac{\partial\Psi}{\partial t} = \vec{X} \cdot \nabla_{\vec{q}} \Psi + Y(\vec{q} \Psi$$
with $X_j(q)$ determined by the reaction rates and stoichiometry. (The absorption
term $Y(q)$ due to immigration is not currently implemented.)

  * Evolve the conjugate coordinates $q_j$ along the flow $X_j$ in Fourier space
    implied by the reaction kinetics using Heun's method, a numerical
    time-stepping procedure.


This latter algorithm was inspired by prior work by
Suresh Moolgavkar in the 1980s and Dennis Quinn 1989, and performs a direct
numerical integration of a transformed version of the chemical master
equation/Kolmogorov forward equation. The previous work on characteristics used
a two-pass approach for non-constant coefficients, and used Euler integration.
This meant the asymptotic complexity was $O(\Delta t^{-2}) = O(\epsilon^{-2})$
for error tolerance $\epsilon$. 

This new algorithm is optimised for constant coefficients and uses Heun's
method, and it has an asymptotic complexity of $O(\Delta t^{-1}) = O(\epsilon^{-1/2})$ .
It shares a similar approach to Quinn's algorithm, by integrating characteristic
curves in the conjugate coordinates using a time-stepping procedure, but does so
in only one pass. 
As it is based on Kolmogorov forward equations rather than Kolmogorov backward
equations, and is noticeably faster than Gillespie, I have named it the ''fast
forward'' method. Other related ''fast forward'' methods may be possible, e.g.
with multiple passes, or higher order Runge-Kutta steps. A ''fast forward''
method is one that 

1. solves the Komogorov forward equations by numerically integrating the characteristics, and 
2. is strictly faster than $O(\epsilon^{-2})$ in the global error $\epsilon$ (as
$\epsilon$ decreases).

Georg Luebeck and Suresh Moolgavkar previously developed a numerical integration
approach based on Gaussian quadrature of Kolmogorov backward equations. I do not
know if this can be extended to arbitrary directed graphs -- it may be possible,
see Sanyi Tang et al. 2023. Dennis Quinn's
algorithm is much more closely related to the experimental algorithm here, but
we have found some radical improvements.

What the test cases should be of:
---------------------------------

Simulate models on graphs with both methods (Gillespie and fast forward). Need a convenient way to implement Kronecker products
(TODO).

Outputs:

  * Survival probabilities for Gillespie algorithm. Generate these using Kaplan-Meier
    plots
  * Corresponding survival curves from the generating function

These are both in CSV format, with a layout mirroring Ruibo Zhang et al. 2022.

Prerequisites:
--------------

Whole project is in C++ (with some Shell wrappers).

Compiles under both G++ and Clang. 

Requires [GSL](https://www.gnu.org/software/gsl/), which the Mac subsection of
the Makefile assumes depends on Homebrew.
