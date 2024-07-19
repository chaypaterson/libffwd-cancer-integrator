cancer integrator
=================

Simulates survival curves in generalised multi-stage clonal expansion models of
Moolgavkar-Vernon-Knudsen type. This project allows us to
define a birth-death-mutation process on an arbitrary directed graph
(representing mutations), and generate survival curves using both Gillespie algorithm simulations and
fast direct integration of the Kolmogorov forward equations (a "fast forward" method).

Integration methods
-------------------

This compares different numerical simulation methods for birth-death-mutation
processes: the classic Gillespie algorithm, and an experimental algorithm based on the method of
characteristics. 

The two methods:

1. Use the Gillespie algorithm for many stochastic, "exact" samples of the
master equation. e.g. for the two-hit model:

$$\Gamma = \sum_{events} rate(event)$$

$$x = Unif(0,\Gamma)$$

    choose event from $x$ so that

$$Pr(event) = rate(event) / \Gamma$$

$$t += Exp(0, 1/\Gamma)$$

Repeat while $t < $ maximum time.

Suppose we want to estimate the probability of an outcome using random sampling.
The error $\epsilon$ in our estimate will scale with the square root of number of simulations
$N$ we run, so $\epsilon = O(N^{-1/2})$. Since the total run time is
proportional to N, even if we parallelise and optimise the algorithm further,
the time complexity of random sampling will be $T = O(\epsilon^{-2})$.

2. Use the method of characteristics to numerically integrate the generating function:
  * Generating function $\Psi =$ the Fourier transform of the probability
    distribution $P$, a function of the conjugate variables $q_j$
  * The $q_j$ are conjugate variables to the populations $N_j$
  * The Kolmogorov forward equation is of the form

$$\frac{\partial\Psi}{\partial t} = \vec{X} \cdot \nabla_{\vec{q}} \Psi + Y(\vec{q}) \Psi$$
with $X_j(q)$ determined by the reaction rates and stoichiometry. (The absorption
term $Y(q)$ due to immigration is not currently implemented.)

  * Evolve the characteristics $q_j$ along the flow $X_j$ in Fourier space
    implied by the reaction kinetics using Heun's method, a numerical
    time-stepping procedure.


This latter algorithm was inspired by prior work by
Suresh Moolgavkar in the 1980s and Dennis Quinn 1989, and performs a direct
numerical integration of a transformed version of the chemical master
equation/Kolmogorov forward equation. The previous work on characteristics used
a two-pass approach for non-constant coefficients, and used Euler integration.
The time complexity of these earlier methods was $T = O(\epsilon^{-2})$, no
better than random sampling.


This new algorithm is optimised for constant coefficients, and uses Heun's method/improved Euler instead of Euler
integration. 
It shares a similar approach to Quinn's algorithm, by integrating characteristic
curves using a time-stepping procedure, but does so
in only one pass, eliminating half of the operations, and uses a better
time-stepping scheme. These optimisations make the method much more efficient
than earlier approaches.  The error scales as $\epsilon = O(\Delta t^{2})$
for error tolerance $\epsilon$, and this new method therefore has a time complexity of $T = O(\Delta t^{-1}) = O(\epsilon^{-1/2})$.

As it is based on Kolmogorov forward equations rather than Kolmogorov backward
equations, and is noticeably faster than Gillespie, we have named it the ''fast
forward'' method. Other related ''fast forward'' methods may be possible, e.g.
with multiple passes, or higher order Runge-Kutta steppers. 

Another ''fast forward'' method is one that 

1. solves the Komogorov forward equations by numerically integrating the characteristics, and 
2. is strictly faster than $O(\epsilon^{-2})$ in the global error $\epsilon$ (as
$\epsilon$ decreases).

Georg Luebeck and Suresh Moolgavkar previously developed a numerical integration
approach based on Gaussian quadrature of Kolmogorov backward equations. I do not
know if this can be extended to arbitrary directed graphs -- it may be possible,
see Sanyi Tang et al. 2023. Dennis Quinn's
algorithm is much more closely related to the experimental algorithm here, but
we have found some radical improvements.

Programs this project builds:
-----------------------------

Various programs simulate models on graphs with both methods (Gillespie and fast forward).

  * Unit tests for the core library
  * Programs to generate survival probability curves in two ways:
    * Survival probabilities for Gillespie algorithm. Generate these using Kaplan-Meier plots
    * Corresponding survival curves from the generating function
  * Programs to estimate numerical errors in both the Gillespie algorithm and
    the new method (under src/errors)
  * A program that demonstrates the use of fast forward integration for
    statistical inference: simulating a clinical study with random sampling, and
    then learns the parameters used to generate this study with maximum
    likelihood estimation. All under src/inference.

Project structure:
------------------

 * `src`: contains source code
 * `src/core`: the core library
 * `src/errors`: programs used to measure numerical error in the new method
 * `src/inference`: the statistical inference harness
 * `src/tests`: unit tests
 * `include`: headers defining the API for the core library

Building:
---------

From a tarball:

    mkdir build
    cd build
    ../configure --prefix=[...your install location here...] --with-eigen="/path/to/eigen3"
    make
    make install

`configure` will throw errors if the system dependencies are not met.

Mac users will also need to pass flags for their installation of GSL. e.g.

	../configure [...] CPPFLAGS="-I/opt/homebrew/include" LDLIBS="-L/opt/homebrew/lib"

Without a tarball:

Alternatively, a hand-written Makefile is provided as Makefile.old. This is the
original but is not machine agnostic, and only targetted 4 machines: my
workstation, my development box at the WMIC, my MacBook, and the [CSF](https://ri.itservices.manchester.ac.uk/csf3/).


Build & run Python binding:
------------------------

*Ensure pybind11 is available in your environment. Install it via pip if necessary:
    
    pip install pybind11

Navigate to the cancer-integrator/src/core directory:

Create and activate a Python virtual environment:
    
    python3 -m venv venv
    
    source venv/bin/activate  

Build and install the Python bindings:
    
    python pyffwdsetup.py build
    
    python pyffwdsetup.py install

Run the Python test script:
    
    python pyffwdtest.py


Requirements:
-------------

Whole project is in C++ (with some Shell wrappers). Compiles under both G++ and Clang. 

Requires 
 * C++17 or newer (was 14 but newer versions of Eigen expect 17)
 * [GSL 2.7](https://www.gnu.org/software/gsl/)
 * [Eigen 3](https://eigen.tuxfamily.org/index.php?title=Main_Page) 
 * make

The Mac subsection of Makefile.old assumes GSL and Eigen are installed under Homebrew.

Sources: 
--------

 * R. Meza PNAS 2008; 
 * K.S. Crump, Risk Analysis, Vol. 25, No. 4, 2005
 * E.G. Luebeck et al 2012, doi: 10.1158/0008-5472.CAN-12-2198
 * D. Quinn, Risk Analysis, Vol. 9, Issue 3, 1989, doi: 10.1111/j.1539-6924.1989.tb01006.x


