# BdG solver for a quasiperiodic superconductor

A self-consistent Bogoliubov-de Gennes (BdG) mean-field solver for the one-dimensional attractive Hubbard chain in an Aubry-André quasiperiodic potential, written in Python with NumPy and SciPy.

## Model

$$
H = -t\sum_{\langle ij\rangle\sigma}\left(c^\dagger_{i\sigma}c_{j\sigma} + \text{h.c.}\right) + \sum_i (V_i - \mu)\, n_i - U\sum_i n_{i\uparrow} n_{i\downarrow}, \qquad V_i = \lambda\cos(2\pi\beta i + \phi)
$$

with $U > 0$ attractive. The interaction is decoupled in both the pairing and Hartree channels, with $\Delta_i = -U\langle c_{i\downarrow}c_{i\uparrow}\rangle$, and the problem is written in the two-component Nambu basis $\Psi = (c_{1\uparrow},\dots,c_{L\uparrow},c^\dagger_{1\downarrow},\dots,c^\dagger_{L\downarrow})^T$ as

$$
H_{\mathrm{BdG}} = \begin{pmatrix} h & D \\ D^\dagger & -h^* \end{pmatrix}, \qquad D = \mathrm{diag}(\Delta_i).
$$

The quasiperiodic potential uses Fibonacci rational approximants, $\beta = F_{n-1}/F_n$ on a ring of $L = F_n$ sites, so periodic boundary conditions are exact. Because translations map $\phi \to \phi + 2\pi k/L$, distinct realisations are sampled on $[0, 2\pi/L)$.

## What the code does

- Solves the self-consistency in $(\Delta_i, n_i)$ at fixed $\mu$ with Anderson mixing, and bisects on $\mu$ to hit a target filling
- Removes the global U(1) phase drift by aligning each iterate to the previous one, leaving the local phase profile untouched
- Tries several seeds at each point and keeps the solution with the lowest grand potential $\Omega$, so the $\Delta = 0$ normal branch is never mistaken for the answer
- Computes the superfluid stiffness $D_s = L^{-1}\,\partial^2\Omega/\partial A^2$ by a Peierls twist, re-converging the order parameter at every twist, with Richardson extrapolation in the step size
- Computes the diamagnetic weight, the staggered CDW order $W$, pair correlations and the coherence length, and the zero-temperature quantum Fisher information of the current operator
- Computes the generalised IPR of every BdG eigenstate and the multifractal dimension $D_2$ from $P_2 \sim L^{-D_2}$ across Fibonacci sizes
- Includes a sweep driver that writes HDF5 files with the full parameter record and stored fields, sweeping $\lambda$ in both directions

## Validation

Running the notebook runs the full validation suite first. All 11 tests pass.

| Test | Checks | Result |
|---|---|---|
| $U = 0$ | $\Delta$ vanishes | $\max\lvert\Delta\rvert = 0$ |
| Particle-hole symmetry | spectrum symmetric about zero | $5.3\times10^{-15}$ |
| Clean BCS, $T = 0$ | agreement with an independent k-space gap equation | $\delta\Delta = 3.0\times10^{-12}$ |
| Clean BCS, $T > 0$ | same, at finite temperature | $\delta\Delta = 1.3\times10^{-13}$ |
| Clean BCS stiffness | $D_s = D_{\mathrm{dia}}$ at $T = 0$ | relative difference $6\times10^{-10}$ |
| Filling | bisection on $\mu$ hits the target density | $3.2\times10^{-10}$ |
| QFI | Wick formula, vectorised vs loop | agree to all printed digits |
| Aubry-André transition | $D_2$ at $\lambda = t, 2t, 3t$ | 0.967, 0.543, 0.001 |
| Pseudospin SU(2) | $\lvert\Delta\rvert^2 + W^2$ invariant at half filling | spread $1.8\times10^{-12}$ |
| $\phi$ sampling | $[0, 2\pi/L)$ grid gives distinct realisations | 12 of 12 distinct |
| Sweep direction | forward and backward sweeps in $\lambda$ agree | $1.7\times10^{-10}$ |

## Running it

Open `BdG_Project.ipynb` in Colab or Jupyter and run the cell. It needs `numpy` and `scipy`. `h5py` is only needed for `sweep()`.

A parameter sweep at fixed filling looks like

```python
sweep("scan.h5", n_approx=12, lam_grid=(0.0, 4.0, 21), U=2.0, T=0.0,
      filling=0.9, phis=32)
```

## Status

The solver and the validation suite are complete. Production sweeps of the phase diagram in $\lambda$ at fixed filling are the next step.

## Author

Adhrit Das. Released under the MIT License.
