# Bosonic Projector QMC

This note summarises the rank-1 projector algorithm used in the bosonic PQMC codebase.  The simulation propagates only the left/right trial vectors and their overlap, so all local updates are matrix-vector operations ($O(N^2)$); equal-time Green's functions are reconstructed on demand at measurement times.

The current code stores one bosonic orbital per propagated vector, but the statistical weight corresponds to two physical flavors.  In the implementation this appears through the modulus square of the bosonic amplitude and through the factor `2 Re[...]` in the HMC force.

## Rank-1 trial state

All $N_b$ bosons occupy the same orbital:

$$|\Phi_T\rangle = \frac{1}{\sqrt{N_b!}} \left( \sum_i a_i^\dagger P_i \right)^{N_b} |0\rangle \;\;\equiv\;\; \frac{1}{\sqrt{N_b!}} \left( a^\dagger P \right)^{N_b} |0\rangle.$$

The form is closed under Gaussian operators:

$$e^{-a^\dagger T a}  \frac{1}{\sqrt{N_b!}} \left( a^\dagger P \right)^{N_b} |0\rangle= \frac{1}{\sqrt{N_b!}} \left( a^\dagger e^{-T} P \right)^{N_b} |0\rangle.$$

Consequently the imaginary-time propagator

$$U(\tau_2,\tau_1) = e^{-a^\dagger h(\tau_2) a} \cdots e^{-a^\dagger h(\tau_1)a}, \qquad
B(\tau_2,\tau_1) = e^{-h(\tau_2)} \cdots e^{-h(\tau_1)}$$

acts on $P$ as a single column vector.

## Configuration weight and code variables

Given a Hubbard-Stratonovich field configuration $\phi$,

$$P\left[\phi\right] = e^{-\tfrac{1}{2}\left|\vec{\phi}\right|^2}
\left(\left[ P^\dagger B(2\theta,0) P \right]^{N_b}\right)
\left(\left[ P^\dagger B(2\theta,0) P \right]^{N_b}\right)^*.$$

This is the two-flavor weight implemented in the local update through

$$|e^{-\tfrac{1}{2}\Delta \phi^2} \, r_b^{N_b} \, (r_b^{N_b})^*|,$$

which matches the code line

```fortran
ratio_abs = abs(ratio_exp * ratio_Pfa * dconjg(ratio_Pfa))
```

The main variable mapping in the code is

- `Nbos` $\leftrightarrow N_b$
- `Beta` $\leftrightarrow 2\theta$
- `Prop%UUR(:,1)` $\leftrightarrow P_R(\tau)=B(\tau,0)P$
- `Prop%UUL(1,:)` $\leftrightarrow P_L^\dagger(\tau)=P^\dagger B(2\theta,\tau)$
- `Prop%overlap` $\leftrightarrow P_L^\dagger(\tau) P_R(\tau)$ for the currently normalized vectors
- `OperatorHubbard%alpha` $\leftrightarrow \alpha=\sqrt{-2 \Delta\tau U}$ for $U<0$ and $\alpha=i\sqrt{2 \Delta\tau U}$ for $U>0$
- `Conf%phi_list(nf,ii,nt)` stores the auxiliary fields for both Hubbard channels `nf=1,2`

During the sweep we store only the right vector $P_R(\tau)=B(\tau,0)P$ and the left vector $P_L^\dagger(\tau)=P^\dagger B(2\theta,\tau)$.  Their norms are tracked separately for stability, so no $N\times N$ matrices are needed until a measurement.

## Local update ratio

For a local update $\phi_i(\tau) \to \phi_i'(\tau)$, let $\Delta = \mathrm{diag}(0,\dots,\Delta_i,\dots,0)$ with $\Delta_i = e^{-(h_i' - h_i)} - 1$.  The Metropolis ratio reads

$$\frac{P[\phi']}{P[\phi]}= e^{-\tfrac{1}{2}(\phi_i'^2 - \phi_i^2)}\left[\frac{P^\dagger B(2\theta,\tau)(1+\Delta)B(\tau,0)P}{P^\dagger B(2\theta,\tau)B(\tau,0)P}\right]^{N_b}$$

$$=\exp\left[ -\tfrac{1}{2}(\phi_i'^2 - \phi_i^2) + N_b \log \left( 1 + \frac{\Delta_i \big[ P_R(\tau) \big]_i  \big[ P_L^\dagger(\tau) \big]_i}{P^\dagger B(2\theta,0) P} \right) \right],$$

where the numerator only requires the two propagated vectors.  This keeps the propagation cost at $O(N^2)$.

## HMC effective action

For HMC sampling the auxiliary field and its conjugate momentum are evolved with

$$H[\phi,\pi] = K[\pi] + S[\phi], \qquad K[\pi] = \frac12 \sum_{\tau,i,n_f} \pi_{n_f,i}^2(\tau),$$

$$S[\phi] = \frac12 \sum_{\tau,i,n_f} \phi_{n_f,i}^2(\tau) - N_b \ln \left| P^\dagger B(2\theta,0) P \right|^2.$$

The HMC force for each auxiliary-field component is

$$F_{n_f,i}(\tau) = -\phi_{n_f,i}(\tau) + 2 \operatorname{Re}\!\left[\alpha_{n_f}\,\bar{G}_{ii}^{(n_f)}(\tau)\right],$$

with

$$\bar{G}_{ii}(\tau) = N_b \frac{\big[P_R(\tau)\big]_i \big[P_L^\dagger(\tau)\big]_i}{P_L^\dagger(\tau) P_R(\tau)}.$$

The factor $2 \operatorname{Re}[\cdots]$ is again the manifestation of the two physical flavors.  The plus sign follows directly from the local-update ratio implemented in `localU.f90`: the HMC force must match the linearized change of `log P[\phi]`.

## Equal-time Green’s function for measurements

The full Green’s function is reconstructed only when observables are evaluated (typically at $\tau=\theta$):

$$G_{ij} \equiv \langle a_i a_j^\dagger \rangle = \delta_{ij} + N_b  \frac{\big[ B(\theta,0)P \big]_i  \big[ P^\dagger B(2\theta,\theta) \big]_j}{P^\dagger B(2\theta,0) P},$$

$$\bar{G}_{ij} \equiv (G - I)_{ij} = N_b  \frac{\big[ B(\theta,0)P \big]_i  \big[ P^\dagger B(2\theta,\theta) \big]_j}{P^\dagger B(2\theta,0) P}.$$

With normalized vectors (see below), this simplifies to

$$\bar{G} = N_b  \frac{P_R P_L^\dagger}{P_L^\dagger P_R},$$

which is a rank-1 matrix that can be formed on demand for measurements.

## Numerical stabilization

To control the norms of $P_R$ and $P_L^\dagger$, periodically rescale them:

$$P_R = \frac{B(\tau,0)P}{|B(\tau,0)P|} = \frac{B(\tau,0)P}{Z_R}, \qquad Z_R^2 = \sum_i \big| \big[ B(\tau,0)P \big]_i \big|^2,$$

$$P_L^\dagger = \frac{P^\dagger B(2\theta,\tau)}{|P^\dagger B(2\theta,\tau)|}= \frac{P^\dagger B(2\theta,\tau)}{Z_L}, \qquad Z_L^2 = \sum_i \big| \big[ P^\dagger B(2\theta,\tau) \big]_i \big|^2.$$

For HMC the cumulative logarithms of the discarded normalization factors must also be tracked.  If the stored vectors are normalized by factors $Z_R$ and $Z_L$, then

$$\ln \left| P^\dagger B(2\theta,0) P \right|^2
= 2 \ln Z_R + 2 \ln Z_L + \ln \left| P_L^\dagger P_R \right|^2,$$

where the last overlap is computed from the normalized vectors kept in memory.

## Runtime input

The sampler is selected through the extra line in `test/paramC_sets.txt`:

```text
is_global   Nfrog   hmc_dt
```

- `is_global = .false.` keeps the local-update sampler.
- `is_global = .true.` enables the HMC sampler.
- `Nfrog` is the leapfrog step count.
- `hmc_dt` is the leapfrog step size.

The local-update path still uses `shiftLoc`.  The HMC path ignores `shiftLoc` during production sweeps and uses the full HMC warm-up/sampling trajectory instead.

## Benchmarks and Tuning

The repository now includes two helper scripts under `test/`:

- `python test/tune_hmc.py ...` scans `Nfrog` and `hmc_dt`, reports acceptance and autocorrelation estimates, and recommends candidates inside the target acceptance window when available.
- `python test/benchmark_hmc.py ...` runs local-update and HMC jobs on the same parameter sets and compares the main scalar observables plus several representative momentum-point observables.

Representative HMC candidates obtained from the offline scans in this repository are

- `weak_u2`: `Nfrog = 8`, `hmc_dt = 0.010`
- `mixed_u1_u2`: `Nfrog = 12`, `hmc_dt = 0.400`
- `strong_u2`: `Nfrog = 10`, `hmc_dt = 0.004`

For the mixed-coupling case, longer trajectories were needed before HMC and local-update runs entered the same high-double-occupancy region.  In practice the benchmark scripts should be run with sufficiently long thermal cuts before comparing means.
