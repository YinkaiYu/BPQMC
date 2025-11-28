# Bosonic Projector QMC

This note summarises the rank-1 projector algorithm used in the bosonic PQMC codebase.  The simulation propagates only the left/right trial vectors and their overlap, so all local updates are matrix–vector operations ($O(N^2)$); equal-time Green’s functions are reconstructed on demand at measurement times.

## Rank-1 trial state

All $N_b$ bosons occupy the same orbital:

$$|\Phi_T\rangle = \frac{1}{\sqrt{N_b!}} \left( \sum_i a_i^\dagger P_i \right)^{N_b} |0\rangle \;\;\equiv\;\; \frac{1}{\sqrt{N_b!}} \left( a^\dagger P \right)^{N_b} |0\rangle.$$

The form is closed under Gaussian operators:

$$e^{-a^\dagger T a}  \frac{1}{\sqrt{N_b!}} \left( a^\dagger P \right)^{N_b} |0\rangle= \frac{1}{\sqrt{N_b!}} \left( a^\dagger e^{-T} P \right)^{N_b} |0\rangle.$$

Consequently the imaginary-time propagator

$$U(\tau_2,\tau_1) = e^{-a^\dagger h(\tau_2) a} \cdots e^{-a^\dagger h(\tau_1)a}, \qquad
B(\tau_2,\tau_1) = e^{-h(\tau_2)} \cdots e^{-h(\tau_1)}$$

acts on $P$ as a single column vector.

## Configuration weight

Given a Hubbard-Stratonovich field configuration $\phi$,

$$P\left[\phi\right] = e^{-\tfrac{1}{2}\left|\vec{\phi}\right|^2}  \langle \Phi_T | U(2\theta,0) | \Phi_T \rangle
= e^{-\tfrac{1}{2}\left|\vec{\phi}\right|^2} \left[ P^\dagger B(2\theta,0) P \right]^{N_b}.$$

During the sweep we store only the right vector $P_R(\tau)=B(\tau,0)P$ and the left vector $P_L^\dagger(\tau)=P^\dagger B(2\theta,\tau)$.  Their norms are tracked separately for stability, so no $N\times N$ matrices are needed until a measurement.

## Local update ratio

For a local update $\phi_i(\tau) \to \phi_i'(\tau)$, let $\Delta = \mathrm{diag}(0,\dots,\Delta_i,\dots,0)$ with $\Delta_i = e^{-(h_i' - h_i)} - 1$.  The Metropolis ratio reads

$$\frac{P[\phi']}{P[\phi]}= e^{-\tfrac{1}{2}(\phi_i'^2 - \phi_i^2)}\left[\frac{P^\dagger B(2\theta,\tau)(1+\Delta)B(\tau,0)P}{P^\dagger B(2\theta,\tau)B(\tau,0)P}\right]^{N_b}$$

$$=\exp\left[ -\tfrac{1}{2}(\phi_i'^2 - \phi_i^2) + N_b \log \left( 1 + \frac{\Delta_i \big[ P_R(\tau) \big]_i  \big[ P_L^\dagger(\tau) \big]_i}{P^\dagger B(2\theta,0) P} \right) \right],$$

where the numerator only requires the two propagated vectors.  This keeps the propagation cost at $O(N^2)$.

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
