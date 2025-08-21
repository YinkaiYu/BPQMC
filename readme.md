# Bosonic Projector QMC

试探波函数（ $N_b$ 个粒子处在同一个态）

$$|\Phi_T\rangle = \frac{1}{\sqrt{N_b!}} \left( \sum_i a_i^\dagger P_i \right)^{N_b} |0\rangle \;\;\equiv\;\; \frac{1}{\sqrt{N_b!}} \left( a^\dagger P \right)^{N_b} |0\rangle$$

该形式对高斯型算符封闭

$$e^{-a^\dagger T a} \, \frac{1}{\sqrt{N_b!}} \left( a^\dagger P \right)^{N_b} |0\rangle= \frac{1}{\sqrt{N_b!}} \left( a^\dagger e^{-T} P \right)^{N_b} |0\rangle$$

虚时传播子

$$U(\tau_2,\tau_1) = e^{-a^\dagger h(\tau_2) a} \cdots e^{-a^\dagger h(\tau_1)a}$$

$$B(\tau_2,\tau_1) = e^{-h(\tau_2)} \cdots e^{-h(\tau_1)}$$

构型权重

$$P\left[\phi\right] = e^{-\tfrac{1}{2}\left|\vec{\phi}\right|^2} \, \langle \Phi_T | U(2\theta,0) | \Phi_T \rangle$$

$$= e^{-\tfrac{1}{2}\left|\vec{\phi}\right|^2} \left[ P^\dagger B(2\theta,0) P \right]^{N_b}$$

等时格林函数

$$G_{ij}(\tau) = \delta_{ij} + N_b \, \frac{\big[ B(\tau,0)P \big]_i \, \big[ P^\dagger B(2\theta,\tau) \big]_j}{P^\dagger B(2\theta,\tau) B(\tau,0) P}$$

构型更新权重比值 (local update:  $\phi_i(\tau) \to \phi_i'(\tau)$ )

$$\frac{P[\phi']}{P[\phi]}= e^{-\tfrac{1}{2}(\phi_i'^2 - \phi_i^2)}\left[\frac{P^\dagger B(2\theta,\tau)(1+\Delta)B(\tau,0)P}{P^\dagger B(2\theta,\tau)B(\tau,0)P}\right]^{N_b}$$

$$= e^{-\tfrac{1}{2}(\phi_i'^2 - \phi_i^2)}\left[1 + \frac{\Delta_i}{N_b}(G_{ii}-1)\right]^{N_b}$$

可见 $\bar{G} \equiv G - I$ 储存了态中所有有用信息，快速更新：

$$\bar{G}' = \frac{1+\Delta_i}{r_b}\bar{G}$$

其虚时演化为

$$\bar{G}(\tau+1) = B(\tau+1,\tau) \, \bar{G} \, B(\tau+1,\tau)^{-1}$$

## 附：数值稳定性

归一化

$$P_R = \frac{B(\tau,0)P}{|B(\tau,0)P|} = \frac{B(\tau,0)P}{Z_R}$$

其中 $B(\tau,0)P$ 是 $N$ 行 1 列向量， $Z_R^2 = \sum_i \big[ B(\tau,0)P \big]_i \, \big[ B(\tau,0)P \big]_i^*$

同理

$$P_L^\dagger = \frac{P^\dagger B(2\theta,\tau)}{|P^\dagger B(2\theta,\tau)|}= \frac{P^\dagger B(2\theta,\tau)}{Z_L}$$


于是

$$\bar{G} = N_b \, \frac{B(\tau,0)P \, P^\dagger B(2\theta,\tau)}{P^\dagger B(2\theta,\tau) B(\tau,0)P}$$

$$= N_b \, \frac{Z_R P_R \, Z_L P_L^\dagger}{Z_L P_L^\dagger \, Z_R P_R}$$

$$= N_b \, \frac{P_R P_L^\dagger}{P_L^\dagger P_R}$$

每隔一段时间用上述式子正规化一下 $\bar{G}$ 即可。
