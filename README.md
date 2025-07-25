# A Multi-layer Model of Notch-Delta Signalling

We introduce a toolkit for modelling and analysing 3D Notch–Delta signalling in multi-layered epithelial tissues, aimed at exploring how cell geometry and depth affect lateral inhibition.

## Mathematical model

The dynamics of Notch-Delta signalling within each cell $i$ may be represented by the following system (Collier et al., [1996](https://www.sciencedirect.com/science/article/pii/S0022519396902337), Binshtok & Sprinzak, [2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6879322/))

$$
\begin{align}
\frac{d}{dt}n_i = f(\langle d_i \rangle)-n_i\\
\frac{d}{dt}d_i = \nu(g(n_i) - d_i)
\end{align}
$$

for $1\leq i \leq N$, where $N$ is the total number of cells. We define $f$ and $g$ as Hill functions

$$
\begin{align}
f(x) = \frac{x^k}{K_a^k + x^k},\quad
g(x) = \frac{K_r^h}{K_r^h + x^h}
\end{align}
$$

where $k$ and $h$ (Hill coefficients) determine the steepness of the response (cooperativity), and $K_a$ and $K_r$ set the half-maximal activation and repression thresholds, respectively. $\nu>0$ is the ratio between Notch and Delta decay rates. To simulate three-dimensional interactions, we introduce

$$
\begin{align}
\langle d_i \rangle=\sum_{k=0}^{n-1} \omega_k \left(\sum_{j\in  \mathbf{nn}(i)}\frac{\ell_{ij,k}}{P_{j,k}}d_{j}\right)
\end{align}
$$

for a total number of signalling layers $n$ (layer range), where, at each layer $k$ ($0\leq k\leq n-1$), $\ell_{ij,k}$ is the length of the shared edge between cells $i$ and neighbouring cell $j$, and $P_{j,k}$ is the cross-sectional perimeter of cell $j$ at that layer. $\mathbf{nn}(i)$ is the set of nearest neighbours of cell $i$, and $\omega_k$ is the signalling weight of layer $k$. The total number of signalling layers can be defined by $n=L/\mathrm{\Delta}L$, where $L$ is the actual apical-to-basal length, determined experimentally, and $\mathrm{\Delta}L$ is the width of each layer.

## Usage
