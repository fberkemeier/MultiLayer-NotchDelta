# A Multi-layer Signalling Model of Notch-Delta

Toolkit for modelling and analysing 3D Notch–Delta signalling in multi-layered epithelial tissues, aimed at exploring how cell geometry and depth affect lateral inhibition.

## Mathematical model

The dynamics of Notch and Delta concentrations within each cell may be represented by the following equations (Collier et al., [1996](https://www.sciencedirect.com/science/article/pii/S0022519396902337), Binshtok & Sprinzak, [2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6879322/))

$$
\begin{align}
\frac{d}{dt}n_i = f(\langle d_i \rangle) - n_i - k_c n_id_i\\
\frac{d}{dt}d_i = \nu(g(n_i) - d_i) - k_c n_id_i
\end{align}
$$

where $f,g:[0,\infty)\to [0,\infty)$ are continuous increasing and decreasing functions, respectively, often taken to be Hill functions

$$
 f(x)=\frac{x^k}{a+x^k}, \quad
 g(x)=\frac{1}{1+bx^h}
$$

for $x\geq 0$ and $h,k\geq 1$. $r_t\equiv 1/a$ and $b$ are the trans-interactions strength and ligand inhibition strength parameters, respectively. $\nu>0$ is the ratio between Notch and Delta decay rates, determining the strength of decay. Finally, $\langle d_i\rangle$ is the average level of Delta activity in the cells adjacent to cell $i$. The model accounts for cis-inhibition effects, which reduce the signalling potential of Delta in a cell by the co-present Notch, modelling more complex interactions within the same cell. Parameter $k_c$ represents the strength of cis-inhibition.
