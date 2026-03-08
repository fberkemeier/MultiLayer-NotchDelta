# A Multi-layer Model of Notch-Delta Signalling

We present a computational toolkit for simulating and analysing Notchâ€“Delta signalling in three-dimensional epithelial tissues using a Multi-layer Signalling Model (MSM), which accounts for depth-resolved cellâ€“cell contacts across apical and lateral surfaces. This framework enables systematic exploration of how tissue geometry and signalling range influence lateral inhibition and pattern formation. For a detailed description of the MSM, see Paci et al. (2025).

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

<p align="center">
  <img src="https://github.com/user-attachments/assets/01b68b32-499a-47fc-aaf1-4430b51d2bd7" width="500" alt="Multi-layer Signalling Model">
</p>

<p align="center"><em>
Multi-layer Signalling Model (MSM) overview. Segmented 3D cellular data across successive tissue layers (left) are used to construct a depth-resolved contact network. The MSM simulates lateral inhibition on this layered structure to predict SOP fate decisions (right), incorporating both apical and lateral cellâ€“cell interactions. Adapted from Paci et al. (2025).
</em></p>

## Usage

The core implementation is in `src/msm_model.py`, and the interactive workflows are provided in `notebooks/msm_notebook.ipynb`.

### Quick setup

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

### Run the notebook

```powershell
jupyter notebook notebooks/msm_notebook.ipynb
```

The first notebook cell is configured to import from `src/`:

```python
import sys
from pathlib import Path

repo_root = Path.cwd()
if not (repo_root / 'src').exists():
    repo_root = repo_root.parent

sys.path.insert(0, str(repo_root / 'src'))
from msm_model import *
```

### Run from a Python script/session

```powershell
$env:PYTHONPATH="src"
python -c "import msm_model; print('msm_model import: OK')"
```

### Data and outputs

- Input datasets are expected under `data/` at repository root.
- Figures are written to `figures/`.
- Sensitivity-analysis outputs are written under `data/sensitivity_analysis/`.

### Structure and functionality

The notebook is divided into the following sections:

- **Imports and model access** - Loads dependencies and imports MSM functions from `src/msm_model.py`.
- **Main simulations** - Runs Notch-Delta simulations for selected wing disc datasets.
- **Spacing and sensitivity analyses** - Computes SOP spacing and parameter heatmaps.
- **Visualisation tools** - Produces graph-based and summary plots, saved to `figures/` when enabled.

## System requirements

This codebase was developed and tested on Pythonâ€¯3.12.3 under both Windowsâ€¯10 and Windowsâ€¯11. No installation procedure is required beyond installing standard Pythonâ€¯3 and the key dependencies. All scripts should remain compatible with standard Pythonâ€¯3 distributions on other operating systems.

## License

This project is openly distributed under the MIT License. This license allows unrestricted use, redistribution, and modification, provided that proper attribution to the original creators is maintained.

## Contact information

For further information, contributions, or queries, please contact:

- **Email**: [fp409@cam.ac.uk](mailto:fp409@cam.ac.uk)
- **GitHub**: [fberkemeier](https://github.com/fberkemeier)

We welcome discussions via GitHub to improve the model or address potential issues.

## References

Paci, G., Berkemeier, F., Baum, B., Page, K.M., & Mao, Y. 3D epithelial cell topology tunes signalling range to promote precise patterning. _bioRxiv_ (2025).

