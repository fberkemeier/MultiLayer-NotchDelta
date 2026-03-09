# A Multi-layer Model of Notch-Delta Signalling

[![Documentation Status](https://readthedocs.org/projects/multilayer-notchdelta/badge/?version=latest)](https://multilayer-notchdelta.readthedocs.io/en/latest/?badge=latest)

We present a computational toolkit for simulating and analysing Notch-Delta signalling in three-dimensional epithelial tissues using a Multi-layer Signalling Model (MSM), which accounts for depth-resolved cell-cell contacts across apical and lateral surfaces. This framework enables systematic exploration of how tissue geometry and signalling range influence lateral inhibition and pattern formation. For a detailed description of the MSM, see [Paci et al. (2025)](https://www.biorxiv.org/content/10.1101/2025.08.08.668674v1).

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
Multi-layer Signalling Model (MSM) overview. Segmented 3D cellular data across successive tissue layers (left) are used to construct a depth-resolved contact network. The MSM simulates lateral inhibition on this layered structure to predict SOP fate decisions (right), incorporating both apical and lateral cell-cell interactions. Adapted from Paci et al. (2025).
</em></p>

## Usage

The core implementation is in `src/msm_model.py`, and interactive analyses are in `notebooks/msm_notebook.ipynb`. For more details, see the [Documentation](https://multilayer-notchdelta.readthedocs.io/en/latest/).

### Quick setup

Install the required packages:

```powershell
pip install -r requirements.txt
```

Open the notebook:

```powershell
jupyter notebook notebooks/msm_notebook.ipynb
```

### Data layout

Provide data at repository root under `data/` using this structure:

- Adjacency matrices: `data/adjacency_matrices/adjacency_matrices_<wing_region>.xlsx`
: Layer-by-layer cell-cell adjacency matrices for one region (`<wing_region>`). Each sheet corresponds to one depth layer.
- Cell geometry: `data/cell_geometry/cell_geometry_<wing_region>.csv`
: Per-cell geometric information (at least frame, centroid coordinates, and area) used for distances, straightening, and plotting.
- Wing disc metadata: `data/wing_region_metadata.csv`
: Region-level configuration table used to build model dictionaries and parse signalling-competent cells.

The metadata file should include at least:

- `wing_region`
: Region key used in filenames and throughout the notebook/model (for example `wd_1`).
- `wing_label`
: Human-readable label used in plot annotations.
- `gap`
: Vertical layer thickness used for integrating signalling weights across depth.
- `n_layers`
: Number of depth layers used for the region.
- `default_height`
: Default analysis band height used in SOP spacing calculations.
- `is_wing_disc`
: Boolean flag to mark canonical wing disc datasets (used for grouping/selection).
- `signalling_labels` (list of cell indices, e.g. `"[1, 5, 9, 12]"`)
: Cell indices considered signalling-competent in that region.

### Region selection

In the notebook, `wing_regions` is set explicitly as a list (manual selection).
For example:

```python
# Example set from Paci et al. (2026)
wing_regions = ['wd_1', 'wd_2', 'wd_3', 'wd_1_mbs', 'wd_2_mbs', 'wd_3_mbs']

# Or a custom subset
wing_regions = ['wd_1', 'wd_2', 'wd_3']
```

Then build dictionaries from these selected regions:

```python
wing_discs = list_wing_discs(wing_regions)
signalling_labels_dict = load_signalling_labels_dict(wing_regions)
wd_dict = build_wd_label_dict(wing_regions)
gap_dict = build_gap_dict(wing_regions)
n_dict = build_n_layers_dict(wing_regions)
heights_dict = build_default_height_dict(wing_regions)

A_dict = build_adjacency_dict(wing_regions)
centroids_dict = build_centroids_dict(wing_regions)
area_apical_dict = build_area_apical_dict(wing_regions)
diam_apical_dict = build_diam_apical_dict(wing_regions)
```

### Notes

- `notch_data` is notebook-owned and can be edited per dataset/experiment.
- Figures are saved to `figures/` when `saveQ=True`.
- Intermediate sensitivity outputs are saved under `data/sensitivity_analysis/`.


## System requirements

This codebase was developed and tested on Python 3.12.3 under both Windows 10 and Windows 11. No installation procedure is required beyond installing standard Python 3 and the key dependencies. All scripts should remain compatible with standard Python 3 distributions on other operating systems.

## License

This project is openly distributed under the MIT License. This license allows unrestricted use, redistribution, and modification, provided that proper attribution to the original creators is maintained.

## Contact information

For further information, contributions, or queries, please contact:

- **Email**: [fp409@cam.ac.uk](mailto:fp409@cam.ac.uk)
- **GitHub**: [fberkemeier](https://github.com/fberkemeier)

We welcome discussions via GitHub to improve the model or address potential issues.

## References

Paci, G., Berkemeier, F., Baum, B., Page, K.M., & Mao, Y. 3D epithelial cell topology tunes signalling range to promote precise patterning. _bioRxiv_ (2026). [doi.org/10.1101/2025.08.08.668674](https://www.biorxiv.org/content/10.1101/2025.08.08.668674v1)

