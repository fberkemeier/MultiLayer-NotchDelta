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

All functionality is provided in a single Jupyter Notebook, designed for use with the standard Python scientific stack. The code can be run in any environment supporting Jupyter (e.g. [Anaconda](https://www.anaconda.com/)). All relevant data is stored in the `data/` folder and includes adjacency matrices, edge data, centroid positions, and SOP labels for three wing discs.

### Structure and functionality

The notebook is divided into the following sections:

- **Library imports** – Imports all necessary Python libraries, including numerical computation (`numpy`, `scipy`), plotting (`matplotlib`), file I/O (`pandas`, `openpyxl`), and optimisation tools.
- **Main functions** – Contains the key computational routines for simulating Notch–Delta dynamics across a 3D tissue structure, including neighbour detection, signalling calculations, and SOP identification.
- **Parameter setup and data loading** – Defines model parameters (e.g. Notch–Delta interaction constants, simulation time, layer configuration) and loads input data such as adjacency matrices, centroids, and cell-specific signalling labels from the `data/` folder for three wing discs.
- **Data fitting** – Processes and fits the experimental Notch intensity profile along the apico-basal axis to guide the construction of depth-dependent signalling weight profiles.
- **Custom plotting and visualisation** – Provides tools to visualise simulation outcomes, such as SOP spacing distributions, network graphs, and comparative plots across experimental conditions or parameter sets.

Each simulation produces a 3D lateral inhibition pattern and allows for comparison between conditions (e.g. different depths or straightening factors). The number of simulations (`sim_number`) used in SOP spacing plots is set to 20 by default, but increasing this value can improve accuracy at the cost of longer computation times.

### Model parameters, wing disc datasets and metadata

The Notch–Delta signalling system is defined by:

- `k`, `h`: Hill coefficients determining cooperativity
- `Ka`, `Kr`: thresholds for Notch activation and Delta repression
- `ν`: relative decay rate of Delta vs Notch
- `t_final`: total simulation time
- `dt`: time step size

The three discs are identified by:

- `wing_regions`: dataset names – `wd_1`, `wd_2`, `wd_3`
- `wd_dict`: maps dataset name to a numeric label for plotting
- `gap_dict`: disc-specific vertical layer heights
- `n_dict`: number of depth layers per dataset

Signalling labels and intensity are stored in:

- `signalling_labels_dict`: maps each dataset to the set of cells competent for Notch–Delta signalling
- `notch_data`: experimental Notch intensity profile across the apico-basal axis (used to define the depth-dependent weights)

The following data structures are created from the Excel files in the `data/` folder:

- `path_dict`: maps dataset name to the corresponding `.xlsx` file path
- `A_dict`: stores normalised adjacency matrices for each wing disc, with each entry returning a list of matrices corresponding to successive tissue layers
- `centroids_dict`: stores 2D centroid coordinates for each cell, returned as a list (one per layer) for each wing disc

Neighbours are categorized according to:

- `signalling_labels_apical_dict`: for each dataset, lists signalling cells that are connected at the apical layer (layer 0)
- `apical_neighbours_dict`: maps each signalling cell to its apical neighbours
- `nonapical_neighbours_dict`: maps each signalling cell to its lateral (non-apical) neighbours

These structures allow the model to distinguish between interactions occurring at the surface and those mediated by deeper 3D contacts.

### Examples

The notebook includes built-in visualisation tools to aid interpretation of simulation results. Below are a few examples:

#### Graph representation

Cells are represented as nodes in a multi-layered graph, with edges indicating contact-based signalling potential. The apical and sub-apical contacts are extracted from segmented imaging data and rendered using `matplotlib`. Edges may be colour-coded to distinguish contact types (e.g. apical vs. non-apical).

Example output:
- Red nodes: SOP cells (selected based on Delta threshold).
- Black edges: apical contacts.
- Light blue edges: exclusively non-apical contacts.

```
wing_region = 'wd_1'
def omega_exp(z):
    return 40 * np.exp(-0.4 * z)
Lmax = 0.5
height_list = [85., 170., 105.]
for Lmax in [Lmax]:
    result = compute_band_distance(wing_region, omega_func=omega_exp, Lmax=Lmax,
        sim_number=1, quad_method='simpson', heights=height_list, plotQ=True,
        normalQ=False, alpha=0, degen_T=1., y_shift_steps=20
    )
```

<img width="537" height="530" alt="image" src="https://github.com/user-attachments/assets/0625fd71-70bb-4003-b5fb-5b1e8216fa80" />


#### SOP spacing

To quantify the effect of 3D connectivity on pattern formation, SOP spacing is computed as the shortest-path distance between SOP cells in the apical contact graph, restricted to a defined region of interest. This metric is plotted across different simulations and experimental conditions.

Example output:
- SOP spacing vs. signalling depth.
- SOP spacing vs. straightening percentage (geometry manipulation).
- Each curve corresponds to a specific wing disc.
- Degenerate patterns (e.g. excessive crowding, irregular spacing) are flagged in red.

```
# exponential signalling
def omega_exp(z):
    return 10 * np.exp(-0.2 * z) + 2

wing_regions = ['wd_1', 'wd_2', 'wd_3']
sim_number = 20
threshold = 0.8
Lmax_list = np.linspace(0.5, 25, 15)
height_list = [85., 170., 105.]
degen_T = 1.
normalQ = False

spacing_dict_exp = {region: [] for region in wing_regions}
it=0
for Lmax in Lmax_list:
    for region in wing_regions:
        d, vr, degenQ = compute_band_distance(
            region,
            omega_func=omega_exp,
            Lmax=Lmax,
            sim_number=sim_number,
            quad_method='simpson',
            heights=height_list,
            degen_T=degen_T, normalQ=normalQ,
            y_shift_steps=20
        )
        d = d[threshold]
        vr = vr[threshold]
        spacing_dict_exp[region].append([d,vr,degenQ])
        it=it+1
        print(f'{it}/{len(Lmax_list)*3}',end='\r')

fancy_plot(spacing_dict_exp, Lmax_list, 'exp', wing_regions, degenplotQ=True, ylim=(1.1,2.65), errorbarQ=False, saveQ=False)
```

<img width="527" height="390" alt="image" src="https://github.com/user-attachments/assets/8b6e0c46-51f3-4338-8a47-7dab86d1bf0e" />

These plots provide intuitive insight into how geometry and signalling range influence the emergent spacing of SOP cells.

