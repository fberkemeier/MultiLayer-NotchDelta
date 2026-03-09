Overview
========

Main components
---------------

The MSM is structured as follows:

- ``src/msm_model.py``: core implementation of the model equations, data-loading utilities, and helper routines used across simulations and analyses.
- ``notebooks/msm_notebook.ipynb``: primary user workflow to select datasets, assemble region dictionaries, run simulations, and evaluate outputs step by step.
- ``data/``: input datasets organized by region, including layered adjacency matrices, cell geometry tables, and region-level metadata.
- ``figures/``: output directory for generated plots, such as single-simulation summaries and SOP spacing analyses.

Mathematical model
------------------

The MSM extends lateral inhibition dynamics to a layered 3D epithelial topology.
For each cell :math:`i`, with :math:`1 \leq i \leq N`, the model follows:

.. math::

   \frac{d n_i}{d t} = f(\langle d_i \rangle) - n_i

.. math::

   \frac{d d_i}{d t} = \nu \left(g(n_i) - d_i\right)

where :math:`n_i` and :math:`d_i` are Notch and Delta levels, and :math:`\nu` is the Delta/Notch timescale ratio.
Activation and repression are represented with Hill functions:

.. math::

   f(x) = \frac{x^k}{K_a^k + x^k}, \qquad g(x) = \frac{K_r^h}{K_r^h + x^h}

Here, :math:`k` and :math:`h` are Hill coefficients, and :math:`K_a`, :math:`K_r` are activation/repression thresholds.

The effective Delta input integrates contacts over signalling layers:

.. math::

   \langle d_i \rangle = \sum_{k=0}^{n-1} \omega_k \left(\sum_{j \in \mathbf{nn}(i,k)} \frac{\ell_{ij,k}}{P_{j,k}} d_j\right)

where :math:`k=0` corresponds to the apical layer, :math:`\ell_{ij,k}` is shared interface length,
:math:`P_{j,k}` is perimeter normalization, and :math:`\omega_k` is the depth-dependent layer weight.
The number of active layers is set by :math:`n = L / \Delta L`, where :math:`L` is tissue depth and :math:`\Delta L` is layer thickness.

Layer weights are obtained by integrating a fitted depth profile :math:`\omega(z)` over each layer bin:

.. math::

   \omega_k = \int_{k\Delta L}^{(k+1)\Delta L} \omega(z)\,dz

Straightening model
-------------------

To quantify how apico-basal tortuosity alters non-apical contacts, we also include a straightening model parameterised by
:math:`\alpha \in [0,1]` that interpolates between measured 3D contacts and
apical topology.

For each layer :math:`k`, the adjacency is binarized:

.. math::

   b_{ij}^{(k)} =
   \begin{cases}
   1, & A_{ij}^{(k)} > 0 \\
   0, & \text{otherwise}
   \end{cases}

Cell displacement from the apical plane is defined as:

.. math::

   d_i^{(k)} = \left\| p_i^{(k)} - p_i^{(0)} \right\|

Two discrepancy sets are then computed relative to apical contacts:

.. math::

   E_{\mathrm{ex}}^{(k)} = \{(i,j)\mid b_{ij}^{(0)}=0,\; b_{ij}^{(k)}=1\}

.. math::

   E_{\mathrm{mi}}^{(k)} = \{(i,j)\mid b_{ij}^{(0)}=1,\; b_{ij}^{(k)}=0\}

with edge severity score:

.. math::

   \delta_{ij}^{(k)} = \max\left(d_i^{(k)}, d_j^{(k)}\right)

and quantile thresholds:

.. math::

   \tau_{\mathrm{ex}}^{(k)}(\alpha)=Q_{1-\alpha}\!\left(\{\delta_{ij}^{(k)}:(i,j)\in E_{\mathrm{ex}}^{(k)}\}\right)

.. math::

   \tau_{\mathrm{mi}}^{(k)}(\alpha)=Q_{\alpha}\!\left(\{\delta_{ij}^{(k)}:(i,j)\in E_{\mathrm{mi}}^{(k)}\}\right)

The straightened adjacency :math:`\widetilde{A}_{ij}^{(k)}(\alpha)` is:

.. math::

   \widetilde{A}_{ij}^{(k)}(\alpha)=
   \begin{cases}
   0, & b_{ij}^{(0)}=0,\; b_{ij}^{(k)}=1,\; \delta_{ij}^{(k)} > \tau_{\mathrm{ex}}^{(k)}(\alpha) \\
   A_{ij}^{(0)}, & b_{ij}^{(0)}=1,\; b_{ij}^{(k)}=0,\; \delta_{ij}^{(k)} < \tau_{\mathrm{mi}}^{(k)}(\alpha) \\
   A_{ij}^{(k)}, & \text{otherwise}
   \end{cases}

So :math:`\alpha=0` leaves measured layers unchanged, while :math:`\alpha=1` maximally restores apical-like
connectivity by pruning extraneous deep contacts and reintroducing apical-missing ones.
