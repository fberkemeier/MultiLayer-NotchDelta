Python Source: msm_model.py
===========================

Core responsibilities
---------------------

``src/msm_model.py`` provides:

- data loaders/builders from metadata and region files,
- Notch-Delta simulation routines,
- spacing and geometry analysis utilities,
- plotting helpers used by the notebook.

Frequently used functions
-------------------------

Data and metadata helpers:

- ``list_wing_discs(...)``
- ``load_signalling_labels_dict(...)``
- ``build_*_dict(...)`` helpers for labels, gaps, layers, heights, adjacency, centroids, and neighbours.

Simulation and analysis:

- ``simulate(...)``
- ``compute_omega_k(...)``
- ``compute_band_distance(...)``

Plotting:

- ``plot_layer_graph(...)``
- ``fancy_plot(...)``
- ``plot_straightening_nonapical(...)``

Signalling profiles
-------------------

Available profile functions include:

- ``omega_exp``
- ``omega_cnt``
- ``omega_lin``
- ``omega_exp0``

These are exposed via ``omega_map``.
