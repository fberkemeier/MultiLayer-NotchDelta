Notebook Workflow
=================

Notebook
--------

Main workflow notebook: ``notebooks/msm_notebook.ipynb``.

1. Import and setup
-------------------

The notebook imports from ``src/msm_model.py`` and manually defines ``wing_regions`` (for example a subset or a single region).
Then it builds all dictionaries required for simulations (adjacency, centroids, labels, metadata-derived parameters).

2. Notch profile setup
----------------------

``notch_data`` is defined in the notebook, so users can replace intensity profiles per dataset.

3. Single simulation
--------------------

The first simulation block calls ``compute_band_distance(...)`` with one region and a chosen parameter set.
This is the fastest section for sanity checks.

4. Spacing analysis
-------------------

Subsequent sections run sweeps over depth (``Lmax``) and/or straightening (``alpha``), storing spacing metrics and plotting with ``fancy_plot(...)``.

5. Neighbour and geometry analyses
----------------------------------

The notebook evaluates neighbour changes under straightening using ``straight_adjacency(...)`` and ``plot_straightening_nonapical(...)``.

6. Plot saving
--------------

Figures are saved to ``figures/`` when the corresponding ``saveQ=True`` flags are enabled.

Practical run order
-------------------

- Run setup cells.
- Run one single-simulation block.
- Run reduced sweeps (small ``sim_number``).
- Increase sweep size only after confirming outputs.
