Function Reference
==================

This page documents the main notebook-facing functions in ``src/msm_model.py``.
It focuses on the APIs used throughout ``msm_notebook.ipynb``.

Region and metadata utilities
-----------------------------

.. code-block:: python

   list_wing_regions()
   list_wing_discs(wing_regions=None)
   get_region_metadata(wing_region)
   build_gap_dict(wing_regions=None)
   build_n_layers_dict(wing_regions=None)
   build_default_height_dict(wing_regions=None)
   build_wd_label_dict(wing_regions=None)

Argument details:

- ``wing_regions``: optional list of region IDs. If ``None``, functions infer available regions from ``data/``.
- ``wing_region``: one region ID (for example ``'wd_1'``) that must exist in metadata/files.

What these functions provide:

- discover available regions,
- filter canonical wing-disc datasets,
- map each region to metadata-derived values (gap, layer count, default height, plot label).

Data-loading functions
----------------------

.. code-block:: python

   get_adjacency_layers(wing_region)
   build_adjacency_dict(wing_regions=None)
   get_centroids_layers(wing_region)
   build_centroids_dict(wing_regions=None)
   get_area_apical(wing_region)
   build_area_apical_dict(wing_regions=None)
   get_diam_apical(wing_region)
   build_diam_apical_dict(wing_regions=None)
   load_signalling_labels_dict(wing_regions=None, labels_dir=SIGNALLING_LABELS_DIR)

Argument details:

- ``wing_region``: region key used to resolve input filenames.
- ``wing_regions``: optional list of regions to batch-load.
- ``labels_dir``: optional legacy path to JSON signalling-label files (metadata column is preferred).

Returns:

- adjacency layers: list of 2D adjacency matrices (one per layer),
- centroids: list of ``(N_layer, 2)`` arrays,
- apical area/diameter vectors,
- signalling-label dictionaries keyed by region.

Neighbour-set builders
----------------------

.. code-block:: python

   build_signalling_labels_apical_dict(wing_regions=None, signalling_labels_dict=None, adjacency_dict=None)
   build_apical_neighbours_dict(wing_regions=None, signalling_labels_dict=None, adjacency_dict=None)
   build_nonapical_neighbours_dict(wing_regions=None, signalling_labels_dict=None, adjacency_dict=None)

Argument details:

- ``wing_regions``: regions to process.
- ``signalling_labels_dict``: optional preloaded dictionary (avoids reloading labels).
- ``adjacency_dict``: optional preloaded adjacency dictionary (avoids reloading adjacency files).

Returns dictionaries keyed by region that define:

- apical signalling-capable labels,
- apical neighbours,
- non-apical neighbours.

Simulation core
---------------

.. code-block:: python

   simulate(M, t_final, dt, signalling_labels=None, k=2, h=2, Ka=0.1, Kr=0.001, nu=1, randomQ=True)

Arguments:

- ``M``: effective weighted interaction matrix.
- ``t_final``: simulation end time.
- ``dt``: integration time step.
- ``signalling_labels``: optional index set; when provided, non-signalling cells are masked.
- ``k``: Hill coefficient for activation.
- ``h``: Hill coefficient for repression.
- ``Ka``: activation half-saturation threshold.
- ``Kr``: repression half-saturation threshold.
- ``nu``: Delta timescale factor.
- ``randomQ``: if ``False``, uses a fixed seed for reproducibility.

Returns ``(notch, delta)`` arrays.

Layer-weight integration
------------------------

.. code-block:: python

   compute_omega_k(omega_func, wing_region, Lmax=None, method='simpson', num_simpson=101, gap=None, n_layers=None)

Arguments:

- ``omega_func``: callable weight function over depth ``z``.
- ``wing_region``: region used for default ``gap`` and ``n_layers``.
- ``Lmax``: optional depth cutoff; weights beyond this depth are set to zero.
- ``method``: integration scheme (``'midpoint'``, ``'trapezoid'``, ``'simpson'``).
- ``num_simpson``: odd number of points per layer interval for Simpson integration.
- ``gap``: optional layer spacing override.
- ``n_layers``: optional layer count override.

Returns one integrated weight per layer.

Main notebook driver
--------------------

.. code-block:: python

   compute_band_distance(
       wing_region,
       weight_list=None,
       omega_func=None,
       t_final=1000.0,
       dt=0.1,
       Lmax=None,
       height=80,
       y_shift_steps=20,
       sim_number=1,
       quad_method='simpson',
       alpha=0.0,
       degen_T=1.0,
       plotQ=False,
       graphsaveQ=False,
       normalQ=False,
       str_type='centroid',
       epsmodelQ=False,
       eps=0.0,
       prot_len=0.0,
       marker_type='delta',
       k=2,
       h=2,
       Ka=0.1,
       Kr=0.001,
       nu=1,
       show_labels=True,
       show_other_layers=True,
       randomQ=True,
       signalling_labels=None,
       adjacency_layers=None,
       centroid_layers=None,
       **quad_kwargs,
   )

Argument groups:

- Region/input selection:
  ``wing_region``, ``signalling_labels``, ``adjacency_layers``, ``centroid_layers``.
- Weight construction:
  ``weight_list`` (direct), or ``omega_func`` + ``Lmax`` + ``quad_method`` + ``quad_kwargs``.
- Simulation control:
  ``t_final``, ``dt``, ``sim_number``, ``k``, ``h``, ``Ka``, ``Kr``, ``nu``, ``randomQ``.
- Spacing statistics:
  ``height``, ``y_shift_steps``.
- Structure/geometry options:
  ``alpha``, ``str_type``, ``epsmodelQ``, ``eps``, ``prot_len``.
- Plotting/output options:
  ``plotQ``, ``graphsaveQ``, ``marker_type``, ``show_labels``, ``show_other_layers``.
- Degeneracy handling:
  ``degen_T`` and ``normalQ``.

Returns:

- ``avg_results``: mean spacing by SOP threshold,
- ``avg_results_v``: associated spread by threshold,
- ``degenQ``: degeneracy flag,
- ``avg_bimodal``: bimodality score summary,
- ``notch``, ``delta``: final simulated states.

Straightening analysis
----------------------

.. code-block:: python

   straight_adjacency(A_dict, centroids_dict, alpha)
   plot_straightening_nonapical(A_str_dict, region, alphas, saveQ=False)

Arguments:

- ``A_dict``: dictionary of adjacency-layer lists by region.
- ``centroids_dict``: centroid-layer dictionary by region.
- ``alpha``: straightening percentage/control between reference and non-apical layers.
- ``A_str_dict``: precomputed straightened adjacency outputs.
- ``region``: region key to plot.
- ``alphas``: sampled straightening values.
- ``saveQ``: save plot to ``figures/``.

Plot helpers and misc
---------------------

.. code-block:: python

   plot_layer_graph(delta, centroids, A_layers, ...)
   fancy_plot(spacing_dict, Lmax_list, weight_type, wing_regions, ...)
   SOP_bimodality(delta, notch)
   height_set(diam_list, hs=0.04)

Use these to:

- draw single-layer network graphs,
- plot spacing trends over depth,
- compute a Delta-Notch contrast metric,
- derive default analysis band heights from apical diameters.
