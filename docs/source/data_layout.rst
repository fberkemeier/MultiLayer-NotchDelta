Data Layout
===========

Expected structure
------------------

Provide data under ``data/`` using the following layout:

- ``data/adjacency_matrices/adjacency_matrices_<wing_region>.xlsx``
- ``data/cell_geometry/cell_geometry_<wing_region>.csv``
- ``data/wing_region_metadata.csv``

Adjacency matrices
------------------

Each Excel file stores one region. Each sheet corresponds to one depth layer.

Cell geometry
-------------

Each CSV should include at least:

- ``frame``
- ``centroid_x``
- ``centroid_y``
- ``area``

Metadata columns
----------------

``wing_region_metadata.csv`` should include:

- ``wing_region``: region key used in filenames/code.
- ``wing_label``: display label for plots.
- ``gap``: layer thickness (Delta L).
- ``n_layers``: number of depth layers.
- ``default_height``: default band height for spacing calculations.
- ``is_wing_disc``: boolean-like selector for canonical wing discs.
- ``signalling_labels``: list of signalling-capable cell indices (for example ``"[1, 5, 9, 12]"``).

Naming consistency
------------------

For each metadata ``wing_region``, matching adjacency and geometry files must exist.
