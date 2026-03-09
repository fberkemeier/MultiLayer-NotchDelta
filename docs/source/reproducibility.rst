Reproducibility
===============

Environment
-----------

.. code-block:: powershell

   python -m venv .venv
   .\.venv\Scripts\Activate.ps1
   pip install -r requirements.txt

Run notebook
------------

.. code-block:: powershell

   jupyter notebook notebooks/msm_notebook.ipynb

Minimal validation run
----------------------

- Select one region in ``wing_regions``.
- Run setup/import cells.
- Run one single-simulation section.
- Confirm spacing output and one figure.

Scaling up
----------

After validation, increase:

- ``sim_number``
- sweep resolution (``Lmax``/``alpha`` ranges)
- runtime parameters (for sensitivity blocks)

Outputs
-------

- Interactive plots in notebook.
- Saved files in ``figures/`` (if ``saveQ=True``).
- Optional intermediate files under ``data/sensitivity_analysis/``.
