Installation
============

Prerequisites
-------------

- Python 3.8+
- A scientific Python environment (recommended: conda or mamba)

Install from GitHub
-------------------

.. code-block:: bash

   python -m pip install "scanpy>=1.9" plotly
   python -m pip install git+https://github.com/haniffalab/scanpy_plus.git

Install for local development
-----------------------------

From the repository root:

.. code-block:: bash

   python -m pip install -e .

Build docs locally
------------------

.. code-block:: bash

   python -m pip install -r docs/requirements.txt
   sphinx-build -b html docs/source docs/build/html

Open ``docs/build/html/index.html`` in a browser to preview.
