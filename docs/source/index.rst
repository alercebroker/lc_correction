.. lc-correction documentation master file, created by
   sphinx-quickstart on Tue Sep 29 14:51:23 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ALeRCE lc_correction's documentation!
================================================

Alert magnitudes are produced after measuring the flux in an image difference, which is produced from subtracting a given observation from a reference image. This means that if the object was present in the reference image, the object’s true magnitude can be corrected. The formulas for the correction and associated error are the following:

.. image:: https://camo.githubusercontent.com/82a7cb4030c6f569c36ff2bc799541fbad1123a7/68747470733a2f2f616c657263652d736369656e63652e73332e616d617a6f6e6177732e636f6d2f696d616765732f636f7272656374696f6e2e6d61782d31363030783930302e706e67

.. image:: https://camo.githubusercontent.com/7f0d7f596f36c6222760258d3fefb7461cf4ec71/68747470733a2f2f616c657263652d736369656e63652e73332e616d617a6f6e6177732e636f6d2f696d616765732f636f7272656374696f6e5f6572726f722e6d61782d31363030783930302e706e67

Where mcorr and δmcorr are the corrected magnitude and error, mref and δmref are the reference magnitude and error, mdiff and δmdiff are the difference magnitude and error, and sgn is the sign of the alert (ifdiffpos). Note that these formulas can diverge if the reference and difference magnitude are the same and sgn is -1, but this should never happen as no alerts should be triggered in that case.

It is important to note that only if the reference object’s flux is known these formulas can be applied, which is not always the case. Moreover, if the reference image changes, it is possible that the object changes from being possible to correct to not being possible to correct, and vice versa.

We approach this problem by always providing both the uncorrected and corrected photometries, and flagging data where we detect inconsistent corrections through time, e.g., if the object changes from not being possible to correct to being possible to correct. We also provide a flag which tells whether we believe the object is unresolved or not, for users to decide whether to use the corrected photometry or not (see discussion on the database).


Installing lc_correction
========================

Clone the repository and install from there

.. code-block:: bash

    git clone https://github.com/alercebroker/lc_correction
    cd lc_correction
    python setup.py install


Main features
=============

- Do correction to a light curve.
- Get magnitude statistics of a light curve.
- Get main statistics of an object.


How to use
==========
.. code-block:: python

   import pandas as pd
   from correction.compute import *
   from correction.helpers import *


The lc_correction use a pandas.DataFrame for all calculations, you must have a dafaframe for detections and non detections.

.. code-block:: python

   detections = pd.read_parquet(path_to_data)
   non_detections = pd.read_parquet(path_to_data)


First, we need a modified julian dates (in case of you have julian dates). We use a function in for apply function of pandas. So we use apply_correction_df function for all unique [objectId, fid] pairs. When you have corrected detections, you can get the magnitude statistics. If you want to get a dm/dy information, only use magstats and non detections.

.. code-block:: python

   detections["mjd"] = detections.jd - 2400000.5
   non_detections["mjd"] = non_detections.jd - 2400000.5

   corrected = detections.groupby(["objectId", "fid"]).apply(apply_correction_df)
   corrected.reset_index(inplace=True)

   magstats = corrected.groupby(["objectId", "fid"]).apply(apply_mag_stats)
   magstats.reset_index(inplace=True)

   dmdt = do_dmdt_df(magstats, non_detections)


.. toctree::
   :maxdepth: 3
   :caption: Contents:

   modules

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
