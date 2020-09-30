![unittest](https://github.com/alercebroker/lc_correction/workflows/unittest/badge.svg?branch=main&event=push)
[![codecov](https://codecov.io/gh/alercebroker/lc_correction/branch/master/graph/badge.svg?token=5C8D7F627W)](https://codecov.io/gh/alercebroker/lc_correction)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/alercebroker/lc_correction/blob/master/LICENSE)

# Light curve correction library

Alert magnitudes are produced after measuring the flux in an image difference, which is produced from subtracting a given observation from a reference image. This means that if the object was present in the reference image, the object’s true magnitude can be corrected. The formulas for the correction and associated error are the following:

![corr](https://alerce-science.s3.amazonaws.com/images/correction.max-1600x900.png)

![corr_error](https://alerce-science.s3.amazonaws.com/images/correction_error.max-1600x900.png)

Where *mcorr* and *δmcorr* are the corrected magnitude and error, mref and δmref are the reference magnitude and error, mdiff and δmdiff are the difference magnitude and error, and sgn is the sign of the alert (ifdiffpos). Note that these formulas can diverge if the reference and difference magnitude are the same and sgn is -1, but this should never happen as no alerts should be triggered in that case.

It is important to note that only if the reference object’s flux is known these formulas can be applied, which is not always the case. Moreover, if the reference image changes, it is possible that the object changes from being possible to correct to not being possible to correct, and vice versa.

We approach this problem by always providing both the uncorrected and corrected photometries, and flagging data where we detect inconsistent corrections through time, e.g., if the object changes from not being possible to correct to being possible to correct. We also provide a flag which tells whether we believe the object is unresolved or not, for users to decide whether to use the corrected photometry or not (see discussion on the database).

## Installing *lc_correction*
For development:

```bash
pip install -e .
```

## Main features

- Do correction to a lightcurve.
- Get magnitude statistics of a lightcurve.
- Get main statistics of an object.

## Dependencies
- [NumPy](https://numpy.org/)
- [Pandas](https://pandas.pydata.org/)


## How to use

First you must import some libraries:

```python
import pandas as pd
from correction.compute import *
from correction.helpers import *
```

The *lc_correction* use a `pandas.DataFrame` for all calculations, you must have a dafaframe for detections and non detections.

```python
detections = pd.read_parquet(path_to_data)
non_detections = pd.read_parquet(path_to_data) 
``` 

First, we need a modified julian dates (in case of you have julian dates).

```
detections["mjd"] = detections.jd - 2400000.5
non_detections["mjd"] = non_detections.jd - 2400000.5
```


### Apply correction to detections
We use a function in for apply function of pandas. So we use `apply_correction_df` function for all unique `[objectId, fid]` pairs.

```
corrected = detections.groupby(["objectId", "fid"]).apply(apply_correction_df)
corrected.reset_index(inplace=True)
```

### Get magnitude statistics
When you have corrected detections, you can get the magnitude statistics

```
magstats = corrected.groupby(["objectId", "fid"]).apply(apply_mag_stats)
magstats.reset_index(inplace=True)
```

### Get dm/dt
If you want to get a dm/dy information, only use magstats and non detections.
```
dmdt = do_dmdt_df(magstats, non_detections)
```


