[![Build Status](https://travis-ci.com/alercebroker/lc_correction.svg?token=ky96CzpxxqojpJq8cck6&branch=master)](https://travis-ci.com/alercebroker/lc_correction)
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

```
pip install -e .
```

