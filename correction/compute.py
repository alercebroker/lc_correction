import numpy as np
import pandas as pd
import logging
from joblib import Parallel, delayed

DISTANCE_THRESHOLD = 1.4
SCORE_THRESHOLD = 0.4
CHINR_THRESHOLD = 2
SHARPNR_MAX = 0.1
SHARPNR_MIN = -0.13
ZERO_MAG = 100.

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(name)s.%(funcName)s: %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',)


def validate_object(candidate, is_first_detection, stellar_object=False): #Algorithm 1 
    stellar_magstats = False
    if is_first_detection:
        if candidate["distpsnr1"] < DISTANCE_THRESHOLD and candidate["distpsnr1"] < DISTANCE_THRESHOLD and candidate["sgscore1"] > SCORE_THRESHOLD:
            stellar_object = True
            stellar_magstats = True
        else:
            if candidate["distnr"] < DISTANCE_THRESHOLD and candidate["chinr"] < CHINR_THRESHOLD and candidate["sharpnr"] > SHARPNR_MIN and candidate["sharpnr"] < SHARPNR_MAX:
                stellar_magstats = True
    else:
        if candidate["distnr"] < DISTANCE_THRESHOLD:
            if stellar_object:
                stellar_magstats = True
            else:
                if candidate["chinr"] < CHINR_THRESHOLD and candidate["sharpnr"] > SHARPNR_MIN and candidate["sharpnr"] < SHARPNR_MAX:
                    stellar_magstats = True
    return stellar_object, stellar_magstats


def validate_magnitudes(candidate, corr_detection=None, flag=None, corr_magstats=None, is_first_detection=False): #Algorithm 2
    if candidate["distnr"] < DISTANCE_THRESHOLD:
        corr_detection = True
        flag = False
        if is_first_detection:
            corr_magstats = True
        else: 
            if corr_magstats == False:
                flag = True
    else: 
        corr_detection = False
        if is_first_detection:
            corr_magstats = False
            isdiffpos = 1 if (candidate["isdiffpos"] in ["t", "1"]) else -1
            if isdiffpos > 0:
                flag = False
            else:
                flag = True
        else: 
            if corr_magstats == True:
                flag = True
            else:
                isdiffpos = 1 if (candidate["isdiffpos"] in ["t", "1"]) else -1
                if isdiffpos > 0:
                    flag = False
                else:
                    flag = True

    return corr_detection, corr_magstats, flag


def correction(magnr, magpsf, sigmagnr, sigmapsf, isdiffpos, oid): #Correction Algorithm
    try:
        aux1 = 10**(-0.4 * magnr)
        aux2 = 10**(-0.4 * magpsf)
        aux3 = aux1 + isdiffpos * aux2
        if aux3 > 0:
            magpsf_corr = -2.5 * np.log10(aux3)
            aux4 = aux2**2 * sigmapsf**2 - aux1**2 * sigmagnr**2

            if aux4 >= 0:
                sigmapsf_corr = np.sqrt(aux4) / aux3
            else:
                sigmapsf_corr = ZERO_MAG

            sigmapsf_corr_ext = aux2 * sigmapsf / aux3
        else:
            magpsf_corr = ZERO_MAG
            sigmapsf_corr = ZERO_MAG
            sigmapsf_corr_ext = ZERO_MAG

        return magpsf_corr, sigmapsf_corr, sigmapsf_corr_ext
    except Exception as e:
        logging.error('Object {}: {}'.format(oid, e))
        return np.nan, np.nan, np.nan


def apply_correction(candidate):
    isdiffpos = 1 if (candidate["isdiffpos"] in ["t", "1"]) else -1
    magnr = candidate["magnr"]
    magpsf = candidate['magpsf']
    sigmagnr = candidate['sigmagnr']
    sigmapsf = candidate['sigmapsf']

    magpsf_corr, sigmapsf_corr, sigmapsf_corr_ext = correction(magnr, magpsf, sigmagnr, sigmapsf, isdiffpos)

    return magpsf_corr, sigmapsf_corr, sigmapsf_corr_ext


def apply_correction_df(data):
    df = data.copy()
    df['isdiffpos'] = df['isdiffpos'].map({'t': 1., 'f': -1.})
    df["corr"] = df["distnr"] < DISTANCE_THRESHOLD
    correction_results = df.apply(
        lambda x: correction(x.magnr, x.magpsf, x.sigmagnr, x.sigmapsf, x.isdiffpos, x.objectId)
        if x["corr"]
        else (np.nan, np.nan, np.nan), axis=1, result_type="expand")

    df["magpsf_corr"], df["sigmapsf_corr"], df["sigmapsf_corr_ext"] = correction_results[0], correction_results[1], correction_results[2]

    idxmin = df.candid.idxmin()
    df["corr_magstats"] = df.loc[idxmin]["corr"]

    mask = ((df["corr"] == False) & (df.isdiffpos == -1)) | (df.corr_magstats & (df["corr"] == False)) | ((df.corr_magstats == False) & df["corr"])
    df["dubious"] = mask

    return df


def apply_parallel(df_groups, func, n=1):
    response = Parallel(n_jobs=n)(delayed(func)(group) for name, group in df_groups)
    return pd.DataFrame.from_dict(response)
