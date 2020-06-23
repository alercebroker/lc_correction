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
TRIPLE_NAN = (np.nan, np.nan, np.nan)

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(name)s.%(funcName)s: %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',)


def validate_object(candidate, is_first_detection, stellar_object=False):
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
            if not corr_magstats:
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


def correction(magnr, magpsf, sigmagnr, sigmapsf, isdiffpos, oid=None, candid=None):
    """Do correction in psf magnitude

    Parameters
    ----------
    magnr : :py:class:`float`
            Descrption
    magpsf : :py:class:`float`
            Descrption
    sigmagnr : :py:class:`float`
            Descrption
    sigmapsf : :py:class:`float`
            Descrption
    isdiffpos : :py:class:`int`
            Descrption
    oid : :py:class:`str`
            Descrption
    candid : :py:class:`str`
            Descrption

    Returns
    -------
    :py:class:`string`
        Tuple -> magpsf_corr, sigmapsf_corr, sigmapsf_corr_ext
    """
    if magnr < 0 or magpsf < 0:
        return TRIPLE_NAN

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
        return TRIPLE_NAN


def apply_correction(candidate):
    isdiffpos = 1 if (candidate["isdiffpos"] in ["t", "1"]) else -1
    magnr = candidate["magnr"]
    magpsf = candidate['magpsf']
    sigmagnr = candidate['sigmagnr']
    sigmapsf = candidate['sigmapsf']

    magpsf_corr, sigmapsf_corr, sigmapsf_corr_ext = correction(magnr, magpsf, sigmagnr, sigmapsf, isdiffpos)

    return magpsf_corr, sigmapsf_corr, sigmapsf_corr_ext


def near_stellar(first_distnr, first_distpsnr1, first_sgscore1, first_chinr, first_sharpnr):
    nearZTF = 0 <= first_distnr < DISTANCE_THRESHOLD
    nearPS1 = 0 <= first_distpsnr1 < DISTANCE_THRESHOLD
    stellarPS1 = first_sgscore1 > SCORE_THRESHOLD
    stellarZTF = first_chinr < CHINR_THRESHOLD and SHARPNR_MIN < first_sharpnr < SHARPNR_MAX
    return nearZTF, nearPS1, stellarPS1, stellarZTF


def is_stellar(nearZTF, nearPS1, stellarPS1, stellarZTF):
    return (nearZTF & nearPS1 & stellarPS1) | (nearZTF & ~nearPS1 & stellarZTF)


def is_dubious(corrected, isdiffpos, corr_magstats):
    return (~corrected & (isdiffpos == -1)) | (corr_magstats & ~corrected) | (~corr_magstats & corrected)


def dmdt(magpsf_first, sigmapsf_first, nd_diffmaglim, mjd_first, nd_mjd):
    dm_sigma = magpsf_first + sigmapsf_first - nd_diffmaglim
    dt = mjd_first - nd_mjd
    dmsigdt = (dm_sigma / dt)
    return dm_sigma, dt, dmsigdt


def apply_correction_df(df, parallel=False):
    # create copy of dataframe
    # df = data.copy()
    df.set_index("candid", inplace=True)

    df['isdiffpos'] = df['isdiffpos'].map({'t': 1., 'f': -1., '1': 1., '0': -1.})
    df["corrected"] = df["distnr"] < DISTANCE_THRESHOLD
    correction_results = df.apply(
        lambda x: correction(x.magnr, x.magpsf, x.sigmagnr, x.sigmapsf, x.isdiffpos, x.objectId)
        if x["corrected"]
        else (np.nan, np.nan, np.nan), axis=1, result_type="expand")

    df["magpsf_corr"], df["sigmapsf_corr"], df["sigmapsf_corr_ext"] = correction_results[0], correction_results[1], correction_results[2]

    corr_magstats = df.loc[df.index.min()]["corrected"]
    df["dubious"] = is_dubious(df["corrected"], df['isdiffpos'], corr_magstats)
    if parallel:
        return df
    return df.drop(["objectId", "fid"], axis=1)


def apply_mag_stats(df):
    response = {}
    # minimum and maximum candid
    idxmin = df.index.min()
    idxmax = df.index.max()

    # corrected at the first detection?
    response['corrected'] = df.loc[idxmin]["corrected"]
    response["nearZTF"], response["nearPS1"], response["stellarZTF"], response["stellarPS1"] = near_stellar(df.loc[idxmin].distnr,
                                                                                                            df.loc[idxmin].distpsnr1,
                                                                                                            df.loc[idxmin].sgscore1,
                                                                                                            df.loc[idxmin].chinr,
                                                                                                            df.loc[idxmin].sharpnr)
    response["stellar"] = is_stellar(response["nearZTF"], response["nearPS1"], response["stellarZTF"], response["stellarPS1"])
    # number of detections and dubious detections
    response["ndet"] = df.shape[0]
    response["ndubious"] = df.dubious.sum()

    # reference id
    response["nrfid"] = len(df.rfid.dropna().unique())

    # psf magnitude statatistics
    response["magpsf_mean"] = df.magpsf.mean()
    response["magpsf_median"] = df.magpsf.median()
    response["magpsf_max"] = df.magpsf.max()
    response["magpsf_min"] = df.magpsf.min()
    response["magpsf_first"] = df.loc[idxmin].magpsf
    response["sigmapsf_first"] = df.loc[idxmin].sigmapsf
    response["magpsf_last"] = df.loc[idxmax].magpsf

    # psf corrected magnitude statatistics
    response["magpsf_corr_mean"] = df.magpsf_corr.mean()
    response["magpsf_corr_median"] = df.magpsf_corr.median()
    response["magpsf_corr_max"] = df.magpsf_corr.max()
    response["magpsf_corr_min"] = df.magpsf_corr.min()
    response["magpsf_corr_first"] = df.loc[idxmin].magpsf_corr
    response["magpsf_corr_last"] = df.loc[idxmax].magpsf_corr

    # corrected psf magnitude statistics
    response["magap_mean"] = df.magap.mean()
    response["magap_median"] = df.magap.median()
    response["magap_max"] = df.magap.max()
    response["magap_min"] = df.magap.min()
    response["magap_first"] = df.loc[idxmin].magap
    response["magap_last"] = df.loc[idxmax].magap

    # time statistics
    response["first_mjd"] = df.loc[idxmin].mjd
    response["last_mjd"] = df.loc[idxmax].mjd
    return pd.Series(response)


def apply_objstats_from_correction(df):
    response = {}
    idxmax = df.mjd.idxmax()
    response["ndethist"] = df.loc[idxmax].ndethist
    response["ncovhist"] = df.loc[idxmax].ncovhist
    response["mjdstarthist"] = df.loc[idxmax].jdstarthist - 2400000.5
    response["mjdendhist"] = df.loc[idxmax].jdendhist - 2400000.5
    response["meanra"] = df.ra.mean()
    response["meandec"] = df.dec.mean()
    response["sigmara"] = df.ra.std()
    response["sigmadec"] = df.dec.std()
    response["firstmjd"] = df.mjd.min()
    response["lastmjd"] = df.loc[idxmax].mjd
    response["deltamjd"] = response["lastmjd"] - response["firstmjd"]
    return pd.Series(response)


def apply_objstats_from_magstats(df):
    response = {}
    response["nearZTF"] = df.nearZTF.all()
    response["nearPS1"] = df.nearPS1.all()
    response["stellar"] = df.stellar.all()
    response["corrected"] = df.corrected.all()
    response["ndet"] = df.ndet.sum()  # sum of detections in all bands
    response["ndubious"] = df.ndubious.sum()  # sum of dubious corrections in all bands
    if len(df.index) == 2:
        oid = df.index[0][0]
        response["g-r_max"] = df.loc[oid, 1]["magpsf_min"] - df.loc[oid, 2]["magpsf_min"]  # 1=g ; 2=r
        response["g-r_max_corr"] = df.loc[oid, 1]["magpsf_corr_min"] - df.loc[oid, 2]["magpsf_corr_min"]
        response["g-r_mean"] = df.loc[oid, 1]["magpsf_mean"] - df.loc[oid, 2]["magpsf_mean"]
        response["g-r_mean_corr"] = df.loc[oid, 1]["magpsf_corr_mean"] - df.loc[oid, 2]["magpsf_corr_mean"]
    else:
        response["g-r_max"] = np.nan
        response["g-r_max_corr"] = np.nan
        response["g-r_mean"] = np.nan
        response["g-r_mean_corr"] = np.nan
    return pd.Series(response)


def apply_object_stats_df(corrected, magstats):
    basic_stats = corrected.groupby("objectId").apply(apply_objstats_from_correction)
    obj_magstats = magstats.groupby("objectId").apply(apply_objstats_from_magstats)
    return basic_stats.join(obj_magstats)


def do_dmdt(nd, magstats, dt_min=0.5):
    response = {}
    nd.reset_index(inplace=True)
    mjd_first = magstats.first_mjd.iloc[0]
    mask = nd.mjd < mjd_first - dt_min
    response["close_nondet"] = nd.loc[mask].mjd.max() < nd.loc[nd.mjd < mjd_first].mjd.max()
    # is there some non-detection before the first detection
    if mask.sum() > 0:
        magpsf_first = magstats.magpsf_first.iloc[0]
        sigmapsf_first = magstats.sigmapsf_first.iloc[0]
        # assume the worst case
        dm_sigma, dt, dmsigdt = dmdt(magpsf_first,
                                     sigmapsf_first,
                                     nd.loc[mask].diffmaglim,
                                     mjd_first,
                                     nd.loc[mask].mjd)
        idxmin = dmsigdt.idxmin()
        response["dmdt_first"] = dmsigdt.loc[idxmin]
        response["dm_first"] = magpsf_first - nd.diffmaglim.loc[idxmin]
        response["sigmadm_first"] = sigmapsf_first - nd.diffmaglim.loc[idxmin]
        response["dt_first"] = dt.loc[idxmin]
    else:
        response["dmdt_first"] = np.nan
        response["dm_first"] = np.nan
        response["sigmadm_first"] = np.nan
        response["dt_first"] = np.nan

    return pd.Series(response)


def do_dmdt_df(magstats, non_dets):
    g_mags = magstats.groupby(["objectId", "fid"])
    g_nd = non_dets.groupby(["objectId", "fid"])
    idxs = g_mags.indices.keys()
    result = []
    for idx in idxs:
        if idx in g_nd.indices.keys():
            non_dets_g = g_nd.get_group(idx)
            magstats_g = g_mags.get_group(idx)
            resp = do_dmdt(non_dets_g, magstats_g)
            resp["objectId"], resp["fid"] = idx[0], idx[1]
            result.append(resp)
    return pd.DataFrame.from_records(result, index=["objectId", "fid"])


def apply_parallel(df_groups, func, n=1):
    response = Parallel(n_jobs=n)(delayed(func)(group) for name, group in df_groups)
    return pd.concat(response)
