import numpy as np
import pandas as pd
from astropy.time import Time
#import logging

#level = logging.INFO
#logging.basicConfig(level=level,
#                    format='%(asctime)s %(levelname)s %(name)s.%(funcName)s: %(message)s',
#                    datefmt='%Y-%m-%d %H:%M:%S',)--


DISTANCE_THRESHOLD = 1.4
SCORE_THRESHOLD = 0.4
CHINR_THRESHOLD = 2
SHARPNR_MAX = 0.1
SHARPNR_MIN = -0.13
WEIRD = 100


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


def correction(magnr, magpsf, sigmagnr, sigmapsf, isdiffpos): #Correction Algorithm 
    aux1 = np.power(10, -0.4 * magnr)
    aux2 = isdiffpos * np.power(10, -0.4 * magpsf)
    aux3 = aux1 + aux2
    if aux3 > 0:
        magpsf_corr = -2.5 * np.log10(aux3)
        aux4 = np.square(aux2) * np.square(sigmapsf) - np.square(aux1) * np.square(sigmagnr)

        if aux4 >= 0:
            sigmapsf_corr = np.sqrt(aux4) / aux3
        else:
            sigmapsf_corr = WEIRD
                
        sigmapsf_corr_ext = aux2 * sigmapsf / aux3
    else:
        magpsf_corr = WEIRD
        sigmapsf_corr = WEIRD
        sigmapsf_corr_ext = WEIRD

    return magpsf_corr, sigmapsf_corr, sigmapsf_corr_ext


def apply_correction(candidate): 
    isdiffpos = 1 if (candidate["isdiffpos"] in ["t", "1"]) else -1
    magnr = candidate["magnr"]
    magpsf = candidate['magpsf']
    sigmagnr = candidate['sigmagnr']
    sigmapsf = candidate['sigmapsf']

    magpsf_corr, sigmapsf_corr, sigmapsf_corr_ext = correction(magnr, magpsf, sigmagnr, sigmapsf, isdiffpos)

    return magpsf_corr, sigmapsf_corr, sigmapsf_corr_ext


def get_prv_candidates(message, light_curve):
    oid = message["objectId"]
    prv_cands = []
    non_dets = []
    if message["prv_candidates"]:
        for prv_cand in message["prv_candidates"]:
            mjd = prv_cand["jd"] - 2400000.5
            if prv_cand["diffmaglim"] is not None:
                non_detection_args = {
                    "diffmaglim": prv_cand["diffmaglim"],
                    "oid": oid,
                    "mjd": mjd
                }

                dt = Time(mjd, format="mjd")
                filters = {"datetime": dt.datetime, "fid": prv_cand["fid"], "oid": message["objectId"]}
                found = list(filter(lambda non_det: ((non_det["datetime"] == filters["datetime"]) and
                                                     (non_det["fid"] == filters["fid"]) and
                                                     (non_det["oid"] == filters["oid"])),
                                    light_curve["non_detections"]))
                if len(found) == 0:
                    non_detection_args.update(filters)
                    if non_detection_args not in non_dets:
                        non_dets.append(non_detection_args)
                        light_curve["non_detections"].append(non_detection_args)
            else:
                found = list(filter(lambda det: det['candid'] == prv_cand["candid"], light_curve["detections"]))
                if len(found) == 0:
                    # prv_cand.update(self.correct_message(prv_cand))
                    detection_args = {
                        "mjd": prv_cand["jd"] - 2400000.5,
                        "fid": prv_cand["fid"],
                        "ra": prv_cand["ra"],
                        "dec": prv_cand["dec"],
                        "rb": prv_cand["rb"],
                        "magap": prv_cand["magap"],
                        "magap_corr": prv_cand["magap_corr"],
                        "magpsf": prv_cand["magpsf"],
                        "magpsf_corr": prv_cand["magpsf_corr"],
                        "sigmagap": prv_cand["sigmagap"],
                        "sigmagap_corr": prv_cand["sigmagap_corr"],
                        "sigmapsf": prv_cand["sigmapsf"],
                        "sigmapsf_corr": prv_cand["sigmapsf_corr"],
                        "oid": message["objectId"],
                        "alert": prv_cand,
                        "candid": str(prv_cand["candid"]),
                        "parent_candidate": str(message["candid"])
                    }
                    prv_cands.append(detection_args) 
                    light_curve["detections"].append(detection_args)
        return light_curve, prv_cands
