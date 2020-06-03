import pandas as pd
import numpy as np
from astropy.time import Time
#import logging

#level = logging.INFO
#logging.basicConfig(level=level,
#                    format='%(asctime)s %(levelname)s %(name)s.%(funcName)s: %(message)s',
#                    datefmt='%Y-%m-%d %H:%M:%S',) test


def validate_object():
    return


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


def correct_magnitude(magref, sign, magdiff):
    result = np.nan
    try:
        aux = np.power(10, (-0.4 * magref)) + sign * np.power(10, (-0.4 * magdiff))
        result = -2.5 * np.log10(aux)
    except Exception as e:
        #logging.exception("Correct magnitude failed: {}".format(e))
        return
    return result


def correct_sigma_mag(magref, sigmagref, sign, magdiff, sigmagdiff):
    result = np.nan
    try:
        auxref = np.power(10, (-0.4 * magref))
        auxdiff = np.power(10, (-0.4 * magdiff))
        aux = auxref + sign * auxdiff

        result = np.sqrt(np.power((auxref * sigmagref), 2) +
                         np.power((auxdiff * sigmagdiff), 2)) / aux

    except Exception as e:
        #logging.exception("Correct sigma magnitude failed: {}".format(e))
        return
    return result