import os
import fastavro
import unittest
import pandas as pd

from correction.compute import *

AVRO_PATH = "data_examples/avros"
CSV_PATH = "data_examples/csv"
PARQUET_PATH = "data_examples/parquets"


def read_avro(avro_path):
    with open(avro_path, "rb") as file:
        try:
            reader = fastavro.reader(file)
            data = reader.next()
        except Exception:
           return None
    return data


def get_avros(oid):
    path = os.path.join(AVRO_PATH, oid)
    avros = sorted(os.listdir(path))
    data = []
    for a in avros:
        content = read_avro(os.path.join(path, a))
        if content:
            data.append(content)
    return data


class TestZTF18aazxcwf(unittest.TestCase):
    def setUp(self) -> None:
        self.oid = "ZTF18aazxcwf"
        self.avros = get_avros(self.oid)
        self.is_first_detection = {"object": True, 1: True, 2: True}
        self.first_magnr = {}
        self.stellar_object = False

    # Validate object to a RRL object
    def test_validate_object(self):
        data = self.avros[0]
        fid = data["candidate"]["fid"]

        if self.is_first_detection[fid]:
            self.first_magnr[fid] = data["candidate"]["magnr"]

        stellar_object, stellar_magstats = validate_object(data["candidate"], self.is_first_detection[fid], self.stellar_object)
        self.stellar_object = stellar_object
        self.assertTrue(stellar_object)
        self.assertTrue(stellar_magstats)

    def test_correction(self):
        candidate = self.avros[0]["candidate"]
        correct_candidate = apply_correction(candidate)
        self.assertEqual(len(correct_candidate), 3)
        self.assertAlmostEqual(correct_candidate[0], 16.387620303061944)
        self.assertAlmostEqual(correct_candidate[1], 100.)
        self.assertAlmostEqual(correct_candidate[2], 0.02641869927986706)


class TestDataframeCorrection(unittest.TestCase):
    def setUp(self) -> None:
        self.data = pd.read_csv(os.path.join(CSV_PATH, "raw_detections.csv"))

    def test_apply_correction_df_c1(self):
        dflarge = self.data.groupby(["objectId", "fid"]).apply(apply_correction_df)
        dubious = dflarge.loc[dflarge.dubious]
        self.assertEqual(len(dubious), 1)
        self.assertTrue(dubious.dubious.values[0])

    def test_apply_correction_df_c2(self):
        data = self.data.copy()
        data.at[180, "isdiffpos"] = 'f'
        dflarge = data.groupby(["objectId", "fid"]).apply(apply_correction_df)
        dubious = dflarge.loc[dflarge.dubious]
        example = dubious.iloc[1]
        self.assertEqual(len(dubious), 2)
        self.assertFalse(example["corrected"])
        self.assertTrue(example["dubious"])

    def test_apply_correction_df_c3(self):
        data = self.data.copy()
        data.at[150, "distnr"] = 0.1
        dflarge = data.groupby(["objectId", "fid"]).apply(apply_correction_df)
        dubious = dflarge.loc[dflarge.dubious]
        example = dubious.iloc[1]
        self.assertEqual(len(dubious), 2)
        self.assertTrue(example["corrected"])
        self.assertTrue(example["dubious"])


class TestDataframeCorrectionChain(unittest.TestCase):
    def setUp(self) -> None:
        self.detections = pd.read_parquet(os.path.join(PARQUET_PATH, "10_objects_detections.parquet"))
        self.non_detections = pd.read_parquet(os.path.join(PARQUET_PATH, "10_objects_non_detections.parquet"))
        self.corrected = pd.read_parquet(os.path.join(PARQUET_PATH, "10_objects_corrected.parquet"))
        self.magstats = pd.read_parquet(os.path.join(PARQUET_PATH, "10_objects_magstats.parquet"))
        self.corrected.reset_index(inplace=True)
        self.magstats.reset_index(inplace=True)
        self.unique_objects = 10
        self.corrected_cols = ["corrected", "magpsf_corr", "sigmapsf_corr", "sigmapsf_corr_ext", "dubious"]
        self.magstats_cols = ['corrected', 'nearZTF', 'nearPS1', 'stellarZTF', 'stellarPS1', 'stellar', 'ndet',
                              'ndubious', 'nrfid', 'magpsf_mean', 'magpsf_median', 'magpsf_max', 'magpsf_min',
                              'magpsf_first', 'sigmapsf_first', 'magpsf_last', 'magpsf_corr_mean', 'magpsf_corr_median',
                              'magpsf_corr_max', 'magpsf_corr_min', 'magpsf_corr_first', 'magpsf_corr_last', 'magap_mean',
                              'magap_median', 'magap_max', 'magap_min', 'magap_first', 'magap_last', 'first_mjd', 'last_mjd']
        self.objstats_cols = ['ndethist', 'ncovhist', 'mjdstarthist', 'mjdendhist', 'meanra', 'meandec', 'sigmara',
                              'sigmadec', 'firstmjd', 'lastmjd', 'deltamjd', 'nearZTF', 'nearPS1', 'stellar', 'corrected',
                              'ndet', 'ndubious']
        self.detections["mjd"] = self.detections.jd - 2400000.5
        del self.detections["jd"]
        self.non_detections["mjd"] = self.non_detections.jd - 2400000.5
        del self.non_detections["jd"]

    def test_check_dataset(self):
        unique_detections = len(self.detections.objectId.unique())
        unique_non_detections = len(self.non_detections.objectId.unique())
        self.assertEqual(unique_detections, self.unique_objects)
        self.assertEqual(unique_non_detections, self.unique_objects)

    def test_apply_correction_df(self):
        corrected = self.detections.groupby(["objectId", "fid"]).apply(apply_correction_df)
        for col in self.corrected_cols:
            self.assertIn(col, corrected.columns)
        self.assertEqual(len(corrected.index.levels[0]), self.unique_objects)
        self.assertEqual(corrected.corrected.values.sum(), 1373)

    def test_apply_mag_stats(self):
        magstats = self.corrected.groupby(["objectId", "fid"]).apply(apply_mag_stats)
        self.assertEqual(len(magstats.index.levels[0]), 10)
        self.assertEqual(len(magstats), 16)
        for col in self.magstats_cols:
            self.assertIn(col, magstats.columns)

    def test_apply_object_stats(self):
        objstats = apply_object_stats_df(self.corrected, self.magstats)
        self.assertEqual(len(objstats), 10)
        for col in self.objstats_cols:
            self.assertIn(col, objstats.columns)
        self.assertEqual(objstats.ndubious.sum(), 4)

    def test_apply_do_dmdt_df(self):
        dmdt = do_dmdt_df(self.magstats, self.non_detections)
        self.assertEqual(len(dmdt.index.levels[0]), 10)
        self.assertEqual(len(dmdt), 16)
        self.assertEqual(dmdt.close_nondet.sum(), 1)
