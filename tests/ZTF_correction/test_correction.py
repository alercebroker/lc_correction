import os
import fastavro
import unittest
import pandas as pd

from correction.compute import *

AVRO_PATH = "data_examples/avros"
CSV_PATH = "data_examples/csv"


def read_avro(avro_path):
    with open(avro_path, "rb") as file:
        try:
            reader = fastavro.reader(file)
            data = reader.next()
        except Exception as e:
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
