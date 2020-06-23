import os
import unittest

from correction.compute import *
from correction.helpers import *

PARQUET_PATH = "data_examples/parquets"


class TestHelpersCorrectionChain(unittest.TestCase):
    def setUp(self) -> None:
        self.detections = pd.read_parquet(os.path.join(PARQUET_PATH, "10_objects_detections.parquet"))
        self.detections["mjd"] = self.detections.jd - 2400000.5
        del self.detections["jd"]

    def test_get_data_quality(self):
        result = get_data_quality(self.detections)
        self.assertEqual(len(result), 1381)

    def test_get_ss_ztf(self):
        result = get_ss_ztf(self.detections)
        self.assertEqual(len(result), 1381)

    def test_get_ps1_ztf(self):
        result = get_ps1_ztf(self.detections)
        self.assertEqual(len(result), 16)

    def test_helpers_clean_corrected(self):
        corrected = self.detections.groupby(["objectId", "fid"]).apply(apply_correction_df)
        corrected = get_clean_corrected(corrected)
        self.assertEqual(len(corrected), 1381)
