import os
import fastavro
import unittest

from correction.compute import *

AVRO_PATH = "tests/data_examples/"


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
        self.assertTrue(stellar_object)
        self.assertTrue(stellar_magstats)
