# after delete this file
from correction.compute import *
import fastavro
import pickle

avro = 'data_examples/avros/1184159740015015000.avro'
with open(avro, "rb") as f:
    reader = fastavro.reader(f)
    data = reader.next()
a = apply_correction(data["candidate"], first_magnr=20)

f = open("data_examples/pickles/ZTF20aatvpww_lc.pkl", "rb")
lc = pickle.load(f)

i = 0
for a in lc["detections"]:
    if a["candid"] == "1184159740015015000":
        lc["detections"].pop(i)
        i+=1
i = 0
print(lc["non_detections"])
print("\n\n\n")
print(data["prv_candidates"])

b = get_prv_candidates(data, lc)
