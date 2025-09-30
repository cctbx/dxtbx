from __future__ import annotations

import os

import dials_data.download

df = dials_data.download.DataFetcher()
data_dir = df("image_examples", pathlib=True)
image_examples = []
for f in data_dir.iterdir():
    ext = os.path.splitext(f)[1]
    if ext == ".h5" and "_master." not in str(f):
        continue
    if ext in (".expt", ".json"):
        continue
    image_examples.append(str(f))
