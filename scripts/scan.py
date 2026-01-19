import firefly as fly
import numpy as np

Ts = [0.25, 0.1, 0.01, 0.001]

mb_cfg = fly.unpack("2_mb.cfg")

outdir = './data/responses'
mb_cfg["CONTROL"]["outdir"] = outdir

for i in range(len(Ts)):
    mb_cfg["CONTROL"]["prefix"] = Ts[i]
    mb_cfg["SYSTEM"]["Temperature"] = Ts[i]
    output, error = fly.launch(mb_cfg)
    match = fly.grep(output, "Max χ:")
    print(f"T = {Ts[i]}, X = {match}")

