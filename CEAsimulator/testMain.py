from TomSim import TomSim
from tools.utils import linear_interpolation
import numpy as np

sim = TomSim(co2=234, temperature=23, fruit_per_truss=7, lon=-83.363174, lat=33.976097, debug=False)

sim.start_simulation()