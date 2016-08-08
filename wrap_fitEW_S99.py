import jrr
import os
import sys
import re
from matplotlib import pyplot as plt
import pandas
import numpy as np


sorted_by_age = ('S0033+0242', 'rcs0327-knotE', 'rcs0327-knotG', 'S0108+0624', 'S0957+0509', 'Horseshoe', 'rcs0327-knotU', 'S2111-0114', 'S1429+1202', 'S0004-0103','S0900+2234', 'stack-A', 'S1527+0652', 'S1226+2152') # 'S1458-0023', 

for rootname in sorted_by_age[0:1] :
#for rootname in sorted_by_age :
    df = jrr.mage.open_S99_spectrum(rootname)
    # Fit some EWs here...
