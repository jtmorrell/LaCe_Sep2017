import os
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from npat.spectroscopy import spectrum

for fnm in list(os.walk('Spe_files'))[0][2]:
	sp = spectrum(fnm,'Spe_files')	
	sp.plot(wfit=True,printout_log='printlogs/'+fnm.replace('.Spe','.txt'),saveas='plots/'+fnm.replace('.Spe','.png'))

