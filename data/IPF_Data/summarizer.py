import os
import numpy as np

foils = {}
for logfile in list(os.walk('printlogs'))[0][2]:
	with open('printlogs/'+logfile,'r') as f:
		foil = logfile.split('_')[1].split('.')[0]
		lns = f.read().split('\n')[:-1]
		isotopes = map(lambda l:l.split(':')[0],lns)
		activities = map(lambda l:float(l.split(',')[3].split('=')[1].split('(')[0]),lns)
		ac = {i:[] for i in list(set(isotopes))}
		for n,a in enumerate(activities):
			ac[isotopes[n]].append(a)
		ac = {a:np.average(ac[a])/3.7E4 for a in ac}
	if foil in foils:
		foils[foil].append(ac)
	else:
		foils[foil] = [ac]

results = {}
for foil in foils:
	dat = foils[foil]
	results[foil] = {}
	for istp in list(set([istp for ac in dat for istp in ac])):
		if istp=='64CO':
			continue
		results[foil][istp] = np.average([d[istp] for d in dat if istp in d])

ss = 'Foil, Isotope, Activity (uCi)\n'
for foil in results:
	for istp in results[foil]:
		ss += ','.join([foil,istp,str(round(results[foil][istp],3))])+'\n'
	ss += ',,\n'
with open('activities.csv','w+') as f:
	f.write(ss)