import os
import numpy as np
from scipy.optimize import curve_fit
from ast import literal_eval
import matplotlib.pyplot as plt
import re
import datetime as dtm
from npat.spectroscopy import spectrum

def map_efficiency(energy,*effcal):
	return np.exp(effcal[0]*np.log(energy)**2+effcal[1]*np.log(energy)+effcal[2])

def scrape_efficiency_data():
	eff_fit = {}
	eff_unc = {}
	for shelf in [5,12,16,18]:
		eff = []
		with open('summaries/ctr76_shelf'+str(shelf)+'.html','r') as f:
			dat = f.read().split('\n')
			for ln in dat:
				if re.findall(r'<font color="#FFFFFF" size="4">',ln):
					eff.append(float(ln.split('<font color="#FFFFFF" size="4">')[1].split('<')[0].replace(',','')))
		eff = np.array(eff).reshape((len(eff)/3,3))
		eng,eff = eff[:,0],eff[:,1]
		fit,unc = curve_fit(map_efficiency,eng,eff,p0=[0.0,0.0,0.0])
		eff_fit[shelf] = fit.tolist()
		eff_unc[shelf] = unc.tolist()
		# plt.plot(eng,eff,ls='None',marker='o')
		# plt.plot(eng,map_efficiency(eng,*fit))
		# plt.xlabel('Energy (keV)')
		# plt.ylabel('Efficiency')
		# plt.show()
	with open('efficiency_fits.txt','w+') as f:
		f.write('\n'.join([str(eff_fit),str(eff_unc)]))




def convert_to_Spe():
	with open('efficiency_fits.txt','r') as f:
		lns = f.read().split('\n')
		eff_fit = literal_eval(lns[0])
		eff_unc = literal_eval(lns[1])
	for fnm in list(os.walk('data_text_files'))[0][2]:
		with open('data_text_files/'+fnm,'r') as f:
			dat = f.read().split('\n')
			ctr_fnm = str(int(dat[0].split(' ')[-1].split('.')[0].split('_')[1]))
			foil_id = 'La'+dat[1].split(' ')[4]
			start_time = dtm.datetime.strftime(dtm.datetime.strptime(' '.join(dat[1].split(' ')[6:8]).split('.')[0],'%Y-%m-%d %H:%M:%S'),'%m/%d/%Y %H:%M:%S')
			spec = [i.split(' ')[1] for i in dat[4:-2]]
		with open('summaries/s'+ctr_fnm+'.txt','r') as f:
			dat = ' '.join(f.read().split('\n')[3:7])+' EOL:'
			fields = [i.split(' ')[-1] for i in dat.split(':')[:-1]]
			ts_idx = fields.index('Timestamp')
			fields = fields[:ts_idx+1]+fields[ts_idx+3:]
			fields[fields.index('acc/rej')] = 'Count acc/rej'
			fields[fields.index('Prog/Yr/Seq')] = 'LP Prog/Yr/Seq'
			header = {f:dat.split(f+':')[1].split(fields[n+1])[0].strip() for n,f in enumerate(fields[:-1])}
			LT,RT = str(int(60.0*float(header['Live']))),str(int(60.0*float(header['True'])))
			shelf_N = int(header['Shelf'])
			m,b = header['Gain'].split(' ')[0],header['Zero'].split(' ')[0]
		ss = '\n'.join(['$SPEC_ID:','No sample description was entered.','$SPEC_REM:','DET# 9480','DETDESC# HPGE','AP# Maestro Version 7.01','$DATE_MEA:',''])
		ss += '\n'.join([start_time,'$MEAS_TIM:',LT+' '+RT,'$DATA:','0 '+str(len(spec)),''])
		ss += '\n'.join([''.join([' ']*(8-len(i)))+i for i in spec])
		ss += '\n'.join(['','$ROI:','0','$PRESETS:','None','0','0','$ENER_FIT:',b+' '+m,'$MCA_CAL:','3',b+' '+m+' 0.00 keV','$SHAPE_CAL:','3',' '.join(map(str,eff_fit[shelf_N])),''])
		with open('Spe_files/'+ctr_fnm+'_'+foil_id+'.Spe','w+') as f:
			f.write(ss)
		sp = spectrum(ctr_fnm+'_'+foil_id+'.Spe','Spe_files/')
		sp.update_params(meta={'res':0.05,'R':0.2,'alpha':0.9,'unc_eff':eff_unc[shelf_N],'istp':['133CEm','133CEg','133LA','132CE','132LA','139CE','137CEg','137CEm','137LA','135CE','135LA','134LA','132CS','7BE','133BAm','133BAg']})
		sp.write_maestro()
		
def convert_monitors_to_Spe():
	eff = []
	with open('summaries/ctr76_shelf5.html','r') as f:
		dat = f.read().split('\n')
		for ln in dat:
			if re.findall(r'<font color="#FFFFFF" size="4">',ln):
				eff.append(float(ln.split('<font color="#FFFFFF" size="4">')[1].split('<')[0].replace(',','')))
	eff = np.array(eff).reshape((len(eff)/3,3))
	eng,eff = eff[:,0],eff[:,1]
	fit,unc = curve_fit(map_efficiency,eng,eff,p0=[0.0,0.0,0.0])
	eff_fit = fit.tolist()
	eff_unc = unc.tolist()
	for fnm in list(os.walk('data_text_files/copper_monitors'))[0][2]:
		with open('data_text_files/copper_monitors/'+fnm,'r') as f:
			dat = f.read().split('\n')
			ctr_fnm = str(int(dat[0].split(' ')[-1].split('.')[0].split('_')[1]))
			foil_id = 'Cu'+dat[1].split(' ')[4]
			start_time = dtm.datetime.strftime(dtm.datetime.strptime(' '.join(dat[1].split(' ')[6:8]).split('.')[0],'%Y-%m-%d %H:%M:%S'),'%m/%d/%Y %H:%M:%S')
			spec = [i.split(' ')[1] for i in dat[4:-2]]
		with open('summaries/copper_monitors/s'+ctr_fnm+'.txt','r') as f:
			dat = ' '.join(f.read().split('\n')[3:7])+' EOL:'
			fields = [i.split(' ')[-1] for i in dat.split(':')[:-1]]
			ts_idx = fields.index('Timestamp')
			fields = fields[:ts_idx+1]+fields[ts_idx+3:]
			fields[fields.index('acc/rej')] = 'Count acc/rej'
			fields[fields.index('Prog/Yr/Seq')] = 'LP Prog/Yr/Seq'
			header = {f:dat.split(f+':')[1].split(fields[n+1])[0].strip() for n,f in enumerate(fields[:-1])}
			LT,RT = str(int(60.0*float(header['Live']))),str(int(60.0*float(header['True'])))
			shelf_N = int(header['Shelf'])
			m,b = header['Gain'].split(' ')[0],header['Zero'].split(' ')[0]
		ss = '\n'.join(['$SPEC_ID:','No sample description was entered.','$SPEC_REM:','DET# 9480','DETDESC# HPGE','AP# Maestro Version 7.01','$DATE_MEA:',''])
		ss += '\n'.join([start_time,'$MEAS_TIM:',LT+' '+RT,'$DATA:','0 '+str(len(spec)),''])
		ss += '\n'.join([''.join([' ']*(8-len(i)))+i for i in spec])
		ss += '\n'.join(['','$ROI:','0','$PRESETS:','None','0','0','$ENER_FIT:',b+' '+m,'$MCA_CAL:','3',b+' '+m+' 0.00 keV','$SHAPE_CAL:','3',' '.join(map(str,eff_fit)),''])
		with open('Spe_files/'+ctr_fnm+'_'+foil_id+'.Spe','w+') as f:
			f.write(ss)
		sp = spectrum(ctr_fnm+'_'+foil_id+'.Spe','Spe_files/')
		sp.update_params(meta={'res':0.05,'R':0.2,'alpha':0.9,'unc_eff':eff_unc,'istp':['62ZN','63ZN','65ZN','7BE','58CO','61CU','57CO','64CO','56CO']})
		sp.write_maestro()
		



convert_to_Spe()
# convert_monitors_to_Spe()