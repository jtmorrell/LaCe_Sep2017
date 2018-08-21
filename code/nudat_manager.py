import re,urllib2,os,numpy as np
from scipy.interpolate import interp1d
import sqlite3
_talys_path = '/home/jmorrell/bin'
_empire_path = '/home/jmorrell/empire'
_db_connection = sqlite3.connect('../data/peak_data.db')
_db = _db_connection.cursor()


class manager(object):
	def __init__(self):
		self.db_connection = _db_connection
		self.db = _db
	def create_plots_dir(self):
		os.chdir('..')
		os.system('mkdir plots')
		os.chdir('plots')
		for dr in ['calibration','decay_curves','peak_fits','cross_sections','monitors']:
			os.system('mkdir '+dr)
		os.chdir('../data')
		os.system('mkdir cross_sections')
		os.chdir('../code')
	def run_talys(self):
		os.chdir('talys')
		os.system(_talys_path+'/talys < input > output')
		os.chdir('..')
		print 'talys output finished'
		self.update_talys()
	def update_talys(self):
		self.db.execute('DELETE FROM model_xs WHERE model=?',('talys',))
		files = {'132CS':'rp055132.tot','133BAg':'rp056133.L00','133BAm':'rp056133.L02','135LA':'rp057135.tot','136LA':'rp057136.tot','134CE':'rp058134.tot','135CE':'rp058135.tot','137CEg':'rp058137.L00','137CEm':'rp058137.L02','139CE':'rp058139.tot'}
		for istp in files:
			self.db.executemany('INSERT INTO model_xs VALUES(?,?,?,?)',[(istp,float(i.split(' ')[1]),float(i.split(' ')[2]),'talys') for i in open('talys/'+files[istp],'r').read().split('\n')[:-1] if not i.startswith('#')])
		self.db_connection.commit()
	def run_empire(self):
		fnm = 'la139'
		E_range = np.arange(32,60.5,0.5)
		lines = open('empire/'+fnm+'.inp','r').read().split('\n')[1:]
		for eng in E_range:
			f = open('empire/'+fnm+'.inp','w+')
			f.write('\n'.join([str(round(eng,1))+' ;INCIDENT ENERGY (IN LAB)']+lines))
			f.close()
			os.chdir('empire')
			os.system(_empire_path+'/scripts/runE '+fnm)
			os.chdir('..')
			try:
				for istp in ['137CEg','137CEm','132CS','133BAm','133BAg','134CE','135CE','135LA','136LA','139CE']:
					if istp in ['137CEg','137CEm','133BAg','133BAm']:
						n = {'137CEm':'  58-Ce-137 isomer','137CEg':'  58-Ce-137 ground','133BAm':'  56-Ba-133 isomer','133BAg':'  56-Ba-133 ground'}[istp]
						xs = [float([k for k in i.split(' ') if k!=''][4]) for i in open('empire/'+fnm+'.out','r').read().split('\n') if i.startswith(n)][0]
					else:
						n = {'132CS':33,'134CE':13,'135CE':12,'135LA':18,'136LA':17,'139CE':8}[istp]
						xs = [float(i.split(' ')[n+1]) for i in open('empire/'+fnm+'.xsc','r').read().split('\n')[2:-1]][0]
					self.db.execute('DELETE FROM model_xs WHERE model=? AND E=? AND isotope=?',('empire',eng,istp))
					self.db.execute('INSERT INTO model_xs VALUES(?,?,?,?)',(istp,eng,xs,'empire'))
				print 'EMPIRE finished successfully for E=',eng
				self.db_connection.commit()
			except:
				print 'EMPIRE threw an error for E=',eng
	def update_exfor(self):
		files = {'58CO':'../data/exfor/58CO_exfor.txt','61CU':'../data/exfor/61CU_exfor.txt'}
		for istp in files:
			self.db.execute('DELETE FROM exfor WHERE isotope=?',(istp,))
			for ln in open(files[istp],'r').read().split('\n'):
				if not ln.startswith('#') and not ln.startswith('//'):
					d = [i.strip() for i in ln.split(' ') if i.strip()!='']
					self.db.execute('INSERT INTO exfor VALUES(?,?,?,?,?,?)',(istp,float(d[0]),float(d[1]),1e3*float(d[2]),1e3*float(d[3]),d[5]))
			self.db_connection.commit()
		print 'EXFOR updated'
	def exp_smooth(self,ls,alpha=0.3):
		R,RR,b = [ls[0]],[ls[-1]],1.0-alpha
		for i,ii in zip(ls[1:],reversed(ls[:-1])):
			R.append(alpha*i+b*R[-1])
			RR.append(alpha*ii+b*RR[-1])
		return [0.5*(R[n]+r) for n,r in enumerate(reversed(RR))]
	def update_xs_prediction(self):
		E_range = np.arange(35,60.5,0.5)
		self.db.execute('DELETE FROM xs_prediction')
		for istp in ['137CEg','137CEm','132CS','133BAm','133BAg','134CE','135CE','135LA','136LA','139CE']:
			for eng in E_range:
				xs = [float(i[2]) for i in self.db.execute('SELECT * FROM model_xs WHERE isotope=? AND E=?',(istp,round(eng,1)))]
				self.db.execute('INSERT INTO xs_prediction VALUES(?,?,?)',(istp,eng,np.average([xs])))
		for istp in ['22NA','24NA','56CO','62ZN','63ZN','65ZN']:
			xs = [(istp,float(i[2]),float(i[3])) for i in self.db.execute('SELECT * FROM monitor_xs WHERE product=?',(istp,))]
			self.db.executemany('INSERT INTO xs_prediction VALUES(?,?,?)',xs)
		for istp in ['58CO','61CU']:
			xs = sorted([(float(i[1]),float(i[3])) for i in self.db.execute('SELECT * FROM exfor WHERE isotope=?',(istp,))],key=lambda h:h[0])
			xs = [x for n,x in enumerate(xs[:-1]) if xs[n+1][0]!=x[0]]
			xs_func = interp1d([x[0] for x in xs],self.exp_smooth([x[1] for x in xs],alpha=0.1))
			self.db.executemany('INSERT INTO xs_prediction VALUES(?,?,?)',[(istp,eng,round(xs_func(eng),2)) for eng in E_range])
		self.db_connection.commit()
		print 'XS prediction up to date'
	def update_files(self):
		self.db.execute('DELETE FROM all_files')
		for p in ['../data/Calibration/','../data/Experiment/']:
			for (dirpath, dirnames, filenames) in os.walk(p):
				self.db.executemany('INSERT INTO all_files VALUES(?,?)',[(i,p) for i in filenames if i.endswith('.Spe')])
				break
		self.db_connection.commit()
		print 'All files up to date'
	def update_ensdf(self):
		all_isotopes = [p.replace('g','') for p in list(set(','.join(list([str(i[0]) for i in self.db.execute('SELECT isotopes FROM sample_map')])).split(','))) if not p.endswith('m')]
		for nucl in all_isotopes:
			print 'Updating ',nucl
			self.search_ensdf(nucl)
		print 'ENSDF decay radiation up to date'
	def ensdf_multiplier(self,s):
		if 'E' in s:
			if '.' in s:
				return 10**(float(s.split('E')[1])-len(s.split('E')[0].split('.')[1]))
			return 10**(float(s.split('E')[1]))
		elif '.' in s:
			return 10**(-1*len(s.split('.')[1]))
		return 1.0
	def search_ensdf(self,nucl='133BA'):
		link = 'http://www.nndc.bnl.gov/chart/decaysearchdirect.jsp?nuc={0}&unc=nds'.format(nucl)
		response = urllib2.urlopen(link)
		self.db.execute('DELETE FROM ensdf WHERE isotope=?',(nucl,))
		for datum in response.read().split('Dataset')[1:]:
			data = datum.split('\n')
			eng,t_half = (i.split('white>')[1].split('<i>')[0].replace('&nbsp;','').strip() for n,i in enumerate(data[5].split('<td')) if n==4 or n==6)
			unc_t_half = [i.split('white>')[1].split('<i>')[1].split('</i>')[0] for n,i in enumerate(data[5].split('<td')) if n==6][0]
			unc_t_half = float(unc_t_half)*self.ensdf_multiplier(t_half.split(' ')[0].replace('+',''))*{'ms':0.001,'s':1.0,'m':60.0,'h':3600.0,'d':86400.0,'y':31557600.0}[t_half.split(' ')[1]]
			eng,t_half = float(eng),round(float(t_half.split(' ')[0])*{'ms':0.001,'s':1.0,'m':60.0,'h':3600.0,'d':86400.0,'y':31557600.0}[t_half.split(' ')[1]],2)
			if nucl=='152EU' and eng>0:
				continue
			gammas = []
			start,N = False,0
			for ln in data:
				if ln=='<p><u>Gamma and X-ray radiation</u>:':
					start = True
				elif ln=='</table>' and start:
					start = False
					break
				if start:
					if ln.startswith('<td nowrap'):
						if N==0:
							N+=1
						elif N==1:
							gammas.append([float(ln.split('align>')[1].split(' ')[0].replace('&nbsp;',''))])
							N+=1
						elif N==2:
							gammas[-1].append(float(ln.split('align>')[1].split(' ')[0].replace('&nbsp;','')))
							if ln.split('%')[1]!='</td>':
								gammas[-1].append(float(ln.split('<i>')[1].split('</i>')[0])*self.ensdf_multiplier(ln.split('align>')[1].split(' ')[0].replace('&nbsp;','')))
							else:
								gammas[-1].append(0.1*gammas[-1][-1])
							N+=1
						else:
							N=0
			self.db.executemany('INSERT INTO ensdf VALUES(?,?,?,?,?,?,?)',[(nucl,eng,t_half,g[0],g[1],g[2],unc_t_half) for g in gammas])
		self.db_connection.commit()
	def update_alice(self):
		elements = {'58':'CE','57':'LA','56':'BA','55':'CS'}
		alice = [[l for l in i.strip().split(' ') if l!=''][:7][:-1] for i in open('../data/cross_sections/p139La_1MeV_bins.txt','r').read().split('\n')[48:1072]]
		for ln in alice:
			if ln[1]=='54':
				continue
			isomer = len(ln)>4
			istp,E = ln[2]+elements[ln[1]],float(ln[0])
			if isomer and istp not in ['135CE','139CE']:
				self.db.execute('INSERT INTO model_xs VALUES(?,?,?,?)',(istp+'g',E,float(ln[4]),'alice'))
				self.db.execute('INSERT INTO model_xs VALUES(?,?,?,?)',(istp+'m',E,float(ln[5]),'alice'))
			else:
				self.db.execute('INSERT INTO model_xs VALUES(?,?,?,?)',(istp,E,float(ln[3]),'alice'))
		self.db_connection.commit()

if __name__ == "__main__":
	mn = manager()
	# mn.create_plots_dir()
	# mn.update_files()
	# mn.update_ensdf()
	# mn.run_empire()
	# mn.run_talys()
	# # mn.update_talys()
	# mn.update_exfor()
	# mn.update_xs_prediction()
	mn.update_alice()
