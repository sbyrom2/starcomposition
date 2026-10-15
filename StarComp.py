import numpy as np

# atomic masses, accounting for isotopes
# using isotope ratios from Lodders 2021
L21_solar_atomic_masses = np.array([1.007844856, 4.002439456, 6.940046604, 9.01218291,10.81175541,12.01038885,14.00532465,15.99904151,
						18.99840317, 20.12763505, 22.98976928, 24.30498176, 26.98153859, 28.08538367, 30.973762, 32.06387933, 35.45284372,
						36.27537813, 39.09946879, 40.07802253, 44.9559119, 47.8668887, 50.94147049, 51.9961191, 54.93804512,
						55.84514458, 58.93319506, 58.69335196, 63.54556016, 65.37776549, 69.72306808, 72.62958875, 74.9215965,
						78.97168088, 79.90360789, 83.79283539, 85.4941539, 87.61750433, 88.9058483, 91.22353759, 92.9063781, 95.95978968,
						101.0650881, 102.9055043, 106.4153292, 107.8681467, 112.4121436, 114.8182664, 118.71035, 121.7597261,
						127.5855887, 126.9044728, 131.1826916, 132.905452, 137.3269159, 138.905438, 140.1157106, 140.9076525,
						144.2446602, 150.3632757, 151.9643705, 157.2520476, 158.9253468, 162.4997627, 162.4997627, 167.2590764,
						168.9342133, 173.0544616, 174.9690303, 178.4865658, 180.9478756, 183.8416979, 186.2415239, 190.2476945, 192.2166171,
						195.0839451, 196.9665687, 200.5923999, 204.3833324, 207.3188698, 208.9803987, 232.0380553, 237.3199066]) 
L21_elems = ['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg',
		'Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr',
		'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br',
		'Kr','Rb','Sr','Y','Zr','Nb','Mo','Ru','Rh','Pd',
		'Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La',
		'Ce','Pr','Nd','Sm','Eu','Gd','Tb','Dy','Ho','Er',
		'Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Au','Pt',
		'Hg','Tl','Pb','Bi','Th','U']

# in the format A(elem) = log10(N_elem/N_H) + 12
GAS07 = np.array([
	 12.00,  10.93,  1.050,  1.380,  2.700,  8.390,  7.780,  8.660,
	 4.560,  7.840,  6.170,  7.530,  6.370,  7.510,  5.360,  7.140,
	 5.500,  6.180,  5.080,  6.310,  3.170,  4.900,  4.000,  5.640,
	 5.390,  7.450,  4.920,  6.230,  4.210,  4.600,  2.880,  3.580,
	 2.290,  3.330,  2.560,  3.250,  2.600,  2.920,  2.210,  2.580,
	 1.420,  1.920, -99.00,  1.840,  1.120,  1.660,  0.940,  1.770,
	 1.600,  2.000,  1.000,  2.190,  1.510,  2.240,  1.070,  2.170,
	 1.130,  1.700,  0.580,  1.450, -99.00,  1.000,  0.520,  1.110,
	 0.280,  1.140,  0.510,  0.930,  0.000,  1.080,  0.060,  0.880,
	-0.170,  1.110,  0.230,  1.250,  1.380,  1.640,  1.010,  1.130,
	 0.900,  2.000,  0.650, -99.00, -99.00, -99.00, -99.00, -99.00,
	-99.00,  0.060, -99.00, -0.520])

# L21 and GAS07 have the same elements they track
# so remove all values from 
solar = GAS07[GAS07 > -10]



# define a class
class StarComp:
	''' 
	A class to represent a "star" of a certain composition
	Class methods allow users to modify X,Y,Z, [Fe/H], [C/Fe], etc
	And the rest of the abundances will scale to maintain 
	a certain solar mixture.
	When you instantiate a StarComp object, it is created with GAS07 solar mixture
	and can be modified from there. So far there is not an option to switch
	to different solar mixtures.
	
    Attributes
    ----------
	Ntot : float 
        Number of atoms in the "star" (default 1e12)
    N : array(float), len(83)
        Number of atoms of each element that is tracked
    idx_dict : dict(str,int)
        mapping of element abbreviation to index within N
    m : array(float)
        the average atomic mass for each element (default is from L21)
        May add ability to change isotope fractions, so this might not be constant
    mij : Uneven np.array(float)
        atomic mass of each isotope j for element i
	N_iso : uneven np.array(float)
        Number of each isotope. Same shape as mij.
	Nsun : array(float), len(83)
        Number of atoms of each element for GAS07 solar mix

    Methods
    -------
    getXYZ(): 
        returns the mass fraction of H, He, and metals
	get_Riso(elem_abbr,iso1,iso2):
        returns mass ratio of isotopes 1 and 2 of element
        e.g. get_Riso('C',12,13) returns mass(C12)/mass(C13)
	FRAC_C (N,O):
        returns the mass fraction of metals contained in C (N,O)

	changeY_fixedZ(Ynew):
		change compositio to match an input Y value, keeping Z fixed
		Other functions are variations on this
	
	There are more methods. Most are of the form: get/set_[value] 
	Read further if you wish to use them.

   '''
	
	# m is the atomic mass for each element
	m = L21_solar_atomic_masses

	# set the solar counts for each element
	Ntot = 1e12
	N_NH = 10**(solar-12)
	N_Ntot = N_NH / np.sum(N_NH)
	# track number of each element
	Nsun = N_Ntot * Ntot

	# list of idxs/elems that have tracked isotopes
	# right now this is H,He,Li,Be,B,C,N,O
	iso_elems = np.arange(8)

	# list of atomic numbers for isotopes
	atomic_number_iso = [[1,2],[3,4],[6,7],[9],[10,11],[12,13],[14,15],[16,17,18]]
	mij = [[1.007825032,2.014101778], # H1, H2
	                [3.01602932,4.002603254], # He3, He4
					[6.015122885,7.016003428], # Li6, Li7
					[9.01218291], # Be9
					[10.01293696,11.00930537], # B10, B11
					[12,13.00335484], # C12, C13
					[14.003074,15.0001089], # N14, N15
					[15.99491462,16.99913176,17.99915961]] # O16, O17, O18
	
	# mapping of element abbreviation to atomic mass
	elem_dict = dict(zip(L21_elems,L21_solar_atomic_masses))
	# mapping of element abbreviation to index
	idx_dict = dict(zip(L21_elems,np.arange(83)))


	def __init__(self,Ntot=1e12): # default starts at GAS07 solar composition
		# track number of each element
		self.N = self.Nsun * Ntot / self.Ntot
		# reset Ntot to the input value
		self.Ntot = Ntot 

		# now track isotopes
		# start with N_iso/N_elem:
		Niso_Nelem = np.array([
			np.array([.999980, .0000197]), # H1, H2
	   	 	np.array([.000166, .999834 ]), # He3, He4
			np.array([.07589, .92411 ]), # Li6, Li7
			np.array([1]), # Be9
			np.array([.19827,.80173]), # B10, B11
			np.array([.98965,.01035 ]), # C12, C13
			np.array([.99774, .00226]), # N14, N15
			np.array([.99777,.00035,.00188])],dtype=object) # O16, O17, O18
		# multiply by the count for each element
		# to get the count for each isotope
		self.N_iso = self.N[:8] * Niso_Nelem

	def X_i(self,i):
		# returns mass fraction of element I
		return self.M_i(i) / self.Mtot()

	def M(self,elem_abbr):
		idx = self.idx_dict[elem_abbr]
		return self.M_i(idx)
	
	def MZ(self):
		# return mass in metals of the StarComp object
		return np.sum(self.m[2:] * self.N[2:])
	
	def N_Z(self):
		# return atom count of metals
		return np.sum(self.N[2:])
	
	def muZ(self):
		# return average atomic mass of metals
		return self.MZ() / self.N_Z()
	
	def Mtot(self):
		# return total mass of the StarComp object
		return np.sum(self.m * self.N)
	
	def M_i(self,i):
		# i is the element's atomic number
		# return the mass of that element
		return self.m[i] * self.N[i]

	def M_ij(self,i,j):
		return self.mij[i][j] * self.N_iso[i][j]
	
	def getX(self):
		# return mass fraction of hydrogen
		return self.M_i(0) / self.Mtot()

	def getY(self):
		# return mass fraction of helium
		return self.M_i(1) / self.Mtot()

	def getZ(self):
		# return mass fraction of metals
		return self.MZ() / self.Mtot()

	def getXYZ(self):
		return (self.getX(),self.getY(),self.getZ())

	def getFeH(self):
		# returns [Fe/H]
		return np.log10(self.N[26] / self.N[0]) - np.log10(self.Nsun[26]/self.Nsun[0])

	def get_bracket_elemFe(self,elem_abbr):
		# returns [elem/Fe]
		idx = self.idx_dict[elem_abbr]
		return np.log10(self.N[idx]/self.N[26]) -  np.log10(self.Nsun[idx]/self.Nsun[26])

	def FRAC_C(self):
		# carbon mass fraction of metals
		return self.M_i(5) / self.MZ()

	def FRAC_N(self):
		# nitrogen mass fraction of metals
		return self.M_i(6) / self.MZ()

	def FRAC_O(self):
		# oxygen mass fraction of metals
		return self.M_i(7) / self.MZ()

	def get_A(self,elem_abbr):
		# return A = log10(Nelem/NH) + 12
		idx = self.idx_dict[elem_abbr]
		return np.log10(self.N[idx]/self.N[0]) + 12

	def get_Riso(self,elem_abbr,iso1,iso2):
		# return mass ratio of isotopes e.g. C12/C13
		idx = self.idx_dict[elem_abbr]
		# TODO: add try, except for incorrect isotopes
		iso1 = self.atomic_number_iso[idx].index(iso1)
		iso2 = self.atomic_number_iso[idx].index(iso2)

		return self.M_ij(idx,iso1) / self.M_ij(idx,iso2)

	def changeY_fixedZ(self,Ynew):
		# change compositio to match an input Y value
		# keeping Z fixed
		Xnew = 1-self.getZ()-Ynew
		Mtot = self.Mtot()
		self.N[0] = Xnew * Mtot / self.m[0]
		self.N[1] = Ynew * Mtot / self.m[1]
		self.rescale_to_Ntot()

	def changeX_fixedZ(self,Xnew):
		# change compositio to match an input X value
		# keeping Z fixed
		Ynew = 1-self.getZ()-Xnew
		Mtot = self.Mtot()
		self.N[0] = Xnew * self.Mtot() / self.m[0]
		self.N[1] = Ynew * self.Mtot() / self.m[1]
		self.rescale_to_Ntot()

	def changeY_fixedX(self,Ynew):
		# change compositio to match an input Y value
		# keeping X fixed
		Znew = 1-self.getX()-Ynew
		Mtot = self.Mtot()
		self.N[1] = Ynew * Mtot / self.m[1]
		NZ_new = Znew * Mtot / self.muZ()
		self.N[2:] = self.N[2:] * NZ_new / self.N_Z()

		self.rescale_to_Ntot()

	def changeX_fixedY(self,Xnew):
		# change compositio to match an input X value
		# keeping Y fixed
		Znew = 1-self.getY()-Xnew
		Mtot = self.Mtot()
		self.N[0] = Xnew * Mtot / self.m[0]
		NZ_new = Znew * Mtot / self.muZ()
		self.N[2:] = self.N[2:] * NZ_new / self.N_Z()
		self.rescale_to_Ntot()

	def changeZ_fixedY(self,Znew):
		# change compositio to match an input Z value
		# keeping Y fixed
		Xnew = 1-self.getY()-Znew
		Mtot = self.Mtot()
		self.N[0] = Xnew * Mtot / self.m[0]
		NZ_new = Znew * Mtot / self.muZ()
		self.N[2:] = self.N[2:] * NZ_new / self.N_Z()
		self.rescale_to_Ntot()

	def changeZ_fixedX(self,Znew):
		# change compositio to match an input Z value
		# keeping X fixed
		Ynew = 1-self.getX()-Znew
		Mtot = self.Mtot()
		self.N[1] = Ynew * Mtot / self.m[1]
		NZ_new = Znew * Mtot / self.muZ()
		self.N[2:] = self.N[2:] * NZ_new / self.N_Z()
		self.rescale_to_Ntot()

	def set_XYZ(self,X,Y,Z):
		if abs(X+Y+Z - 1) > 1e-9:
			raise Exception(f'input values {X}, {Y}, {Z} must sum to 1\n These values sum to ({X+Y+Z})')
		self.changeX_fixedY(X)
		self.changeY_fixedX(Y)
		return self.getXYZ()

	def rescale_to_Ntot(self):
		Ntot_new = np.sum(self.N)
		self.N = self.N * (Ntot_new / self.Ntot)
		# also rescale isotopes
		self.N_iso = self.N_iso * (Ntot_new / self.Ntot)
		return None
	
	def set_Ntot(self,Ntot,i_know_what_im_doing=False):
		if Ntot < 1e6 and not i_know_what_im_doing:
			raise Exception("Don't make Ntot small or you risk precision errors")
		self.N = self.N * Ntot/self.Ntot 
		self.Ntot =  Ntot
		return 0

	def set_bracket_elemFe(self,elem_abbr,A):
		# set [elem/Fe]
		idx = self.idx_dict[elem_abbr]
		self.N[idx] = self.N[26] * (10**A) * self.Nsun[idx] / self.Nsun[26]
		self.rescale_to_Ntot()
		return self.get_bracket_elemFe(elem_abbr=elem_abbr)

	def set_alphaFe(self,A):
		# change all alpha elements: [alpha/Fe]
		alpha_elems = ['O','Ne','Mg','Si','S','Ar','Ca','Ti']
		for alpha in alpha_elems:
			self.set_bracket_elemFe(alpha,A)

	def set_FeH(self,A):
		# scale metals such that [Fe/H] = A
		# and [metal/Fe] remains unchanged 
		self.N[2:] = self.N[0] * (10**A) * self.Nsun[2:] / self.Nsun[0]
		self.rescale_to_Ntot()

	def get_iso(self,elem_abbr,iso):
		# return mass fraction of a given isotope
		idx = self.idx_dict[elem_abbr]
		isoidx = self.atomic_number_iso[idx].index(iso)
		Niso = self.N_iso[idx][isoidx]
		return Niso * self.mij[idx][isoidx] / self.Mtot()
	
	# TODO later
	def rescale_avg_atomic_mass(self):
		# called after changing isotope ratios
		# I'm not doing this yet, so I don't need to implement it
		return 0 

s = StarComp()
