

from epitopes.download import fetch_data
import pandas as pd


class LazyLoader(object):
	def __init__(self, 
			filename, 
			header, 
			server = "ftp://ftp.ensembl.org", 
			dir_template = "/pub/release-%(RELEASE)s/mysql/%(SPECIES)s_core_%(RELEASE)s_37/"):
		self.filename  = filename 
		self.header = header
		self.server = server
		self.dir_template = dir_template
		
		# maps each species and release  to a dataframe
		self._df_cache = {}
		self._path_cache = {}

	def _normalize_species(self, species):
		return species.lower().replace(" ", "_")
	def get_path(self, release = "75", species = "homo sapiens"):
		species_no_spaces = self._normalize_species(species)
		# convert to string in case given int
		release = str(release)
		key = (species_no_spaces, release)
		if key in self._path_cache:
			return self._path_cache[key]
		
		server_dir = self.dir_template % {"RELEASE" : release, "SPECIES" : species_no_spaces}
		url = self.server + server_dir + self.filename 
		local_filename = self.filename.replace(".gz", "")
		local_path = fetch_data(local_filename, url, "ensembl")
		
		self._path_cache[key] = local_path
		return local_path

	def get_dataframe(self, release = "75", species = "homo sapiens"):
		# conver to string in case given int
		relase = str(release)
		species_no_spaces = self._normalize_species(species)
		key = (species_no_spaces, release)
		if key in self._df_cache:
			return self._df_cache[key]
		
		local_path = self.get_path(release, species)
		df = pd.read_csv(local_path, sep='\t')
		
		self._df_cache[key] = df
		return df

	def get_db(self, release="75", species="homo sapiens"):
		df = self.get_dataframe(release, species)
		assert False, "Conversion to sqlite3 not implemented"


