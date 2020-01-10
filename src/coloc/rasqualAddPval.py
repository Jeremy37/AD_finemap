import gzip
import argparse
from scipy import stats

def main():
	parser = argparse.ArgumentParser(description = "Add column with p value, calculated from the chisq column of the RASQUAL output.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("--rasqualfile", type=file, help = "Path to the merged RASQUAL output file")
	args = parser.parse_args()

	rasqual_file = openg(args.rasqualfile)
	for line in rasqual_file:
		line = line.rstrip()
		fields = line.split("\t")
		chi_stat = float(fields[10])
		p_value = stats.chi2.sf(chi_stat, 1)
		print("\t".join(fields + [str(p_value)]))


def isgzfile(fname):
	namelen = len(fname)
	return fname[namelen-3:namelen] == ".gz"

def openg(f, mode=None):
	if type(f) is file:
		if isgzfile(f.name):
			return gzip.GzipFile(fileobj=f)
		return f
	elif type(f) is str:
		if isgzfile(f):
			if mode is None:
				mode = 'rb'
			return gzip.open(f, mode)
		if mode is None:
			mode = 'r'
		return open(f, mode)
	else:
		die("openg: unrecognized type for parameter 'f': {0}".format(type(f)))
 

main()
