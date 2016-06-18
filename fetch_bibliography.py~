# Python script for fetching bibtex entries from NASA ADS.
# To run on the command line, use:

#	python fetch_bibliography.py $Paper.tex

# where $Paper.tex contains a list of citations in the form

#	%citation{vkercj10}{2010ApJ...722L.157V}

# below \end{document} ("vkercj10" is a citation label and 
# "2010ApJ...722L.157V" is the ADS extension to the paper).
# Output is "AutoBibliography.bib" (which can be changed in
# the extension)

# If you wish to include additional citations not on ADS,
# use the -addition flag, and specify the file they're
# in.

import urllib, urllib2
import argparse

def list_split(big, n):
	big_split = []
	for i in range(0, len(big), n):
		big_split.append(big[i:i+n])
	return big_split

def get_ads_entries(wantstuff):
	"""For a list of ADS handles, returns (unsorted) list of bibtex entries from ADS website."""

	ads_handles = ""
	for item in wantstuff:
		ads_handles = ads_handles + item + "\n"

	# Doing this single query is much faster than doing multiple queries!
	ads_query = "http://adsabs.harvard.edu/cgi-bin/nph-abs_connect?db_key=ALL&warnings=YES&version=1&bibcode={0:s}&nr_to_return=$nr&start_nr=1&data_type=BIBTEX".format(urllib.quote_plus(ads_handles))

	response = urllib2.urlopen(ads_query)
	return (response.read()).split("@")[1:]	# Retrieve ADS entries (zeroth entry is just preamble)


parser = argparse.ArgumentParser(description = "Python script for fetching bibtex entries from NASA ADS.")
parser.add_argument('-filename', metavar='N1', type=str, help="name of TeX file", default="Thesis.tex")
parser.add_argument('-bibname', metavar='N2', type=str, help="name of BibTeX output file", default="AutoBibliography.bib")
parser.add_argument('-missed', help="Report handles not retrived by ADS", action='store_true', default=True)
parser.add_argument('-addition', metavar='N3', help="Include additional bibtex file of sources not on ADS", type=str)

args = parser.parse_args()

f = open(args.filename, 'r')
#f = open("Thesis.tex", 'r')

wantstuff = []
wantlabels = []

for line in f.readlines():
	if line.startswith("%citation{"):
		wantlabels.append(line.split("{")[1][:-1])
		wantstuff.append(line.split("{")[2].split("}")[0])

f.close()

# If we only wanted to return unique values (while preserving the order; otherwise we could just use set()), we'd do below

#seen = set()
#seen_add = seen.add
#[ x for x in seq if not (x in seen or seen_add(x))]

# See http://www.secnetix.de/olli/Python/list_comprehensions.hawk and 
# http://stackoverflow.com/questions/18212574/why-does-checking-a-variable-against-multiple-values-with-or-only-check-the-fi
# Note that seen_add(x) returns None, which is treated as Boolean false

# BUT we actually want to the intersection of unique stuff and labels, so we have to do something like this:

uniquestuff = []
uniquelabels = []
seenstuff = set()
seenlabels = set()
for i in range(len(wantstuff)):
	if wantstuff[i] not in seenstuff:
		seenstuff.add(wantstuff[i])
		stuff_trig = False
	else:
		stuff_trig = True
	if wantlabels[i] not in seenlabels:
		seenlabels.add(wantlabels[i])
		labels_trig = False
	else:
		labels_trig = True
	if (not stuff_trig) and (not labels_trig):
		uniquestuff.append(wantstuff[i])
		uniquelabels.append(wantlabels[i])
	elif stuff_trig and (not labels_trig):
		print "DUPLICATE ADS HANDLE FOUND:", wantstuff[i], "; SKIPPING!"
	elif stuff_trig and (not labels_trig):
		print "DUPLICATE LABEL FOUND:", wantlabels[i], "; SKIPPING!"
	else:
		print "DUPLICATE ENTRY FOUND:", wantlabels[i], wantstuff[i], "; SKIPPING!"

wantstuff = uniquestuff
wantlabels = uniquelabels

#wantstuff = ["2010ApJ...725..296F", "2015ApJ...788..BIG", "2014ApJ...788...75R", "2012ApJ...748...35S", "2012MNRAS.425.3024V", "2002A&A...381..923S"]
#wantlabels = ["frye+10", "BIG", "rask+14", "something", "shen+12", "spru02"]

print "fetch_bibliography found {0:d} unique entries".format(len(wantstuff))
#if len(wantstuff) > 200:
#	print "WARNING: UNIQUE CITATIONS EXCEED 200 ENTRIES!  NASA ADS MAY TRUNCATE RESULTS!"

wantstuff_split = list_split(wantstuff, 200)
html_split = []
for item in wantstuff_split:
	html_split.extend(get_ads_entries(item))

# There's discussion on how to do a better sort here: http://stackoverflow.com/questions/8251541/numpy-for-every-element-in-one-array-find-the-index-in-another-array

f = open(args.bibname, 'w')
#f = open("Autobibligraphy.bib", 'w')
foundstuff = []
for item in html_split:
	item_split = item.split("\n")											# Split bib entry by line
	foundstuff.append(item_split[0].split("{")[1][:-1])						# Append ADS handle of bib entry
	newname = wantlabels[wantstuff.index(foundstuff[-1])]					# Find corresponding label
	item_split[0] = "@" + item_split[0].split("{")[0] + "{" + newname + ","	# Reconstitute first line
	newitem = "\n".join(item_split)											# Reconstitute all lines
	f.write(newitem)														# Write out

if args.addition:
	fappd = open(args.addition, 'r')
	fappd_lines = fappd.readlines()
	fappd.close()

	f.write("\n")
	for item in fappd_lines:
		f.write(item)

f.close()

missingstuff = list(set(wantstuff) - set(foundstuff))
if len(missingstuff) > 0 and args.missed:
	print "WARNING: the following ADS handles were NOT found using an ADS search!"
	for item in missingstuff:
		print item
