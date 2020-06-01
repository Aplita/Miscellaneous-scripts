################################################################################
# Transitions and transversions           #
###########################################
# Practice - script to find transitions and transversions.
#
# Written by Ana Paula Vargas.
# Bioinformatics lab - UPCH. Lima, Peru.
# Last updated: 10-sep-20.
###############################################################################

import xlrd
import xlwt

# Load Excel:
file_location = "/home/vargas/Desktop/HRM/ALL_FILTERED.xlsx"
workbook = xlrd.open_workbook(file_location)
sheet = workbook.sheet_by_index(0)

# Define inputs:
# RefBases -> column 5
# MutBases -> column 6

RefBases = []
MutBases = []

for h in range(sheet.nrows):
	RefBases.append(sheet.cell_value(h, 5))
	MutBases.append(sheet.cell_value(h, 6))


A = len(RefBases)

# Append TRANSITIONS to X
X = ["Transitions"]
for i in range (1, A):
	if (RefBases[i] == "A"):
		if (MutBases[i] == "G"):
			X.append(1)
		else:
			X.append(0)
	elif (RefBases[i] == "G"):
		if (MutBases[i] == "A"):
			X.append(1)
		else:
			X.append(0)
	elif (RefBases[i] == "C"):
		if (MutBases[i] == "T"):
			X.append(1)
		else:
			X.append(0)
	elif (MutBases[i] == "C"):
		X.append(1)
	else:
		X.append(0)


# Append TRANSVERSIONS to Y
# Important:
#	A-T & C-G
Y = ["Total transversions"]
AT = ["A <-> T transversions"]
CG = ["C <-> G transversions"]
for j in range (1, A):
	if (RefBases[j] == "A"):
		CG.append(0)
		if (MutBases[j] == "T"):
			Y.append(1)
			AT.append(1)
		elif (MutBases[j] == "C"):
			Y.append(1)
			AT.append(0)
		else:
			Y.append(0)
			AT.append(0)
	elif (RefBases[j] == "C"):
		AT.append(0)
		if (MutBases[j] == "G"):
			Y.append(1)
			CG.append(1)
		elif (MutBases[j] == "A"):
			Y.append(1)
			CG.append(0)
		else:
			Y.append(0)
			CG.append(0)
	elif (RefBases[j] == "G"):
		AT.append(0)
		if (MutBases[j] == "C"):
			Y.append(1)
			CG.append(1)
		elif (MutBases[j] == "T"):
			Y.append(1)
			CG.append(0)
		else:
			Y.append(0)
			CG.append(0)
	elif (MutBases[j] == "A"):
		Y.append(1)
		AT.append(1)
		CG.append(0)
	elif (MutBases [j] == "G"):
		Y.append(1)
		AT.append(0)
		CG.append(0)
	else:
		Y.append(0)
		AT.append(0)
		CG.append(0)

# Export data to Excel SpreadSheet
NewBook = xlwt.Workbook(encoding="utf-8")
s1 = NewBook.add_sheet("Sheet 1")


# 	Row, Column
for k in range(0, len(RefBases)):
	s1.write(k, 0, RefBases[k])
	s1.write(k, 1, MutBases[k])
	s1.write(k, 2, X[k])
	s1.write(k, 3, Y[k])
	s1.write(k, 4, AT[k])
	s1.write(k, 5, CG[k])

NewBook.save("470-genome_muts.xlsx")
