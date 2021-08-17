PDBCUTOFF = 33.0
DSSPCUTOFF = 0.55
TEST=False

three2oneAA={
	"ALA": "A",
	"ARG": "R",
	"ASN": "N",
	"ASP": "D",
	"CYS": "C",
	"GLN": "Q",
	"GLU": "E",
	"GLY": "G",
	"HIS": "H",
	"ILE": "I",
	"LEU": "L",
	"LYS": "K",
	"MET": "M",
	"PHE": "F",
	"PRO": "P",
	"SER": "S",
	"THR": "T",
	"TRP": "W",
	"TYR": "Y",
	"VAL": "V"
}

"""
Maximum Allowed Solvent Accessibilites of Residues in Proteins

    Matthew Z. Tien,
    Austin G. Meyer,
    Dariya K. Sydykova,
    Stephanie J. Spielman,
    Claus O. Wilke

2013, PLOS ONE

https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0080635
"""
RSA_norm={
	"A": 129,
	"R": 274,
	"N": 195,
	"D": 193,
	"C": 167,
	"Q": 225,
	"E": 223,
	"G": 104,
	"H": 224,
	"I": 197,
	"L": 201,
	"K": 236,
	"M": 224,
	"F": 240,
	"P": 159,
	"S": 155,
	"T": 172,
	"W": 285,
	"Y": 263,
	"V": 174
}