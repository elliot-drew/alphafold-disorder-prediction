# alphafold-disorder-prediction

This is a simple implementation of the disorder prediction approach described [here](https://github.com/normandavey/ProcessedAlphafold) by Norman Davey [github](https://github.com/normandavey) and Balint Meszaros [github](https://github.com/BalintMeszaros).

Both the pLDDT method and the RSA method to predict disorder are implemented. The suggested window widths and cutoffs are used. It is possible to supply your own cutoffs as optional arguments.

An alphafold model in PDB format and/or a DSSP file produced using the model are required to perform the predictions. You can supply one or both using the arguments below.

## CAUTION

This program is fairly thick so will happily predict disorder with information that isn't pLDDT values - e.g. actual crystallographic temperature factors from an experimental PDB file. So only use Alphafold derived models.

Also! I make no claims as to the accuracy of predictions made using this method (or even whether this code gives the correct answers). Test it thoroughly yourself before using in your work/research.

## Usage:

You need python 3.X and numpy installed. It would probably be fairly easy to change it to work with python 2.x by modifying the print statements in the files.

```
usage: make_prediction.py [-h] [--pdb PDB] [--dssp DSSP]
                          [--pdbcutoff PDBCUTOFF] [--dsspcutoff DSSPCUTOFF]

optional arguments:
  -h, --help            show this help message and exit
  --pdb PDB             path to AlphaFold PDB file
  --dssp DSSP           path to DSSP file produced from AlphaFold PDB file
  --pdbcutoff PDBCUTOFF
                        lPDDT method cutoff value (0.0 - 100.0) default = 33
  --dsspcutoff DSSPCUTOFF
                        RSA method cutoff value (0.0 - 1.0) default = 0.55
```

The program prints to Stdout, so you could do the following in a terminal to write the results to a file called output.txt.

```
python make_prediction.py --pdb XXXX.pdb > output.txt
```

Batch predictions can also be made easily using some simple scripting. I may add that soon.


