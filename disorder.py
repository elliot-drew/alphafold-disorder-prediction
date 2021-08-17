from config import *
from dssp import DSSPData

import numpy as np 
import os

"""

Predict disorder from AF models.

Two different methods:

1. from pLDDT in PDB files
    - window 15, cutoff value 0.33
2. from ACC in DSSP
    - normalised to RSA using max poss. acc.
    - window 15, cutoff value 0.55

Both methods from https://github.com/normandavey/ProcessedAlphafold


"""

def getWindow(vec, width=15):
    """
    Return a list as long as input where each value is average
    of window of width centred on ith input value.
    """    
    offset=(width-1)/2

    window_av=[]
    for i in range(0, len(vec)):
        first=int(i-offset)
        last=int(i+offset+1)
        if(first<0):
            first=0
        if(last>(len(vec)-1)):
            last=len(vec)-1
        window=vec[first:last]
        window_av.append(np.mean(window))
    
    return(window_av)


def getRSA(ddob):
    """
    Get ACC value from DSSP file, then normalise using config.RSA_norm values
    """
    seq = ddob.getAAs()
    resnums = ddob.getResnums()
    sasa = ddob.getACC()
    ss = ddob.getSecStruc()
    
    ss_ind={
        "H":0,
        "G":1,
        "I":2,
        "E":3,
        "B":4,
        "T":5,
        "S":6,
        "O":7
        }
    ss_counts=[]
    ss_mod=[]
    rsa=[]
    for aa, acc, s in zip(seq, sasa, ss):
        norm=RSA_norm[aa]
        rsa.append(float(acc)/norm)

        ss_count_tmp=[0]*8
        if s=="":
            ss_mod.append("O")
            ss_count_tmp[7]=1
        else:
            ss_mod.append(s)
            ss_count_tmp[ss_ind[s]]=1
        
        ss_counts.append(ss_count_tmp)
    
    ss_av_counts=np.mean(ss_counts, axis=0)

    return(rsa, resnums, seq, ss_mod, ss_av_counts)

def getpLDDT(pdb):
    """
    From alphafold PDB file get backbone atom pLDDT values then average for each residue.

    pLDDT is stored in B factor column - [60-66]

    Atom name is stored - [12-16]
    """
    pLDDT=[]
    pLDDT_tmp=[]
    resnums=[]
    resns=[]

    backbone=["C", "CA", "CB", "N", "O"]
    curr_resnum=1

    for l in pdb.readlines():
        if(l.startswith("ATOM")):
            rn = l[17:20].strip()
            rnum = int(l[22:26])
            n = l[12:16].strip()
            b = float(l[60:66])
            
            if(n in backbone):
                if(rnum == curr_resnum):
                    pLDDT_tmp.append(b)
                else:
                    pLDDT.append(np.mean(pLDDT_tmp))
                    pLDDT_tmp=[b]
                    resnums.append(rnum)
                    resns.append(rn)
                    curr_resnum=rnum

    pLDDT.append(np.mean(pLDDT_tmp))
    
    return(pLDDT, resnums, resns)

def disorderPred(vec, cutoff, lessthan):
    pred=[]
    for v in vec:
        if lessthan:
            if(v<=cutoff):
                pred.append(1)
            else:
                pred.append(0)
        else:
            if(v>=cutoff):
                pred.append(1)
            else:
                pred.append(0)
    return(pred)


    

