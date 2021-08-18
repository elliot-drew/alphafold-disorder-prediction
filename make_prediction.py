from disorder import *
from dssp import DSSPData

import os
import argparse

def runPrediction(pdbFlag, pdbPath,  pdbCutoff, dsspFlag, dsspPath, dsspCutoff, errors):
    
    variables={}

    if pdbFlag:
        with open(pdbPath, "r") as PDBIN:
            confidence, resnums, resnames=getpLDDT(PDBIN)
        confWindow = getWindow(confidence, width=15)
        pdbPred = disorderPred(confWindow, pdbCutoff, True)
        variables["resnums"] = resnums
        variables["resnames"] = resnames
        variables["confWindow"] = confWindow
        variables["pdbPred"] = pdbPred
    

    if dsspFlag:
        with open(dsspPath, "r") as DSSPIN:
            ddob=DSSPData()
            ddob.parseDSSP(DSSPIN)
            rsa, resnums, resnames, ss, ss_av_counts=getRSA(ddob)
        rsaWindow = getWindow(rsa, width=15)
        dsspPred = disorderPred(rsaWindow, dsspCutoff, False)
        variables["resnums"] = resnums
        variables["resnames"] = resnames
        variables["rsa"] = rsa
        variables["ss"] = ss
        variables["dsspPred"] = dsspPred

        if pdbFlag:
            if len(variables["dsspPred"]) != len(variables["pdbPred"]):
                return(
                    "error",
                    "Number of residues in DSSP file ({0}) does not match number in PDB file ({1})".format(
                        len(variables["dsspPred"]),
                        len(variables["pdbPred"])
                        ) +
                    "\nProgram Exiting."
                )
    output=[]
    for i in range(0, len(variables["resnums"])):
        output_tmp=[]
        output_tmp.append(str(variables["resnums"][i]))
        output_tmp.append(str(variables["resnames"][i]))
        if pdbFlag and dsspFlag:
            if i==0:
                output.append("\t".join(["ResNum", "ResName","pLDDT_wind", "RSA_wind", "SS8", "pLDDT_pred", "RSA_pred"]))
            output_tmp.append(str(variables["confWindow"][i]))
            output_tmp.append(str(variables["rsa"][i]))
            output_tmp.append(str(variables["ss"][i]))
            output_tmp.append(str(variables["pdbPred"][i]))
            output_tmp.append(str(variables["dsspPred"][i]))
        elif pdbFlag:
            if i==0:
                output.append("\t".join(["ResNum", "ResName","pLDDT_wind", "pLDDT_pred"]))
            output_tmp.append(str(variables["confWindow"][i]))
            output_tmp.append(str(variables["pdbPred"][i]))
        elif dsspFlag:
            if i==0:
                output.append("\t".join(["ResNum", "ResName","RSA_wind", "SS8", "RSA_pred"]))
            output_tmp.append(str(variables["rsa"][i]))
            output_tmp.append(str(variables["ss"][i]))
            output_tmp.append(str(variables["dsspPred"][i]))
        
        else:
            if len(errors)>0:
                return(
                    "error",
                    "\n".join(errors) +
                    "\nProgram Exiting."
                    )
            else:
                return(
                    "error",
                    "Neither PDB nor DSSP file paths provided." +
                    "\nProgram Exiting."
                )
        
        output.append("\t".join(output_tmp))

    outputText="\n".join(output)
    if(pdbFlag):
        outputText+="\n"
        outputText+="\nAverage pLDDT_pred: {0}".format(str(round(np.mean(variables["pdbPred"])*1000)/1000))
        outputText+="\n"
    if(dsspFlag):
        outputText+="\nAverage 8 state DSSP SS: {0}".format("\t".join([str(e) for e in ss_av_counts]))
        outputText+="\nAverage RSA_pred: {0}".format(str(round(np.mean(variables["dsspPred"])*1000)/1000))
        outputText+="\n"
    
    return(
        "Success!",
        outputText
    )

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", help="path to AlphaFold PDB file",
                        type=str)
    parser.add_argument("--dssp", help="path to DSSP file produced from AlphaFold PDB file",
                        type=str)
    parser.add_argument("--pdbcutoff", help="lPDDT method cutoff value (0.0 - 100.0) default = 33",
                        type=float)
    parser.add_argument("--dsspcutoff", help="RSA method cutoff value (0.0 - 1.0) default = 0.55",
                        type=float)                                     
    args = parser.parse_args()

    pdbFlag=False
    dsspFlag=False

    pdbPath = ""
    dsspPath = ""

    errors=[]

    if args.pdbcutoff is not None:
        if(args.pdbcutoff>0.0 and args.pdbcutoff<=100.0):
            pdbCutoff=args.pdbcutoff
        else:
            errors.append("Error: pLDDT cutoff value supplied ({0}) is not between 0 and 100".format(args.pdbcutoff))
    else:
        pdbCutoff=PDBCUTOFF
    
    if args.dsspcutoff is not None:
        if(args.dsspcutoff>0.0 and args.dsspcutoff<=1):
            dsspCutoff=args.dsspcutoff
        else:
            errors.append("Error: RSA cutoff value supplied ({0}) is not between 0 and 1.".format(args.dsspcutoff))
    else:
        dsspCutoff=DSSPCUTOFF

    if args.pdb is not None:
        if os.path.isfile(args.pdb):
            pdbFlag = True
            pdbPath=args.pdb
        else:
            errors.append("Error: Path {0} for PDB file does not exist.".format(pdbPath))

    if args.dssp is not None:
        if os.path.isfile(args.dssp):
            dsspFlag = True
            dsspPath=args.dssp
        else:
            errors.append("Error: Path {0} for DSSP file does not exist.".format(dsspPath))
    
    if(len(errors)>0):
        print("\n".join(errors)+"\nExiting program.")
    else:
        log, text = runPrediction(pdbFlag, pdbPath, pdbCutoff, dsspFlag, dsspPath, dsspCutoff, errors)

        print(text)
    
