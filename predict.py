import numpy as np
from propy.PyPro import GetProDes
import argparse

def EstimatePAACNumber(estimator,paacvector):
    """applies estimator to a paacvector"""

    number = estimator[0]
    for i in range(len(paacvector)):
        number += estimator[i+1]*paacvector[i]

    return number

def main():
    """predict.py must be used after train.py, and uses the trained model to make predictions.
    The script needs three mandatory arguments:
    -e or --estimator, which is the name of the file containing the parameters from the ridge regression obtained from train.py
    -i or --inputfile, which is the name of the file with a sequence of new peptides
    -o or --outputfile, which is the name of the file where the numbers corresponding to the sequences from -i will be stored
    """

    #Get file name with estimator parameters obtained using train.py
    parser = argparse.ArgumentParser()
    parser.add_argument("-e","--estimator",required=True,help="estimator obtained as output file from train.py")
    parser.add_argument("-i","--inputfile",required=True,help="new peptide sequences which numbers are unknown")
    parser.add_argument("-o","--outputfile",required=True,help="number corresponding to sequences obtained using estimator")
    arguments = vars(parser.parse_args())
    
    #Read estimator
    estimatorparameters = np.loadtxt(arguments["estimator"])

    #Read new sequences
    datafromfile = np.loadtxt(arguments["inputfile"],dtype='S4')
    
    sequences = []
    if datafromfile.size == 1:
        sequences.append(str(datafromfile))
    else:
        for i in range(datafromfile.size):
            sequences.append(datafromfile[i])
    
    #Convert new peptide sequence into pseudo amino acid composition (paac) vector and store in list
    paaclist = []
    for i in range(len(sequences)):
        paac = GetProDes(sequences[i]).GetPAAC(lamda=3,weight=0.05)
        vector = [value for value in paac.itervalues()]
        paaclist.append(vector)
    
    X = np.array(paaclist)
    
    #Do predictions
    pepnumber = np.zeros(len(X))
    for i in range(len(X)):
        pepnumber[i] = EstimatePAACNumber(estimatorparameters,X[i])

    #Save estimated numbers in output file
    np.savetxt(arguments["outputfile"],pepnumber)

if __name__ == "__main__":
    main()
