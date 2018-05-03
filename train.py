import numpy as np
from propy.PyPro import GetProDes
from sklearn import linear_model
import argparse

def main():
    """train.py uses a training set, a file containing a peptide sequence and a number, and generates a predictor using
    ridge regression. The peptide sequence is converted to pseudo aminoacid composition to create a vector of
    samples using the propy module.
    The script needs two mandatory arguments:
    -i or --inputfile, which is the name of the file contaning the training set
    -o or --outputfile, which is the name of the file where the parameters of the estimator will be saved
    """
    #Get file name with training data
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--inputfile",required=True,help="input file with two columns: tetra peptide sequences and number")
    parser.add_argument("-o","--outputfile",required=True,help="output file containing ridge regression parameters")
    arguments = vars(parser.parse_args())

    #Read data and store in lists
    datafromfile = np.loadtxt(arguments["inputfile"],dtype='S4,f4')

    sequences = []
    numbers = []
    for i in range (len(datafromfile)):
        sequences.append(datafromfile[i][0])
        numbers.append(datafromfile[i][1])

    #Convert peptide sequence into pseudo amino acid composition (paac) vector and store in a list
    paaclist = []
    for i in range(len(sequences)):
        paac = GetProDes(sequences[i]).GetPAAC(lamda=3,weight=0.05)
        vector = [value for value in paac.itervalues()]
        paaclist.append(vector)

    #Convert lists into numpy arrays
    X = np.array(paaclist)
    y = np.array(numbers)

    #Make a ridge regression with default values (alpha=1, solver='auto', tol=0.001)
    rreg = linear_model.Ridge()
    rreg.fit(X,y)

    #print results of regression
    print "coefficients of regression"
    print rreg.coef_
    print "intercept"
    print rreg.intercept_
    print "R2"
    print rreg.score(X,y)

    #Save parameters of estimator in output file
    estimator = rreg.intercept_
    estimator = np.append(estimator,rreg.coef_)
    np.savetxt(arguments["outputfile"],estimator)
    print "estimator generated, use predict.py next"

if __name__ == "__main__":
    main()

