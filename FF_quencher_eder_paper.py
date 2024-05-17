# In here goes the line by line confirmation that the FF_quencher_eder function is workingimport numpy as np
import numpy as np
import pandas as pd
from ederMatrixGenerator_paper import ederMatrix

# The following function reads in an array of counts from a dukecevns spectrum, for 1 kg and 1 GWhr exposure
# as the detector number and the desired form factor model (klein-nystrand or helm)
# it also takes in desired mass and exposure
def FF_quencher_eder(dataIn,detectorNumber,FF_model,mass,exposure):
    
    # split the data in to energy and dukecevnsSpec
    energy = dataIn[:, 0]
    dukecevnsSpec = dataIn[:, 1]

    # ***** UNIT ADJUSTMENTS *****

    # Convert counts down to counts per keV
    dukecevnsSpec = dukecevnsSpec / 1000

    # Account for the 0.1 keV binning
    dukecevnsSpec = dukecevnsSpec / 10

    # Convert energy from MeV to keV
    energy = energy * 1000

    # round the energy to 1 decimal place
    energy = np.round(energy, 1)

    # ***** FORM FACTORS *****

    m = 72.63 * 931.5 # MeV/c^2
    def TtoQ(T):
        #return np.sqrt(2*(74 * 931.5)*T)
        return np.sqrt((2*(72.63 * 931.5) * T) + (T**2))
    
    def FF_KleinNystrand(Q):
        Q = Q / 197.3  # to convert to fm^-1
        Rn = 1.2 * (72.63)**(1/3)
        a = 0.7
        return (3 / (Q**3 * Rn**3)) * (np.sin(Q * Rn) - (Q * Rn * np.cos(Q * Rn))) * (1 / (1 + (Q**2 * a**2)))
    
    # def FF_Helm(Q):
    #     Q = Q / 197.3 # to convert to fm^-1
    #     R = 1.2 * (42 + 32)**(1/3)
    #     sval = 0.9
    #     return (3 / (Q * R)) * ((np.sin(Q * R)/(Q**2 * R**2)) - (np.cos(Q * R)/(Q * R))) * np.exp((-1/2) * Q**2 * sval**2)

    # Generate an array of Erecoil values from 0.1 to 85.6 keV in steps of 0.1 keV
    Erecoil_arr = np.arange(0.1, 85.6, 0.1)
    # Convert Erecoil_arr to MeV
    Erecoil_arr = Erecoil_arr / 1000
    # Round Erecoil_arr to 5 decimal places
    Erecoil_arr = np.round(Erecoil_arr, 5)
    # Calculate the corresponding Q values
    Q_arr = TtoQ(Erecoil_arr)

    # Determine the FF model to use
    if FF_model == 'KN':
        FF_arr = FF_KleinNystrand(Q_arr)
        # print length of FF_arr
        print("Length of FF_arr: ", len(FF_arr))
    elif FF_model == 'Helm':
        FF_arr = FF_Helm(Q_arr)

    # ***** END OF FORM FACTORS *****



    # ***** QUENCHING MATRIX *****

    # Import the quenchingMatrix.npy file as a numpy array
    quenchingMatrix = np.load('quenchingMatrix.npy')

    # ***** END OF QUENCHING MATRIX *****



    # ***** EDER MATRIX *****

    # Define the EDER matrix based off the detector number

    ederSmearingMatrix = ederMatrix(detectorNumber)

    # ***** END OF EDER MATRIX *****

   

    # ***** TIMING EFFICIENCY *****
    # Import the timing efficiency for the detector
    timingEfficiency = np.load('timing_efficiencies/eff_interp_Ge' + str(detectorNumber) + '.npy')



    # ***** DUKE CEVNS SPECTRUM *****
    # print length of dukecevnsSpec
    print("Length of dukecevnsSpec: ", len(dukecevnsSpec))
    countsFF = dukecevnsSpec * FF_arr * FF_arr
    # print length of countsFF
    print("Length of countsFF: ", len(countsFF))
    # print length of quenching matrix
    print("Shape of quenchingMatrix: ", quenchingMatrix.shape)
    countsFFquenched = np.matmul(quenchingMatrix, countsFF)
    countsFFquenchedSmeared = np.matmul(ederSmearingMatrix, countsFFquenched)
    countsFFquenchedSmearedTiming = countsFFquenchedSmeared * timingEfficiency

    # scaling the counts to desired mass and exposure
    countsFFquenchedSmearedTiming = countsFFquenchedSmearedTiming * mass * exposure

    return countsFFquenchedSmeared, energy