import numpy as np
import pandas as pd
import scipy.stats as stats

# Function for energy dependent energy resolution

def eder(x):
    if detectorNum == 21:
        a = 0.132
        b = 0.036
        c = 0.003
    elif detectorNum == 23:
        a = 0.114
        b = 0.040
        c = 0.002
    elif detectorNum == 25:
        a = 0.157
        b = 0.035
        c = 0.002
    elif detectorNum == 26:
        a = 0.14
        b = 0.037
        c = 0.001
    elif detectorNum == 28:
        a = 0.126
        b = 0.034
        c = 0.003
    return ((a + b * np.sqrt(x + c * x**2)))/(2.355)

def ederMatrix(detectorNumber):
    # define a global variable "detectorNum" to be used in the eder function
    global detectorNum
    detectorNum = detectorNumber
    # Generate a dataframe df_eder with columns 'Erecoil' and 'sigma'
    df_eder = pd.DataFrame(columns=['Erecoil', 'sigma'])
    # populate the Erecoil column from 0.05 to 85.55 keV in steps of 0.1 keV
    df_eder['Erecoil'] = np.arange(0.1, 85.6, 0.1)
    # Round the Erecoil values to 2 decimal places
    df_eder['Erecoil'] = df_eder['Erecoil'].round(2)
    # Calculate the sigma value for each Erecoil value
    df_eder['sigma'] = df_eder['Erecoil'].apply(eder)

    pdf_dictEDER = {}  # To store the gaussian columns per Erecoil

    for Erecoil in df_eder['Erecoil']:
        mean = Erecoil
        sigma = df_eder.loc[df_eder['Erecoil'] == Erecoil, 'sigma'].values[0]
        y = stats.norm.pdf(df_eder['Erecoil'], mean, sigma)
        y = y / sum(y)
        pdf_dictEDER[Erecoil] = y

    df_ederMatrix = pd.DataFrame(pdf_dictEDER) # generating the dataframe from the dictionary
    EDER_Matrix = df_ederMatrix.values # converting the dataframe to a numpy array

    return EDER_Matrix