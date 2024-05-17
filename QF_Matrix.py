import matplotlib.pyplot as plt
import numpy as np
from quenchingFactors import QF_StandLind
import pandas as pd
import scipy.stats as stats
import math
from scipy.optimize import curve_fit
import seaborn as sns
import random
import h5py
import os
import csv


# The following few code blocks are using data from Long's thesis, specifically around the Quenching Smearing Factor from pages 92-93
# The error bars of the Erecoil values are not given, so those are actually taken off the plot and put into webplot digitizer to get the values

df_long = pd.DataFrame(columns=['Erecoil', 'Erecoil_eV','Eee','Eee_err', 'sigSq_nr','sigSq_nr_err', 'sigSq_noise','sigSq_noise_err','sigSq_E','sigSq_E_err','Erec_errMin','Erec_errMax'])

# Manually insert data on a per row basis
df_long.loc[0] = [1.75, 1750, 292, 21, 458, 66, 1038, 15, 3739, 216, 0.1321, 0.1234]
df_long.loc[1] = [1.97, 1970, 332, 20, 595, 86, 1050, 17, 4290, 231, 0.1361, 0.152]
df_long.loc[2] = [2.31, 2310, 403, 20, 882, 128, 1071, 21, 5282, 256, 0.1648, 0.1766]
df_long.loc[3] = [1.99, 1990, 350, 20, 1310, 190, 1692, 28, 6820, 439, 0.2049, 0.2086]
df_long.loc[4] = [4.20, 4200, 827, 28, 3422, 496, 1816, 70, 15372, 600, 0.2894, 0.2894]
df_long.loc[5] = [2.37, 2370, 425, 28, 1980, 287, 1715, 27, 8730, 579, 0.248, 0.249]
df_long.loc[6] = [4.89, 4890, 999, 35, 5402, 783, 1871, 69, 19552, 924, 0.3546, 0.3655]

# Now we want to reorganize the rows by increasing Erecoil
df_long.sort_values(by=['Erecoil'], inplace=True)

# Now define the function for the sigma "spread" of the quenching smearing factor
def QSF_sigma(Eee, sigSq_E, sigSq_noise, sigSq_nr):
    return (sigSq_E - sigSq_noise - sigSq_nr) / (Eee**2)

# Add a new column to df_long for the quenching smearing factor
df_long['QSF'] = QSF_sigma(df_long['Eee'], df_long['sigSq_E'], df_long['sigSq_noise'], df_long['sigSq_nr'])

# Since we're going to do a fit to get a functional form, we need to calculate and propagate the error of the quenching smearing factor
# Define a function to calculate the error for a single row
def error_QSF_sigma_row(row):
    Eee = row['Eee']
    sigSq_E = row['sigSq_E']
    sigSq_noise = row['sigSq_noise']
    sigSq_nr = row['sigSq_nr']
    Eee_err = row['Eee_err']
    sigSq_E_err = row['sigSq_E_err']
    sigSq_noise_err = row['sigSq_noise_err']
    sigSq_nr_err = row['sigSq_nr_err']

    # Calculate the partial derivatives
    dQ_dEee = -2 * (sigSq_E - sigSq_noise - sigSq_nr) / (Eee**3)
    dQ_dsigSq_E = 1 / (Eee**2)
    dQ_dsigSq_noise = -1 / (Eee**2)
    dQ_dsigSq_nr = -1 / (Eee**2)

    # Calculate the error using error propagation formula
    delta_QSF_sigma = math.sqrt(
        (dQ_dEee * Eee_err)**2 +
        (dQ_dsigSq_E * sigSq_E_err)**2 +
        (dQ_dsigSq_noise * sigSq_noise_err)**2 +
        (dQ_dsigSq_nr * sigSq_nr_err)**2
    )

    return delta_QSF_sigma

# Apply the error_QSF_sigma_row function to each row of your DataFrame
df_long['QSF_err'] = df_long.apply(error_QSF_sigma_row, axis=1)

# Now we want to fit the data to a linear function like 1/x
# Define the function
def QSF_func(x, a):
    return a * (1/x)

# Fit the data to the function
popt, pcov = curve_fit(QSF_func, df_long['Erecoil'], df_long['QSF'], sigma=df_long['QSF_err'])

# Calculate the standard deviations (errors) for the parameters
perr = np.sqrt(np.diag(pcov))

# plot the fit curve from 0 to 6 keV
x = np.arange(0.92, 6, 0.01)
# Plot the curve with the label in equation form, with parameter errors
plt.figure(figsize=(10, 6))
plt.plot(x, QSF_func(x, *popt), label='QSF = a/Erecoil (keVnr) \na = {:.5f} $\pm$ {:.5f}'.format(popt[0], perr[0]))

# Plot the QSF vs Erecoil with error bars
plt.errorbar(df_long['Erecoil'], df_long['QSF'], yerr=df_long['QSF_err'],xerr=(df_long['Erec_errMin'], df_long['Erec_errMax']), 
            fmt='o', markersize=4, capsize=2, ecolor='black', color='red', elinewidth=1)

# Plot a vertical line at 2.577 keV
plt.axvline(x=2.577, color='red', linestyle='--', label='example: 0.5 keVee threshold')

# Plot a horizontal line at 0.056 which ends at 2.6 keV
x_values = [0, 1.2]
y_value = 0.056
plt.plot(x_values, [y_value, y_value], color='C0', linestyle='--', label='QSC Linhard prediction\n at low energy (0.056)')

plt.xlabel('Recoil Energy (keVnr)')
plt.ylabel('Quenching Smearing Factor')
plt.title('Leading order $1/E$ fit to Quenching Smearing Factor\n (adapted from L. Li PhD dissertation)')
plt.xlim(0,6)
plt.ylim(0,0.065)
plt.legend()
plt.show()

# Generate a dataframe with column Erecoil and variance. These will get populated as we go
df_data = pd.DataFrame(columns=['Erecoil', 'Eee', 'variance'])
# Starting at 0 keV and going to 85 keV, populate the dataframe with Erecoil values in steps of 0.1 keV
df_data['Erecoil'] = np.arange(0.1, 85.6, 0.1)
# Round the Erecoil values to 2 decimal places
df_data['Erecoil'] = df_data['Erecoil'].round(2)
# Calculate the corresponding Eee values by using the corresponding Erecoil value in QF_StandLind
df_data['Eee'] = df_data['Erecoil'].apply(QF_StandLind)
# Add a column that gives the QSF value for each Erecoil value
df_data['QSF'] = QSF_func(df_data['Erecoil'], *popt)
# Calculate the variance (refer to eq. 4.10 in Lindhard's paper)
df_data['variance'] = df_data['QSF'] * (df_data['Eee'])**2 
# Add a column for the standard deviation
df_data['std'] = np.sqrt(df_data['variance'])



# Plots to show the effect of quenching (with a spread! Instead of a one-to-one mapping of E_nr to E_ee)
# Plot a gaussian distribution. The mean is the Eee value, and the standard deviation is the corresponding std in df_data
# Generate x values (range of values where you want to plot the distribution)
mean1 = df_data.loc[df_data['Erecoil'] == 0.5, 'Eee']
mean2 = df_data.loc[df_data['Erecoil'] == 2.0, 'Eee']
mean3 = df_data.loc[df_data['Erecoil'] == 10.0, 'Eee']
mean4 = df_data.loc[df_data['Erecoil'] == 80.0, 'Eee']
std_dev1 = df_data.loc[df_data['Erecoil'] == 0.5, 'std']
std_dev2 = df_data.loc[df_data['Erecoil'] == 2.0, 'std']
std_dev3 = df_data.loc[df_data['Erecoil'] == 10.0, 'std']
std_dev4 = df_data.loc[df_data['Erecoil'] == 80.0, 'std']
# Generate the x values from 0 to 85.05 keV in steps of 0.1 keV
en_arr = np.arange(0.1, 85.6, 0.1) #OFFSET THIS BY 0.05
# Generate the y values (probability density function) pulling from the df_data mean = Eee and std_dev = std, for Erecoil = 0.5 keV
y1 = stats.norm.pdf(en_arr, mean1, std_dev1)
y2 = stats.norm.pdf(en_arr, mean2, std_dev2)
y3 = stats.norm.pdf(en_arr, mean3, std_dev3)
y4 = stats.norm.pdf(en_arr, mean4, std_dev4)
# normalize the y values to 1
y1 = y1 / sum(y1)
y2 = y2 / sum(y2)
y3 = y3 / sum(y3)
y4 = y4 / sum(y4)
# Plot the distribution
plt.figure(figsize=(10, 6))
plt.plot(en_arr, y1, label='E_nr = 0.5 keV')
plt.plot(en_arr, y2, label='E_nr = 2.0 keV')
plt.plot(en_arr, y3, label='E_nr = 10.0 keV')
plt.plot(en_arr, y4, label='E_nr = 80.0 keV')
# plt.xlim(20,30)
# plt.ylim(0,0.2)
plt.xlabel('Eee (keV)')
plt.ylabel('Probability')
plt.legend()
plt.title('Gaussian Distribution of Eee for Different E_(nuclear recoil) Values')
# plt.show()


# Now to generate a dataframe for important values, that'll be saved as our "detector effects" matrix
# This is sort of redundant as some of these calculations were done earlier, but I'm copying from an old jupyter notebook
# Generate a dataframe with column Erecoil and variance. These will get populated as we go
df_QFMatrix = pd.DataFrame(columns=['Erecoil', 'Eee', 'variance'])
# Starting at 0 keV and going to 85 keV, populate the dataframe with Erecoil values in steps of 0.1 keV
df_QFMatrix['Erecoil'] = np.arange(0.1, 85.6, 0.1)
# round the Erecoil values to 2 decimal places
df_QFMatrix['Erecoil'] = df_QFMatrix['Erecoil'].round(2)
# Calculate the corresponding Eee values by using the corresponding Erecoil value in QF_StandLind
df_QFMatrix['Eee'] = df_QFMatrix['Erecoil'].apply(QF_StandLind)
# Add a column that gives the QSF value for each Erecoil value
df_QFMatrix['QSF'] = QSF_func(df_QFMatrix['Erecoil'], *popt)
# Calculate the variance
df_QFMatrix['variance'] = df_QFMatrix['QSF'] * (df_QFMatrix['Eee'])**2
# Add a column for the standard deviation
df_QFMatrix['std'] = np.sqrt(df_QFMatrix['variance'])

# Generate a new dataframe called df_GaussianPDF with a column called Eee which goes from 0.1 to 85.6 keV in steps of 0.1 keV
df_GaussianPDF = pd.DataFrame(columns=['Eee'])
df_GaussianPDF['Eee'] = np.arange(0.1, 85.6, 0.1)
# Round the Eee values to 2 decimal places
df_GaussianPDF['Eee'] = df_GaussianPDF['Eee'].round(2)

pdf_dict = {}

for Erecoil in df_QFMatrix['Erecoil']:
    # Calculate the mean and std_dev for the Erecoil value
    mean = df_QFMatrix.loc[df_QFMatrix['Erecoil'] == Erecoil, 'Eee'].values[0]
    std_dev = df_QFMatrix.loc[df_QFMatrix['Erecoil'] == Erecoil, 'std'].values[0]
    # Calculate the y values (probability density function)
    y = stats.norm.pdf(df_GaussianPDF['Eee'], mean, std_dev)
    # Normalize the y values to 1
    y = y / sum(y)
    # Store the PDF in the dictionary
    pdf_dict[Erecoil] = y

# Create the df_GaussianPDF DataFrame from the dictionary
df_GaussianPDF = pd.DataFrame(pdf_dict)

# It seems that when you initialize the DF with the dictionary, it overwrites what was there before. So we need to add the Eee column back in
df_GaussianPDF['Eee'] = np.arange(0.1, 85.6, 0.1)
# Round the Eee values to 2 decimal places
df_GaussianPDF['Eee'] = df_GaussianPDF['Eee'].round(2)

# let's do a plot to check this came out correctly
# Plot the Gaussian PDF vs Eee 
plt.figure(figsize=(10, 6))
plt.plot(df_GaussianPDF['Eee'], df_GaussianPDF[0.5], label='Erecoil = 0.5 keV')
plt.plot(df_GaussianPDF['Eee'], df_GaussianPDF[2.0], label='Erecoil = 2.0 keV')
plt.plot(df_GaussianPDF['Eee'], df_GaussianPDF[10.0], label='Erecoil = 10.0 keV')
plt.plot(df_GaussianPDF['Eee'], df_GaussianPDF[80.0], label='Erecoil = 80.0 keV')
plt.plot(df_GaussianPDF['Eee'], df_GaussianPDF[5.5], label='Erecoil = 5.5 keV')
plt.plot(df_GaussianPDF['Eee'], df_GaussianPDF[20.0], label='Erecoil = 20.0 keV')
plt.plot(df_GaussianPDF['Eee'], df_GaussianPDF[30.0], label='Erecoil = 30.0 keV')
plt.plot(df_GaussianPDF['Eee'], df_GaussianPDF[40.0], label='Erecoil = 40.0 keV')
plt.plot(df_GaussianPDF['Eee'], df_GaussianPDF[50.0], label='Erecoil = 50.0 keV')
plt.plot(df_GaussianPDF['Eee'], df_GaussianPDF[60.0], label='Erecoil = 60.0 keV')
plt.plot(df_GaussianPDF['Eee'], df_GaussianPDF[70.0], label='Erecoil = 70.0 keV')
plt.xlabel('Eee (keV)')
plt.ylabel('Probability')
plt.legend()
plt.title('Gaussian Distribution of Eee for Different E_(nuclear recoil) Values')
plt.show()


# Convert the df_GaussianPDF into a numpy array, without the row and column labels
gaussianPDF_arr = df_GaussianPDF.values

# print the shape of the array
print('Shape of the array:', gaussianPDF_arr.shape)

# remove the last column which is empty
gaussianPDF_arr = gaussianPDF_arr[:, :-1]

# print the shape of the array
print('Shape of the array (adjusted):', gaussianPDF_arr.shape)

# save the array as a .npy file
np.save('quenchingMatrix.npy', gaussianPDF_arr)

# save the array as a .txt file
np.savetxt('quenchingMatrix_test.txt', gaussianPDF_arr)

