import numpy as np
import matplotlib.pyplot as plt
from ederMatrixGenerator_paper import ederMatrix
from FF_quencher_eder_paper import FF_quencher_eder

# Generate the EDER matrix for detector 21
detectorNumber = 21
ederMatrix_21 = ederMatrix(detectorNumber)

#print the values in index 500 of the EDER matrix
print(ederMatrix_21[500])
print(ederMatrix_21[854])

# generate an energy array for plotting purposes (0.1 to 85.5 keV in steps of 0.1 keV)
energy = np.arange(0.1, 85.6, 0.1)

# Plot the gaussian values for a variety of Erecoil values
plt.figure(figsize=(10, 6))
plt.plot(energy, ederMatrix_21[0], label='Erecoil = 0.1 keV')
plt.plot(energy, ederMatrix_21[100], label='Erecoil = 10.1 keV')
plt.plot(energy, ederMatrix_21[500], label='Erecoil = 50.1 keV')
plt.plot(energy, ederMatrix_21[800], label='Erecoil = 80.1 keV')
plt.plot(energy, ederMatrix_21[854], label='Erecoil = 85.5 keV')
plt.xlabel('Energy (keV)')
plt.ylabel('Probability Density')
plt.title('Energy dependent energy resolution "smearing" (Detector 21)')
plt.legend()
plt.show()



# *************** Check of unit conversion ***************
# read in the sns_diff_rates file
dataIn = np.loadtxt('sns_diff_rates_alliso-gemini_campaign2_crosscheck-Ge-unity.out')

energyIn = dataIn[:, 0]
dukecevnsSpec = dataIn[:, 1]

    # ***** UNIT ADJUSTMENTS *****

    # Convert counts down to counts per keV
dukecevnsSpec = dukecevnsSpec / 1000

    # Account for the 0.1 keV binning
dukecevnsSpec = dukecevnsSpec / 10

    # Convert energy from MeV to keV
energyIn = energyIn * 1000

    # round the energy to 1 decimal place
energyIn = np.round(energyIn, 1)

# print the sum of the counts in the spectrum
print('Sum of counts in the spectrum:', np.sum(dukecevnsSpec))

# plot the spectrum to check the unit conversion
plt.figure(figsize=(10, 6))
plt.plot(energyIn, dukecevnsSpec)
plt.xlabel('Energy (keV_nr)')
plt.ylabel('Counts / 0.1 keV_nr')
plt.title('Raw DukeCEvNS spectrum')
# show a grid
plt.grid()
plt.show()


# *************** Check of form factor calculation ***************

# Here are the functions used for form factors and TtoQ in FF_quencher_eder_paper.py
m = 72.63 * 931.5 # MeV/c^2
def TtoQ(T):
    #return np.sqrt(2*(74 * 931.5)*T)
    return np.sqrt((2*(72.63 * 931.5) * T) + (T**2))
    
def FF_KleinNystrand(Q):
    Q = Q / 197.3  # to convert to fm^-1
    Rn = 1.2 * (72.63)**(1/3)
    a = 0.7
    return (3 / (Q**3 * Rn**3)) * (np.sin(Q * Rn) - (Q * Rn * np.cos(Q * Rn))) * (1 / (1 + (Q**2 * a**2)))

# Generate an array of Erecoil values from 0.1 to 85.6 keV in steps of 0.1 keV
energy_ff = np.arange(0.1, 85.6, 0.1)
# Convert Erecoil_arr to MeV
energy_ff = energy_ff / 1000
# Convert to Q
Q_ff = TtoQ(energy_ff)
# Calculate the FF values for klein nystrand 
FF_KN = FF_KleinNystrand(Q_ff)

# apply FF to the dukecevns spectrum
FF_dukecevnsSpec = FF_KN * FF_KN * dukecevnsSpec

# plot the spectrum to check the form factor calculation
plt.figure(figsize=(10, 6))
plt.plot(energyIn, FF_dukecevnsSpec, label="Klein-Nystrand FF")
plt.plot(energyIn, dukecevnsSpec, label="Unity FF")
plt.xlabel('Energy (keV_nr)')
plt.ylabel('Counts / 0.1 keV_nr')
plt.title('DukeCEvNS spectrum FF comparison')
plt.legend()
plt.show()