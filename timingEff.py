import numpy as np
import matplotlib.pyplot as plt

# This serves to interpolate a series of timing efficiencies on a per-detector basis. 
# Will output a .npy file for each detector that has timing efficiencies (as a function of keVee)
# The output is intended to be used as the last matrix in FF_quencher_eder.py

# Read in each detector efficiency file from the timing_efficiencies directory
# The files are named in format eff_Ge**.txt

eff21_in = np.loadtxt('timing_efficiencies/eff_Ge21.txt')
eff23_in = np.loadtxt('timing_efficiencies/eff_Ge23.txt')
eff25_in = np.loadtxt('timing_efficiencies/eff_Ge25.txt')
eff26_in = np.loadtxt('timing_efficiencies/eff_Ge26.txt')
eff28_in = np.loadtxt('timing_efficiencies/eff_Ge28.txt')

# separate data into energy and efficiency
energy21 = eff21_in[:, 0]
eff21 = eff21_in[:, 1]
energy23 = eff23_in[:, 0]
eff23 = eff23_in[:, 1]
energy25 = eff25_in[:, 0]
eff25 = eff25_in[:, 1]
energy26 = eff26_in[:, 0]
eff26 = eff26_in[:, 1]
energy28 = eff28_in[:, 0]
eff28 = eff28_in[:, 1]

# Define an energy array from 0.1 to 85.5 keV in steps of 0.1 keV
energy = np.arange(0.1, 85.6, 0.1)
# Round the energy to 1 decimal place
energy = np.round(energy, 1)

# prevent the interpolated data from being in scientific notation
np.set_printoptions(suppress=True)

# Interpolate the efficiency data to the energy array
eff21_interp = np.interp(energy, energy21, eff21)
eff23_interp = np.interp(energy, energy23, eff23)
eff25_interp = np.interp(energy, energy25, eff25)
eff26_interp = np.interp(energy, energy26, eff26)
eff28_interp = np.interp(energy, energy28, eff28)


# plot the interpolated data
plt.figure(figsize=(10, 6))
plt.plot(energy, eff21_interp, label='Ge21')
plt.plot(energy, eff23_interp, label='Ge23')
plt.plot(energy, eff25_interp, label='Ge25')
plt.plot(energy, eff26_interp, label='Ge26')
plt.plot(energy, eff28_interp, label='Ge28')
plt.xlim(0, 15)
plt.xlabel('Energy (keVee)')
plt.ylabel('Timing Efficiency')
plt.title('Timing Efficiency Interpolation')
plt.legend()
plt.show()

# Save the interpolated data as .txt files
np.savetxt('timing_efficiencies/eff_Ge21_interp.txt', np.column_stack((energy, eff21_interp)))
np.savetxt('timing_efficiencies/eff_Ge23_interp.txt', np.column_stack((energy, eff23_interp)))
np.savetxt('timing_efficiencies/eff_Ge25_interp.txt', np.column_stack((energy, eff25_interp)))
np.savetxt('timing_efficiencies/eff_Ge26_interp.txt', np.column_stack((energy, eff26_interp)))
np.savetxt('timing_efficiencies/eff_Ge28_interp.txt', np.column_stack((energy, eff28_interp)))


