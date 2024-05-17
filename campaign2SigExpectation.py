import numpy as np
import matplotlib.pyplot as plt
from FF_quencher_eder_paper import FF_quencher_eder
import os


cwd = os.getcwd()

# Import the standard dukeCEvNS output, scaled in 1 kg and 1 GWHr

data = np.loadtxt('sns_diff_rates_alliso-gemini_campaign2_crosscheck-Ge-unity.out')

# define the masses and exposures
# below are the nominal masses from mirion
# ge21mass = 2.191
# ge23mass = 2.201
# ge25mass = 2.200
# ge26mass = 2.270
# ge28mass = 2.290

# below are the active masses estimates
ge21mass = 2.13
ge23mass = 2.13
ge25mass = 2.13 
ge26mass = 2.13
ge28mass = 2.14

ge21expo = 0.65
ge23expo = 0.77
ge25expo = 0.74
ge26expo = 0.36
ge28expo = 0.85

# Apply FF, quench, smear, and timing efficiency

ge21Counts_KN, energy = FF_quencher_eder(data, 21, "KN", ge21mass, ge21expo)
ge23Counts_KN, _ = FF_quencher_eder(data, 23, "KN", ge23mass, ge23expo)
ge25Counts_KN, _ = FF_quencher_eder(data, 25, "KN", ge25mass, ge25expo)
ge26Counts_KN, _ = FF_quencher_eder(data, 26, "KN", ge26mass, ge26expo)
ge28Counts_KN, _ = FF_quencher_eder(data, 28, "KN", ge28mass, ge28expo)

energy_bins = energy[1] - energy[0]

# Define the ROI
ROI_low = 14        # 1.5 keV
ROI_high = 85       # 8.6 keV (because non-inclusive)
highROI_low = 105   # 10.6 keV
highROI_high = 200  # 20.1 keV (because non-inclusive)

print("the energy value at index ", ROI_low, " is ", energy[ROI_low])
print("the energy value at index ", ROI_high, " is ", energy[ROI_high])
print("the energy value at index ", highROI_low, " is ", energy[highROI_low])
print("the energy value at index ", highROI_high, " is ", energy[highROI_high])

# Construct a total counts for KN and Helm

totCounts_KN = ge21Counts_KN + ge23Counts_KN + ge25Counts_KN + ge26Counts_KN + ge28Counts_KN

# Integrate the regions
ge21ROI_counts_KN = np.sum(ge21Counts_KN[ROI_low:ROI_high])
ge23ROI_counts_KN = np.sum(ge23Counts_KN[ROI_low:ROI_high])
ge25ROI_counts_KN = np.sum(ge25Counts_KN[ROI_low:ROI_high])
ge26ROI_counts_KN = np.sum(ge26Counts_KN[ROI_low:ROI_high])
ge28ROI_counts_KN = np.sum(ge28Counts_KN[ROI_low:ROI_high])
totCountsROI_KN = np.sum(totCounts_KN[ROI_low:ROI_high])


ge21Wide_counts_KN = np.sum(ge21Counts_KN[ROI_low:highROI_high])
ge23Wide_counts_KN = np.sum(ge23Counts_KN[ROI_low:highROI_high])
ge25Wide_counts_KN = np.sum(ge25Counts_KN[ROI_low:highROI_high])
ge26Wide_counts_KN = np.sum(ge26Counts_KN[ROI_low:highROI_high])
ge28Wide_counts_KN = np.sum(ge28Counts_KN[ROI_low:highROI_high])
totCountsWide_KN = np.sum(totCounts_KN[ROI_low:highROI_high])


ge21High_counts_KN = np.sum(ge21Counts_KN[highROI_low:highROI_high])
ge23High_counts_KN = np.sum(ge23Counts_KN[highROI_low:highROI_high])
ge25High_counts_KN = np.sum(ge25Counts_KN[highROI_low:highROI_high])
ge26High_counts_KN = np.sum(ge26Counts_KN[highROI_low:highROI_high])
ge28High_counts_KN = np.sum(ge28Counts_KN[highROI_low:highROI_high])
totCountsHigh_KN = np.sum(totCounts_KN[highROI_low:highROI_high])

totcountsAnalysisROI = totCountsROI_KN + totCountsHigh_KN

print(f"Total counts in the ROI 1.5-8.5 + 10.6-20 keVee for KN: {totcountsAnalysisROI:.5f}")


# Final plot
plt.figure(figsize=(8,6))
plt.hist(energy, bins=energy, weights=ge21Counts_KN, histtype = 'step', label=f'Ge21, {(ge21ROI_counts_KN + ge21High_counts_KN):.5f} counts')
plt.hist(energy, bins=energy, weights=ge23Counts_KN, histtype = 'step', label=f'Ge23, {(ge23ROI_counts_KN + ge23High_counts_KN):.5f} counts')
plt.hist(energy, bins=energy, weights=ge25Counts_KN, histtype = 'step', label=f'Ge25, {(ge25ROI_counts_KN + ge25High_counts_KN):.5f} counts')
plt.hist(energy, bins=energy, weights=ge26Counts_KN, histtype = 'step', label=f'Ge26, {(ge26ROI_counts_KN + ge26High_counts_KN):.5f} counts')
plt.hist(energy, bins=energy, weights=ge28Counts_KN, histtype = 'step', label=f'Ge28, {(ge28ROI_counts_KN + ge28High_counts_KN):.5f} counts')
plt.xlim(0,20)
plt.xlabel('Energy (keVee)')
plt.ylabel(f'Counts / {energy_bins} keVee')
plt.title(f'Ge-Mini Campaign 2 exposure-scaled CEvNS recoil spectra\n (k=0.157, Klein-Nystrand FF, Quenching-smearing, EDER, timing eff)\n ROI: 1.5-8.5 + 10.6-20.0 keVee. Total ROI counts: {totcountsAnalysisROI:.5f}')
plt.legend()

# Plot the spectrum 
# plt.figure(figsize=(8,6))
# plt.hist(energy, bins=energy, weights=ge21Counts_KN, histtype = 'step', label=f'Ge21, {ge21ROI_counts_KN:.5f} counts')
# plt.hist(energy, bins=energy, weights=ge23Counts_KN, histtype = 'step', label=f'Ge23, {ge23ROI_counts_KN:.5f} counts')
# plt.hist(energy, bins=energy, weights=ge25Counts_KN, histtype = 'step', label=f'Ge25, {ge25ROI_counts_KN:.5f} counts')
# plt.hist(energy, bins=energy, weights=ge26Counts_KN, histtype = 'step', label=f'Ge26, {ge26ROI_counts_KN:.5f} counts')
# plt.hist(energy, bins=energy, weights=ge28Counts_KN, histtype = 'step', label=f'Ge28, {ge28ROI_counts_KN:.5f} counts')
# plt.xlim(0,20)
# plt.xlabel('Energy (keVee)')
# plt.ylabel(f'Counts / {energy_bins} keVee')
# plt.title('Campaign 2 exposure-scaled CEvNS recoil spectra\n (k=0.157, Klein-Nystrand FF, counts from 1.5-8.5 keVee)')
# plt.legend()
# # plt.savefig(cwd + '/campaign2_recoilPlots/campaign2_CEvNS_KN_ROI_individual.png', dpi=600)

# plt.figure(figsize=(8,6))
# plt.hist(energy, bins=energy, weights=ge21Counts_KN, histtype='step', label=f'Ge21, {ge21High_counts_KN:.5f} counts')
# plt.hist(energy, bins=energy, weights=ge23Counts_KN, histtype='step', label=f'Ge23, {ge23High_counts_KN:.5f} counts')
# plt.hist(energy, bins=energy, weights=ge25Counts_KN, histtype='step', label=f'Ge25, {ge25High_counts_KN:.5f} counts')
# plt.hist(energy, bins=energy, weights=ge26Counts_KN, histtype='step', label=f'Ge26, {ge26High_counts_KN:.5f} counts')
# plt.hist(energy, bins=energy, weights=ge28Counts_KN, histtype='step', label=f'Ge28, {ge28High_counts_KN:.5f} counts')
# plt.xlim(0,20)
# plt.xlabel('Energy (keVee)')
# plt.ylabel(f'Counts / {energy_bins} keVee')
# plt.title('Campaign 2 exposure-scaled CEvNS recoil spectra\n (k=0.157, Klein-Nystrand FF, counts from 10.6-20.0 keVee)')
# plt.legend()
# plt.savefig(cwd + '/campaign2_recoilPlots/campaign2_CEvNS_KN_HighE_individual.png', dpi=600)


# plt.figure(figsize=(8,6))
# plt.hist(energy, bins=energy, weights=totCounts_KN, histtype='step', label=f'Klein-Nystrand FF, {totCountsROI_KN:.8f} counts')
# plt.hist(energy, bins=energy, weights=totCounts_Helm, histtype='step', label=f'Helm FF, {totCountsROI_Helm:.8f} counts')
# plt.xlim(0,20)
# plt.xlabel('Energy (keVee)')
# plt.ylabel(f'Counts / {energy_bins} keVee')
# plt.title('Campaign 2 exposure-scaled CEvNS recoil spectra\n (k=0.157, FF comparison, counts from 1.5-8.5 keVee)')
# plt.legend()
# plt.savefig(cwd + '/campaign2_recoilPlots/campaign2_CEvNS_KNvHelm_ROI.png', dpi=600)


# plt.figure(figsize=(8,6))
# plt.hist(energy, bins=energy, weights=totCounts_KN, histtype='step', label=f'Klein-Nystrand FF, {totCountsHigh_KN:.8f} counts')
# plt.hist(energy, bins=energy, weights=totCounts_Helm, histtype='step', label=f'Helm FF, {totCountsHigh_Helm:.8f} counts')
# plt.xlim(0,20)
# plt.xlabel('Energy (keVee)')
# plt.ylabel(f'Counts / {energy_bins} keVee')
# plt.title('Campaign 2 exposure-scaled CEvNS recoil spectra\n (k=0.157, FF comparison, counts from 10.6-20.0 keVee)')
# plt.legend()
# plt.savefig(cwd + '/campaign2_recoilPlots/campaign2_CEvNS_KNvHelm_HighE.png', dpi=600)

# plt.figure(figsize=(8,6))
# plt.hist(energy, bins=energy, weights=totCounts_KN, histtype='step', label=f'Klein-Nystrand FF, {totCountsWide_KN:.8f} counts')
# plt.hist(energy, bins=energy, weights=totCounts_Helm, histtype='step', label=f'Helm FF, {totCountsWide_Helm:.8f} counts')
# plt.xlim(0,20)
# plt.xlabel('Energy (keVee)')
# plt.ylabel(f'Counts / {energy_bins} keVee')
# plt.title('Campaign 2 exposure-scaled CEvNS recoil spectra\n (k=0.157, FF comparison, counts from 1.5-20.0 keVee)')
# plt.legend()
# plt.savefig(cwd + '/campaign2_recoilPlots/campaign2_CEvNS_KNvHelm_WideE.png', dpi=600)

# klein-nystrand minus helm
# plt.figure(figsize=(8,6))
# plt.hist(energy, bins=energy, weights=totCounts_KN-totCounts_Helm, histtype='step', label=f'KN-Helm (1.5-8.5 keV), {totCountsROI_KN-totCountsROI_Helm:.8f} counts')
# plt.hist(energy, bins=energy, weights=totCounts_KN-totCounts_Helm, histtype='step', label=f'KN-Helm (10.6-20 keV), {totCountsHigh_KN-totCountsHigh_Helm:.8f} counts')
# plt.xlim(0,20)
# plt.xlabel('Energy (keVee)')
# plt.ylabel(f'Counts difference / {energy_bins} keVee')
# plt.title('Campaign 2 exposure-scaled FF difference for ROI and HighE regions\n (k=0.157, counts from 1.5-8.5 keVee and 10.6-20.0 keVee)')
# plt.legend()


# Show the plots
plt.show()