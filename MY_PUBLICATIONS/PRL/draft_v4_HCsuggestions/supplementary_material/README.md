#SUPPLEMENTAL MATERIAL
This directory contains relevant work for supplemental material of the PRL paper.

Directory Contents:

/supp_tables: containes Xsec, statistical and systematic uncertainties as
well as correction factors. These should be organized in Tables as supplemental
material.


/latex_files: contains the .tex code to produce the supplemental_material.pdf
file


/root_files: contains root files for all_thnq, thnq=35, 45, 75 deg settings
to make plots of relevant analysis plots and decide which setting is
best to show.

/supp_plots: contains plots produced from root_files that can potentially
be in the supplementary materials

/codes: contains relevant codes to read root files and produce plots,
as well as tables.


#The tables contain momentum distribution data, correction factors and
relative errors (sys_norm + sys_kin + stats)

#Plots that need to be made:
--------
spectrometer acceptance plots that show the HMS determined the SHMS acceptance.
DATA/SIMS comparisons for both HMS and SHMS of: delta vs. (xptar), delta vs. (yptar) for the 80 MeV setting only, integrated over all th_nq,
and all standard cuts applied.
--------

--------
Data Analysis Cuts (Data/SIMC), showing the cut region. Apply all cuts except the cut on the plot that is being shown. 

Emiss, HMS_Collimator, coin_time, Z_diff, calorimeter, th_nq, Q2
Show these for the 80, 580 and 750 with the greatest statistics.
-------
