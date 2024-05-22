#!/bin/bash
# Bash strict mode:
set -euo pipefail
IFS=$'\n\t'
# End of strict mode

energies=(5 10 20 50 100) # 5 10 20 50 100

for energy in "${energies[@]}"; do
    # Calculate the form factors F1 - F5 using the GENIE library
    #  This adds five new columns: F1, F2, F3, F4, F5
    #  NOTE: the output of these are not used, as they give unphysical values 
    genie -l -b -q "analytic_tau_pol_dis_int.C(\"../data/test_genie_NuTau_${energy}.0_GeV_event_info_dis.csv\",\"../data/test_genie_NuTau_${energy}.0_GeV_event_info_f.csv\",false)"
    # Calculate the GRV98LO PDF values for each DIS event from x and Q^2
    #  This adds 6 new columns: fuv, fus, fdv, fds, fs, fc
    #  Each column corresponds to the PDF value for the event for the valence up quark, sea up quark, valence down quark, sea down quark, strange quark, and charm quark respectively.
    #  These are used for the polarization calculations
    #  When the last input argument is true, the charm-corrected x value is used, as described in Hagiwara et al. (2003). This is the default.
    genie -l -b -q "grv98lo_pdf.C(\"../data/test_genie_NuTau_${energy}.0_GeV_event_info_f.csv\",\"../data/test_genie_NuTau_${energy}.0_GeV_event_info_pdf.csv\",true)"
    # Calculate the Sigma values of for RES events using the Berger-Sehgal model
    #  This adds two new columns: sigmm, sigpp
    #  The longitudinal polarization is then given by (sigmm - sigpp) / (sigmm + sigpp)
    #  The goal is to also add a way of calculating the transverse polarization, but this seems a bit more complicated
    #  as it is not relevant for the GENIE cross section calculations, and therefore not implemented in the library.
    genie -l -b -q "analytic_tau_pol_res_int.C(\"../data/test_genie_NuTau_${energy}.0_GeV_event_info_pdf.csv\",\"../data/test_genie_NuTau_${energy}.0_GeV_event_info_sig.csv\")"
done
