#!/bin/bash
# Bash strict mode:
set -euo pipefail
IFS=$'\n\t'
# End of strict mode
energies=(20 50 100)

for energy in "${energies[@]}"; do
    genie -l -b -q "analytic_tau_pol_dis_int.C(\"../data/test_genie_NuTau_${energy}.0_GeV_event_info_dis.csv\",\"../data/test_genie_NuTau_${energy}.0_GeV_event_info_f.csv\")"
    genie -l -b -q "crv98lo_pdf.C(\"../data/test_genie_NuTau_${energy}.0_GeV_event_info_f.csv\",\"../data/test_genie_NuTau_${energy}.0_GeV_event_info_pdf.csv\")"
    genie -l -b -q "analytic_tau_pol_res_int.C(\"../data/test_genie_NuTau_${energy}.0_GeV_event_info_pdf.csv\",\"../data/test_genie_NuTau_${energy}.0_GeV_event_info_sig.csv\")"
done