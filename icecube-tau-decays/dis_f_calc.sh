#!/bin/bash
# Bash strict mode:
set -euo pipefail
IFS=$'\n\t'
# End of strict mode

genie -l -b -q 'analytic_tau_pol_dis.cc("../data/test_genie_NuTau_5.0_GeV_event_info_dis.csv","../data/test_genie_NuTau_5.0_GeV_event_info_f.csv")'
genie -l -b -q 'crv98lo_pdf.C("../data/test_genie_NuTau_5.0_GeV_event_info_f.csv","../data/test_genie_NuTau_5.0_GeV_event_info_pdf.csv")'

genie -l -b -q 'analytic_tau_pol_dis.cc("../data/test_genie_NuTau_10.0_GeV_event_info_dis.csv","../data/test_genie_NuTau_10.0_GeV_event_info_f.csv")'
genie -l -b -q 'crv98lo_pdf.C("../data/test_genie_NuTau_10.0_GeV_event_info_f.csv","../data/test_genie_NuTau_10.0_GeV_event_info_pdf.csv")'

genie -l -b -q 'analytic_tau_pol_dis.cc("../data/test_genie_NuTau_20.0_GeV_event_info_dis.csv","../data/test_genie_NuTau_20.0_GeV_event_info_f.csv")'
genie -l -b -q 'crv98lo_pdf.C("../data/test_genie_NuTau_20.0_GeV_event_info_f.csv","../data/test_genie_NuTau_20.0_GeV_event_info_pdf.csv")'

genie -l -b -q 'analytic_tau_pol_dis.cc("../data/test_genie_NuTau_50.0_GeV_event_info_dis.csv","../data/test_genie_NuTau_50.0_GeV_event_info_f.csv")'
genie -l -b -q 'crv98lo_pdf.C("../data/test_genie_NuTau_50.0_GeV_event_info_f.csv","../data/test_genie_NuTau_50.0_GeV_event_info_pdf.csv")'

genie -l -b -q 'analytic_tau_pol_dis.cc("../data/test_genie_NuTau_100.0_GeV_event_info_dis.csv","../data/test_genie_NuTau_100.0_GeV_event_info_f.csv")'
genie -l -b -q 'crv98lo_pdf.C("../data/test_genie_NuTau_100.0_GeV_event_info_f.csv","../data/test_genie_NuTau_100.0_GeV_event_info_pdf.csv")'
