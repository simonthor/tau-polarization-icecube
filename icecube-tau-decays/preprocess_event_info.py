import vector
import numpy as np
import pandas as pd
from argparse import ArgumentParser


def mass_corrected_xi(x, Q2, charm):
    xi = x.copy()
    m_charm = 1.27 # GeV
    # If a charm quark is produced, replace x with xi
    xi[charm] = x[charm] / (Q2[charm] / (Q2[charm] + m_charm**2))
    # Cap the x to 1
    xi[xi > 1] = 1
    # print(xi)
    return xi


def main():
    argparser = ArgumentParser()
    argparser.add_argument("-ip", type=str, required=True)  # Input particle info file name
    argparser.add_argument("-ie", type=str, required=True)  # Output particle info file name
    argparser.add_argument("-o", type=str, required=True)  # Output file name for csv file
    
    args = argparser.parse_args()

    particle_info = pd.read_csv(args.ip)
    event_info = pd.read_csv(args.ie)
    
    nutau = particle_info.groupby("event_num").nth(1)
    nutau4m = vector.array({"E": nutau["E"], "px": nutau["px"], "py": nutau["py"], "pz": nutau["pz"]})
    
    assert set(particle_info.groupby("event_num").first()["pdg"].unique().tolist()) == {2212, 1000080160}

    event_info["atom"] = 16
    event_info.loc[(particle_info.groupby("event_num").first()["pdg"] == 2212).values, "atom"] = 1
    
    event_info["Mnuc"] = np.sqrt(event_info["En"]**2 - event_info["pxn"]**2 - event_info["pyn"]**2 - event_info["pzn"]**2)
    event_info["Enu"] = nutau4m["E"]
    event_info["pxnu"] = nutau4m["px"]
    event_info["pynu"] = nutau4m["py"]
    event_info["pznu"] = nutau4m["pz"]

    event_info["xi"] = mass_corrected_xi(event_info["xs"], event_info["Q2"], event_info["charm"])

    # Write to file
    event_info.to_csv(args.o, index=False)


if __name__ == "__main__":
    main()
