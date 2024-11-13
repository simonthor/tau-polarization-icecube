import numpy as np
import vector
import pandas as pd
import itertools
import warnings
from argparse import ArgumentParser

# Tau mass
mtau = 1.77682 # GeV


def polarization_vector(
        tau4m: vector.MomentumNumpy4D, 
        nu4m: vector.MomentumNumpy4D, 
        nucleon4m: vector.MomentumNumpy4D, 
        int_type: str, 
        x: np.ndarray = None, 
        W: np.ndarray = None, 
        charm: np.ndarray = None, 
        structure_function_values: pd.DataFrame = None,
        res_sigs: pd.DataFrame = None):
    """Calculate the polarization vector of the tau lepton that is produced in neutrino-nucleon scattering.
    
    The formulas are based on Hagiwara et al. (https://arxiv.org/abs/hep-ph/0305324) for most of the equations
    and Kuzmin et al. (https://arxiv.org/abs/hep-ph/0312107v2) for the QEL structure functions.

    Args:
        tau4m: The 4-momentum of the outgoing tau lepton. Can be an array of 4-vectors.
        nu4m: The 4-momentum of the incoming tau neutrino. Can be an array of 4-vectors. Must be pointing in the +z direction.
        nucleon4m: The 4-momentum of the nucleon. Can be an array of 4-vectors. Must be at rest.
        int_type: The type of interaction that is happening. Can be one of the following:
            - "qel": Quasi-elastic scattering
            - "res": Resonance scattering (not implemented yet, will raise an error)
            - "dis": Deep inelastic scattering
        x: The Björken scaling variable. If not passed, it is calculated from the 4-momenta
        W: The hadronic invariant mass. If not passed, it is calculated from the 4-momenta
        charm: Whether a charm quark is produced in the interaction or not.
            Only used for DIS interactions.
            Must be an array of bools that is the same size as the 4-vectors. 
        structure_function_values: The PDF values for the specific neutrino-nucleus interaction.
            Only used for DIS interactions.
            Must be some kind of map from the quark PDG ID (all values from -6 to 6, excluding 0). 
            Can e.g. be a pandas dataframe with the columns being the PDG ID and the rows being the values for each event, 
            or a dict with keys being the PDG IDs and the values being arrays containing the PDF values for each event.
        res_sigs: The resonance partial differential cross sections for the specific neutrino-nucleus interaction.
            Only used for RES interactions.
            Must be some kind of map from the keys "sigmm", "sigpp", "sigmp" to arrays.
            Can e.g. be a pandas dataframe with the columns being "sigmm", "sigpp", "sigmp" and the rows being the values for each event, 
            or a dict with keys being "sigmm", "sigpp", "sigmp" and the values being arrays containing the partial cross sections for each event.

    Returns:
        A tuple of two arrays, the first being the transverse polarization component, and the second being the longitudinal polarization component.

    Raises:
        ValueError: If the interaction type is not supported.
        NotImplementedError: If the interaction type is "res".

    Note:
        The formulas only work when
            - The nucleon is at rest
            - The neutrino is pointing in the +z direction
        If they are not, one must use boost_rotated_4m to change the reference frame.
        See main() for an example of how this is done (it is relatively simple).
        
        Currently, only tau leptons are supported, not anti-tau leptons. In theory, this change should be relatively simple to implement, 
        with only some changes to this function and not the rest of the code.

        Also note that the resonance scattering calculations do not use the structure functions at all. 
        Instead, almost all calculations must have been done before-hand using GENIE.
    """
    # Define variables that are used and based on the 4-momenta
    # The same naming convention as in Hagiwara et al. is used
    M = nucleon4m.m
    ptau = tau4m.p
    Etau = tau4m.E
    theta = nu4m.deltaangle(tau4m)
    costheta = np.cos(theta)
    Enu = nu4m.E
    p = nucleon4m
    k = nu4m
    kprime = tau4m
    
    q = k - kprime
    Q2 = -q**2
    
    if x is None:
        # Raise a warning
        warnings.warn("x is not passed as a parameter, calculating it manually")
        x = Q2 / (2*p.dot(q)) # Björken scaling variable
    if W is None:
        # Raise a warning
        warnings.warn("W is not passed as a parameter, calculating it manually")
        W = (p + q)**2 # Hadronic invariant mass
    
    if int_type == "qel":
        W1 = W1qel(x, Q2, p, q, M)
        W2 = W2qel(x, Q2, p, q, M)
        W3 = W3qel(x, Q2, p, q, M)
        W4 = W4qel(x, Q2, p, q, M)
        W5 = W5qel(x, Q2, p, q, M)
    
    elif int_type == "res":
        return (
            # transverse polarization component
            (res_sigs["sigmp"] + res_sigs["sigmp"]) / (res_sigs["sigpp"] + res_sigs["sigmm"]).values, 
            # longitudinal polarization
            (res_sigs["sigpp"] - res_sigs["sigmm"]) / (res_sigs["sigpp"] + res_sigs["sigmm"]).values
        )
        
    elif int_type == "dis":
        if charm is None:
            raise ValueError("charm parameter must be set to an array of bools of the same size as the 4-vectors")
        if structure_function_values is None:
            raise ValueError("structure_function_values parameter must be set to a map from the quark PDG ID to the PDF values")
        
        W1 = Wndis_f(x, Q2, p, q, M, charm, structure_function_values, 1)
        W2 = Wndis_f(x, Q2, p, q, M, charm, structure_function_values, 2)
        W3 = Wndis_f(x, Q2, p, q, M, charm, structure_function_values, 3)
        W4 = Wndis_f(x, Q2, p, q, M, charm, structure_function_values, 4)
        W5 = Wndis_f(x, Q2, p, q, M, charm, structure_function_values, 5)

    else:
        raise ValueError(f"Unsupported interaction type {int_type = }")
    
    F = (
        (2*W1 + mtau**2 / M**2 * W4) * (Etau - ptau * costheta)
        + W2 * (Etau + ptau * costheta)
        + W3 / M * (Enu * Etau + ptau**2 - (Enu + Etau) * ptau * costheta)
        - mtau**2 / M * W5
    )

    return (
        # The transverse polarization component, in the tau-nu plane
        -mtau*np.sin(theta) *
        (2*W1 - W2 + Enu / M * W3 - mtau**2/M**2 * W4 + Etau/M * W5)
        / F,
        # The longitudinal polarization component
        -(
            (2*W1 - mtau**2 / M**2 * W4) * (ptau - Etau * costheta)
            + W2 * (ptau + Etau * costheta)
            + W3 / M * ((Enu + Etau) * ptau - (Enu * Etau + ptau**2) * costheta)
            - mtau**2 / M * W5 * costheta
        ) / F
    )


def boost_rotated_4m(nutau4m: vector.MomentumNumpy4D, nucleon4m: vector.MomentumNumpy4D, tau4m: vector.MomentumNumpy4D, rotate_tau: bool = False):
    """Change reference frame to one where the nucleon is at rest, 
    the neutrino has a pure momentum in the +z direction, and the tau only has a momentum in the x and z directions.

    Args:
        nutau4m: The 4-momentum of the tau neutrino.
        nucleon4m: The 4-momentum of the nucleon.
        tau4m: The 4-momentum of the tau lepton.
        rotate_tau: Whether to rotate the tau lepton so that it only has 4-momentum components along the z and x axes.
            Strictly speaking, the tau polarization formula assumes that the tau 4-momentum is only in the z and x directions.
            However, based on the tests I have run, this does not seem to affect the polarization results.
            This is likely because I take the delta angle between the neutrino and the tau lepton, which should be invariant under rotations.
            Therefore, this parameter is set to False by default.

    Returns:
        A tuple of the 4-momenta of the tau neutrino, nucleon, and tau lepton in the new reference frame.

    Raises:
        ValueError: If the nucleon 4-momentum is 0. In this case, the boost does not work.
    """
    # Check if the nucleon 4-momentum is 0
    if np.any(nucleon4m.E == 0):
        raise ValueError("Nucleon 4-momentum is 0, cannot boost to that reference frame")
    
    # Boost all vectors such that the nucleon is at rest
    tau4m_boosted = tau4m.boostCM_of(nucleon4m)
    nucleon4m_boosted = nucleon4m.boostCM_of(nucleon4m)
    nutau4m_boosted = nutau4m.boostCM_of(nucleon4m)

    # Get the angles to rotate nutau4m_boosted so that it is aligned with the z axis, pointing towards positive z
    phi = nutau4m_boosted.phi
    theta = nutau4m_boosted.theta

    # Rotate nutau4m_boosted so that it is aligned with the z axis, pointing towards positive z
    nutau4m_rotated = nutau4m_boosted.rotateZ(-phi).rotateY(-theta)
    # Rotate all other momentum vectors, thereby preserving the 4-momentum conservation
    tau4m_rotated = tau4m_boosted.rotateZ(-phi).rotateY(-theta)
    nucleon4m_rotated = nucleon4m_boosted.rotateZ(-phi).rotateY(-theta)
    
    if not rotate_tau:
        return nutau4m_rotated, nucleon4m_rotated, tau4m_rotated
    
    tau_phi = tau4m_rotated.phi
    # Rotate nutau4m_boosted so that it is aligned with the z axis, pointing towards positive z
    nutau4m_rotated2 = nutau4m_rotated.rotateZ(-tau_phi)
    # Rotate all other momentum vectors, thereby preserving the 4-momentum conservation
    tau4m_rotated2 = tau4m_rotated.rotateZ(-tau_phi)
    nucleon4m_rotated2 = nucleon4m_rotated.rotateZ(-tau_phi)

    return nutau4m_rotated2, nucleon4m_rotated2, tau4m_rotated2


# Cabbibo angle, from https://en.wikipedia.org/wiki/Cabibbo%E2%80%93Kobayashi%E2%80%93Maskawa_matrix
theta_c = 13.02 * np.pi/180 # radians
# Vector mass, from Hagiwara et al.
M_V_qel = 0.84 # GeV
# difference between the anomalous magnetic moments of the proton anf neutron, from Hagiwara et al.
xi_qel = 3.706
# Axial vector mass, from Hagiwara et al.
M_A_qel = 1 # GeV
# Axial vector form factor at q^2 = 0, from Hagiwara et al.
F_A_0_qel = -1.23
# Pion mass
m_pi = 0.139 # GeV


def w(p, q, M):
    return p.dot(q)/M**2

def G_V_E(q):
    """From eq. 34 in Hagiwara et al. (2003)"""
    return 1 / (1-q**2/M_V_qel**2)**2

def G_V_M(q):
    """From eq. 34 in Hagiwara et al. (2003)"""
    return (1+xi_qel) / (1-q**2/M_V_qel**2)**2

def F_V(q, M):
    """From eq. 33 in Hagiwara et al. (2003)"""
    return (G_V_E(q) - q**2 / (4*M**2) * G_V_M(q)) / (1 - q**2 / (4*M**2))

def F_M(q, M):
    """From eq. 33 in Hagiwara et al. (2003)"""
    return (G_V_M(q) - G_V_E(q)) / (xi_qel * (1 - q**2 / (4*M**2)))

def F_A(q):
    """From eq. 35 in Hagiwara et al. (2003)"""
    return F_A_0_qel / (1 - q**2 / M_A_qel**2)**2

def F_T(q):
    """Equivalent to F_3^A in Hagiwara et al. (2003). See first paragraph, page 7."""
    return 0

def F_S(q):
    """Equivalent to F_3^V in Hagiwara et al. (2003). See first paragraph, page 7."""
    return 0

def F_p(q, M):
    """From eq. 36 in Hagiwara et al. (2003)"""
    return 2*M**2 * F_A(q) / (m_pi**2 - q**2)

def Wqel_coefficient(p, q, M):
    """The coefficient that coes in front of all the QEL structure functions. 
    From the last equation on page 4 in Kuzmin et al. (2003)"""
    return np.cos(theta_c)**2 * 1/w(p, q, M) 

def W1qel(x, Q2, p, q, M):
    """The first structure function for QEL interactions.
    Based on the last equation on p. 4 and the first equation on p. 5 in Kuzmin et al. (2003)"""
    xprime = Q2 / (4*M**2)
    return (
        Wqel_coefficient(p, q, M)
        * (F_A(q)**2 + xprime * (F_A(q)**2 + (F_V(q, M) + F_M(q, M))**2))
    )

def W2qel(x, Q2, p, q, M):
    """The second structure function for QEL interactions.
    Based on the last equation on p. 4 and the second equation on p. 5 in Kuzmin et al. (2003)"""
    xprime = Q2 / (4*M**2)
    return (
        Wqel_coefficient(p, q, M)
        * (F_V(q, M)**2 + F_A(q)**2 + xprime * (F_M(q, M)**2 + 4*F_T(q)**2))
    )

def W3qel(x, Q2, p, q, M):
    """The third structure function for QEL interactions.
    Based on the last equation on p. 4 and the third equation on p. 5 in Kuzmin et al. (2003)"""
    return (
        Wqel_coefficient(p, q, M)
        * -2 * np.real(np.conj(F_A(q)) * (F_V(q, M) + F_M(q, M)))
    )

def W4qel(x, Q2, p, q, M):
    """The fourth structure function for QEL interactions.
    Based on the last equation on p. 4 and the fourth equation on p. 5 in Kuzmin et al. (2003)"""
    xprime = Q2 / (4*M**2)
    return (
        Wqel_coefficient(p, q, M)
        # omega_4
        * (
            np.real(np.conj(F_V(q, M)) * (F_S(q) - 1/2 * F_M(q, M)) - np.conj(F_A(q)) * (F_T(q) + F_p(q, M)))
            + xprime * (1/2 * (F_M(q, M) - F_S(q))**2 + (F_T(q) + F_p(q, M))**2)
            - 1/4* (1+xprime) * F_M(q, M)**2 + (1+1/2 * xprime) * F_S(q)**2
        )
    )

def W5qel(x, Q2, p, q, M):
    """The fifth structure function for QEL interactions.
    Based on the last equation on p. 4 and the fifth equation on p. 5 in Kuzmin et al. (2003)"""
    xprime = Q2 / (4*M**2)
    return (
        Wqel_coefficient(p, q, M)
        * (2 * np.real(np.conj(F_S(q)) * (F_V(q, M) - xprime*F_M(q, M)) - np.conj(F_T(q)) * (F_A(q) - 2*xprime*F_p(q, M))))
        + W2qel(x, Q2, p, q, M)
    )

# Deep inelastic scattering
m_charm = 1.27 # GeV

def mass_corrected_xi(x, Q2, charm):
    """Correction that is applied to the Björken x variable when a charm quark is produced.
    From the last paragraph of p. 11 in Hagawira et al. (2003)."""
    xi = x.copy()
    # If a charm quark is produced, replace x with xi
    xi[charm] = x[charm] / (Q2[charm] / (Q2[charm] + m_charm**2))

    return xi

def Wndis_f(x, Q2, p, q, M, charm, fvalues, n):
    """Calculate the nth structure function using the nth form factor.
    
    Args:
        x (np.array): The Bjorken x values
        Q2 (np.array): The four-momentum transfer squared
        p (np.array): The four-momentum of the nucleon
        q (np.array): The four-momentum of the W boson (Q2 = -q^2)
        M (np.array): The mass of the nucleon
        charm (np.array): Ignored. Only included for backwards-compatibility
        fvalues (dict, pd.DataFrame): A mappable object containing the form factor values F1, F2, F3, F4, F5
        n (int): The structure function to calculate. Should be a number between 1 and 5.
    
    Returns:
        np.array: The values of the structure function W
    """
    if n == 1:
        xi = mass_corrected_xi(x, Q2, charm)
        return (1+xi/w(p, q, M)) * fvalues[f"F{n}"].values

    return 1/w(p, q, M) * fvalues[f"F{n}"].values


def calculate_pol_all():
    """Calculate the polarization vector for all simulated IceCube events."""
    neutrino_energies = (5, 10, 20, 50, 100)
    
    all_pols = {}
    all_event_infos = {}

    for e in neutrino_energies:
        # Load files
        particle_info = pd.read_csv(f"../data/test_genie_NuTau_{e}.0_GeV_particles.csv")
        event_info = pd.read_csv(f"../data/test_genie_NuTau_{e}.0_GeV_event_info_sig.csv")
    
        pols = calc_pol(particle_info, event_info)
        
        all_pols[e] = pols
        all_event_infos[e] = event_info

    return all_pols, all_event_infos
    

def calc_pol(particle_info, event_info):
    """Calculate the tau polarization for the given event information and particle information."""
    pols = pd.DataFrame() # columns=["polx", "poly", "polz", "event_num"]
        
    for int_type in ("qel", "res", "dis"):
        int_type_col = int_type

            # Select qel, res and dis particles
        events = event_info[event_info[int_type_col]]

        particles = particle_info[
                particle_info["event_num"]
                .isin(events.loc[events[int_type_col], "event_num"].values)
            ]

        selected_taus = particles[particles["pdg"] == 15]
        selected_nus = particles.groupby("event_num").nth(1)

        tau4m = vector.array({"E": selected_taus["E"], "px": selected_taus["px"], "py": selected_taus["py"], "pz": selected_taus["pz"]})
        nutau4m = vector.array({"E": selected_nus["E"], "px": selected_nus["px"], "py": selected_nus["py"], "pz": selected_nus["pz"]})

        x = events["xs"].values
        W = events["Ws"].values
        nucleon4m = vector.array({"E": events["En"], "px": events["pxn"], "py": events["pyn"], "pz": events["pzn"]})

        df = None
        res_sigs = None

        if int_type == "dis":
            df = events.loc[:, "F1":"F5"]
        elif int_type == "res":
            res_sigs = events.loc[:, ["sigmm", "sigpp", "sigmp"]]
            
            # Rotate all vectors such that the nucleon is at rest
        nutau4m_rotated, nucleon4m_rotated, tau4m_rotated = boost_rotated_4m(nutau4m, nucleon4m, tau4m)

        s = np.array(polarization_vector(
                tau4m_rotated,
                nutau4m_rotated,
                nucleon4m_rotated,
                int_type,
                x=x,
                W=W,
                charm=events["charm"].values,
                structure_function_values=df,
                res_sigs=res_sigs
            ))

            # Transfer the polarization vector to the lab frame
            # s[1] should be along the tau momentum direction.
            # Project the first component of the polarization vector onto the plane of the tau lepton and the neutrino
            # The second component should be in the tau neutrino-tau lepton plane, orthogonal to the tau and the tau-tau neutrino plane normal
        pol_l = vector.MomentumNumpy3D(tau4m) * s[1] / tau4m.p
        transverse_direction = vector.MomentumNumpy3D(nutau4m).cross(vector.MomentumNumpy3D(tau4m)).cross(vector.MomentumNumpy3D(tau4m))
        pol_t = transverse_direction * s[0] / transverse_direction.p
        pol_vec = pol_l + pol_t
            # Check that the longitudinal and transverse components are orthogonal
        assert np.allclose(pol_l.dot(pol_t), 0)
            # Check that the angle between the tau momentum and the polarization vector is the same as theta_p
        close_angles = np.isclose(tau4m.deltaangle(pol_vec), np.arccos(s[1] / np.linalg.norm(s, axis=0)))
        assert np.all(close_angles | (s[0] == 0)), (tau4m.deltaangle(pol_vec)[~close_angles], np.arccos(s[1] / np.linalg.norm(s, axis=0))[~close_angles])

        new_pols = pd.DataFrame({"event_num": selected_taus["event_num"], "polx": pol_vec.px, "poly": pol_vec.py, "polz": pol_vec.pz})
            
        # Concatenate the data with the pols dataframe
        pols = pd.concat([pols, new_pols])
        
        # Fill in the missing polarizations with 0, i.e. assume that they are unpolarized.
    events_not_in_pols = event_info.loc[~event_info["event_num"].isin(pols["event_num"].values), "event_num"].unique()
        
    # NOTE for other interactions (which in practice is just coherent scattering), we assume unpolarized (or fully left-handed is a better approximation perhaps?)
    pols = pd.concat([pols, pd.DataFrame({"event_num": events_not_in_pols, "polx": np.zeros(events_not_in_pols.shape[0]), "poly": np.zeros(events_not_in_pols.shape[0]), "polz": np.zeros(events_not_in_pols.shape[0])})])
    
    if events_not_in_pols.shape[0] > 0:
        print(events_not_in_pols.shape[0], "events are not qel, res or dis")
        
    print("Fraction of events with norm > 1:", (np.linalg.norm(pols[["polx", "poly", "polz"]], axis=1) > 1).sum() / pols.shape[0])
        
        # Fill up nan values. In practice, this should never happen
    if pols["polz"].isnull().any():
        print("Number of NaN values:", pols["polz"].isnull().sum())
        pols = pols.fillna(0)
            
    pols = pols.sort_values("event_num")
    
    return pols


def main():
    """Calculate tau polarization for one set of events, with the event and particle info file names given as input via the command line.
    Then create a new file containing the tau polarization information and tau leptons.
    """
    argparser = ArgumentParser()
    argparser.add_argument("-ip", type=str, required=True)  # Input particle info file name
    argparser.add_argument("-ie", type=str, required=True)  # Output particle info file name
    argparser.add_argument("-o", type=str, required=True)  # Output file name for csv file
    
    args = argparser.parse_args()

    particle_info = pd.read_csv(args.ip)
    event_info = pd.read_csv(args.ie)
    calc_pol(particle_info, event_info).to_csv(args.o, index=False)


if __name__ == "__main__":
    main()
