import numpy as np

# definition of the subdomain for multi processing
def domain_facture(mesh, Ncore: int) -> list:
    """
    Divide the computational mesh into subdomains for multiprocessing.

    Args:
        mesh (Mesh): The computational mesh containing coordinate arrays.
        Ncore (int): The number of cores to divide the domain.

    Returns:
        list: A list of dictionaries, each containing a subdomain with 
              'x_coord', 'z_coord', and 'tau_coord' arrays.
    """
    domain_length = mesh.x_coord.shape[0] * mesh.x_coord.shape[1]
    delta = domain_length // Ncore
    final_coord = []
    for temp in range(Ncore - 1):
        temp_dict = {}
        temp_dict['x_coord'] = mesh.x_coord[temp * delta : (temp + 1) * delta]
        temp_dict['z_coord'] = mesh.z_coord[temp * delta : (temp + 1) * delta]
        temp_dict['tau_coord'] = mesh.tau_coord[temp * delta : (temp + 1) * delta]
        final_coord.append(temp_dict)

    temp_dict = {}
    temp_dict['x_coord'] = mesh.x_coord[(Ncore-1) * delta : ]
    temp_dict['z_coord'] = mesh.z_coord[(Ncore-1) * delta : ]
    temp_dict['tau_coord'] = mesh.tau_coord[(Ncore-1) * delta : ]
    final_coord.append(temp_dict)

    return final_coord

#function for reconstructing the output data from the mp modules
def output_construct(output: list, mesh, Ncore: int, Nfreq: int, Nblade: int, Nseg: int):
    """
    Reconstruct output data from multiple processing subdomains into a unified dataset.

    Args:
        output (list): List of dictionaries containing computed data for each subdomain.
        mesh (Mesh): The computational mesh used in the calculations.
        Ncore (int): Number of processing cores used.
        Nfreq (int): Number of frequency bands.
        Nblade (int): Number of wind turbine blades.
        Nseg (int): Number of blade segments.

    Returns:
        tuple: A tuple containing reconstructed arrays:
            - SPL_tot (np.array): Total Sound Pressure Level.
            - SPL_TIN (np.array): Tonal noise SPL.
            - SPL_TEN (np.array): Broadband noise SPL.
            - SWL_tot (np.array): Sound power levels.
            - AoA (np.array): Angle of Attack values.
            - U_inf (np.array): Freestream wind velocity.
            - epsilon (np.array): Turbulence dissipation rate.
            - U_rot (np.array): Rotational wind velocity.
            - U_rel (np.array): Relative wind velocity.
            - a (np.array): Axial induction factor.
            - adash (np.array): Tangential induction factor.
    """
    x1 = mesh.x_coord.shape[0]
    x2 = mesh.x_coord.shape[1]
    print(x1,x2)
    delta = x1 * x2 // Ncore
    SPL_tot = np.zeros((x1 * x2, Nseg, Nblade, Nfreq))
    SPL_TIN = np.zeros((x1 * x2, Nseg, Nblade, Nfreq))
    SPL_TEN = np.zeros((x1 * x2, Nseg, Nblade, Nfreq))
    SWL_tot = np.zeros((x1 * x2, Nseg, Nblade, Nfreq))
    for temp in range(Ncore - 1):
        SPL_tot[temp * delta : (temp + 1) * delta, :, :, :] = output[temp]['SPL_tot']
        SPL_TIN[temp * delta : (temp + 1) * delta, :, :, :] = output[temp]['SPL_TIN']
        SPL_TEN[temp * delta : (temp + 1) * delta, :, :, :] = output[temp]['SPL_TEN']
        SWL_tot[temp * delta : (temp + 1) * delta, :, :, :] = output[temp]['SWL_tot']
        

    SPL_tot[(temp + 1) * delta : , :, :, :] = output[temp + 1]['SPL_tot']
    SPL_TIN[(temp + 1) * delta : , :, :, :] = output[temp + 1]['SPL_TIN']
    SPL_TEN[(temp + 1) * delta : , :, :, :] = output[temp + 1]['SPL_TEN']
    SWL_tot[(temp + 1) * delta : , :, :, :] = output[temp + 1]['SWL_tot']
    SPL_tot = SPL_tot.reshape(x1 * x2, Nseg*Nblade, Nfreq)
    SPL_TIN = SPL_TIN.reshape(x1 * x2, Nseg*Nblade, Nfreq)
    SPL_TEN = SPL_TEN.reshape(x1 * x2, Nseg*Nblade, Nfreq)
    SWL_tot = SWL_tot.reshape(x1 * x2, Nseg*Nblade, Nfreq)
    AoA = output[0]['AoA'].reshape(Nseg*Nblade)
    U_inf = output[0]['U_inf'].reshape(Nseg*Nblade)

    epsilon = output[0]['epsilon'].reshape(Nseg*Nblade)
    U_rot = output[0]['U_rot'].reshape(Nseg*Nblade)
    U_rel = output[0]['U_rel'].reshape(Nseg*Nblade)
    a = output[0]['a'].reshape(Nseg*Nblade)
    adash = output[0]['adash'].reshape(Nseg*Nblade)

    return SPL_tot, SPL_TIN,SPL_TEN,SWL_tot,AoA,U_inf, epsilon, U_rot,U_rel,a,adash


def reshuffle_output(output_total: list, Ncore: int) -> list:
    """
    Reorder the output from parallel processing based on reconstruction indices.

    Args:
        output_total (list): List of dictionaries containing output data from processing.
        Ncore (int): Number of processing cores.

    Returns:
        list: The reordered list of output data sorted by their reconstruction indices.
    """
    order = []
    for temp in range(Ncore):
        order.append(output_total[temp]['recon_index'])
    order_index = np.argsort(order)
    final_output = []
    for temp in range(Ncore):
        final_output.append(output_total[order_index[temp]])
    return final_output
