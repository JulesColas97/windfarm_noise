import numpy as np

# definition of the subdomain for multi processing
def domain_facture(mesh, Ncore):
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
def output_construct(output, mesh, Ncore, Nfreq, Nblade, Nseg):
    x1 = mesh.x_coord.shape[0]
    x2 = mesh.x_coord.shape[1]
    print(x1,x2)
    delta = x1 * x2 // Ncore
    #OASPL = np.zeros((x1 * x2, 1))
    #OASPL_beta = np.zeros((x1 * x2, Nblade))
    #AM_rec = np.zeros((x1 * x2, 1))
    #SPL_freq = np.zeros((x1 * x2, Nfreq))
    #SWL_freq = np.zeros((x1 * x2, Nfreq))
    SPL_tot = np.zeros((x1 * x2, Nseg, Nblade, Nfreq))
    SPL_TIN = np.zeros((x1 * x2, Nseg, Nblade, Nfreq))
    SPL_TEN = np.zeros((x1 * x2, Nseg, Nblade, Nfreq))
    SWL_tot = np.zeros((x1 * x2, Nseg, Nblade, Nfreq))
    #Spp_tot = np.zeros((x1 * x2, Nseg, Nblade, Nfreq))
    #Spp_TIN = np.zeros((x1 * x2, Nseg, Nblade, Nfreq))
    #Spp_TEN = np.zeros((x1 * x2, Nseg, Nblade, Nfreq))
    for temp in range(Ncore - 1):
        #OASPL[temp * delta : (temp + 1) * delta, :] = output[temp]['OASPL']
        #OASPL_beta[temp * delta : (temp + 1) * delta, :] = output[temp]['OASPL_beta']
        #AM_rec[temp * delta : (temp + 1) * delta, :] = output[temp]['AM_rec']
        #SPL_freq[temp * delta : (temp + 1) * delta, :] = output[temp]['SPL_freq']
        #SWL_freq[temp * delta : (temp + 1) * delta, :] = output[temp]['SWL_freq']
        SPL_tot[temp * delta : (temp + 1) * delta, :, :, :] = output[temp]['SPL_tot']
        SPL_TIN[temp * delta : (temp + 1) * delta, :, :, :] = output[temp]['SPL_TIN']
        SPL_TEN[temp * delta : (temp + 1) * delta, :, :, :] = output[temp]['SPL_TEN']
        SWL_tot[temp * delta : (temp + 1) * delta, :, :, :] = output[temp]['SWL_tot']
        
        #Spp_tot[temp * delta : (temp + 1) * delta, :, :, :] = output[temp]['Spp_tot'] 
        #Spp_TIN[temp * delta : (temp + 1) * delta, :, :, :] = output[temp]['Spp_TIN']
        #Spp_TEN[temp * delta : (temp + 1) * delta, :, :, :] = output[temp]['Spp_TEN']

    #OASPL[(temp + 1) * delta : , :] = output[temp + 1]['OASPL']
    #OASPL_beta[(temp + 1) * delta : , :] = output[temp + 1]['OASPL_beta']
    #AM_rec[(temp + 1) * delta : , :] = output[temp + 1]['AM_rec']
    #SPL_freq[(temp + 1) * delta : , :] = output[temp + 1]['SPL_freq']
    #SWL_freq[(temp + 1) * delta : , :] = output[temp + 1]['SWL_freq']
    SPL_tot[(temp + 1) * delta : , :, :, :] = output[temp + 1]['SPL_tot']
    SPL_TIN[(temp + 1) * delta : , :, :, :] = output[temp + 1]['SPL_TIN']
    SPL_TEN[(temp + 1) * delta : , :, :, :] = output[temp + 1]['SPL_TEN']
    SWL_tot[(temp + 1) * delta : , :, :, :] = output[temp + 1]['SWL_tot']
    #Spp_tot[(temp + 1) * delta : , :, :, :] = output[temp + 1]['Spp_tot']
    #Spp_TIN[(temp + 1) * delta : , :, :, :] = output[temp + 1]['Spp_TIN']
    #Spp_TEN[(temp + 1) * delta : , :, :, :] = output[temp + 1]['Spp_TEN']

    # OASPL = OASPL.reshape(x1, x2, 1)
    # OASPL_beta = OASPL_beta.reshape(x1, x2, Nblade)
    # AM_rec = AM_rec.reshape(x1, x2)
    # SPL_freq = SPL_freq.reshape(x1, x2, Nfreq)
    # SWL_freq = SWL_freq.reshape(x1, x2, Nfreq)
    SPL_tot = SPL_tot.reshape(x1 * x2, Nseg*Nblade, Nfreq)
    SPL_TIN = SPL_TIN.reshape(x1 * x2, Nseg*Nblade, Nfreq)
    SPL_TEN = SPL_TEN.reshape(x1 * x2, Nseg*Nblade, Nfreq)
    SWL_tot = SWL_tot.reshape(x1 * x2, Nseg*Nblade, Nfreq)
    # Spp_tot = Spp_tot.reshape(x1 * x2, Nseg*Nblade, Nfreq)
    # Spp_TIN = Spp_TIN.reshape(x1 * x2, Nseg*Nblade, Nfreq)
    # Spp_TEN = Spp_TEN.reshape(x1 * x2, Nseg*Nblade, Nfreq)
    AoA = output[0]['AoA'].reshape(Nseg*Nblade)
    U_inf = output[0]['U_inf'].reshape(Nseg*Nblade)

    epsilon = output[0]['epsilon'].reshape(Nseg*Nblade)
    U_rot = output[0]['U_rot'].reshape(Nseg*Nblade)
    U_rel = output[0]['U_rel'].reshape(Nseg*Nblade)
    a = output[0]['a'].reshape(Nseg*Nblade)
    adash = output[0]['adash'].reshape(Nseg*Nblade)

    # return OASPL, OASPL_beta, AM_rec, SPL_freq, SWL_freq, SPL_tot,Spp_tot,Spp_TIN,Spp_TEN,beta
    return SPL_tot, SPL_TIN,SPL_TEN,SWL_tot,AoA,U_inf, epsilon, U_rot,U_rel,a,adash


# reshuffle output from multiprocessing acording to recon index
def reshuffle_output(output_total, Ncore):
    order = []
    for temp in range(Ncore):
        order.append(output_total[temp]['recon_index'])
    order_index = np.argsort(order)
    final_output = []
    for temp in range(Ncore):
        final_output.append(output_total[order_index[temp]])
    return final_output
