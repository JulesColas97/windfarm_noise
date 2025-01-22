import numpy as np
import matplotlib.pyplot as plt
import os
import os.path
import scipy.interpolate as interp
import scipy.special as spc
from func_timeout import func_timeout, FunctionTimedOut


def c_round(temp):
    return np.round(temp * 10 ** 6) / 10 ** 6


# write xfoil input file
def XFOIL_pre(inf:dict, foil:dict, path:str):
    """Create thye Xoil prompt to be run in order to obtain data for a given airfoile profile, alpha

    Args:
        inf (dict): dictionary containing the flow input
        foil (dict): dictionary containing thea airfoil input 
        path (str): ath to write and read the metadata
    """
    rho_inf = inf['rho']; nu_inf = inf['nu']; T_inf = inf['T']
    foil_name = foil['name']; foil_chord = foil['chord']; foil_re = foil['re']; foil_mach = foil['mach']; foil_aoa = foil['AoA']; foil_trip_top = foil['trip_top']; foil_trip_bottom = foil['trip_bottom']
    U_inf = foil_mach * (1.4 * 287 * T_inf) ** 0.5


     
    if os.path.isfile(path + "metadata1.dat"):
        os.remove(path + "metadata1.dat")
    if os.path.isfile(path + "metadata2.dat"):
        os.remove(path + "metadata2.dat")
    if os.path.isfile(path + "metadata3.dat"):
        os.remove(path + "metadata3.dat")
    if os.path.isfile("xfoil_driver.txt"):
        os.remove("xfoil_driver.txt")

    #for the development of the xfoil_driver file
    outfile  = open("xfoil_driver.txt", mode = "w")
    # outfile.write("PLOP\nG F\n\n")
    outfile.write("LOAD " + str(foil_name) + "\n")
    # outfile.write("63415\n")
    outfile.write("GDES\nGSET\nSCAL\n" + str(foil_chord) + "\nEXEC\nGSET\n\n")
    # outfile.write("PPAR\nN 200\n t 0.8\n\n\n")
    outfile.write("PPAR\nN 200\n\n\n")

    # 
    outfile.write("OPER\n")
    outfile.write("VPAR\n")
    
    # outfile.write("N 1\n\n")
    
    outfile.write("XTR\n%s\n%s\n\n"%(foil_trip_top,foil_trip_bottom))
    outfile.write("ITER\n%s\n"%(400))
    outfile.write("VISC\n%s\n"%(foil_re))
    outfile.write("M\n%s\n"%(foil_mach))

    #outfile.write("OPER\nRE " + str(foil_re) + "\nM " + str(foil_mach) + "\nITER 400\nVISC\nVPAR\nXTR\n" + str(foil_trip_top) + "\n" + str(foil_trip_bottom) + "\n\n")
    #outfile.write("OPER\nRE " + str(foil_re) + "\nM " + str(foil_mach) + "\nITER 400\nVISC\nVPAR\n\n")

    outfile.write("pacc\n")
    outfile.write(path+"metadata3.dat\n\n")
    # outfile.write("ALFA 1 \n")
    outfile.write("ALFA  1\n")
    outfile.write("ALFA " + str(foil_aoa) + "\n")
    outfile.write("pacc\n")
    outfile.write("DUMP " +path+"metadata1.dat\n")
    outfile.write("CPWR " + path+"metadata2.dat\n\n")
    outfile.write("QUIT\n")
    outfile.close()


def run_Xvfb(n=5):
    """  create a X virtual buffer to run xfoil wthout the frame popping

    Args:
        n (int, optional): token of the fram buffer (could be set to anything I think). Defaults to 5.
    """
    os.system("Xvfb :"+str(n)+" &")


def XFOIL_run(n=5):
    """Ryun xfoil inside the virtual buffer 

    Args:
        n (int, optional): must be set to the same as in  `run_Xvfb`. Defaults to 5.
    """
    os.system("DISPLAY=:"+str(5) +" xfoil < xfoil_driver.txt")


def kill_Xvfb():
    """kill all virtual buffers. 
    """
    os.system("kill $!")


def XFOIL_post(path:str):
    """read and create interpoland from boundary layer quantities writen in metadata1 and metadata2

    Args:
        path (str): path where to read BL quantities 

    Returns:
        tuple[dict,int]: dictionary of interpoland for BL quantities, success flag
    """
    # read metada1.dat
    if os.path.isfile(path + '/metadata1.dat'):
        run_flag = 1
        #now for reading the required data from the metadata1/2.dat files
        final_data = {}
        final_data['x'] = []; final_data['x_cp'] = []; final_data['cf'] = []; final_data['Ue'] = []; final_data['Dstar'] = []; final_data['theta'] = []; final_data['cp'] = []
        file_name = path + "metadata1.dat"
        infile = open(file_name, mode = "r+")
        data_in = infile.read()
        data_in = data_in.split("\n")
        # read data line by line
        for temp in range(len(data_in) - 2):
            temp_array = data_in[temp + 1].split(' ')
            repo = []
            for temp1 in range(len(temp_array)):
                if temp_array[temp1] != '':
                    repo.append(float(temp_array[temp1]))
            final_data['x'].append(repo[1])
            final_data['Ue'].append(repo[3])
            final_data['Dstar'].append(repo[4])
            final_data['theta'].append(repo[5])
            final_data['cf'].append(repo[6])
        infile.close()

        # read metadata 2
        file_name = path + "metadata2.dat"
        infile = open(file_name, mode = "r+")
        data_in = infile.read()
        data_in = data_in.split("\n")
        # read data line by line 
        for temp in range(len(data_in) - 4):
            temp_array = data_in[temp + 3].split(' ')
            repo = []
            for temp1 in range(len(temp_array)):
                if temp_array[temp1] != '':
                    repo.append(float(temp_array[temp1]))
            final_data['cp'].append(repo[-1])
            final_data['x_cp'].append(repo[0])
        infile.close()

        #for the fitting of splines between the stored data
        #for Dstar, theta and edge velocity
        # find leading edge index 
        for temp in range(len(final_data['x'])):
            if final_data['x'][temp] < final_data['x'][temp + 1]:
                break
        LE_index = temp
        tck_dict = {}

        # create BL interpoland 
        tck_dict['Ue_suction'] = interp.splrep(final_data['x'][: LE_index][::-1], final_data['Ue'][: LE_index][::-1], k = 2)
        tck_dict['Dstar_suction'] = interp.splrep(final_data['x'][: LE_index][::-1], final_data['Dstar'][: LE_index][::-1], k = 2)
        tck_dict['theta_suction'] = interp.splrep(final_data['x'][: LE_index][::-1], final_data['theta'][: LE_index][::-1], k = 2)
        tck_dict['Cf_suction'] = interp.splrep(final_data['x'][: LE_index][::-1], final_data['cf'][: LE_index][::-1], k = 2)

        tck_dict['Ue_pressure'] = interp.splrep(final_data['x'][LE_index + 1 : ], final_data['Ue'][LE_index + 1 : ], k = 2)
        tck_dict['Dstar_pressure'] = interp.splrep(final_data['x'][LE_index + 1 : ], final_data['Dstar'][LE_index + 1 : ], k = 2)
        tck_dict['theta_pressure'] = interp.splrep(final_data['x'][LE_index + 1 : ], final_data['theta'][LE_index + 1 : ], k = 2)
        tck_dict['Cf_pressure'] = interp.splrep(final_data['x'][LE_index + 1 : ], final_data['cf'][LE_index + 1 : ], k = 2)

        #for the coefficient of pressure
        # find leading edge index 
        for temp in range(len(final_data['x_cp'])):
            if final_data['x_cp'][temp] < final_data['x_cp'][temp + 1]:
                break
        LE_index_cp = temp
        tck_dict['cp_suction'] = interp.splrep(final_data['x_cp'][: LE_index_cp][::-1], final_data['cp'][: LE_index_cp][::-1], k = 2)
        tck_dict['cp_pressure'] = interp.splrep(final_data['x_cp'][LE_index_cp + 1 : ], final_data['cp'][LE_index_cp + 1 : ], k = 2)
    else:
        run_flag = 0
        print('No result generated!!!!')
        tck_dict = {}

    return tck_dict, run_flag


def XFOIL_post2(path:str)->np.ndarray:
    """read polar from metadata3.dat

    Args:
        path (str): path to metadata files

    Returns:
        np.ndarray: array containing AoA, CL and CD coeff 
    """
    file_name = path + 'metadata3.dat'
    infile = open(file_name, mode = "r+")
    data_in = infile.read()
    data_in = data_in.split("\n")
    
    row_offset = 13
    iline = 0
    aoa = []; Cl = []; Cd=[]
    while data_in[iline + row_offset] != '':
        temp_array = data_in[iline + row_offset].split(' ')
        ichar = 0
        aoa_flag = 0; Cl_flag = 0;Cd_flag=0
        while Cd_flag == 0:
            if temp_array[ichar] != '' and aoa_flag == 0:
                aoa.append(float(temp_array[ichar]))
                aoa_flag = 1
                ichar = ichar + 1
            elif temp_array[ichar] != '' and Cl_flag == 0 and aoa_flag ==1:
                Cl.append(float(temp_array[ichar]))
                Cl_flag = 1
                ichar = ichar + 1
            elif temp_array[ichar] != '' and Cl_flag == 1:
                Cd.append(float(temp_array[ichar]))
                Cd_flag = 1
                ichar = ichar + 1
            elif temp_array[ichar] == '':
                ichar = ichar + 1
        iline = iline + 1
    # tck_polar = interp.splrep(aoa, Cl, k = 2)
    # return tck_polar

    # remove metadata file (needed for newt run )
    return np.array(np.squeeze([aoa,Cl,Cd]))


def XFOIL_mp(inf:dict, foil:dict,path:str):
    """run xfoil, check if time limit is reached. if xfoil dont converge after 10s: break

    Args:
        inf (dict): dictionary containing the flow input
        foil (dict): dictionary containing thea airfoil input 
        path (str): ath to write and read the metadata
    """
    # create prompt to input in xfoil 
    XFOIL_pre(inf, foil,path)
    # run xfoil 
    try:
        doitReturnValue = func_timeout(10, XFOIL_run)
        print("XFOIL has converged satisfactorily!!")
    except FunctionTimedOut:
        os.system("TASKKILL /F /IM xfoil.exe")
        print("XFOIL has timed out!!!")


def XFOIL(path:str, target_frac:float, inf:dict, foil:dict):
    """create prompt, launch Xcode and return values of interest acording at the given target fraction of the chord (closest to the trailing edge)


    Args:
        path (str): path to write and read the metadata
        target_frac (float): fraction of the chord to measure the BL quantities
        inf (dict): dictionary containing the flow input
        foil (dict): dictionary containing thea airfoil input 

    Returns:
        tuple[dict,np.ndarray]: outdir contains boundary layer quantities, polar contains the lift and drag coeeficient 
    """
    # compute dynamic pressur ?? dont know why 
    dyn_press = 0.5 * inf['rho'] * (inf['mach'] * (1.4 * 287 * inf['T']) ** 0.5) ** 2
    out_dict = {}

    # create directory to save data file from xfoil 
    if not os.path.exists(path):
        os.makedirs(path)
    
    # run Xfoil, check if function timed out 
    XFOIL_mp(inf, foil,path)

    # read data and post process data from X foil 
    tck_dict, run_flag = XFOIL_post(path)
    # sample BL quantities at the given target_fraction of the chord 
    if run_flag == 1:
        Dstar_suction = c_round(interp.splev(target_frac * foil['chord'], tck_dict['Dstar_suction'], der = 0))
        theta_suction = c_round(interp.splev(target_frac * foil['chord'], tck_dict['theta_suction'], der = 0))
        Ue_suction = c_round(interp.splev(target_frac * foil['chord'], tck_dict['Ue_suction'], der = 0))
        Cf_suction = c_round(interp.splev(target_frac * foil['chord'], tck_dict['Cf_suction'], der = 0))
        dpdx_suction = c_round(interp.splev(target_frac * foil['chord'], tck_dict['cp_suction'], der = 1) * dyn_press)
        dpdx_suction1 = c_round((interp.splev(target_frac * foil['chord'], tck_dict['cp_suction'], der = 0) - interp.splev((target_frac - 0.05) * foil['chord'], tck_dict['cp_suction'], der = 0)) / (0.05 * foil['chord']) * dyn_press)

        Dstar_pressure = c_round(interp.splev(target_frac * foil['chord'], tck_dict['Dstar_pressure'], der = 0))
        theta_pressure = c_round(interp.splev(target_frac * foil['chord'], tck_dict['theta_pressure'], der = 0))
        Ue_pressure = c_round(interp.splev(target_frac * foil['chord'], tck_dict['Ue_pressure'], der = 0))
        Cf_pressure = c_round(interp.splev(target_frac * foil['chord'], tck_dict['Cf_pressure'], der = 0))
        dpdx_pressure = c_round(interp.splev(target_frac * foil['chord'], tck_dict['cp_pressure'], der = 1) * dyn_press)
        dpdx_pressure1 = c_round((interp.splev(target_frac * foil['chord'], tck_dict['cp_pressure'], der = 0) - interp.splev((target_frac - 0.05) * foil['chord'], tck_dict['cp_pressure'], der = 0)) / (0.05 * foil['chord']) * dyn_press)

    elif run_flag == 0:
        Dstar_suction = float("nan"); 
        theta_suction = float("nan"); 
        Ue_suction = float("nan"); 
        dpdx_suction = float("nan"); 
        Cf_suction = float('nan'); 
        dpdx_suction1 = float('nan')
        Dstar_pressure = float("nan"); 
        theta_pressure = float("nan"); 
        Ue_pressure = float("nan"); 
        dpdx_pressure = float("nan"); 
        Cf_pressure = float('nan'); 
        dpdx_pressure1 = float('nan')

    out_dict['Dstar_suction'] = Dstar_suction; 
    out_dict['theta_suction'] = theta_suction; 
    out_dict['Ue_suction'] = Ue_suction; 
    out_dict['dpdx_suction'] = dpdx_suction; 
    out_dict['Cf_suction'] = Cf_suction; 
    out_dict['dpdx_suction1'] = dpdx_suction1
    out_dict['Dstar_pressure'] = Dstar_pressure; 
    out_dict['theta_pressure'] = theta_pressure; 
    out_dict['Ue_pressure'] = Ue_pressure; 
    out_dict['dpdx_pressure'] = dpdx_pressure; 
    out_dict['Cf_pressure'] = Cf_pressure; 
    out_dict['dpdx_pressure1'] = dpdx_pressure1

    # post process polar files to read Cl, and Cd 
    polar = XFOIL_post2(path)
    
    return out_dict,polar



# use for computring lift and g=drag coef 
#---------------------------------------------------------------------------------
# def XFOIL_pre2(inf, foil, aoa_start, aoa_stop, aoa_delta):

#     rho_inf = inf['rho']; nu_inf = inf['nu']; T_inf = inf['T']
#     foil_name = foil['name']; foil_chord = foil['chord']; foil_re = foil['re']; foil_mach = foil['mach']; foil_trip_top = foil['trip_top']; foil_trip_bottom = foil['trip_bottom']
#     foil['chord'] = 1.0
#     U_inf = foil_mach * (1.4 * 287 * T_inf) ** 0.5
    
#     #for the development of the xfoil_driver file
#     outfile  = open("xfoil_driver2.txt", mode = "w")
#     outfile.write("LOAD " + str(foil_name) + "\n")
#     # outfile.write("GDES\nSCAL\n" + str(foil_chord) + "\nEXEC\n\n")
#     outfile.write("PPAR\nN 450\n\n\n")
#     outfile.write("OPER\nRE " + str(foil_re) + "\nM " + str(foil_mach) + "\nITER 100\nVISC\nVPAR\nXTR\n" + str(foil_trip_top) + "\n" + str(foil_trip_bottom) + "\n\n")
#     outfile.write("pacc\n")
#     outfile.write("polar_data.txt\n\n")
#     outfile.write("aseq " + str(aoa_start) + " " + str(aoa_stop) + " " + str(aoa_delta) + "\n")
#     outfile.write("pacc\n")
#     outfile.write("QUIT\n")
#     outfile.close()

# def XFOIL_run2():
#     os.system("xfoil < xfoil_driver2.txt")
    

# def XFOIL_mp2(inf, foil, aoa_start, aoa_stop, aoa_delta):
#     XFOIL_pre2(inf, foil, aoa_start, aoa_stop, aoa_delta)
#     try:
#         doitReturnValue = func_timeout(50, XFOIL_run2)
#         print("XFOIL has converged satisfactorily!!")
#     except FunctionTimedOut:
#         os.system("TASKKILL /F /IM xfoil.exe")
#         print("XFOIL has timed out!!!")

# def XFOIL2(path, target_frac, inf, foil, aoa_start, aoa_stop, aoa_delta):
#     XFOIL_mp2(inf, foil, aoa_start, aoa_stop, aoa_delta)
#     tck_polar = XFOIL_post2(path)
#     os.remove('xfoil_driver2.txt')
#     os.remove('polar_data.txt')
#     return tck_polar
