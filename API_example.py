from scipy.constants import e
import sys
sys.path.append('/Users/fosterd/Release/bin/mac/gcc_Catalina_10.15.7/python3')
import os
import pyfispact as pf
import pandas as pd 

#
#  This is a FISPACT-II example irradiating a sample of Co and Li with 20MeV
#  protons. Identifying how much Ni57 is produced.
#
    
# callback for compute
def computefunc(p, index, total):
    print(" [{}/{}]  processing {}".format(index, total, p))

# Load nuclear data 
def load_nuclear_data(m):
    print("Reading nuclear data")
    nd_base_path = '/Users/fosterd/nuclear_data'
    nd = pf.NuclearData(m)
    ndr = pf.io.NuclearDataReader(m)
    nd.setprojectile(pf.PROJECTILE_PROTON() )
    ndr.setpath(pf.io.ND_IND_NUC_KEY(),  os.path.join(nd_base_path, 'TENDL2017data', 'tendl17_decay12_index'))
    ndr.setpath(pf.io.ND_XS_ENDF_KEY(),  os.path.join(nd_base_path, 'TENDL2017data', 'tal2017-p', 'gxs-162'))
    ndr.setpath(pf.io.ND_DK_ENDF_KEY(),  os.path.join(nd_base_path, 'decay', 'decay_2012'))
    ndr.setpath(pf.io.ND_ABSORP_KEY(),   os.path.join(nd_base_path, 'decay', 'abs_2012'))
    ndr.load(nd)

    return nd

# Set input
def set_input(m, flux_data, Proton_Flux,  ip, cooling_times_hours):
    ip.setname("Ra226")
    # Change from FISPACT input to API
    flux_data.reverse()
    ip.setflux(pf.groups.G162(), flux_data)
    # Name the flux
    ip.setfluxname("162 Proton flux")

    ip.setdensity(5.5)

    hours_to_seconds  = 3600

    # Heating time
    ip.appendschedule(14*24*hours_to_seconds, Proton_Flux)
    
    # Cooling time
    
    for t in cooling_times_hours:
        ip.appendschedule(t*hours_to_seconds, 0.0)



    # extract activity = generic routine
def extract_activity_from_nuc(o, inv_index , nuclide_zai ):
    activity = -1.0
    # check if nuclide present at a given step
    if (o.findinventoryexists(inv_index , nuclide_zai )): # get index of nuclide in inventory nuclide l i s t
        nuclide_index = o.findinventoryindex(inv_index , nuclide_zai)
        # get nuclide list at given step
        nuclide_list = o.getinventorynuclides(inv_index)
        # get heating value for that nuclide from inventory
        activity = nuclide_list[ nuclide_index ].activity 
        return activity
    
def extract_mass_from_nuc(o, inv_index , nuclide_zai ):
    mass = -1.0
    # check if nuclide present at a given step
    if (o.findinventoryexists(inv_index , nuclide_zai )): # get index of nuclide in inventory nuclide l i s t
        nuclide_index = o.findinventoryindex(inv_index , nuclide_zai)
        # get nuclide list at given step
        nuclide_list = o.getinventorynuclides(inv_index)
        # get heating value for that nuclide from inventory
        mass = nuclide_list[ nuclide_index ].grams 
        return mass

if __name__ == "__main__":
    runname = "Example_py"
    # Proton current 750 u Amps
    Proton_Current = 400.0e-6
    Proton_Flux = Proton_Current / e
    print("A proton current of ", Proton_Current,
          " Amps, is the same as ",Proton_Flux," protons per second")

    # Create a 19MeV flux
    # set flux to zero array
    flux = [0.0] * 162
    # find 19MeV bin and set to 1
    for i in pf.groups.g162():
        if i == 20E6:
            flux[pf.groups.g162().index(i)] = 1.0
            print('Got bin', i, pf.groups.g162().index(i))
            break

    m = pf.Monitor()
    # Initialise and instance of FISPACT
    pf.initialise(m)

    print("End initialise")

    print("Setting input")
    # Set FISPACT input data object
    ip = pf.InputData(m)
    cooling_times_hours = [24, 25, 26, 27, 24, 24, 24, 24, 24, 24]
    set_input(m, flux, Proton_Flux, ip, cooling_times_hours)
    nd = load_nuclear_data(m)
    print("FISPACT-II running")
    #  Define the output data object
    o = pf.OutputData(m)
    target_isotope_name = 'Ga69'
    nuc_zai_target_iso = pf.util.zai_from_name(m, target_isotope_name)

    ip.setmasstotal(1.0)
    elements = ['Co', 'Li']
    percentages = [50, 50]
    elements_z = [pf.util.z_from_element(m, element) for element in elements]
    print(elements, elements_z)
    ip.setmass(elements_z, percentages)

    # Call FISPACT
    pf.process(ip, nd, o, m, op=computefunc)

    name = 'Ni57'
    nuc_zai = pf.util.zai_from_name(m, name)
    final_mass_g = extract_mass_from_nuc(o, 1, nuc_zai)
    print('mass of ',name, ' after irradiation', final_mass_g, ' g')

    pf.io.to_file(o, m, "{}.json".format(runname))
    nuclei_name = []
    original_mass = []
    n_mass = []
    p_mass = []
    zai, sort_mass = o.getsortedinventory(1, pf.INVENTORY_TOTAL_MASS())
    for nuc in zai:
        name = pf.util.nuclide_from_zai(m, nuc)
        mass0 = extract_mass_from_nuc(o, 0, nuc)
        mass_p = extract_mass_from_nuc(o, 1, nuc)
        if mass_p > 1e-6:
            nuclei_name.append(name)
            original_mass.append(mass0)
            p_mass.append(mass_p)
    material_list = {'nuceli': nuclei_name,
                    'original_mass': original_mass,
                    'p_mass': p_mass,
                    'p activity' : o.getinventoryvalue(1, pf.INVENTORY_TOTAL_ACTIVITY())
                    }
    df = pd.DataFrame(material_list)
    print(df.to_latex(index=False))
    df.to_excel("masses2.xlsx")
