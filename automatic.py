# Influence of parameters on building simulations of
# Building Information Modelling (BIM) models
# using Reduced Order Model (ROM) approach

import os
import pickle
from morris_sobol import SAMorris, SASobol


def run(comp_path):
    # Function for selecting method of sensitivity analysis
    while True:
        sa_method = input("Which method for Sensitivity Analysis: \n1   Morris\n2   Sobol\n")
        if int(sa_method) == 1:
            print("Morris Selected")
            morris(comp_path)
            break
        elif int(sa_method) == 2:
            print("Sobol selected")
            sobol(comp_path)
            break
        else:
            print("Wrong Entry")


def problem_morris(complete_path):
    # Problem defined for Morris method
    pickle_read = open(complete_path, "rb")
    prj = pickle.load(pickle_read)
    trajectories = 10
    no_variables = 17
    num_models = trajectories*(no_variables + 1)
    prob = {
        'num_vars': no_variables,
        'names': ['infiltration_rate',
                  'person_per_m2',
                  'lighting_power',
                  'machines',
                  'Heating_Set_Point',
                  'Cooling_Set_Point',
                  'inner_walls_thickness',
                  'inner_walls_thermal_conduc',
                  'inner_walls_heat_capac',
                  'inner_walls_density',
                  'outer_walls_thickness',
                  'outer_walls_thermal_conduc',
                  'outer_walls_solar_absorp',
                  'outer_walls_heat_capac',
                  'outer_walls_density',
                  'window_SHGC',
                  'window_Uvalue'],
        'bounds': [[0.050, 0.500], [0.05, 0.2], [2.00, 10.00], [2.00, 10.00], [293.15, 298.15],
                   [303.15, 308.15], [0.015, 0.240], [0.03, 0.045], [0.8, 0.9], [40.0, 100.0],
                   [0.025, 0.120], [0.025, 0.04], [0.4, 0.6], [1.3, 1.7], [28.0, 45.0],
                   [0.33, 0.86], [0.75, 5.4]]
    }
    return prj, trajectories, prob, num_models


def problem_sobol(complete_path):
    # Problem defined for Sobol method
    pickle_read = open(complete_path, "rb")
    prj = pickle.load(pickle_read)
    no_variables = 8
    num_models = (no_variables + 2) * 500
    prob = {
        'num_vars': no_variables,
        'names': ['infiltration_rate',
                  'person_per_m2',
                  'lighting_power',
                  'Heating_Set_Point',
                  'Cooling_Set_Point',
                  'outer_walls_thickness',
                  'outer_walls_thermal_conduc',
                  'window_SHGC'],
        'bounds': [[0.050, 0.500], [0.05, 0.2], [2.00, 10.00], [293.15, 298.15],
                   [303.15, 308.15],
                   [0.025, 0.120], [0.025, 0.04],
                   [0.33, 0.86]]
    }
    return prj, prob, num_models


def morris(complete_path):
    # Methodology for Morris sensitivity analysis
    prj, num_trajectories, prob, no_models = problem_morris(complete_path)
    obj = SAMorris()
    samples = obj.sampling(problem=prob, trajectories=num_trajectories)
    export_path = input('Enter path where Dymola models need to be exported : ')
    obj.export(prj=prj, x=samples, ex_path=export_path)
    obj.simulation(export_path, no_models)
    obj.analysis(export_path, no_models, samples, prob)


def sobol(complete_path):
    # Methodology for Sobol analysis
    prj, prob, no_models = problem_sobol(complete_path)
    obj = SASobol()
    samples = obj.sampling(problem=prob)
    pkl_sample = open("problem_sample.pkl", "wb")
    pickle.dump(samples, pkl_sample)
    pkl_sample.close()
    export_path = input('Enter path where Dymola models need to be exported : ')
    obj.export(prj=prj, x=samples, ex_path=export_path)
    obj.simulation(export_path, no_models)
    obj.analysis(export_path, no_models, prob)


# Initiation of script
if __name__ == '__main__':
    pkl_path = input("Enter full path of pickle model :")
    pkl_name = input("Enter name of pickle model with extension :")
    complet_path = pkl_path + '\\' + pkl_name
    if os.path.exists(complet_path):
        print("File exists")
        run(complet_path)
    else:
        print("File not exists")
