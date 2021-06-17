import os
from os import path
from dymola.dymola_interface import DymolaInterface
import multiprocessing as mp
import time
import re


def simulation1(case_no, pa, mo_name):
    dymola = DymolaInterface()
    dir_aixlib = 'E:\\Aixlib'
    dymola.openModel(path=os.path.join(dir_aixlib, 'package.mo'))
    dir_result = pa + '\\case' + str(case_no)
    folder_name = re.findall('^[^.]+', mo_name)
    dir_cases = pa + '\\case' + str(case_no) + '\\' + folder_name[0]
    dymola.openModel(path=os.path.join(dir_cases, 'package.mo'))
    print('Translating model', case_no)
    dymola.translateModel(mo_name)
    print('Translating Done, Running Simulation of model', case_no)
    output = dymola.simulateExtendedModel(
        problem=mo_name,
        startTime=0.0,
        stopTime=3.1536e+07,
        outputInterval=3600,
        method="Dassl",
        tolerance=0.0001,
        resultFile=os.path.join(dir_result, 'demo_results'),
        finalNames=['multizone.VAir'],
    )
    print('Simulation complete of', case_no)
    dymola.clear()
    dymola.close()


def simulation_exist():
    non_simulated_cases = []
    case_no = 1
    while case_no <= len(os.listdir('E:\\yash\\MYTOOLCHAIN\\new_export_institute_merged')):
        dir_result = 'E:\\yash\\MYTOOLCHAIN\\new_export_institute_merged\\case' + str(case_no)
        path_model = os.path.join(dir_result, 'demo_results.mat')
        boolean = path.exists(path_model)
        if boolean is False:
            non_simulated_cases.append(case_no)
        case_no = case_no + 1
    print(non_simulated_cases)


def run_simulation(path1, sim_no_sim, mod_name):
    start_time = time.time()
    i = 1
    num = int(sim_no_sim)
    k = []
    while i <= len(os.listdir(path1)):
        while i <= num:
            p = mp.Process(target=simulation1, args=(i, path1, mod_name,))
            k.append(p)
            p.start()
            if i < num:
                time.sleep(5)
            i = i + 1
        for l in k:
            l.join()
        print(time.time() - start_time)
        start_time = time.time()
        num = num + int(sim_no_sim)
        if num > len(os.listdir(path1)):
            num = len(os.listdir(path1))

    print('Finish')


if __name__ == '__main__':
    path2 = 'E:/yash/MYTOOLCHAIN/high_level_institute'
    model_name = 'ProjektBuerogebaeude.Buerogebaeude.Buerogebaeude'
    sim_no = 5
    run_simulation(path2, sim_no, model_name)

    # Check which simulations are not done
    # simulation_exist()
