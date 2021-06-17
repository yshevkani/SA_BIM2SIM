import os
from SALib.sample.morris import morris as morris_sample
from SALib.sample import saltelli as salt
from SALib.analyze import morris as morris_analyze
from SALib.analyze import sobol as sobol_analyze
from multiprocessing_simulation import run_simulation
from dymola.dymola_interface import DymolaInterface
import numpy as np
from modelicares import SimRes
from ebcpy_modelica_simres import to_pandas
import re


class SAMorris:
    # Class defined for Morris analysis
    def sampling(self, problem, trajectories):
        # Sampling of parameters
        samples = morris_sample.sample(problem, trajectories, num_levels=4)
        return samples

    def export(self, prj, x, ex_path):
        # Export of models for all the samples in modelica format using TEASER
        self.prj = prj
        ref_no = 1
        self.x = x

        for [infiltration, per_m2, l_p, mach, set_temp_heat, set_temp_cool,
             iwal_thick, iwal_cond, iwal_heatcap, iwal_density,
             owal_thick, owal_cond, owal_solarabs, owal_heatcap, owal_density,
             shgc_win, u_window] in self.x:
            path = ex_path + '\\case' + str(ref_no)
            os.mkdir(path=path)
            zone = 0
            while zone < len(self.prj.buildings[0].thermal_zones):
                    self.prj.buildings[0].thermal_zones[zone].use_conditions.infiltration_rate = infiltration
                    self.prj.buildings[0].thermal_zones[zone].use_conditions.persons = per_m2
                    self.prj.buildings[0].thermal_zones[zone].use_conditions.lighting_power = l_p
                    self.prj.buildings[0].thermal_zones[zone].use_conditions.machines = mach
                    self.prj.buildings[0].thermal_zones[zone].use_conditions.CoolerOn = True

                    hp = []
                    hp_index = 0
                    while hp_index <= 24:
                        hp.append(set_temp_heat)
                        hp_index = hp_index + 1
                    self.prj.buildings[0].thermal_zones[zone].use_conditions.heating_profile = hp

                    cp = []
                    cp_index = 0
                    while cp_index <= 24:
                        cp.append(set_temp_cool)
                        cp_index = cp_index + 1
                    self.prj.buildings[0].thermal_zones[zone].use_conditions.cooling_profile = cp

                    iwal = 0
                    while iwal < len(self.prj.buildings[0].thermal_zones[zone].inner_walls):
                        self.prj.buildings[0].thermal_zones[zone].inner_walls[iwal].layer[1].thickness = iwal_thick
                        self.prj.buildings[0].thermal_zones[zone].inner_walls[iwal].layer[1].material.thermal_conduc\
                            = iwal_cond
                        self.prj.buildings[0].thermal_zones[zone].inner_walls[iwal].layer[1].material.heat_capac\
                            = iwal_heatcap
                        self.prj.buildings[0].thermal_zones[zone].inner_walls[iwal].layer[1].material.density\
                            = iwal_density
                        iwal = iwal + 1

                    owal = 0
                    while owal < len(self.prj.buildings[0].thermal_zones[zone].outer_walls):
                        if len(self.prj.buildings[0].thermal_zones[zone].outer_walls[owal].layer) == 2:
                            self.prj.buildings[0].thermal_zones[zone].outer_walls[owal].layer[-1].thickness = owal_thick
                            self.prj.buildings[0].thermal_zones[zone].outer_walls[owal].layer[-1].material.thermal_conduc \
                                = owal_cond
                            self.prj.buildings[0].thermal_zones[zone].outer_walls[owal].layer[-1].material.heat_capac \
                                = owal_heatcap
                            self.prj.buildings[0].thermal_zones[zone].outer_walls[owal].layer[-1].material.density \
                                = owal_density
                        else:
                            self.prj.buildings[0].thermal_zones[zone].outer_walls[owal].layer[2].thickness = owal_thick
                            self.prj.buildings[0].thermal_zones[zone].outer_walls[owal].layer[2].material.thermal_conduc \
                                = owal_cond
                            self.prj.buildings[0].thermal_zones[zone].outer_walls[owal].layer[2].material.heat_capac\
                                = owal_heatcap
                            self.prj.buildings[0].thermal_zones[zone].outer_walls[owal].layer[2].material.density\
                                = owal_density
                        owal = owal + 1

                    win = 0
                    while win < len(self.prj.buildings[0].thermal_zones[zone].windows):
                        self.prj.buildings[0].thermal_zones[zone].windows[win].g_value = \
                            shgc_win
                        win = win + 1

                    self.prj.buildings[0].thermal_zones[zone].use_conditions.use_constant_infiltration = True

                    self.prj.buildings[0].thermal_zones[zone].calc_zone_parameters()
                    zone = zone + 1
            ref_no = ref_no + 1
            # prj.buildings[0].internal_gains_mode = 2
            self.prj.calc_all_buildings()
            zone = 0
            while zone < len(self.prj.buildings[0].thermal_zones):
                self.prj.buildings[0].thermal_zones[zone].model_attr.u_value_win = u_window
                self.prj.buildings[0].thermal_zones[zone].model_attr.solar_absorp_ow = owal_solarabs
                self.prj.buildings[0].thermal_zones[zone].model_attr.cool_load = -100000
                zone = zone + 1
            self.prj.weather_file_path = 'E:\\Bim2Sim\\bim2sim-coding\\teaser\\data\\input\\inputdata\\' \
                                    'weatherdata\\DEU_BW_Mannheim_107290_TRY2010_12_Jahr_BBSR.mos'

            self.prj.export_aixlib(path=path)

    def simulation(self, exp_path, num_models):
        # Function to choose one at a time or multiprocessing simulation
        multi_proc = input("If you want multiprocessing for simulations, Enter 1 else 0 : ")
        model_name = input(
            "Enter model name in dymola (which has to be simulated) : \neg. ProjektBuerogebaeude.Buerogebaeude.Buerogebaeude\n")
        if int(multi_proc) == 1:
            num = input("No. of simulations to run simultaneously : ")
            run_simulation(exp_path, num, model_name)
        else:
            print("No multiprocessing")
            obj = SAMorris()
            obj.fzk_dymola(num_models, exp_path, model_name)

    def analysis(self, exp_path, num_models, samples, problem):
        # Analysis of results to get influential parameters
        sim_obj = SAMorris()
        analysis_path = exp_path + '\\case'
        result_tp = sim_obj.total_power(num_models, analysis_path)
        morris_analyze.analyze(problem, samples, result_tp, conf_level=0.95, print_to_console=True, num_levels=4)
        result_tc = sim_obj.total_cooling(num_models, analysis_path)
        morris_analyze.analyze(problem, samples, result_tc, conf_level=0.95, print_to_console=True, num_levels=4)
        result_ma = sim_obj.mean_air_temp(num_models, analysis_path)
        morris_analyze.analyze(problem, samples, result_ma, conf_level=0.95, print_to_console=True, num_levels=4)
        result_ph = sim_obj.peak_heat_load(num_models, analysis_path)
        morris_analyze.analyze(problem, samples, result_ph, conf_level=0.95, print_to_console=True, num_levels=4)
        result_pc = sim_obj.peak_cool_load(num_models, analysis_path)
        morris_analyze.analyze(problem, samples, result_pc, conf_level=0.95, print_to_console=True, num_levels=4)
        result_td = sim_obj.thermal_discomfort(num_models, analysis_path)
        morris_analyze.analyze(problem, samples, result_td, conf_level=0.95, print_to_console=True, num_levels=4)

    def fzk_dymola(self, num_models, pa, mo_name):
        # Simulations using dymola-python interface
        case_no = 1
        while case_no <= int(num_models):
            dymola = DymolaInterface()
            dymola.openModel(path=os.path.join('E:\\Aixlib', 'package.mo'))
            dir_result = pa + '\\case' + str(case_no)
            folder_name = re.findall('^[^.]+', mo_name)
            dir_cases = pa + '\\case' + str(case_no) + '\\' + folder_name[0]
            dymola.openModel(path=os.path.join(dir_cases, 'package.mo'))
            dymola.translateModel(mo_name)
            print('Simulation start...')
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
            dymola.clear()
            print('Successful')
            print(case_no)
            case_no = case_no + 1

        dymola.close()

        # dymola.experimentSetupOutput(textual=True)
        # dymola.plot(['multizone.VAir'])
        # dymola.ExportPlotAsImage("E:/yash/plot.png")
        # values = output[1]
        # print(output)

    def total_power(self, num_models, pa):
        # Calculate total heating load
        def zerolistmaker(n):
            listofzeros = [0] * n
            return listofzeros

        case_no = 1
        Y = zerolistmaker(num_models)

        while case_no <= num_models:
            mat = SimRes(pa + str(case_no) + '/demo_results.mat')
            P_Heater = mat.names('^multizone.PHeater', re=True)
            res = to_pandas(mat, P_Heater)
            m = sum(res.sum().tolist())
            Y[case_no-1] = m
            case_no = case_no + 1
        Y = np.array(Y,dtype='float64')
        return Y

    def total_cooling(self, num_models, pa):
        # Calculate total cooling load
        def zerolistmaker(n):
            listofzeros = [0] * n
            return listofzeros

        case_no = 1
        Y = zerolistmaker(num_models)

        while case_no <= num_models:
            mat = SimRes(pa + str(case_no) + '/demo_results.mat')
            P_Cooler = mat.names('^multizone.PCooler', re=True)
            res = to_pandas(mat, P_Cooler)
            m = sum(res.sum().tolist())
            Y[case_no-1] = m
            case_no = case_no + 1
        Y = np.array(Y, dtype='float64')
        return Y

    def mean_air_temp(self, num_models, pa):
        # Calculate mean air temperature
        def zerolistmaker(n):
            listofzeros = [0] * n
            return listofzeros

        case_no = 1
        Y = zerolistmaker(num_models)

        while case_no <= num_models:
            mat = SimRes(pa + str(case_no) + '/demo_results.mat')
            T_Air = mat.names('^multizone.TAir', re=True)
            res = to_pandas(mat, T_Air)
            n = sum(res.mean().tolist())/len(res.mean().tolist())
            Y[case_no - 1] = n
            case_no = case_no + 1
        Y = np.array(Y, dtype='float64')
        return Y

    def peak_heat_load(self, num_models, pa):
        # Calculate peak heating load
        def zerolistmaker(n):
            listofzeros = [0] * n
            return listofzeros

        case_no = 1
        Y = zerolistmaker(num_models)

        while case_no <= num_models:
            mat = SimRes(pa + str(case_no) + '/demo_results.mat')
            P_Heater = mat.names('^multizone.PHeater', re=True)
            res = to_pandas(mat, P_Heater)
            p = max(res.sum(axis=1).tolist())
            Y[case_no - 1] = p
            case_no = case_no + 1
        Y = np.array(Y, dtype='float64')
        return Y

    def peak_cool_load(self, num_models, pa):
        # Calculate peak cooling load
        def zerolistmaker(n):
            listofzeros = [0] * n
            return listofzeros

        case_no = 1
        Y = zerolistmaker(num_models)

        while case_no <= num_models:
            mat = SimRes(pa + str(case_no) + '/demo_results.mat')
            P_Cooler = mat.names('^multizone.PCooler', re=True)
            res = to_pandas(mat, P_Cooler)
            p = min(res.sum(axis=1).tolist())
            Y[case_no - 1] = p
            case_no = case_no + 1
        Y = np.array(Y, dtype='float64')
        return Y

    def thermal_discomfort(self, num_models, pa):
        # Calculate thermal discomfort
        def zerolistmaker(n):
            listofzeros = [0] * n
            return listofzeros

        case_no = 1
        Y = zerolistmaker(num_models)

        while case_no <= num_models:
            mat = SimRes(pa + str(case_no) + '/demo_results.mat')
            T_Air = mat.names('^multizone.TAir', re=True)
            T_Set_Heater = mat.names('^multizone.TSetHeat', re=True)
            de = to_pandas(mat, T_Air)
            fg = to_pandas(mat, T_Set_Heater)
            p = de.sort_index(axis=1)
            q = fg.sort_index(axis=1)
            x = 0
            for row in range(len(p)):
                s = 0
                for column in range(len(p.columns)):
                    r = p.iat[row, column] - q.iat[row, column]
                    s = s + abs(r)
                x = x + s / len(p.columns)
            u = x / len(p)
            Y[case_no - 1] = u
            case_no = case_no + 1
        Y = np.array(Y, dtype='float64')
        return Y

class SASobol:
    # Class defined for Sobol analysis
    def sampling(self, problem):
        # Sampling of parameters
        samples = salt.sample(problem, 500, calc_second_order=False)
        return samples

    def export(self, prj, x, ex_path):
        # Export of models for all the samples in modelica format using TEASER
        self.prj = prj
        ref_no = 1
        self.x = x

        for [infiltration, per_m2, l_p, set_temp_heat, set_temp_cool,
             owal_thick, owal_cond,
             shgc_win] in self.x:
            path = ex_path + '\\case' + str(ref_no)
            os.mkdir(path=path)
            zone = 0
            while zone < len(self.prj.buildings[0].thermal_zones):
                    self.prj.buildings[0].thermal_zones[zone].use_conditions.infiltration_rate = infiltration
                    self.prj.buildings[0].thermal_zones[zone].use_conditions.persons = per_m2
                    self.prj.buildings[0].thermal_zones[zone].use_conditions.lighting_power = l_p
                    self.prj.buildings[0].thermal_zones[zone].use_conditions.CoolerOn = True

                    hp = []
                    hp_index = 0
                    while hp_index <= 24:
                        hp.append(set_temp_heat)
                        hp_index = hp_index + 1
                    self.prj.buildings[0].thermal_zones[zone].use_conditions.heating_profile = hp

                    cp = []
                    cp_index = 0
                    while cp_index <= 24:
                        cp.append(set_temp_cool)
                        cp_index = cp_index + 1
                    self.prj.buildings[0].thermal_zones[zone].use_conditions.cooling_profile = cp


                    owal = 0
                    while owal < len(self.prj.buildings[0].thermal_zones[zone].outer_walls):
                        if len(self.prj.buildings[0].thermal_zones[zone].outer_walls[owal].layer) == 2:
                            self.prj.buildings[0].thermal_zones[zone].outer_walls[owal].layer[-1].thickness = owal_thick
                            self.prj.buildings[0].thermal_zones[zone].outer_walls[owal].layer[-1].material.thermal_conduc \
                                = owal_cond
                        else:
                            self.prj.buildings[0].thermal_zones[zone].outer_walls[owal].layer[2].thickness = owal_thick
                            self.prj.buildings[0].thermal_zones[zone].outer_walls[owal].layer[2].material.thermal_conduc \
                                = owal_cond
                        owal = owal + 1

                    win = 0
                    while win < len(self.prj.buildings[0].thermal_zones[zone].windows):
                        self.prj.buildings[0].thermal_zones[zone].windows[win].g_value = \
                            shgc_win
                        win = win + 1

                    self.prj.buildings[0].thermal_zones[zone].use_conditions.use_constant_infiltration = True

                    self.prj.buildings[0].thermal_zones[zone].calc_zone_parameters()
                    zone = zone + 1
            ref_no = ref_no + 1
            # prj.buildings[0].internal_gains_mode = 2
            self.prj.calc_all_buildings()
            zone = 0
            while zone < len(self.prj.buildings[0].thermal_zones):
                self.prj.buildings[0].thermal_zones[zone].model_attr.cool_load = -100000
                zone = zone + 1
            self.prj.weather_file_path = 'E:\\Bim2Sim\\bim2sim-coding\\teaser\\data\\input\\inputdata\\' \
                                    'weatherdata\\DEU_BW_Mannheim_107290_TRY2010_12_Jahr_BBSR.mos'

            self.prj.export_aixlib(path=path)

    def simulation(self, exp_path, num_models):
        # Function to choose one at a time or multiprocessing simulation
        multi_proc = input("If you want multiprocessing for simulations, Enter 1 else 0 : ")
        model_name = input(
            "Enter model name in dymola (which has to be simulated) : \neg. ProjektBuerogebaeude.Buerogebaeude.Buerogebaeude\n")
        if int(multi_proc) == 1:
            num = input("No. of simulations to run simultaneously : ")
            run_simulation(exp_path, num, model_name)
        else:
            print("No multiprocessing")
            obj = SAMorris()
            obj.fzk_dymola(num_models, exp_path, model_name)

    def analysis(self, exp_path, num_models, problem):
        # Analysis of results to get influential parameters
        sim_obj = SASobol()
        analysis_path = exp_path + '\\case'
        # result_tp = sim_obj.total_power(num_models, analysis_path)
        # sobol_analyze.analyze(problem, result_tp, print_to_console=True, calc_second_order=False)
        # result_tc = sim_obj.total_cooling(num_models, analysis_path)
        # sobol_analyze.analyze(problem, result_tc, print_to_console=True, calc_second_order=False)
        # result_ma = sim_obj.mean_air_temp(num_models, analysis_path)
        # sobol_analyze.analyze(problem, result_ma, print_to_console=True, calc_second_order=False)
        # result_ph = sim_obj.peak_heat_load(num_models, analysis_path)
        # sobol_analyze.analyze(problem, result_ph, print_to_console=True, calc_second_order=False)
        result_pc = sim_obj.peak_cool_load(num_models, analysis_path)
        sobol_analyze.analyze(problem, result_pc, print_to_console=True, calc_second_order=False)
        # result_td = sim_obj.thermal_discomfort(num_models, analysis_path)
        # sobol_analyze.analyze(problem, result_td, print_to_console=True, calc_second_order=False)

    def fzk_dymola(self, num_models, pa, mo_name):
        # Simulations using dymola-python interface
        case_no = 1
        while case_no <= int(num_models):
            dymola = DymolaInterface()
            dymola.openModel(path=os.path.join('E:\\Aixlib', 'package.mo'))
            dir_result = pa + '\\case' + str(case_no)
            folder_name = re.findall('^[^.]+', mo_name)
            dir_cases = pa + '\\case' + str(case_no) + '\\' + folder_name[0]
            dymola.openModel(path=os.path.join(dir_cases, 'package.mo'))
            dymola.translateModel(mo_name)
            print('Simulation start...')
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
            dymola.clear()
            print('Successful')
            print(case_no)
            case_no = case_no + 1

        dymola.close()

        # dymola.experimentSetupOutput(textual=True)
        # dymola.plot(['multizone.VAir'])
        # dymola.ExportPlotAsImage("E:/yash/plot.png")
        # values = output[1]
        # print(output)

    def total_power(self, num_models, pa):
        # Calculate total heating load
        def zerolistmaker(n):
            listofzeros = [0] * n
            return listofzeros

        case_no = 1
        Y = zerolistmaker(num_models)

        while case_no <= num_models:
            mat = SimRes(pa + str(case_no) + '/demo_results.mat')
            P_Heater = mat.names('^multizone.PHeater', re=True)
            res = to_pandas(mat, P_Heater)
            m = sum(res.sum().tolist())
            Y[case_no-1] = m
            case_no = case_no + 1
        Y = np.array(Y,dtype='float64')
        return Y

    def total_cooling(self, num_models, pa):
        # Calculate total cooling load
        def zerolistmaker(n):
            listofzeros = [0] * n
            return listofzeros

        case_no = 1
        Y = zerolistmaker(num_models)

        while case_no <= num_models:
            mat = SimRes(pa + str(case_no) + '/demo_results.mat')
            P_Cooler = mat.names('^multizone.PCooler', re=True)
            res = to_pandas(mat, P_Cooler)
            m = sum(res.sum().tolist())
            Y[case_no-1] = m
            case_no = case_no + 1
        Y = np.array(Y, dtype='float64')
        return Y

    def mean_air_temp(self, num_models, pa):
        # Calculate mean air temperature
        def zerolistmaker(n):
            listofzeros = [0] * n
            return listofzeros

        case_no = 1
        Y = zerolistmaker(num_models)

        while case_no <= num_models:
            mat = SimRes(pa + str(case_no) + '/demo_results.mat')
            T_Air = mat.names('^multizone.TAir', re=True)
            res = to_pandas(mat, T_Air)
            n = sum(res.mean().tolist())/len(res.mean().tolist())
            Y[case_no - 1] = n
            case_no = case_no + 1
        Y = np.array(Y, dtype='float64')
        return Y

    def peak_heat_load(self, num_models, pa):
        # Calculate peak heating load
        def zerolistmaker(n):
            listofzeros = [0] * n
            return listofzeros

        case_no = 1
        Y = zerolistmaker(num_models)

        while case_no <= num_models:
            mat = SimRes(pa + str(case_no) + '/demo_results.mat')
            P_Heater = mat.names('^multizone.PHeater', re=True)
            res = to_pandas(mat, P_Heater)
            p = max(res.sum(axis=1).tolist())
            Y[case_no - 1] = p
            case_no = case_no + 1
        Y = np.array(Y, dtype='float64')
        return Y

    def peak_cool_load(self, num_models, pa):
        # Calculate peak cooling load
        def zerolistmaker(n):
            listofzeros = [0] * n
            return listofzeros

        case_no = 1
        Y = zerolistmaker(num_models)

        while case_no <= num_models:
            mat = SimRes(pa + str(case_no) + '/demo_results.mat')
            P_Cooler = mat.names('^multizone.PCooler', re=True)
            res = to_pandas(mat, P_Cooler)
            p = min(res.sum(axis=1).tolist())
            Y[case_no - 1] = p
            case_no = case_no + 1
        Y = np.array(Y, dtype='float64')
        return Y

    def thermal_discomfort(self, num_models, pa):
        # Calculate thermal discomfort
        def zerolistmaker(n):
            listofzeros = [0] * n
            return listofzeros

        case_no = 1
        Y = zerolistmaker(num_models)

        while case_no <= num_models:
            mat = SimRes(pa + str(case_no) + '/demo_results.mat')
            T_Air = mat.names('^multizone.TAir', re=True)
            T_Set_Heater = mat.names('^multizone.TSetHeat', re=True)
            de = to_pandas(mat, T_Air)
            fg = to_pandas(mat, T_Set_Heater)
            p = de.sort_index(axis=1)
            q = fg.sort_index(axis=1)
            x = 0
            for row in range(len(p)):
                s = 0
                for column in range(len(p.columns)):
                    r = p.iat[row, column] - q.iat[row, column]
                    s = s + abs(r)
                x = x + s / len(p.columns)
            u = x / len(p)
            Y[case_no - 1] = u
            case_no = case_no + 1
        Y = np.array(Y, dtype='float64')
        return Y
