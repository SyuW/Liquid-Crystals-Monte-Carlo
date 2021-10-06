import argparse
import numpy as np
import os


class LCSystem:

    def __init__(self, lc_data_path):

        # get parameters from simulation
        self.sim_params = dict()
        with open(os.path.join(lc_data_path, "MonteCarlo_Annulus_SimNotes.txt")) as file:
            for line in file.readlines()[2:]:
                info = line.strip().split(": ")
                # special case of acceptance rate
                if info[0] == "Acceptance Rate":
                    self.sim_params[info[0]] = float(info[1].strip(" %")) / 100
                else:
                    self.sim_params[info[0]] = float(info[1])

        # get paths of all data files and retrieve system state information for each Monte Carlo step
        self.state_file_paths = [os.path.join(lc_data_path, p) for p in os.listdir(lc_data_path) if p.endswith(".csv")]
        self.state_file_paths = [os.path.normpath(p) for p in self.state_file_paths]
        self.system_state_at_step = dict()
        for p in self.state_file_paths:
            basename = os.path.basename(p)
            # final state of MC sim
            if basename == "FinalPosArray.csv":
                MC_step_no = self.sim_params["Monte Carlo steps"]
            # initial state of MC sim
            elif basename == "PosArray.csv":
                MC_step_no = 0
            # intermediate step during MC sim
            else:
                MC_step_no = int(basename.strip("PosArray").strip(".csv"))
            particle_positions = []
            with open(os.path.join(lc_data_path, p)) as data_file:
                for line in data_file:
                    # get position of particle and store in list
                    particle_pos = tuple(float(x) for x in line.split(","))
                    particle_positions.append(particle_pos)
            self.system_state_at_step[MC_step_no] = particle_positions

        thetas = [pos[-1] for pos in particle_positions]
        print(thetas)
        return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a liquid crystal system dataset")
    parser.add_argument("--data-path", help="path to dataset")

    args = parser.parse_args()

    lc_systems = []
    for _path_ in os.listdir(args.data_path):
        full_path = os.path.join(args.data_path, _path_, "instanceRun")
        lc_systems.append(LCSystem(lc_data_path=full_path))
        break