import hydra
from bioutils.simulator import simulator_factory
from util.path_helpers import get_data_path, get_ref_path


@hydra.main(config_path="../config/simulator", config_name="pbsim2")
def main(cfg):
    simulator = simulator_factory('pbsim2', dict(cfg))
    data_path = get_data_path()
    reference_path = get_ref_path() / 'species_x1'
    simulator.run(reference_path, data_path, chr_dict={'chr1': 3})

if __name__ == "__main__":
    main()