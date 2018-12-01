
NUM_GA_ARGS = 12

def read_ga_input(ga_input_file):
    ga_input_dict = dict()
    num_args = 0
    try:
        with open(ga_input_file, 'r') as f:
            for line in f:
                # remove \n char and separate keys and values
                line = line[:-1].lower().split(' = ')
                # ignore blank lines/long input
                if len(line) != 2:
                    print(f"Warning: line - {line} - was ignored from the ga_input_file.")
                    continue
                ga_input_dict[line[0]] = line[1]
                num_args += 1
                # go through each line and pull data and key
    except FileNotFoundError:
        raise FileNotFoundError("No such ga_input_file.")
    if num_args != NUM_GA_ARGS:
        raise ValueError(f"Incorrect number of GA arguments. Given: {num_args}. Needed: {NUM_GA_ARGS}")
    print(ga_input_dict)
    return ga_input_dict

def verify_ga_input(ga_input_dict):
    # make sure dict is non-empty
    assert len(ga_input_dict)
    # make sure it is actually a dictionary
    assert isinstance(ga_input_dict, dict)
    # key names
    key_names = {"num_mevs", "num_slots", "num_filled", "num_geoms",
                 "num_atoms", "t_size", "num_muts", "num_swaps", "pmem_dist",
                 "fit_form", "coef_energy", "coef_rmsd"}
    for key in key_names:
        try:
            assert key in ga_input_dict
        except AssertionError:
            raise ValueError(f"Misspelled/incorrectly formatted ga input parameter: {key}.")
    # make sure the inputs are of the correct format
    try:
        for key, value in ga_input_dict.items():
            if not key.startswith('coef'):
                ga_input_dict[key] = int(value)
            else:
                ga_input_dict[key] = float(value)
        # num_slots
        assert 0 < ga_input_dict["num_slots"] >= ga_input_dict["num_filled"]
        # num_filled
        assert 0 < ga_input_dict["num_filled"]
        # num_mevs
        assert ga_input_dict["num_mevs"] > 0
        # num_swaps
        assert 0 <= ga_input_dict["num_swaps"] <= ga_input_dict["num_geoms"]
        # num_muts
        assert 0 <= ga_input_dict["num_muts"] <= ga_input_dict["num_atoms"] - 3
        # num_geoms
        assert 0 < ga_input_dict["num_geoms"]
        # num_atoms
        assert 3 < ga_input_dict["num_atoms"]
        # fit_form
        assert ga_input_dict["fit_form"] == 0
        # pmem_dist
        assert 0 <= ga_input_dict["pmem_dist"] < ga_input_dict["num_slots"]/2
        # coef_energy
        assert ga_input_dict["coef_energy"] >= 0
        # coef_rmsd
        assert ga_input_dict["coef_rmsd"] >= 0
        # t_size
        assert ga_input_dict["t_size"] <= ga_input_dict["num_filled"]
#    except AssertionError as e:
#        raise ValueError("GA input values do not satisfy their constraints.")
    except ValueError as e:
        raise ValueError("GA input values should be of integer type.")
