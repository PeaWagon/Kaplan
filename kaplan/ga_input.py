
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
        raise ValueError("Incorrect number of GA arguments.")
    return ga_input_dict

def verify_ga_input(ga_input_dict):
    assert len(ga_input_dict)
    assert isinstance(ga_input_dict, dict)
    # make sure the inputs are of the correct format
    try:
        for key, value in ga_input_dict.items():
            ga_input_dict[key] = int(value)
        assert ga_input_dict["num_slots"] \
               >= ga_input_dict["num_filled"]
    except AssertionError as e:
        raise ValueError("GA input values do not satisfy their constraints.")
        print('hello', e)
    except ValueError as e:
        raise ValueError("GA input values should be of integer type.")
