

def read_ga_input(ga_input_file):
    ga_input_dict = dict()
    try:
        with open(ga_input_file, 'r') as f:
            for line in f:
                # remove \n char and separate keys and values
                line = line[:-1].split(' = ')
                ga_input_dict[line[0]] = line[1]
                # go through each line and pull data and key
    except FileNotFoundError:
        raise FileNotFoundError("No such ga_input_file.")
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
