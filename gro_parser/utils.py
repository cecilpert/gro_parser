def add_to_gro_number(number, to_add):
    """
    Function: add_to_gro_number

    Add a value to a GROMACS-style atom or residue number while handling rollover.

    Parameters:
    - number (int): The current atom or residue number.
    - to_add (int): The value to add to the current number.

    Returns:
    - int: The new atom or residue number after addition, considering rollover.

    This function is used to add a specified value to a GROMACS-style atom or residue number, while taking
    into account the rollover behavior. In GROMACS, atom and residue numbers are typically limited to five digits
    (ranging from 0 to 99999). When adding a value to the current number, if the result exceeds 99999, the
    number rolls over to 0 and continues counting. This function ensures that the rollover behavior is properly
    handled and the new number stays within the valid range.

    Example Usage:
    current_number = 99998
    added_value = 5
    new_number = add_to_gro_number(current_number, added_value)
    # In this case, new_number would be 3, as it rolls over from 99998 to 99999 and then to 1.

    Notes:
    This function assumes that the maximum atom or residue number in GROMACS is 99999.
    """
    new_number = number + to_add
    if new_number >= 99999 :
        new_number = number - 99999
    return new_number

def get_relation_params(header, funct, line, relation_type):
    is_registered = header.get(funct, False)
    if not is_registered:
        raise Exception(f"{funct} funct is not registered for {relation_type}")
    
    other_params = header[funct]["params"]
    other_params_idx = [other_params.index(p) + len(header['default']) for p in other_params]
    params = {}
    for i in range(len(other_params)):
        params[other_params[i]] = line[other_params_idx[i]]

    return params