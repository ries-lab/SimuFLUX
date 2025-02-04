def replace_in_list(sin, *args):
    sstr = args[0::2]  # Extract search strings (odd indices)
    rstr = args[1::2]  # Extract replace strings (even indices)
    
    for j in range(len(sin)):
        for k in range(len(sstr)):
            if not isinstance(sin[j], str) or not isinstance(sstr[k], str):
                continue
            if sin[j] == sstr[k]:
                sin[j] = rstr[k]
    
    return sin
