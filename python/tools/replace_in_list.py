def replace_in_list(sin, *args):
    sstr = args[0::2]  # Extract search strings (odd indices)
    rstr = args[1::2]  # Extract replace values (even indices)

    for j in range(len(sin)):
        if not isinstance(sin[j], str):
            continue
        for k in range(len(sstr)):
            if sin[j] == sstr[k]:
                sin[j] = rstr[k]
                break
    
    return sin
