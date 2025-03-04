
def copyfields(destination, source=None, fields=None):
    missing = []
    if fields is None:
        if source is None:
            return destination, missing
        # if isinstance(source, dict) and isinstance(destination, dict):
        #     fn = source.keys() & destination.keys()
        # else:
        fn = source.keys()
        
        for k in fn:
            destination[k] = source[k]
    else:
        if source is not None:
            fnsource = source.keys()
            for k in fields:
                if k in fnsource:
                    destination[k] = source[k]
                else:
                    missing.append(k)
        else:
            missing = fields
    return destination, missing