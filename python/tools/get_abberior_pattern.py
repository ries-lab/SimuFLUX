import numpy as np

def get_abberior_pattern(itr, seq):
    if "focFldRing" in itr["Mode"]["epsf"]:
        phasemask = "tophat"
        sigma_est_ph = 130
    elif "focFldVortex" in itr["Mode"]["epsf"]:
        phasemask = "vortex"
        sigma_est_ph = 190

    pinholeorbit = ("phl" in itr["Mode"]["modulated"]) and ("hexagon" in itr["Mode"]["pattern"])
    probecenter = itr["ccrLimit"] != -1

    L = np.array(itr["patGeoFactor"]).squeeze()*360.0 # nm
    arg2D = {'makepattern': 'orbitscan'}
    dim = (0, 1)

    if itr["Mode"]["pattern"] == "hexagon":
        arg = arg2D
        patternpoints = 6
    elif itr["Mode"]["pattern"] == "square":
        arg = arg2D
        patternpoints = 4
    elif itr["Mode"]["pattern"] == "triangle":
        arg = arg2D
        patternpoints = 3
    elif itr["Mode"]["pattern"] == "zline":
        dim = 2
        lenL = len(L) if (isinstance(L, list) or isinstance(L, np.ndarray)) else 1
        patternpoints = lenL*2
        patternpos = np.zeros((patternpoints,3))
        patternpos[:,2] = np.array([-L[0], L[0], -L[1], L[1]])/2
        if probecenter:
            # probecenter, pattern points argument ignored if not makepattern
            patternpos = np.vstack([patternpos, 
                                    np.zeros(patternpos.shape[1], 
                                             dtype=patternpos.dtype)])
            
        patternpos = np.vstack([patternpos, np.zeros(patternpos.shape[1], dtype=patternpos.dtype)])
        arg = {'patternpos': patternpos}
    elif itr["Mode"]["pattern"] == "zline2":
        dim = 2
        lenL = len(L) if (isinstance(L, list) or isinstance(L, np.ndarray)) else 1
        patternpoints = lenL*2
        patternpos = np.zeros((patternpoints,3))
        if (isinstance(L, list) or isinstance(L, np.ndarray)):
            patternpos[:,2] = np.array([-L[0], L[0]])/2
        else:
            patternpos[:,2] = np.array([-L, L])/2
        if probecenter:
            # probecenter, pattern points argument ignored if not makepattern
            patternpos = np.vstack([patternpos, 
                                    np.zeros(patternpos.shape[1], 
                                             dtype=patternpos.dtype)])
        arg = {'patternpos': patternpos}
    elif itr["Mode"]["pattern"] == "octahedron":
        dim = (0,1,2)
        patternpoints = 6
        patternpos = np.zeros((6,3))
        patternpos[0,0] = L/2
        patternpos[1,1] = L/2
        patternpos[2,0] = -L/2
        patternpos[3,1] = -L/2
        patternpos[4,2] = -L/2
        patternpos[5,2] = L/2
        if probecenter:
            # probecenter, pattern points argument ignored if not makepattern
            patternpos = np.vstack([patternpos,
                                    np.zeros(patternpos.shape[1], 
                                             dtype=patternpos.dtype)])
        arg = {'patternpos': patternpos}
    else:
        raise ValueError(f"{itr['Mode']['id']} not implemented, SimSequenceFile")
    
    pointdwelltime=itr["patDwellTime"]/itr["patRepeat"]*1e3/patternpoints
    if probecenter:
        pointdwelltime = np.array([pointdwelltime, pointdwelltime*patternpoints*seq["ctrDwellFactor"]])
    laserpower = itr["pwrFactor"]

    arg2 = {"phasemask": phasemask, "orbitpoints": patternpoints, "orbitL": L,
            "probecenter": probecenter, "pointdwelltime": pointdwelltime, 
            "laserpower": laserpower, "repetitions": itr['patRepeat'], "pinholeorbit": pinholeorbit}
    parg = arg | arg2  # actually kwargs

    esth = {}

    esth["dim"] = dim

    if dim == (0,1):
        if "phl" in itr["Mode"]["modulated"]:
            esth["function"] = "est_pinholeorbit"
            esth["par"] = ["patternpos", L, sigma_est_ph, probecenter]
        else:
            esth["function"] = "est_donutLSQ1_2D"
            esth["par"] = ["patternpos", L, 310, 0]
    elif dim == 2:
        if itr["Mode"]["pattern"] == "zline":
            # 5 points: now with 3, but make a 5 point estimaotr
            esth["function"] = "est_zline"
            esth["par"] = [L]
        elif itr["Mode"]["pattern"] == "zline2":
            # 3 points
            esth["function"] = "est_qLSQiter1D"
            esth["par"] = [L]
    elif dim == (0,1,2):
        if itr["Mode"]["pattern"] == "octahedron":
            esth["function"] = "est_octahedron"
            esth["par"] = ["patternpos", L, 310, 0]

    return parg, esth