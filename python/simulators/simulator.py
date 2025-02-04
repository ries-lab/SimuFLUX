from dataclasses import dataclass
from types import FunctionType, SimpleNamespace
import copy

import numpy as np

from ..tools import replace_in_list, appendstruct
from ..fluorophores import FlCollection, FlCollectionBlinking, FlBlinkBleach

@dataclass
class Component:
    functionhandle: FunctionType
    type: str
    dim: tuple
    parameters: list

@dataclass
class Pattern:
    repetitions: int = 1
    par: dict = None
    pos: list = None
    psf: list = None
    phasemask: list = None
    zeropos: list = None
    backgroundfac: list = None
    pointdwelltime: list = None
    laserpower: list = None
    L: float = 75.0
    dim: tuple = (0,1)
    type: str = "pattern"

@dataclass
class PatternScan:
    phot: np.typing.ArrayLike = None
    photbg: float = None
    intensity: np.typing.ArrayLike = None
    flpos: np.typing.ArrayLike = None
    flint = np.typing.ArrayLike = None
    averagetime: float = None
    patterntotaltime: float = None
    counter: int = None
    repetitions: int = 1
    par: dict  = None
    time: float = 0

@dataclass
class Par:
    patternpos: list = None  # list of 3D positions where the PSF is moved to. Use either zeropos or pos
    makepattern: list = None
    orbitorder: list = None  # change the order of measurement points
    phasemask: str = ""  # parameter that defines the shape of the phase mask
    zeropos: float = 0  # position of the zero when calculating PSFs (e.g., PhaseFLUX)
    orbitpoints: int = 4
    orbitL: float = 100
    probecenter: bool = True 
    pointdwelltime: float = .01  # us
    laserpower: float = 1.0  # usually we use relative, but can also be absolute.
    repetitions: int = 1
    dim: tuple = (0,1)

@dataclass
class Summary:
    photch: np.typing.ArrayLike = None
    bg_photons: float = None
    phot: float = None
    phot_signal: float = None
    pos: np.typing.ArrayLike = None
    std: np.typing.ArrayLike = None
    rmse: float = None
    bias: float = None
    locp: np.typing.ArrayLike = None
    sCRB: np.typing.ArrayLike = None
    sCRB1: np.typing.ArrayLike = None

def initlocs(maxlocalizations, fields):
    locs = {}
    for f in fields:
        locs[f] = np.zeros((maxlocalizations,1))
    return SimpleNamespace(**locs)

def removeempty(loc, loccounter):
    loco = {}
    try:
        fields = loc.keys()
    except AttributeError:
        fields = loc.__dict__.keys()
    for f in fields:
        loco[f] = getattr(loc,f)[:loccounter]
    
    return SimpleNamespace(**loco)
    

class Simulator:
    def __init__(self, fluorophores=None):
        self.patterns = {}
        self.sequences = {}
        self.fluorophores = fluorophores
        self.background = 0  # kHz from AF, does not count towards photon budget
        self.posgalvo = np.array([0, 0, 0])  # nm, position of pattern center
        self.posEOD = np.array([0, 0, 0])  # nm, not descanned, with respect to posgalvo
        self.time = 0

    def defineComponent(self, key, type_, functionhandle, parameters=[], dim=(0, 1, 2)):
        self.patterns[key] = Component(functionhandle=functionhandle,
                                       type=type_,
                                       dim=dim,
                                       parameters=parameters)

    def definePattern(self, key, psf, **kwargs):
        pattern = Pattern(repetitions = kwargs.get("repetitions", 1),
                          par = Par(**kwargs),
                          pos = kwargs.get("patternpos", [[0, 0, 0]]),
                          zeropos = kwargs.get("zeropos", [0]),
                          L = kwargs.get("orbitL", 100),
                          dim = kwargs.get("dim", (0, 1)),
                          type = "pattern")

        if "makepattern" in kwargs and kwargs["makepattern"]:
            pattern.pos = self.make_orbit_pattern(kwargs["makepattern"], kwargs.get("orbitpoints", 4), 
                                                  kwargs.get("probecenter", True), kwargs.get("orbitorder", []))
            pattern.pos *= kwargs.get("orbitL", 100) / 2
        
        zeropos = kwargs.get("zeropos", [0])
        if (pattern.pos.shape[0] == 1) and (len(zeropos) > 1):
            pattern.pos = np.tile(pattern.pos, (len(zeropos),1))
        if (pattern.pos.shape[0] > 1) and (isinstance(zeropos, int) or (len(zeropos) <= 1)):
            zeropos = np.tile(zeropos,(1,pattern.pos.shape[0]))
        
        phasemask = kwargs.get("phasemask", "")

        pattern.psf = []
        pattern.backgroundfac = []
        pattern.pointdwelltime = []
        pattern.laserpower = []
        pattern.phasemask = []

        pattern.zeropos = zeropos
        
        for k in range(pattern.pos.shape[0]):
            pattern.psf.append(psf)
            pattern.phasemask.append(phasemask)
            pattern.backgroundfac.append(1 / psf.normfactor(phasemask, pattern.zeropos[:,k]))
            pattern.pointdwelltime.append(kwargs.get("pointdwelltime", 0.01))
            pattern.laserpower.append(kwargs.get("laserpower", 1))
        
        self.patterns[key] = pattern

    def patternscan(self, key):
        pattern = self.patterns[key]
        fluorophores = self.fluorophores
        numpoints = pattern.zeropos.shape[1]
        intall = np.zeros((numpoints,1))
        repetitions = pattern.repetitions 
        timestart = self.time

        timep = 0
        posgalvo = self.posgalvo
        posEOD = self.posEOD
        flpos = np.zeros((fluorophores.numberOfFluorophores, 3))
        flintall = np.zeros((fluorophores.numberOfFluorophores,1))
        flproperties = fluorophores.getproperties()  # for performance
        bgphot = 0 
        time = self.time
        background = self.background

        for _ in range(repetitions):
            for k in range(numpoints):
                timep = timep + time   # for calculating average time point
                flposh, isactive = fluorophores.position(self.time, flproperties)
                flposrel = flposh - posgalvo
                intensityh, pinholehfac = pattern.psf[k].intensity(flposrel,
                                                                   pattern.pos[k,:] + posEOD,
                                                                   pattern.phasemask[k], 
                                                                   pattern.zeropos[:,k])
                intensityh = intensityh * pattern.laserpower[k]
                flint = fluorophores.intensity(intensityh,
                                               pattern.pointdwelltime[k],
                                               time,
                                               pinholehfac,
                                               flproperties)
                intensity = np.sum(flint, axis=0)
                flpos[isactive,:] += flposh
                flintall[isactive,:] += flint
                time = time + pattern.pointdwelltime[k]
                bgphoth = pattern.backgroundfac[k] * background
                bgphot = bgphoth + bgphot
                intall[k] += intensity + bgphoth  # sum over repetitions, fluorophores
            fluorophores.updateonoff(time)
        
        out = PatternScan()
        out.phot = np.random.poisson(intall).squeeze()  # later: fl.tophot(intenall): adds bg, multiplies with brightness, does 
        out.photbg = bgphot
        out.intensity = intall 
        out.flpos = flpos/numpoints/repetitions
        out.flint = flintall
        out.averagetime = timep/numpoints/repetitions
        out.patterntotaltime = time-timestart
        out.counter = 1
        out.repetitions = repetitions
        out.par = copy.deepcopy(pattern.par)
        out.par.L = pattern.L
        out.par.patternpos = pattern.pos
        out.par.zeropos = pattern.zeropos
        out.par.dim = pattern.dim
        self.time = time

        return out
    
    def runSequenceintern(self, seq, maxlocalizations):
        numseq = len(seq)
        numpat = 0 
        ls = np.zeros(numseq, dtype=np.int32)
        for s in reversed(range(numseq)):
            pat = self.patterns[seq[s]]
            if pat.type == "pattern":
                ls[s] = pat.zeropos.shape[1]
                numpat += 1

        photch = np.zeros((maxlocalizations*numpat,max(ls)))-1
        photbg = np.zeros((maxlocalizations*numpat,1))
        loc = initlocs(maxlocalizations*numpat,['xnm','ynm','znm','xfl1',
            'yfl1','zfl1','phot','time','abortcondition', 'xgalvo', 'ygalvo', 'zgalvo',
            'xeod','yeod','zeod'])
        
        bleached = False
        loccounter = -1

        pos = [None]*(maxlocalizations*numpat)
        int = [None]*(maxlocalizations*numpat)
        fl = SimpleNamespace(pos=pos,int=int)

        for _ in range(maxlocalizations):
            if self.fluorophores.remainingphotons < 1:
                bleached = True
                break
            posgalvo_beforecenter = self.posgalvo
            xest = np.array([0,0,0])
            for s in range(numseq):
                curr_seq = seq[s]
                component = self.patterns[curr_seq]
                type_ = component.type
                if type_ == "pattern":

                    scanout = self.patternscan(curr_seq)
                    loccounter += 1  # every localization gets a new entry, as in abberior
                    loc.phot[loccounter] += np.sum(scanout.phot,axis=0)
                    photch[loccounter,:len(scanout.phot)] = scanout.phot
                    photbg[loccounter] = scanout.photbg
                    loc.time[loccounter] += scanout.averagetime
                    flpos = scanout.flpos[0,:]
                    loc.xfl1[loccounter] = flpos[0]
                    loc.yfl1[loccounter] = flpos[1]
                    loc.zfl1[loccounter] = flpos[2]
                    fl.pos[loccounter] = np.zeros((scanout.flpos.shape[0],3))
                    fl.pos[loccounter][:scanout.flpos.shape[0],:] = scanout.flpos
                    fl.int[loccounter] = np.zeros((scanout.flpos.shape[0],3))
                    fl.int[loccounter][:scanout.flpos.shape[0],:] = scanout.flint
                elif type_ == "estimator":
                    # replace placeholder names by values
                    component_par=replace_in_list(component.parameters,
                                                  'patternpos',
                                                  scanout.par.patternpos,
                                                  'L',
                                                  scanout.par.L,
                                                  'probecenter',
                                                  scanout.par.probecenter,
                                                  'bg_phot',
                                                  scanout.photbg)
                    xesth=component.functionhandle(scanout.phot,*component_par)
                    if len(xesth)==3:
                        xest[np.array(component.dim)] = xesth[np.array(component.dim)]
                    else:
                        xest[np.array(component.dim)] = xesth
                    xesttot = xest + posgalvo_beforecenter
                    loc.xnm[loccounter] = xesttot[0]
                    loc.ynm[loccounter] = xesttot[1]
                    loc.znm[loccounter] = xesttot[2]
                    loc.xgalvo[loccounter] = self.posgalvo[0]
                    loc.ygalvo[loccounter] = self.posgalvo[1]
                    loc.zgalvo[loccounter] = self.posgalvo[2]
                    loc.xeod[loccounter] = self.posEOD[0]
                    loc.yeod[loccounter] = self.posEOD[1]
                    loc.zeod[loccounter] = self.posEOD[2]
                elif type_ == "positionupdater":
                    posgalvo, posEOD = component.functionhandle(xest, 
                                                                self.posgalvo, 
                                                                self.posEOD, 
                                                                *component.parameters)
                    self.posgalvo[np.array(component.dim)] = posgalvo[np.array(component.dim)]
                    self.posEOD[np.array(component.dim)] = posEOD[np.array(component.dim)]
                else:
                    raise ValueError(f"Unknown component type {type_}.")
        loc = removeempty(loc, loccounter)
        loc.abortcondition = np.zeros(loc.phot.shape)
        loc.abortcondition[-1] = 1 + 2*bleached
        out = SimpleNamespace(loc=loc,
                              raw=photch[:loccounter,:],
                              fluorophores = fl,
                              bg_photons = photbg[:loccounter])
        return out
    
    def runSequence(self, seq, **kwargs):
        out = None
        for _ in range(kwargs.get('repetitions', 1)):
            self.fluorophores.reset(self.time)
            out2 = self.runSequenceintern(seq,kwargs.get('maxlocs',1000))
            out = self.appendout(out,out2)
        out.sequence = seq

        return out

    def appendout(self, out1=None, out2=None):
        if out2 is None:
            return out1
        if out1 is None:
            out1 = out2
            out1.raw=out2.raw
            out1.loc.rep=1+0*out1.loc.xnm
            out1.fluorophores.pos=out2.fluorophores.pos
            out1.fluorophores.int=out2.fluorophores.int
        else:
            sr = out2.raw.shape
            if len(sr)==2:  # Abberior
                out1.raw = np.concatenate(out1.raw, out2.raw, axis=0)
            else:
                out1.raw = np.concatenate(out1.raw, np.expand_dims(out2.raw, axis=0), axis=0)
            out1.fluorophores.pos = np.concatenate((out1.fluorophores.pos, out2.fluorophores.pos), axis=0)
            out1.fluorophores.int = np.concatenate((out1.fluorophores.int, out2.fluorophores.int), axis=0)
            out2['loc']['rep'] = out1['raw'].shape[0] + 0 * out2['loc']['xnm']
            out1['loc'] = appendstruct(out1.loc, out2.loc)
        
        return out1

    def make_orbit_pattern(self, pattern_type, orbit_points, use_center, orbit_order=None):
        if pattern_type == "orbitscan":
            phi = np.linspace(0, 2 * np.pi, orbit_points, endpoint=False)
            pos_pattern = np.column_stack((np.cos(phi), np.sin(phi), np.zeros_like(phi)))
            if use_center:
                pos_pattern = np.vstack([pos_pattern, [0, 0, 0]])
            if orbit_order:
                pos_pattern = pos_pattern[orbit_order]
            return pos_pattern
        elif pattern_type == "zscan":
            pos_pattern = np.array([[-1, 0, 0], [1, 0, 0]])
            if use_center:
                pos_pattern = np.vstack([pos_pattern, [0, 0, 0]])
            return pos_pattern
        return np.array([])
    
    def locprec(self, photons, L):
        return L/np.sqrt(8*photons)
    
    def calculateCRBdirect(self, patternnames, dim=(0, 1, 2), position=None):
        if position is None:
            position, _ = self.fluorophores.position(0)
        
        flpos = position
        brightness = self.fluorophores.brightness
        bg = self.background / brightness[0]  # Hack, to not have to calculate with all fluorophores

        eps = 1
        IFisher = np.zeros((len(dim), len(dim)))
        
        pattern = self.patterns[patternnames]

        def pi(dpos):
            flposh = flpos.copy()
            flposh[0, :] -= dpos
            bgh = bg * pattern.backgroundfac[k]
            ih = np.sum(pattern.psf[k].intensity(flposh - self.posgalvo, 
                                                pattern.pos[k,:], 
                                                pattern.phasemask[k], 
                                                pattern.zeropos[:,k]) + bgh)

            ihm = 0
            for m in reversed(range(len(pattern.zeropos))):
                bgh = bg * pattern.backgroundfac[m]
                ihm += np.sum(pattern.psf[m].intensity(flposh - self.posgalvo, 
                                                    pattern.pos[m,:], 
                                                    pattern.phasemask[m], 
                                                    pattern.zeropos[:,m]) + bgh)

            return ih / ihm

        for k in reversed(range(len(pattern.zeropos))):
            dpdc = np.zeros(len(dim))
            for coord in range(len(dim)):
                dposa = np.zeros(3)
                dposa[dim[coord] - 1] = eps / 2
                dposa2 = dposa.copy()
                dposa2[dim[coord] - 1] = -eps / 2
                dpdc[dim[coord] - 1] = (pi(dposa) - pi(dposa2)) / eps

            for coord in range(len(dim)):
                for coord2 in range(len(dim)):
                    IFisher[coord, coord2] += dpdc[dim[coord] - 1] * dpdc[dim[coord2] - 1] / (pi(np.zeros(3)) + 1e-5)

        crlb = np.linalg.inv(IFisher)
        locprech = np.sqrt(np.diag(crlb))
        locprec = np.zeros(3)
        locprec[np.array(dim)] = locprech

        return locprec

    def calculateCRB(self, patternnames, dim=(0, 1, 2), position=None):
        if position is None:
            position, _ = self.fluorophores.position(0)

        if isinstance(self.fluorophores, (FlCollectionBlinking, FlBlinkBleach)):
            return self.calculateCRBdirect(patternnames, dim=dim, position=position)

        if isinstance(self.fluorophores, FlCollection):
            fl1 = self.fluorophores.flall[0]
        else:
            fl1 = self.fluorophores

        oldpar = copy.deepcopy(fl1)
        fl1.remainingphotons = np.inf
        eps = 1
        IFisher = np.zeros((len(dim), len(dim)))

        out0 = self.patternscan(patternnames)
        pi0 = (out0.intensity / np.sum(out0.intensity)).squeeze()

        dpdc = np.zeros((len(pi0), len(dim)))

        for coord in range(len(dim)):
            dposa = np.zeros(3)
            dposa[dim[coord]] = eps / 2
            dposa2 = dposa.copy()
            dposa2[dim[coord]] = -eps / 2

            fl1.pos = oldpar.pos + dposa
            out1 = self.patternscan(patternnames)
            fl1.pos = oldpar.pos + dposa2
            out2 = self.patternscan(patternnames)
            dpdc[:, dim[coord]] = ((out1.intensity / np.sum(out1.intensity, axis=0) - out2.intensity / np.sum(out2.intensity, axis=0)) / eps).squeeze()

        for coord in range(len(dim)):
            for coord2 in range(len(dim)):
                fh = dpdc[:, dim[coord]] * dpdc[:, dim[coord2]] / (pi0 + 1e-12)
                IFisher[coord, coord2] = np.sum(fh)

        fl1.pos = oldpar.pos
        fl1.remainingphotons = oldpar.remainingphotons

        try:
            crlb = np.linalg.inv(IFisher)
        except np.linalg.LinAlgError:
            crlb = np.ones_like(IFisher)*np.nan
        locprech = np.sqrt(np.diag(crlb))
        locprec = np.zeros(3)
        locprec[np.array(dim)] = locprech

        return locprec

    def summarize_results(self, out, display=True, filter=None):
        if filter is None:
            filter = np.ones(out.loc.xnm.shape, dtype=bool)
            
        ind = filter.squeeze()
        photraw = out.raw.copy()
        photraw[photraw == -1] = np.nan
        
        if len(photraw.shape) > 3:
            photraw[photraw < 0] = np.nan
            photch = np.nanmean(np.nanmean(photraw[ind, :], axis=0), axis=0).squeeze()
        else:
            photch = np.nanmean(photraw[ind], axis=0).squeeze()
        
        xest = np.column_stack((out.loc.xnm[ind], out.loc.ynm[ind], out.loc.znm[ind]))
        flpos = np.column_stack((out.loc.xfl1[ind], out.loc.yfl1[ind], out.loc.zfl1[ind]))
        phot = out.loc.phot[ind]
        
        seq = out.sequence
        lp = np.zeros(3)
        sigmaCRB = np.zeros(3)
        
        for k in range(len(seq)):
            pattern = self.patterns[seq[k]]
            if pattern.type == "pattern":
                sh = self.calculateCRB(seq[k], dim=pattern.dim)
                sigmaCRB[np.array(pattern.dim)] = sh[np.array(pattern.dim)]
                Lh = self.locprec(np.nanmean(phot), pattern.L)
                lp[np.array(pattern.dim)] = Lh
        
        st = Summary()
        st.photch = photch
        st.bg_photons = np.nanmean(out.bg_photons[ind])
        st.phot = np.nanmean(phot)
        st.phot_signal = st.phot - st.bg_photons
        st.pos = np.nanmean(xest, axis=0)
        st.std = np.nanstd(xest, axis=0)
        st.rmse = np.sqrt(np.nanmean((xest - flpos) ** 2, axis=0))
        st.bias = np.nanmean(xest - flpos, axis=0)
        st.locp = lp
        st.sCRB = sigmaCRB / np.sqrt(st.phot)
        st.sCRB1 = sigmaCRB

        
        if display:
            ff1 = '.1f'
            print(f'photch: {", ".join(list(map("{:.1f}".format, st.photch)))}, mean(phot): {st.phot:{ff1}}, signal phot: {st.phot_signal:{ff1}}')
            ff = '{:.2f}'
            print(f'std:  {", ".join(list(map(ff.format, st.std)))}, rmse: {", ".join(list(map(ff.format, st.rmse)))}, pos: {", ".join(list(map(ff.format, st.pos)))}, bias: {", ".join(list(map(ff.format, st.bias)))}')
            print(f'locp: {", ".join(list(map(ff.format, lp)))}, sCRB: {", ".join(list(map(ff.format, st.sCRB)))}, sCRB*sqrt(phot): {", ".join(list(map("{:.1f}".format, st.sCRB1)))}')
        
        return st
