from dataclasses import dataclass, field
from types import FunctionType, SimpleNamespace
import copy

import numpy as np

from ..tools import replace_in_list, appendstruct
from ..fluorophores import FlCollection, FlCollectionBlinking, FlBlinkBleach, FlMoving

@dataclass
class Component:
    functionhandle: FunctionType  # handle to the function of the component (e.g. an estimator)
    type: str                     # estimator, position updater, background
    dim: tuple                    # for estimators: which dimension is estimated, updated
    parameters: list              # parameters are passed on to the function

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
class Time:
    averagetime: float = 0
    patterntotaltime: float = 0
    patterntime: float = 0

@dataclass
class PatternScan:
    phot: np.typing.ArrayLike = None
    photrate: np.typing.ArrayLike = None
    pointdwelltime: np.typing.ArrayLike = None
    bg_photons_gt: np.typing.ArrayLike = None
    bgphot_est: np.typing.ArrayLike = None
    intensity: np.typing.ArrayLike = None
    flpos: np.typing.ArrayLike = None
    flint = np.typing.ArrayLike = None
    time: Time = field(default_factory=lambda: Time())
    counter: int = None
    repetitions: int = 1
    par: dict  = None
    time: float = 0
    pattern: Pattern = field(default_factory=lambda: Pattern())

@dataclass
class Par:
    patternpos: list = None  # list of 3D positions where the PSF is moved to. Use either zeropos or pos
    makepattern: list = None  # orbitscan, zscan, None (default) no pattern is made
    orbitorder: list = None  # change the order of measurement points. can create problems with some estimators
    phasemask: str = ""  # parameter that defines the shape of the phase mask
    zeropos: float = 0  # position of the zero when calculating PSFs (e.g., PhaseFLUX)
    orbitpoints: int = 4  # orbitpoints
    orbitL: float = 100  # diamter of the scan pattern
    probecenter: bool = True  # If to perform a central measurement (cfr check)
    pointdwelltime: float = .01  # (us) length of a point measurement. Single value, vector with length of the pattern or vector with length 2, then the second value is used for central measurement
    laserpower: float = 1.0  # usually we use relative, but can also be absolute.
    repetitions: int = 1  # repetitions of the patern scan before position estimation
    dim: tuple = (0,1)  # dimensions in which the scan is performed

@dataclass
class Summary:
    photch: np.typing.ArrayLike = None
    bg_photons_gt: np.typing.ArrayLike = None
    bg_photons_est: np.typing.ArrayLike = None
    phot: float = None
    phot_signal: float = None
    pos: np.typing.ArrayLike = None
    std: np.typing.ArrayLike = None
    rmse: float = None
    bias: float = None
    locp: np.typing.ArrayLike = None
    sCRB: np.typing.ArrayLike = None
    sCRB1: np.typing.ArrayLike = None
    duration : float = 0
    stdraw : np.typing.ArrayLike = None
    efo : np.typing.ArrayLike = None
    photrate: np.typing.ArrayLike = None

@dataclass
class Deadtimes:
    point: float = 0
    pattern: float = 0
    estimator: float = 0
    positionupdate: float = 0
    localization: float = 0

@dataclass
class FOVScan:
    std: np.typing.ArrayLike = None
    rmse: np.typing.ArrayLike = None
    bias: np.typing.ArrayLike = None
    sCRB: np.typing.ArrayLike = None
    pos: np.typing.ArrayLike = None
    stdrel: np.typing.ArrayLike = None
    biasrel: np.typing.ArrayLike = None

def initlocs(maxlocalizations, fields):
    locs = {f: np.zeros((maxlocalizations,1)) for f in fields}
    return SimpleNamespace(**locs)

def removeempty(loc, loccounter):
    try:
        fields = loc.keys()
    except AttributeError:
        fields = loc.__dict__.keys()
    
    loco = {f: getattr(loc,f)[:loccounter] for f in fields}
    
    return SimpleNamespace(**loco)
    

class Simulator:
    def __init__(self, 
                 fluorophores: np.typing.ArrayLike = None, 
                 background: float = 0, 
                 background_estimated: float = 0, 
                 loadfile: list = None):
        self.patterns = {}  # patterns that are defined
        # self.sequences = {}
        self.fluorophores = fluorophores  # Fluorophore or FlCollection
        self.background = 0  # kHz from AF, does not count towards photon budget
        self._posgalvo = np.array([0, 0, 0], dtype=float)  # nm, position of pattern center
        self._posEOD = np.array([0, 0, 0], dtype=float)  # nm, not descanned, with respect to 
                                           # posgalvo
        self.time = 0  # That is the master time
        self.background = background  # Constant autofluorescence background
        self.background_estimated = background_estimated  # Estimated background, 
                                                          # used for estimators
        self.loadfile = loadfile  # sequence, only for SimSequencefile
        self.deadtimes = Deadtimes()  # Dead times added to time after each point, 
                                      # pattern, estimation, position update or 
                                      # localization (sequence)

        if loadfile is not None:
            self.loadsequence(*loadfile)
    
    @property
    def posgalvo(self):
        return self._posgalvo
    
    @posgalvo.setter
    def posgalvo(self, val):
        if not isinstance(val, np.ndarray):
            self._posgalvo = np.array(val, dtype=float)
        else:
            self._posgalvo = val

    @property
    def posEOD(self):
        return self._posEOD
    
    @posEOD.setter
    def posEOD(self, val):
        if not isinstance(val, np.ndarray):
            self._posEOD = np.array(val, dtype=float)
        else:
            self._posEOD = val

    def defineComponent(self, key, type_, functionhandle, parameters=None, dim=(0, 1, 2)):
        self.patterns[key] = Component(functionhandle=functionhandle,  
                                       type=type_,  
                                       dim=dim,  
                                       parameters=parameters)

    def definePattern(self, key, psf, **kwargs):
        orbitL = kwargs.get("orbitL", 100)
        pattern = Pattern(repetitions = kwargs.get("repetitions", 1),  
                          par = Par(**kwargs),
                          pos = kwargs.get("patternpos", np.atleast_2d(np.array([0, 0, 0]))),
                          zeropos = kwargs.get("zeropos", np.array([[0]])),
                          L = orbitL,
                          dim = kwargs.get("dim", (0, 1)),
                          type = "pattern")

        if "makepattern" in kwargs and kwargs["makepattern"]:
            pattern.pos = self.make_orbit_pattern(kwargs["makepattern"], kwargs.get("orbitpoints", 4), 
                                                  kwargs.get("probecenter", True), kwargs.get("orbitorder", None))
            pattern.pos *= orbitL / 2
        
        zeropos = kwargs.get("zeropos", np.array([[0]]))

        if (pattern.pos.shape[0] == 1) and (zeropos.shape[1] > 1):
            pattern.pos = np.tile(pattern.pos, (zeropos.shape[1], 1))
        if (pattern.pos.shape[0] > 1) and (isinstance(zeropos, int) or (zeropos.shape[1] <= 1)):
            zeropos = np.tile(zeropos, (1, pattern.pos.shape[0]))
        
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
            pattern.psf[k].calculatePSFs(phasemask, pattern.zeropos[:,k])
            pattern.laserpower.append(kwargs.get("laserpower", 1))
        
        pdt = kwargs.get("pointdwelltime", [0.01])
        if not isinstance(pdt, (int, float)) and len(pdt) == pattern.pos.shape[0]:
            pattern.pointdwelltime = pdt
        else:
            try:
                pattern.pointdwelltime = np.zeros((1,pattern.pos.shape[0]))+pdt[0]
            except (TypeError, IndexError):
                pattern.pointdwelltime = np.zeros((1,pattern.pos.shape[0]))+pdt
        if not isinstance(pdt, (int, float)) and len(pdt) == 2:
            pattern.pointdwelltime[-1] = pdt[1]

        self.patterns[key] = pattern

    def patternscan(self, key):
        pattern = self.patterns[key]
        fluorophores = self.fluorophores
        numpoints = pattern.zeropos.shape[1]
        intall = np.zeros((numpoints,1))
        repetitions = pattern.repetitions 
        # timestart = self.time
        deadtimes = self.deadtimes

        timep = 0
        posgalvo = self.posgalvo
        posEOD = self.posEOD
        flpos = np.zeros((fluorophores.numberOfFluorophores, 3))
        flintall = np.zeros((fluorophores.numberOfFluorophores, 1))
        flproperties = fluorophores.getproperties()  # for performance
        bgphot = 0 
        time = self.time
        background = self.background

        for _ in range(repetitions):
            for k in range(numpoints):
                timep += time   # for calculating average time point
                flposh, isactive = fluorophores.position(self.time, flproperties)
                try:
                    flposrel = flposh - posgalvo
                except TypeError:
                    # print(f"Wrong type flposh: {flposh} and posgalvo: {posgalvo}")
                    flposrel = flposh - posgalvo
                intensityh, pinholehfac = pattern.psf[k].intensity(flposrel[isactive,:],
                                                                   pattern.pos[k,:] + posEOD,
                                                                   pattern.phasemask[k], 
                                                                   pattern.zeropos[:,k])
                intensityh *= pattern.laserpower[k]
                flint = fluorophores.intensity(intensityh,
                                               pattern.pointdwelltime[:,k],
                                               time,
                                               pinholehfac,
                                               flproperties)
                intensity = np.sum(flint, axis=0)
                flpos += flposh
                flintall[isactive,:] += flint
                time = time + pattern.pointdwelltime[:,k] + deadtimes.point
                # bgphoth = pattern.backgroundfac[k] * background * pattern.pointdwelltime[:,k]
                bgphoth = background * pattern.pointdwelltime[:,k] * pattern.laserpower[k]
                bgphot += bgphoth
                intall[k] += intensity + bgphoth  # sum over repetitions, fluorophores
            time = time + deadtimes.pattern
            fluorophores.updateonoff(time)
        
        out = PatternScan()
        out.phot = np.random.poisson(intall)  # later: fl.tophot(intenall): adds bg, multiplies with brightness, does 
        out.photrate = out.phot / pattern.pointdwelltime.T
        out.pointdwelltime = pattern.pointdwelltime.T
        out.bg_photons_gt = np.array([bgphot], dtype=float)
        out.bgphot_est = np.array([0], dtype=float)
        out.intensity = intall 
        out.flpos = flpos/numpoints/repetitions
        out.flint = flintall
        out.time = Time()
        out.time.averagetime = timep/numpoints/repetitions
        out.time.patterntotaltime = pattern.pointdwelltime.sum()*repetitions
        out.time.patterntime = pattern.pointdwelltime.sum()
        out.counter = 1
        out.repetitions = repetitions
        out.par = pattern.par  # copy.deepcopy(pattern.par)
        out.par.L = pattern.L
        out.par.patternpos = pattern.pos
        out.par.zeropos = pattern.zeropos
        out.par.dim = pattern.dim
        out.par.pattern = pattern  # copy.deepcopy(pattern)
        
        self.time = time

        return out
    
    def runSequenceintern(self, seq, maxlocalizations):
        timestart = self.time
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
            'xeod','yeod','zeod','photrate','seq'])
        
        bleached = False
        loccounter = -1

        pos = np.zeros(maxlocalizations*numpat)  # [None]*(maxlocalizations*numpat)
        int = np.zeros(maxlocalizations*numpat)  # [None]*(maxlocalizations*numpat)
        fl = SimpleNamespace(pos=pos,int=int)
        par = []
        deadtimes = self.deadtimes

        for k in range(maxlocalizations):
            if self.fluorophores.remainingphotons < 1:
                bleached = True
                break
            posgalvo_beforecenter = self.posgalvo
            xest = np.array([0,0,0], dtype=float)
            loccounter_seq = loccounter + 1  # XXXX 9.2.25
            for s in range(numseq):
                curr_seq = seq[s]
                component = self.patterns[curr_seq]
                type_ = component.type
                if type_ == "pattern":
                    scanout = self.patternscan(curr_seq)
                    loccounter += 1  # every localization gets a new entry, as in abberior
                    loc.phot[loccounter] += np.sum(scanout.phot)
                    # print(loc.photrate[loccounter].shape, scanout.photrate.shape)
                    loc.photrate[loccounter] += np.mean(scanout.photrate)
                    photch[loccounter,:len(scanout.phot)] = scanout.phot.squeeze()
                    photbg[loccounter] = scanout.bg_photons_gt
                    loc.time[loccounter] += scanout.time.averagetime
                    flpos = scanout.flpos[0,:]
                    loc.xfl1[loccounter] = flpos[0]
                    loc.yfl1[loccounter] = flpos[1]
                    loc.zfl1[loccounter] = flpos[2]
                    fl.pos[loccounter] = np.zeros((scanout.flpos.shape[0],3))
                    fl.pos[loccounter][:scanout.flpos.shape[0],:] = scanout.flpos
                    fl.int[loccounter] = np.zeros((scanout.flpos.shape[0],3))
                    fl.int[loccounter][:scanout.flpos.shape[0],:] = scanout.flint
                    if k == 1:
                        while len(par) < (s+1):
                            par.append([])
                        par[s] = scanout.par
                    loc.seq[loccounter,0] = s
                elif type_ == "estimator":
                    # replace placeholder names by values
                    component_par=replace_in_list(component.parameters, 
                                                  'patternpos', scanout.par.patternpos, 
                                                  'L', scanout.par.L, 
                                                  'probecenter', scanout.par.probecenter, 
                                                  'bg_photons_gt', scanout.bg_photons_gt,
                                                  'background_estimated', self.background_estimated,
                                                  'iteration', s)
                    xesth = component.functionhandle(scanout.photrate,*component_par)
                    if len(xesth)==3:
                        xest[np.array(component.dim)] = xesth[np.array(component.dim)]
                    else:
                        xest[np.array(component.dim)] = xesth
                    xesttot = xest + posgalvo_beforecenter
                    loc.xgalvo[loccounter] = self.posgalvo[0]
                    loc.ygalvo[loccounter] = self.posgalvo[1]
                    loc.zgalvo[loccounter] = self.posgalvo[2]
                    loc.xeod[loccounter] = self.posEOD[0]
                    loc.yeod[loccounter] = self.posEOD[1]
                    loc.zeod[loccounter] = self.posEOD[2]
                    self.time += deadtimes.estimator
                elif type_ == "positionupdater":
                    posgalvo, posEOD = component.functionhandle(xest, 
                                                                self.posgalvo, 
                                                                self.posEOD, 
                                                                *component.parameters)
                    self.posgalvo[np.array(component.dim)] = posgalvo[np.array(component.dim)]
                    self.posEOD[np.array(component.dim)] = posEOD[np.array(component.dim)]
                    self.time += deadtimes.positionupdate
                elif type_ == "background":
                    component_par = replace_in_list(component.parameters,
                                                    'patternpos', scanout.par.patternpos,
                                                    'L', scanout.par.L,
                                                    'probecenter', scanout.par.probecenter,
                                                    'bg_photons_gt', scanout.bg_photons_gt,
                                                    'background_estimated', self.background_estimated,
                                                    'iteration', s)
                else:
                    raise ValueError(f"Unknown component type {type_}.")
            loc.xnm[loccounter_seq:(loccounter+1)] = xesttot[0]
            loc.ynm[loccounter_seq:(loccounter+1)] = xesttot[1]
            loc.znm[loccounter_seq:(loccounter+1)] = xesttot[2]
            self.time += deadtimes.localization
        loc = removeempty(loc, loccounter)
        loc.abortcondition = np.zeros(loc.phot.shape)
        loc.abortcondition[-1] = 1 + 2*bleached
        out = SimpleNamespace(loc=loc,
                              raw=photch[:loccounter,:],
                              fluorophores = fl,
                              bg_photons_gt = photbg[:loccounter].T,
                              par = par,
                              duration = self.time-timestart)
        if hasattr(scanout, 'bgphot_est'):
            out.bg_photons_est = np.sum(scanout.bgphot_est, axis=0)
        return out
    
    def runSequence(self, seq, **kwargs):
        out = None
        timestart = self.time
        for _ in range(kwargs.get('repetitions', 1)):
            self.fluorophores.reset(self.time)
            out2 = self.runSequenceintern(seq, kwargs.get('maxlocs',1000))
            out = self.appendout(out, out2)
        out.sequence = seq
        out.duration = self.time-timestart

        return out

    def appendout(self, out1=None, out2=None):
        if out2 is None:
            return out1
        if out1 is None:
            out1 = out2
            out1.raw = out2.raw
            out1.loc.rep = 1 + 0*out1.loc.xnm
            out1.fluorophores.pos = out2.fluorophores.pos
            out1.fluorophores.int = out2.fluorophores.int
            out1.bg_photons_gt = out2.bg_photons_gt
            out1.bg_photons_est = out2.bg_photons_est
        else:
            sr = out2.raw.shape
            if len(sr)==2:  # Abberior
                out1.raw = np.concatenate(out1.raw, out2.raw, axis=0)
            else:
                out1.raw = np.concatenate(out1.raw, np.expand_dims(out2.raw, axis=0), axis=0)
            out1.fluorophores.pos = np.concatenate((out1.fluorophores.pos, out2.fluorophores.pos), axis=0)
            out1.fluorophores.int = np.concatenate((out1.fluorophores.int, out2.fluorophores.int), axis=0)
            out2.loc.rep = out1.raw.shape[0] + 0*out2.loc.xnm
            out1.loc = appendstruct(out1.loc, out2.loc)
            out1.bg_photons_gt = np.concatenate(out1.bg_photons_gt, out2.bg_photons_gt)
            out1.bg_photons_est = np.concatenate(out1.bg_photons_est, out2.bg_photons_est)

        return out1

    def make_orbit_pattern(self, pattern_type, orbit_points, use_center, orbit_order=None):
        if pattern_type == "orbitscan":
            phi = np.linspace(0, 2 * np.pi, orbit_points, endpoint=False)
            pos_pattern = np.column_stack((np.cos(phi), np.sin(phi), np.zeros_like(phi)))
            if use_center:
                pos_pattern = np.vstack([pos_pattern, [0, 0, 0]], dtype=float)
            if orbit_order:
                pos_pattern = pos_pattern[orbit_order]
            return pos_pattern
        elif pattern_type == "zscan":
            pos_pattern = np.array([[-1, 0, 0], [1, 0, 0]], dtype=float)
            if use_center:
                pos_pattern = np.vstack([pos_pattern, [0, 0, 0]], dtype=float)
            return pos_pattern
        return np.array([])
    
    def locprec(self, photons, L):
        return L/np.sqrt(8*photons)
    
    def calculateCRBdirect(self, patternnames, dim=(0, 1, 2), position=None):
        if position is None:
            position, _ = self.fluorophores.position(0)
        
        flpos = position
        brightness = self.fluorophores.brightness
        try:
            bg = self.background / brightness[0]  # Hack, to not have to calculate with all fluorophores
        except TypeError:
            bg = self.background / brightness
        pattern = self.patterns[patternnames]

        def pi(dpos):
            flposh = flpos.copy()
            try:
                flposh[0, :] -= dpos
            except IndexError:
                flposh -= dpos
            # bgh = bg * pattern.backgroundfac[k]
            ih, _ = pattern.psf[k].intensity(flposh - self.posgalvo, 
                                             pattern.pos[k,:], 
                                             pattern.phasemask[k], 
                                             pattern.zeropos[:,k])
            ih = np.sum(ih + bg)

            ihm = 0
            for m in reversed(range(pattern.zeropos.shape[1])):
                # bgh = bg * pattern.backgroundfac[m]
                ihmp, _ = pattern.psf[m].intensity(flposh - self.posgalvo, 
                                                   pattern.pos[m,:], 
                                                   pattern.phasemask[m], 
                                                   pattern.zeropos[:,m])
                ihm += np.sum(ihmp + bg)

            return ih / ihm
        
        eps = 1
        IFisher = np.zeros((len(dim), len(dim)))

        for k in reversed(range(pattern.zeropos.shape[1])):
            dpdc = np.zeros(len(dim))
            for coord in range(len(dim)):
                dposa = np.zeros(3)
                dposa[dim[coord]] = eps / 2
                dposa2 = dposa.copy()
                dposa2[dim[coord]] = -eps / 2
                dpdc[dim[coord]] = (pi(dposa) - pi(dposa2)) / eps

            for coord in range(len(dim)):
                for coord2 in range(len(dim)):
                    # IFisher[coord, coord2] += dpdc[dim[coord] - 1] * dpdc[dim[coord2] - 1] / (pi(np.zeros(3)) + 1e-5)
                    IFisher[coord, coord2] += dpdc[dim[coord]] * dpdc[dim[coord2]] / (pi(np.zeros(3)) + 1e-4)

        crlb = np.linalg.inv(IFisher)
        locprech = np.sqrt(np.diag(crlb)).T
        locprec = np.zeros(3)
        locprec[np.array(dim)] = locprech

        return locprec

    def calculateCRBpattern(self, patternnames, dim=(0, 1, 2), position=None):
        if position is None:
            position, _ = self.fluorophores.position(0)

        if isinstance(self.fluorophores, (FlCollectionBlinking, FlBlinkBleach, FlMoving)):
            return self.calculateCRBdirect(patternnames, dim=dim, position=position)

        if isinstance(self.fluorophores, FlCollection):
            fl1 = self.fluorophores.flall[0]
        else:
            fl1 = self.fluorophores

        oldpar = copy.deepcopy(fl1)
        fl1.remainingphotons = np.inf
        posold, _ = fl1.position(self.time)
        eps = 1
        IFisher = np.zeros((len(dim), len(dim)))

        out0 = self.patternscan(patternnames)
        pi0 = (out0.intensity / np.sum(out0.intensity)).squeeze()

        # dpdc = np.zeros((len(pi0), len(dim)))
        dpdc = np.zeros((len(pi0), 3))

        for coord in range(len(dim)):
            dposa = np.zeros(3)
            dposa[dim[coord]] = eps / 2
            dposa2 = dposa.copy()
            dposa2[dim[coord]] = -eps / 2

            fl1.pos = posold + dposa
            out1 = self.patternscan(patternnames)
            fl1.pos =  posold + dposa2
            out2 = self.patternscan(patternnames)
            dpdc[:, dim[coord]] = ((out1.intensity / np.sum(out1.intensity, axis=0) - out2.intensity / np.sum(out2.intensity, axis=0)) / eps).squeeze()

        for coord in range(len(dim)):
            for coord2 in range(len(dim)):
                fh = dpdc[:, dim[coord]] * dpdc[:, dim[coord2]] / (pi0 + 1e-4)
                IFisher[coord, coord2] = np.sum(fh)

        fl1.pos = oldpar.pos
        fl1.remainingphotons = oldpar.remainingphotons

        try:
            crlb = np.linalg.inv(IFisher)
        except np.linalg.LinAlgError:
            crlb = np.ones_like(IFisher)*np.nan
        locprech = np.sqrt(np.diag(crlb)).T
        locprec = np.zeros(3)
        locprec[np.array(dim)] = locprech

        return locprec
    
    def calculateCRB(self, out, filter):
        locprecL = np.zeros(3)
        sigmaCRB1 = np.zeros(3)
        sigmaCRB = np.zeros(3)
        phot = 0
        seq = out.sequence

        for k, seqk in enumerate(seq):
            pattern = self.patterns[seqk]
            if pattern.type == "pattern":
                sh = self.calculateCRBpattern(seqk, dim=pattern.dim)
                if hasattr(out.loc, 'seq'):
                    ind = (out.loc.seq.squeeze() == k) & filter
                else:
                    ind = filter
                phot = np.mean(out.loc.phot[ind],axis=0)
                dim = np.array(pattern.dim)
                sigmaCRB1[dim] = sh[dim]
                sigmaCRB[dim] = sigmaCRB1[dim]/np.sqrt(phot)
                Lh = self.locprec(np.mean(phot), pattern.L)
                locprecL[dim] = Lh

        return sigmaCRB, sigmaCRB1, locprecL, phot

    def summarize_results(self, out, display=True, filter=None):
        """
        Parameters
        ----------
        out : object
            Structure created by a sequence
        display : bool
            If true, print the summary to standard out
        filter : np.typing.ArrayLike, optional
            Which localizations to use for statistics
        """
        if filter is None:
            filter = np.ones(out.loc.xnm.shape, dtype=bool)
        
        ind = filter.squeeze()
        photraw = out.raw.copy()
        xest = np.column_stack((out.loc.xnm[ind], out.loc.ynm[ind], out.loc.znm[ind]))
        # print("xnm ynm ", out.loc.xnm[:6], out.loc.ynm[:6], xest.shape)
        flpos = np.column_stack((out.loc.xfl1[ind], out.loc.yfl1[ind], out.loc.zfl1[ind]))
        # print("xest ", xest.shape, xest[:6])
        # print("flpos ", flpos.shape, flpos[:6])

        if hasattr(out.loc, 'seq'):
            sunique = np.unique(out.loc.seq[ind])
            phot = 0
            photch = []
            for k in range(len(sunique)):
                indu = (out.loc.seq == sunique[k]).squeeze()
                phot += np.mean(out.loc.phot[ind&indu])
                photch = np.hstack([photch, np.mean(photraw[ind&indu,:],axis=0)])
            photch = photch[photch!=-1]
        else:
            photch = np.mean(photraw[ind,:],axis=0)
            phot = np.sum(photch)

        sigmaCRB, sigmaCRB1, locprecL, _ = self.calculateCRB(out, ind)
        
        st = Summary()
        st.photch = photch
        st.bg_photons_gt = np.nanmean(out.bg_photons_gt[:,ind])
        st.phot = np.nanmean(phot)
        st.phot_signal = st.phot - st.bg_photons_gt
        st.pos = np.nanmean(xest, axis=0)
        st.std = np.nanstd(xest-flpos, axis=0)
        st.stdraw = np.nanstd(xest, axis=0)
        st.rmse = np.sqrt(np.nanmean((xest - flpos) ** 2, axis=0))
        # print("rmse ", st.rmse, st.rmse[:6])
        st.bias = np.nanmean(xest - flpos, axis=0)

        st.locp = locprecL
        st.sCRB = sigmaCRB #/ np.sqrt(st.phot)
        st.sCRB1 = sigmaCRB1
        st.duration = out.duration
        st.bg_photons_est = np.nanmean(out.bg_photons_est)
        
        if display:
            ff1, ff = "{:.1f}", "{:.2f}"

            if st.bg_photons_est != 0:
                bgtxt = f"bg_est: {ff1.format(st.bg_photons_est)}"
            else:
                bgtxt = ""
            
            if st.bg_photons_gt != 0:
                bgtxtg = f"bg_gt: {ff1.format(st.bg_photons_gt)}"
            else:
                bgtxtg = ""

            if hasattr(out.loc, "efo"):
                st.efo = np.mean(out.loc.efo)
                ratetxt = f"efo: {ff1.format(st.efo)}"
            else:
                st.photrate = np.mean(out.loc.photrate)
                ratetxt = f"photon rate: {ff1.format(st.photrate)}"

            print(f'photch: {", ".join(list(map(ff1.format, st.photch)))}',
                  f'mean(phot): {ff1.format(st.phot)}', 
                  bgtxt, bgtxtg, ratetxt, 
                  f'duration ms: {", ".join(list(map(ff1.format, st.duration)))}')
            print(f'std:  {", ".join(list(map(ff.format, st.std)))}', 
                  f'rmse: {", ".join(list(map(ff.format, st.rmse)))}', 
                  f'pos: {", ".join(list(map(ff.format, st.pos)))}',
                  f'bias: {", ".join(list(map(ff.format, st.bias)))}')
            print(f'locp: {", ".join(list(map(ff.format, locprecL)))}', 
                  f'sCRB: {", ".join(list(map(ff.format, st.sCRB)))}',
                  f'sCRB*sqrt(phot): {", ".join(list(map(ff1.format, st.sCRB1)))}')
        
        return st
    
    def scan_fov(self, 
                 seq,                    # cell of keys of patterns, elements
                 xcoords,                # coordinates to scan
                 maxlocs = 1000,         # how often to localize in each position
                 repetitions = 1,        # repetitions in each position
                 display = True,         # if to plot the results
                 dimplot = 0,            # which dimension is used to calculate statistics
                 dimscan = 0,            # which dimension is scanned 
                 title = "FoV scan",     # title of the output figure
                 fluorophorenumber = 0,  # which fluorophore is scanned. Use 0 for a FoV scan, use 1 to test effect of close-by fluorophore
                 ax1 = "std",            # std, bias, rmse, sCRB, pos: what to displax in the figure, arrray of strings. 
                 clearfigure = False,    # overwrite figure. if false: new plots are added
                 tag = None,             # name of a plot, used in the figure legend
                 linestyle = None):
        
        if isinstance(ax1, (str, int, float)):
            ax1 = [ax1]

        so = FOVScan()
        so.std = np.zeros((len(xcoords),3))
        so.rmse = np.zeros((len(xcoords),3))
        so.bias = np.zeros((len(xcoords),3))
        so.sCRB = np.zeros((len(xcoords),3))
        so.bias = np.zeros((len(xcoords),3))
        so.pos = np.zeros((len(xcoords),3))

        if isinstance(self.fluorophores, FlCollection):
            # print("collection")
            fl = self.fluorophores.flall[fluorophorenumber]
        else:
            fl = self.fluorophores

        # print(f"scan_fov fl.pos: {fl.pos}")

        coords = ["x", "y", "z"]
        posold = fl.pos[dimscan]
        for k in range(len(xcoords)):
            fl.pos[dimscan] = xcoords[k]
            # print(f"scan_fov fl.pos after shift: {fl.pos}")
            out = self.runSequence(seq, maxlocs=maxlocs,repetitions=repetitions)
            stats = self.summarize_results(out, display=False)
            so.std[k,:] = stats.std
            so.bias[k,:] = stats.bias
            so.rmse[k,:] = stats.rmse
            # print(f"scan_fov rmse after shift: {so.rmse[:,dimplot]}")
            so.sCRB[k,:] = stats.sCRB
            so.pos[k,:] = stats.pos
        so.stdrel = so.std/so.sCRB
        so.biasrel = so.bias/so.pos

        if display:
            import matplotlib.pyplot as plt
            axh = plt.gca()

            if axh.get_legend() is not None and not clearfigure:
                ltxt = [x.get_text() for x in axh.get_legend().texts]
            else:
                ltxt = []

            if clearfigure:
                axh.clear()

            # Left axis plotting
            ylab = f"{coords[dimplot]}: "
            for m, ax in enumerate(ax1):
                yval = getattr(so, ax)
                plt.plot(xcoords, yval[:, dimplot], linestyle)
                ylab += ax
                if m < (len(ax1) - 1):
                    ylab += "/"

            plt.ylabel(f"{ylab} (nm)")
            fl.pos[dimscan] = posold  # Assumes `fl_pos` is mutable and linked
            plt.xlabel(f"{coords[dimscan]} position (nm)")
            plt.title(title)

            # Construct legend
            ltxt += [f"{ax}: {tag}" for ax in ax1]
            plt.legend(ltxt)

            # Optional grid for "bias" conditions
            if any("bias" in ax for ax in ax1):
                plt.grid(True)

    def loadsequence(self, *args):
        raise UserWarning("No sequence loader implemented for Simulator class")