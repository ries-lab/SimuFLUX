import os
import pathlib
import json
from dataclasses import dataclass
from types import SimpleNamespace

import numpy as np

from .simulator import Simulator, Deadtimes, PatternScan, initlocs
from ..psfs import PsfVectorial
from ..psfs import PsfGauss2D
from ..psfs import PsfDonut2D
from ..tools import copyfields
from ..tools import psf_sequence
from ..tools import sumstruct
from ..tools import replace_in_list
from ..estimators import backgroundsubtractor

# To ensure they are found in globals()
from ..estimators import est_GaussLSQ1_2D
from ..estimators import est_donutLSQ1_2D
from ..estimators import est_qLSQiter2D

class SimSequencefile(Simulator):
    def __init__(self, fluorophores = None, background = 0, background_estimated = 0, loadfile = None):
        super().__init__(fluorophores, background, background_estimated, loadfile)

        self.sequence = {}
        self.estimators = []
        self.scoutingcoordinates = None
        self.psfvec = PsfVectorial()
        self.deadtimes = Deadtimes(point=0.011,
                                   pattern=0,
                                   estimator=0.015,
                                   positionupdate=0,
                                   localization=0.4)
        
    def loadsequence(self, *args):
        update = True

        for arg in args:
            if arg == "noupdate":
                update = False
                continue
            with open(arg, 'r') as fid:
                data = json.load(fid)
                self.sequence, _ = copyfields(self.sequence, data)

        if update:
            self.makepatterns()

    def makepatterns(self, psfs=None, phasemasks=None):
        if phasemasks is None:
            if 'PSF' in self.sequence.keys() and 'global' in self.sequence['PSF'].keys():
                psfs, phasemasks = psf_sequence(self.sequence['PSF'], self.psfvec, self.sequence)
            else:
                psfs, phasemasks = [], []
                psfs.append(PsfGauss2D())
                psfs.append(PsfDonut2D())
                phasemasks.append({'phasemask': 'gauss2D'})
                phasemasks.append({'phasemask': 'donut2D'})
                base_path = pathlib.Path(__file__).parent.parent
                self.loadsequence(os.path.join(base_path, 'settings', 'defaultestimators.json'), 'noupdate')

        itrs = self.sequence['Itr']
        for k, itr in enumerate(itrs):
            probecenter = False if itr['ccrLimit'] == -1 else True
            L = itr['patGeoFactor']*360
            kmin = np.minimum(k,len(phasemasks)-1)
            psf = psfs[kmin]
            phasemask = phasemasks[kmin]['phasemask']

            estimatorh = self.sequence['PSF']['Itr'][kmin]['estimator']
            while len(self.estimators) < (k+1):
                self.estimators.append([])
            self.estimators[k] = estimatorh

            if itr['Mode']['pattern'] == "hexagon":
                patternpoints = 6
            elif itr['Mode']['pattern'] == "square":
                patternpoints = 4
            elif itr['Mode']['pattern'] == "triangle":
                patternpoints = 3
            else:
                raise UserWarning(f"Pattern {itr['Mode']['id']} not implemented")
            
            # patterntime = (itr['patDwellTime']/itr['patRepeat'](1+probecenter*self.sequence['ctrDwellFactor']))*1e3  # ms
            pointdwelltime = itr['patDwellTime']/itr['patRepeat']*1e3/patternpoints
            if probecenter:
                pointdwelltime = [pointdwelltime, pointdwelltime*patternpoints*self.sequence['ctrDwellFactor']]
            # pointdwelltime=patterntime/(patternpoints+probecenter)

            laserpower = itr['pwrFactor']
            self.definePattern(f"itr{k}", psf, 
                               phasemask=phasemask, 
                               makepattern='orbitscan', 
                               orbitpoints=patternpoints, 
                               orbitL=L,
                               probecenter=probecenter,
                               pointdwelltime=pointdwelltime, 
                               laserpower=laserpower,
                               repetitions=itr['patRepeat'])
            
    def runSequenceintern(self):
        # compatible with MFSimulator: repetitions.
        itrs = self.sequence['Itr']
        maxiter = len(itrs)-1
        starttime = self.time

        deadtimes = self.deadtimes

        stickiness = self.sequence['stickiness']
        loclimit = self.sequence['locLimit']

        if 'maxOffTime' not in self.sequence.keys() or \
            ((not (isinstance(self.sequence['maxOffTime'], (int, float))) \
                and self.sequence['maxOffTime'] == 'unspecified')):
            maxOffTime = 3
        else:
            maxOffTime = self.sequence['maxOffTime']*1e6  # to us
        
        offtimestamp = 0

        if loclimit == -1:
            loclimit = 10000  # instead of inf. safety to avoid running forever

        loc = initlocs(loclimit,['xnm','ynm','znm','xfl1',
                    'yfl1','zfl1','phot','time','abortcondition', 'xgalvo', 'ygalvo', 'zgalvo',
                    'xeod','yeod','zeod','vld', 'ecc', 'eco', 'efc', 'efo', 'measuretime', 
                    'patternrepeat', 'cfr', 'itr', 'numitr', 'loccounter'])

        loccounter = 0

        numitr, itr = 0, 0
        abortphot = False
        stickinesscounter = 0
        xest = np.array([0,0,0])
        par = []
        raw = []
        flpos = []
        flint = []
        bg_photons_gt = []
        bg_photons_est = []
        # sofar = 1

        while numitr<loclimit and stickinesscounter<stickiness and not abortphot:
            itrname = f"itr{itr}"
            stickiness = self.sequence['stickiness']
            photsum, abortccr, abortphot = 0, False, False
            scanout = PatternScan()
            while photsum < itrs[itr]['phtLimit'] and stickinesscounter < stickiness and not abortphot:
                # repeat scanning until sufficient photons collected,
                # or abort condition met
                scanouth = self.patternscan(itrname)
                # print(f"scanouth.photrate: {scanouth.photrate.shape}")
                scanouth = backgroundsubtractor(scanouth, self.background_estimated)
                # print(f"scanouth.photrate no bg: {scanouth.photrate.shape}")
                scanout = sumstruct(scanout, scanouth)
                # print(f"scanout.photrate: {scanout.photrate.shape}")

                photsum = np.sum(scanout.phot)
                photbg = scanout.bg_photons_gt
                bgest = scanout.bgphot_est

                if itrs[itr]['ccrLimit'] > -1:
                    probecenter = True
                    cfr = scanout.phot[-1]/self.sequence['ctrDwellFactor']/np.sum(scanout.phot[:-1])
                    abortccr = cfr > itrs[itr]['ccrLimit']
                else:
                    probecenter = False
                    cfr = -1

                minphot = itrs[itr]['bgcThreshold']*scanouth.time.patterntotaltime*1e-6

                if np.sum(scanouth.phot) < minphot:
                    if scanouth.time.averagetime > (offtimestamp+maxOffTime):
                        abortphot = True
                else:
                    abortphot = False
                    offtimestamp = scanouth.time.averagetime
                
                if abortccr:
                    stickinesscounter += 1
                else:
                    stickinesscounter = 0

            if photsum > 0:
                loc.xgalvo[loccounter] = self.posgalvo[0]
                loc.ygalvo[loccounter] = self.posgalvo[1]
                loc.zgalvo[loccounter] = self.posgalvo[2]
                loc.xeod[loccounter] = self.posEOD[0]
                loc.yeod[loccounter] = self.posEOD[1]
                loc.zeod[loccounter] = self.posEOD[2]

                if not abortphot and not abortccr:  # recenter only for valid
                    # estimate position
                    patternpos = self.patterns[itrname].pos
                    L = itrs[itr]['patGeoFactor']*360  #nm
                    estimator = self.estimators[itr]
                    estf = globals()[estimator['function']]
                    estpar = estimator['par'].copy()
                    lenbe = (len(self.background_estimated)-1) if (isinstance(self.background_estimated, list) or isinstance(self.background_estimated, np.ndarray)) else 1
                    bg_est = self.background_estimated[np.minimum(itr,lenbe)] if lenbe > 1 else self.background_estimated
                    estpar = replace_in_list(estpar,
                                            'patternpos', patternpos,
                                            'L', L,
                                            'probecenter', probecenter,
                                            'background_est', bg_est,
                                            'iteration', itr)
                    xesth = estf(scanout.photrate, *estpar)
                    xest[np.array(estimator['dim'])] = xesth[np.array(estimator['dim'])]
                    self.time = self.time + deadtimes.estimator

                    xesttot = xest + self.posgalvo + self.posEOD
                    
                    # recenter
                    if (itr==maxiter) and not np.any(np.isnan(xesttot)):
                        dampf = 2**(-self.sequence['damping'])
                        xold = self.posgalvo.copy()
                        self.posgalvo = (1-dampf)*self.posgalvo+dampf*(xesttot)
                        self.posEOD = self.posEOD + xold - self.posgalvo + xest
                        self.time += deadtimes.positionupdate
                    elif not np.any(np.isnan(xesttot)):
                        self.posEOD = self.posEOD + xest
                    vldphotccr = True
                else:
                    xesttot = np.array([np.nan,np.nan,np.nan])
                    vldphotccr = False
                
                if itrs[itr]['ccrLimit'] > -1:  # probe center
                    loc.eco[loccounter,0] = np.sum(scanout.phot[:-1])
                    loc.ecc[loccounter,0] = np.sum(scanout.phot[-1])
                    
                    loc.efo=loc.eco/(np.sum(scanout.par.pattern.pointdwelltime[:-1]))/scanout.repetitions
                    loc.efc=loc.ecc/(scanout.par.pattern.pointdwelltime[-1])/scanout.repetitions
                else:
                    loc.eco[loccounter,0] = np.sum(scanout.phot)
                    loc.efo = loc.eco/(np.sum(scanout.par.pattern.pointdwelltime))/scanout.repetitions
                    loc.ecc[loccounter,0] = -1
                    loc.efc[loccounter,0] = -1
                    cfr = -1

                loc.xnm[loccounter,0] = xesttot[0]
                loc.ynm[loccounter,0] = xesttot[1]
                loc.znm[loccounter,0] = xesttot[2]
                loc.xfl1[loccounter,0] = scanout.flpos[0,0]/scanout.counter
                loc.yfl1[loccounter,0] = scanout.flpos[0,1]/scanout.counter
                loc.zfl1[loccounter,0] = scanout.flpos[0,2]/scanout.counter
                loc.time[loccounter,0] = scanout.time.averagetime
                loc.itr[loccounter,0] = itr
                loc.numitr[loccounter,0] = numitr
                loc.loccounter[loccounter,0] = loccounter
                loc.cfr[loccounter,0] = cfr
                loc.phot[loccounter,0] = photsum
                loc.vld[loccounter,0] = (stickinesscounter<stickiness) & vldphotccr
                loc.abortcondition[loccounter,0] = 1*(abortphot) + 2*(abortccr)
                loc.patternrepeat[loccounter,0] = scanout.counter
                loc.measuretime[loccounter,0] = scanout.time.patterntotaltime

                while len(raw) < (loccounter + 1):
                    raw.append([None])
                    flpos.append([None])
                    flint.append([None])
                    bg_photons_gt.append([None])
                    bg_photons_est.append([None])
                raw[loccounter] = scanout.phot
                flpos[loccounter] = scanout.flpos/scanout.counter
                flint[loccounter] = scanout.flint      
                bg_photons_gt[loccounter] = photbg
                bg_photons_est[loccounter] = np.sum(bgest)
                
                loccounter += 1
                try:
                    par[itr] = scanout.par
                except IndexError:
                    while len(par) < (itr+1):
                        par.append([None])
                    par[itr] = scanout.par

            itr += 1
            if itr > maxiter:
                itr += self.sequence['headstart']
            numitr += 1
        
        # Convert from raggedarray to complete array
        if len(flpos) > 0:
            rawsh1 = np.max([len(x) for x in raw])
            flpossh1 = np.max([x.shape[0] for x in flpos])
            flpossh2 = np.max([x.shape[1] for x in flpos])

            rawarr = np.zeros((loccounter, rawsh1))
            flposarr = np.zeros((loccounter, flpossh1, flpossh2))
            flintarr = np.zeros((loccounter, flpossh1))

            for i in range(loccounter):
                rawarr[i, :len(raw[i])] = raw[i].squeeze()
                flposarr[i, :flpos[i].shape[0], :flpos[i].shape[1]] = flpos[i].squeeze()
                flintarr[i, :flpos[i].shape[0]] = flint[i].squeeze()

            # Now delete everything in loc > loccounter
            for k, v in loc.__dict__.items():
                setattr(loc, k, v[:loccounter,:])

            out = SimpleNamespace(loc=loc,
                                  raw = rawarr,
                                  fluorophores = SimpleNamespace(pos=flposarr, int=flintarr),
                                  bg_photons_gt = np.vstack(bg_photons_gt).T,
                                  bg_photons_est = np.vstack(bg_photons_est),
                                  sequence = [f"itr{np.max(loc.itr).astype(int)}"],
                                  par = par,
                                  duration = self.time-starttime)
        else:
            out = None
                
        self.time += deadtimes.localization

        return out
    
    def runSequence(self, **kwargs):
        repetitions = kwargs.get('repetitions', 1)
        resetfluorophores = kwargs.get('resetfluorophores', False)

        starttime = self.time
        out = None

        for k in range(repetitions):
            if resetfluorophores:
                try:
                    self.fluorophores[0].reset()
                except TypeError:
                    self.fluorophores.reset()
            out2 = self.runSequenceintern()
            out = self.appendout(out, out2)
        
        if out is not None:
            out.duration = self.time-starttime

        return out
    
    def makescoutingpattern(self, fov, distance=0, show=False):
        if distance == 0:
            geofactor = self.sequence['field']['fldGeoFactor']
            distance = geofactor*360/640*self.sequence['Itr'][0]['wavelength']*1e9
        self.scoutingcoordinates = makehexgrid(fov, distance)
        if show:
            import matplotlib.pyplot as plt
            plt.gca()
            plt.plot(self.scoutingcoordinates[:,0], self.scoutingcoordinates[:,1],'o')

    def scoutingSequence(self, maxrep=100000, maxtime=1e6):
        timestart = self.time
        allbleached = False
        out = None
        print("scouting, progress: 0%", end="", flush=True)
        for reps in range(maxrep):
            if self.time > (timestart+maxtime):
                break
            prog=reps/maxrep*100
            try:
                prog = np.maximum(prog, (self.time-timestart)/maxtime*100)[0]
            except IndexError:
                prog = np.maximum(prog, (self.time-timestart)/maxtime*100)
            print(f"\rscouting, progress: {prog:3.0f}%", end="", flush=True)
            for pind in range(self.scoutingcoordinates.shape[0]):
                self.posgalvo[:2] = self.scoutingcoordinates[pind,:]
                self.posEOD = np.array([0,0,0])
                out2 = self.runSequence()
                out = self.appendout(out, out2)
                if self.fluorophores.allbleached:
                    allbleached = True
                    break
            if allbleached:
                break
        if out is not None:
            out.duration = self.time-timestart
        print("")  # newline

        return out

    def plotpositions(self, out, figure=None, coordinate=0, axis=None, xvalues="itr"):
        try:
            out.loc.loccounter
        except (NameError, AttributeError):
            raise ValueError("No localizations found")
        
        import matplotlib.pyplot as plt
        
        xnmn = ["xnm", "ynm", "znm"]
        xfln = ["xfl1", "yfl1", "zfl1"]
        xgn = ["xgalvo", "ygalvo", "zgalvo"]
        xen = ["xeod", "yeod", "zeod"]

        if axis is None:
            if figure is not None:
                f = plt.figure(figure)
            ax = plt.gca()
        else:
            ax = axis

        if xvalues == "itr":
            xv = out.loc.locounter
            xtxt = "time (itr)"
        elif xvalues == "time":
            xv = out.loc.time
            xtxt = "time (ms)"

        c = coordinate
        ax.plot(xv, getattr(out.loc, xnmn[c]), color='k')
        ax.plot(xv, getattr(out.loc, xfln[c]), color='r')
        ax.plot(xv, getattr(out.loc, xgn[c]), color=(0,0.7,0))
        ax.plot(xv, getattr(out.loc, xen[c]), color='b')
        ax.set_xlabel(xtxt)
        ax.set_ylabel('x position (nm)')
        ax.legend(['Estimated','Fluorophore','Galvo','EOD'])


def makehexgrid(roi, d):
    h = d*np.cos(np.pi/6)  # 30 degrees
    numpx = (roi[1,0]-roi[0,0])/h
    numpy = (roi[1,1]-roi[0,1])/d
    # TODO: Not totally comfortable with the +1
    pos = np.zeros((np.ceil(numpx*numpy).astype(int)+1,2))
    ind = 0
    numprange = numpy - np.arange(int(numpx+numpy)+1)
    for ll in numprange:
        for k in np.arange(numpx):
            x = k*h+roi[0,0]
            y = ll*d+k*d/2+roi[0,1]
            if x<= roi[1,0] and y<= roi[1,1] and x>=roi[0,0] and y>=roi[0,1]:
                pos[ind,:] = np.array([x,y])
                ind += 1
    pos = pos[:ind,:]

    return pos