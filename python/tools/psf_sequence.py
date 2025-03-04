# For globals
from ..psfs import PsfDonut2D
from ..psfs import PsfGauss2D
from ..psfs import PsfVectorial

def psf_sequence(sq, psfvec=None, seq=None):
    """
    if vectorial in global: use only one PSF object, if also in Itr. other
    parameters are ignored
    """
    if seq is not None:  # abberior iteration: extrac pinhole from iterations, overwrite default.
        sq['global']['addpar']['pinhole'] = {'AU': seq['Itr'][-1]['Mode']['phDiaAU']}
    
    if 'global' in sq.keys() and 'Mode' in sq['global'].keys() and sq['global']['Mode'] == "PSF_vectorial":
        if psfvec is None:
            from ..psfs import PsfVectorial

            psfvec = PsfVectorial()
        
        fng = sq['global'].keys()

        for k in fng:
            ssth = sq['global'][k]
            if isinstance(ssth, dict):
                psfvec.setpar(**ssth)
            else:
                psfvec.setpar(**{k: ssth})
        
        if "pinhole" in sq['global']['addpar'].keys():
            psfvec.setpinhole(**sq['global']['addpar']['pinhole'])
    
    nItr = len(seq['Itr'])
    psfs = [None]*nItr
    phasemasks = [None]*nItr
    for k, itr in enumerate(sq['Itr']):
        if itr['Mode'] == 'PSF_vectorial':
            psfs[k] = psfvec
        else:
            # TODO: itr['Mode'] is supposed to be eval()ed. This means
            # we have to figure out where it's calling functions from
            # import them and then pull from there
            psfs[k] = globals()[itr['Mode']]
        phasemasks[k] = itr['par']

    return psfs, phasemasks