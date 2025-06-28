import numpy as np

def test_backgroundsubtractor():
    from ..estimators import backgroundsubtractor
    from types import SimpleNamespace

    patDwellTime = 0.001
    patRepeat = np.random.randint(1,7)
    patternpoints = np.random.choice([6,4,3])
    pdt = patDwellTime*patRepeat*1e3/patternpoints*np.ones(patternpoints)
    phot = np.random.poisson(np.random.randint(1, 1000, patternpoints))
    reps = np.random.randint(1, 10)
    pr = phot/pdt

    out = SimpleNamespace(repetitions=reps,
                          pointdwelltime=pdt,
                          photrate=pr)
    
    bg = np.random.rand()*100
    
    out2 = backgroundsubtractor(out, bg)

    np.testing.assert_allclose(out2.photrate, pr - bg*reps)
    np.testing.assert_allclose(out2.bgphot_est, pdt*bg*reps)

# def test_est_donutLSQ1_2D():
