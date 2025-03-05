import numpy as np

from .Fluorophore import Fluorophore

class FlMoving(Fluorophore):
    def __init__(self, pos=[0,0,0], brightness=1000):
        super().__init__(pos, brightness)
        self.posind = 0
        self.posparameters = []
        self.posmode = None
        self.posfunction = None

    def position(self, time, props=None):
        posout = np.array([0,0,0])
        posmode = self.posmode
        if posmode == "static":
            posh = self.pos
        elif posmode in ["trace", "diffusion", "steps"]:
            # XXXX add test if ind. already too high (e.g. 2x
            # position asked. If asked 2x same time: should work)
            pos = self.pos
            posind = np.maximum(0, self.posind-1)
            if pos[posind,0] > time:
                result = np.where(pos[:,0]<time)[0]
                if len(result) > 0:
                    posind = result[-1]
                else:
                    posind = 0
            while pos[posind,0] < time:
                if posind >= pos.shape[0]-1:  # reached end
                    self.extendtrace()
                    pos = self.pos
                posind += 1
            posh = pos[posind,1:]
            self.posind = posind
        elif posmode == "function":
             posfunction = self.posfunction
             for k in range(len(posfunction),-1,-1):
                  posh[k] = posfunction[k](time)
        else:
             raise ValueError(f"Posmode {posmode} not implemented.")
        
        posout[:len(posh)] = posh

        return posout, True
    
    def makediffusion(self, 
                      D, # um^2/s
                      dt,  # us
                      startpos=np.array([0,0,0]), 
                      dim=2,  # number of dimensions to simulate
                      numpoints=10000, 
                      boundarybox=None): # nm, half length of box, periodic boundary conditions
        # print("makediffusion")
        kwargs = {'startpos': startpos, 
                  'dim': dim, 
                  'numpoints': numpoints, 
                  'boundarybox': boundarybox}
        self.posparameters = [D, dt, kwargs]

        # to obtain nm^2/ms multiply by 1000;
        Dstep = D*dt*1000
        time = np.atleast_2d((np.arange(numpoints)+1)*dt)
        jumps = np.random.randn(numpoints, dim)*Dstep
        # print("max jump: ", jumps.max(0))
        # print(f"startpos: {startpos[:dim]}")
        pos = np.cumsum(jumps, axis=0) + startpos[:dim]
        # print("max pos: ", pos.max(0))
        if boundarybox is not None:  # periodic boundary conditions
            if not (isinstance(boundarybox, list) or isinstance(boundarybox, np.ndarray)) or len(boundarybox) == 1:
                boundarybox = boundarybox*np.ones((1,dim))
            if not isinstance(boundarybox, np.ndarray):
                boundarybox = np.array(boundarybox)
            boundaryboxpos = boundarybox/2
            pos2 = pos + boundaryboxpos
            pos3 = np.mod(pos2, 2*boundaryboxpos)
            pos = pos3 - boundaryboxpos
        # print("max pos after boundary box: ", pos.max(0))
        
        self.pos = np.hstack([time.T, pos])
        self.posmode = "diffusion"

    def makesteps(self, 
                  stepsize,  # nm
                  dwelltime,  # ms
                  dt, 
                  startpos=np.array([0,0,0]), 
                  dim=1, 
                  numpoints=10000,
                  angle=0):
        # print("makesteps")
        kwargs = {'startpos': startpos, 
                  'dim': dim, 
                  'numpoints': numpoints, 
                  'angle': angle}
        self.posparameters = [stepsize, dwelltime, dt, kwargs]

        time = np.atleast_2d(((np.arange(numpoints)+1)*dt)).T

        numjumps = np.ceil(numpoints*dt/(dwelltime)*2+5).astype(int)
        jmp=np.random.exponential(dwelltime/dt,(numjumps,1))
        jumppos=np.cumsum(jmp)
        jumppos = jumppos[jumppos<=len(time)]
        hist, bin_edges = np.histogram(jumppos, bins=np.arange(len(time) + 1))
        xjh = hist * stepsize
        xpos=np.cumsum(xjh.T)
        angle_rad = angle*np.pi/180
        
        R = np.array([[np.cos(angle_rad) - np.sin(angle_rad)], 
                      [np.sin(angle_rad), np.cos(angle_rad)]])
        pos = (R*np.hstack([xpos, np.zeros(xpos.shape)]).T).T + startpos[:dim]
        self.pos = np.hstack([time, pos])
        self.posmode = "steps"

    def extendtrace(self):
        currentpos = self.pos[-1,:]
        oldpos = self.pos.copy()
        kwargs = self.posparameters[-1]
        if self.posmode == "diffusion":
            self.makediffusion(self.posparameters[0], 
                               self.posparameters[1],
                               dim=kwargs['dim'],
                               numpoints=kwargs['numpoints'],
                               startpos=currentpos[1:],
                               boundarybox=kwargs['boundarybox'])
        elif self.posmode == "steps":
            self.makesteps(self.posparameters[0],
                           self.posparameters[1],
                           self.posparameters[2],
                           dim=kwargs['dim'],
                           numpoints=kwargs['numpoints'],
                           startpos=currentpos[1:],
                           angle=kwargs['angle'])
        
        self.pos[:,0] += currentpos[0]  # time continues
        self.pos = np.vstack([oldpos, self.pos])
