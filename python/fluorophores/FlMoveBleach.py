from .FlMoving import FlMoving
from .FlBleach import FlBleach

class FlMoveBleach(FlMoving, FlBleach):
    def __init__(self, pos=[0,0,0], brightness=1000):
        FlBleach.__init__(self, pos, brightness)
        self.posind = 1
        self.posparameters = []
        self.posmode = None
        self.posfunction = None