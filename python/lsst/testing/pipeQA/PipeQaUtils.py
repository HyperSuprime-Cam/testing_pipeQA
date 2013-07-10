import numpy as num

#
### HELPER CLASSES
#

class SdqaMetric(object):
    @staticmethod
    def MAX(val, maxval): return val <= maxval
    @staticmethod
    def MIN(val, minval): return val >= minval

    def __init__(self, value = 0.0, limits = {}, comment = None):
        self.value   = value
        self.limits  = limits
        self.comment = comment

    def setValue(self, value):
        self.value = value
        
    def evaluate(self):
        for lim in self.limits.keys():
            if not lim(self.value, self.limits[lim]):
                return False
        return True


#
### HELPER FUNCTIONS
#


# Example on how to use the butler
def getAllKeys(inButler):
    allkeys = []
    visits = inButler.queryMetadata('raw', 'visit')
    for visit in visits:
        rafts = inButler.queryMetadata('raw', 'raft', visit=visit)
        for raft in rafts:
            sensors = inButler.queryMetadata('raw', 'sensor', visit=visit, raft=raft)
            for sensor in sensors:
                allkeys.append(dict(visit=visit, raft=raft, sensor=sensor))
    return allkeys

def getAllKeysOpt(opt, butler):
    allkeys = []
    if len(opt.visit) == 1 and len(opt.raft) == 1 and len(opt.sensor) == 1:
        allkeys.append(dict(visit = opt.visit[0], raft = opt.raft[0], sensor = opt.sensor[0]))
    elif len(opt.visit) == 1 and len(opt.raft) == 1:
        sensors = butler.queryMetadata('raw', 'sensor', visit=opt.visit[0], raft=opt.raft[0])
        for sensor in sensors:
            allkeys.append(dict(visit = opt.visit[0], raft = opt.raft[0], sensor = sensor))
    elif len(opt.visit) == 1:
        rafts = butler.queryMetadata('raw', 'raft', visit = opt.visit[0])
        for raft in rafts:
            sensors = butler.queryMetadata('raw', 'sensor', visit = opt.visit[0], raft = raft)
            for sensor in sensors:
                allkeys.append(dict(visit = opt.visit[0], raft = raft, sensor = sensor))
    else:
        allkeys = getAllKeys(butler)
    return allkeys

def pointInsidePolygon(x,y,poly):
    n = len(poly)
    inside = False

    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x,p1y = p2x,p2y
    return inside
    

def sigIQR(data, min = None, max = None):
    data  = num.sort(data)
    if (min != None) and (max != None):
        idx   = num.where( (data > min) & (data < max) )
        data  = data[idx]

    if len(data) == 0:
        return 0.0
    
    d25 = data[int(0.25 * len(data))]
    d75 = data[int(0.75 * len(data))]
    return 0.741 * (d75 - d25)

