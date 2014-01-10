#!/usr/bin/env python
def getRadiation(r):
    rad = {1:"X-Rays", 2:"Neutrons", 3:"Electrons"}
    if rad.has_key(r): 
        return rad[r]
    else:
        return "Invalid"

def getRadiationPar(rad):
    if rad[0].lower()=='n':
        r=2
        lxray=False
    elif rad[0].lower()=='e':
        r=3
        lxray=True
    else:
        r=1
        lxray=True
    return r,lxray
