#!/usr/bin/env python
from discuspy import *
from structure import *
from pdf import *
from powder import *
from copy import copy
import numpy as np
setup()

class MMC(object):
    """
    MMC calculation.
    """
    def __init__(self,structure):
        self.structure=structure

class Fourier(object):
    """
    Fourier calculation.
    """
    def __init__(self,structure):
        self.structure=structure

