import os
import sys
import re
import types
import itertools
import pandas as pd

from CGATReport.Tracker import *
from collections import OrderedDict as odict
from project34Report import *
from CGATReport.Utils import PARAMS as P

class imagesTracker(TrackerImages):

    '''Convience Tracker for globbing images for gallery plot'''
    def __init__(self, *args, **kwargs):
        Tracker.__init__(self, *args, **kwargs)
        if "glob" not in kwargs:
            raise ValueError("TrackerImages requires a:glob: parameter")
        self.glob = kwargs["glob"]

