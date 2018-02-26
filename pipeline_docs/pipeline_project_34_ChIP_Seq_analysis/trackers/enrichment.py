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


def order(frame, var):
    ''' reorder the columns of a dataframe so that the columns in 
    var come at the end'''
    if not type(var) == list:
        var = list(var)
    varlist = [w for w in frame.columns if w not in var]
    varlist.extend(var)
    frame = frame[varlist]

    return frame


class EnrichmentTTest(project34Tracker):

    def __call__(self, track, slice=None):
        table = "enrichment_ttest"
        statement = '''SELECT * FROM %(table)s;''' % locals()

        df = pd.DataFrame(self.getAll(statement))
        return order(df, ["p_value", "FDR"])


class EnrichmentBlockedAnova(project34Tracker):

    def __call__(self, track, slice=None):
        table = "enrichment_blockedANOVA"
        statement = '''SELECT * FROM %(table)s;''' % locals()

        return self.getAll(statement)

        df = pd.DataFrame(self.getAll(statement))
        return order(df, ["p_value", "FDR"])
