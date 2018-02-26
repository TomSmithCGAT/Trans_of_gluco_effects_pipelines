from CGATReport.Tracker import *
from project34Report import *

###########################################################################


class miRNACountsGenome(project34Tracker):
    def __call__(self, track, slice=None):
        statement = '''SELECT * FROM miRNA_feature_counts_bwa_counts'''
        return self.getAll(statement)


class piRNACountsGenome(project34Tracker):
    def __call__(self, track, slice=None):
        statement = '''SELECT * FROM piRBase_rn5_feature_counts_bwa_counts'''
        return self.getAll(statement)


class tRFCountsGenome(project34Tracker):
    def __call__(self, track, slice=None):
        statement = '''
        SELECT * FROM tRNA_no_duplicate_transcripts_feature_counts_bwa_counts'''
        return self.getAll(statement)


class matureMiRNACountsIterative(project34Tracker):
    def __call__(self, track, slice=None):
        statement = '''SELECT * FROM counts_mature_miRNA'''
        return self.getAll(statement)


class piRNACountsIterative(project34Tracker):
    def __call__(self, track, slice=None):
        statement = '''SELECT * FROM counts_piRNA'''
        return self.getAll(statement)


class tRFCountsIterative(project34Tracker):
    def __call__(self, track, slice=None):
        statement = '''SELECT * FROM counts_tRF_sequence '''
        return self.getAll(statement)


class tRFmapSequenceLoci(project34Tracker):
    def __call__(self, track, slice=None):
        statement = '''SELECT * FROM tRNA_sequences'''
        return self.getAll(statement)

