import re
from CGATReport.Tracker import *
from CGATReport.Utils import PARAMS as P
import CGATPipelines.PipelineTracks as PipelineTracks


# get from config file
UCSC_DATABASE = "hg19"
EXPORTDIR = "export"

###################################################################
###################################################################
###################################################################
###################################################################
# Run configuration script

EXPORTDIR = P.get('alleleSpecificExpression_exportdir',
                  P.get('exportdir', 'export'))
DATADIR = P.get('alleleSpecificExpression_datadir',
                P.get('datadir', '.'))
DATABASE = P.get('alleleSpecificExpression_backend',
                 P.get('sql_backend', 'sqlite:///./csvdb'))

TRACKS = PipelineTracks.Tracks(PipelineTracks.Sample).loadFromDirectory(
    glob.glob("%s/*.somatic.snvs.vcf.gz" % DATADIR),
    "(\S+).somatic.snvs.vcf.gz")

###########################################################################


def linkToUCSC(contig, start, end):
    '''build URL for UCSC.'''

    ucsc_database = UCSC_DATABASE
    link = "`%(contig)s:%(start)i..%(end)i <http://genome.ucsc.edu/cgi-bin/hgTracks?db=%(ucsc_database)s&position=%(contig)s:%(start)i..%(end)i>`_" \
        % locals()
    return link

###########################################################################


class project34Tracker(TrackerSQL):

    '''Define convenience tracks for plots'''

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, backend=DATABASE, **kwargs)


class imagesTracker(TrackerImages):

    '''Convience Tracker for globbing images for gallery plot'''
    def __init__(self, *args, **kwargs):
        Tracker.__init__(self, *args, **kwargs)
        if "glob" not in kwargs:
            raise ValueError("TrackerImages requires a:glob: parameter")
        self.glob = kwargs["glob"]
