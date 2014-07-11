################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_snps.py 2870 2010-03-03 10:20:29Z andreas $
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#################################################################################
"""
===========================
Pipeline rrbs
===========================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

A pipeline template.

Overview
========

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file. 

The sphinxreport report requires a :file:`conf.py` and :file:`sphinxreport.ini` file 
(see :ref:`PipelineDocumenation`). To start with, use the files supplied with the
:ref:`Example` data.

Input
-----

Optional inputs
+++++++++++++++

Requirements
------------

The pipeline requires the results from :doc:`pipeline_annotations`. Set the configuration variable 
:py:data:`annotations_database` and :py:data:`annotations_dir`.

On top of the default CGAT setup, the pipeline requires the following software to be in the 
path:

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|                    |                   |                                                |
+--------------------+-------------------+------------------------------------------------+

Pipeline output
===============

The major output is in the database file :file:`csvdb`.

Example
=======

Example data is available at http://www.cgat.org/~andreas/sample_data/pipeline_rrbs.tgz.
To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_rrbs.tgz
   tar -xvzf pipeline_rrbs.tgz
   cd pipeline_rrbs
   python <srcdir>/pipeline_rrbs.py make full

.. note:: 
   For the pipeline to run, install the :doc:`pipeline_annotations` as well.

Glossary
========

.. glossary::


Code
====

"""
from ruffus import *

import sys, glob, gzip, os, itertools, re, math, types, collections, time
import optparse, shutil
import sqlite3

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Database as Database

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################

# load options from the config file
import CGAT.Pipeline as P
P.getParameters( 
    ["%s.ini" % __file__[:-len(".py")],
     "../pipeline.ini",
     "pipeline.ini" ] )

PARAMS = P.PARAMS
PARAMS_ANNOTATIONS = P.peekParameters( PARAMS["annotations_dir"],
                                       "pipeline_annotations.py" )

###################################################################
###################################################################
## Helper functions mapping tracks to conditions, etc
###################################################################
import CGATPipelines.PipelineTracks as PipelineTracks

# define some tracks if needed
TRACKS = PipelineTracks.Tracks( PipelineTracks.Sample ).loadFromDirectory( 
    glob.glob("*.ini" ), "(\S+).ini" )


###################################################################
###################################################################
###################################################################
def connect():
    '''connect to database.

    Use this method to connect to additional databases.

    Returns a database connection.
    '''

    dbh = sqlite3.connect( PARAMS["database"] )
    statement = '''ATTACH DATABASE '%s' as annotations''' % (PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute( statement )
    cc.close()

    return dbh


#########################################################################
#########################################################################
#########################################################################
# Read mapping
#########################################################################


SEQUENCESUFFIXES = ("*.fastq.1.gz",
                    "*.fastq.gz",
                    "*.sra",
                    "*.export.txt.gz",
                    "*.csfasta.gz",
                    "*.csfasta.F3.gz",
                    )

SEQUENCEFILES = tuple([os.path.join(DATADIR, suffix_name)
                      for suffix_name in SEQUENCESUFFIXES])


SEQUENCEFILES_REGEX = regex(
    r".*/(\S+).(fastq.1.gz|fastq.gz|sra|csfasta.gz|csfasta.F3.gz|export.txt.gz)")

###################################################################
###################################################################
###################################################################
## Prepare reads
###################################################################

@follows(mkdir("trimmed_reads.dir"))
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           r"trimmed_reads.dir/\1.fastq.1.gz")
def runTrimGalore(infile,outfile):
    '''trim reads with trim galore'''

    job_options = "-pe dedicated %i -R y -l mem_free=%s" % (
        PARAMS["trim_galore_threads"],
        PARAMS["trim_galore_memory"])

    to_cluster = True
    m = PipelineRRBS.TrimGalore()
        #remove_non_unique=PARAMS["remove_non_unique"],
        #strip_sequence=PARAMS["strip_sequence"],
        #read_group_header=PARAMS["bwa_read_group_header"])

    statement = m.build((infile,), outfile)
    P.run()

###################################################################
###################################################################
###################################################################
## Map reads
###################################################################

@follows(mkdir("bismark.dir"),
         runTrimGalore)
@transform(runTrimGalore,
           SEQUENCEFILE_REGEX,
           r"bismark.dir/\1.bam")
def runBismark(infile,outfile):
    '''map reads with Bismark'''

    job_options = "-pe dedicated %i -R y -l mem_free=%s" % (
        PARAMS["trim_galore_threads"],
        PARAMS["trim_galore_memory"])

    to_cluster = True
    m = PipelineRRBS.TrimGalore()
        #remove_non_unique=PARAMS["remove_non_unique"],
        #strip_sequence=PARAMS["strip_sequence"],
        #read_group_header=PARAMS["bwa_read_group_header"])

    statement = m.build((infile,), outfile)
    P.run()


###################################################################
###################################################################
###################################################################
## primary targets
###################################################################
@follows( loadDummyTask )
def full(): pass

###################################################################
###################################################################
###################################################################
## primary targets
###################################################################

@follows( mkdir( "report" ) )
def build_report():
    '''build report from scratch.'''

    E.info( "starting report build process from scratch" )
    P.run_report( clean = True )

@follows( mkdir( "report" ) )
def update_report():
    '''update report.'''

    E.info( "updating report" )
    P.run_report( clean = False )

@follows( update_report )
def publish():
    '''publish report and data.'''

    E.info( "publishing report" )
    P.publish_report()

if __name__== "__main__":

    # P.checkFiles( ("genome.fasta", "genome.idx" ) )
    sys.exit( P.main(sys.argv) )
