import glob
import os
import re
import sys
import pysam
import rmdup_calc
import fnmatch
import util

class aligned_stats(object):


  def __init__(self, id):
    self._id = id
    self._number_of_files = 0
    self._contig_rmdup_pct = 0
    self._readgroup_rmdup_pct = 0
    self._cov_8x = 0
    self._depth_of_coverage = 0

  @property
  def id(self):
    return self._id

  @property
  def number_of_files(self):
    """determined by counting readgroups in aligned bam"""
    return self._number_of_files

  @number_of_files.setter
  def number_of_files(self, num):
    self._number_of_files = num

  @property
  def contig_rmdup_pct(self):
    """aggregate of contig md.metrics files"""
    return self._contig_rmdup_pct

  @contig_rmdup_pct.setter
  def contig_rmdup_pct(self, pct):
    self._contig_rmdup_pct = pct

  @property
  def readgroup_rmdup_pct(self):
    """aggregate of readgroup md.metrics files"""
    return self._readgroup_rmdup_pct

  @readgroup_rmdup_pct.setter
  def readgroup_rmdup_pct(self, pct):
    self._readgroup_rmdup_pct = pct

  @property
  def cov_8x(self):
    """pct of genomes covered at minimum 8x"""
    return self._cov_8x

  @cov_8x.setter
  def cov_8x(self, cov):
    self._cov_8x = cov

  @property
  def depth_of_coverage(self):
    """x coverage depth of alignments"""
    return self._depth_of_coverage

  @depth_of_coverage.setter
  def depth_of_coverage(self, depth):
    self._depth_of_coverage = depth

  @property
  def rg_file(self):
    """RGfiles.txt location, abs path"""
    return self._rg_file

  @property
  def rg_file_path(self):
    return self._rg_file_path

  @property
  def rg_file_name(self):
    return self._rg_file_name
  
  @property
  def total_reads(self):
    return self._total_reads

  @property
  def mapped_reads(self):
    return self._mapped_reads

  @property
  def pct_mapped(self):
    return self._pct_mapped

  def process_readgroups(self):
    readgoup_file = self._rg_file_path
    num_lines = sum(1 for line in open(self._rg_file))
    self._number_of_files = num_lines * 2

  def calc_rmdup_stats(self):
    #try:
    self.contig_rmdup_pct = rmdup_calc.combine(glob.glob(
      self._rg_file_path + '/' + self._id + '.aln*metrics'))
    #except:
    #  print >>sys.stderr, "Can't open contig metrics files for " + self._id
    #  self.contig_rmdup_pct = 0.0
    #try:
    self.readgroup_rmdup_pct = rmdup_calc.combine(glob.glob(
      self._rg_file_path + '/' + self._id + '.' + self._id + '*metrics'))
    #except:
    #  print >>sys.stderr, "Can't open readgroup metrics files for " + self._id
    #  self.readgroup_rmdup_pct = 0.0

  def calc_pct_8x(self):
    try:
      filename = self._rg_file_path + "/" + self._id + ".coverage"
      cov_file = open(filename, 'r')
      lines = cov_file.readlines()
      self._cov_8x = (float(sum(int(line.split()[2]) for line in lines[8:])) / util.GENOMES['hg19']['size']) * 100.0
    except IOError as ioe:
      print >>sys.stderr, "No coverage file: " + filename
      self._cov_8x = 0.0

  def calc_depth_of_coverage(self):
    try:
      filename = self._rg_file_path + "/" + self._id + ".DoC"
      with open(filename, 'r') as f:
        self._depth_of_coverage = float(f.readline().strip())
    except IOError as ioe:
      print >>sys.stderr, "Unable to open .DoC file: " + filename
      self._depth_of_coverage = 0.0

  def calc_mapability(self):
    try:
      '''Grab the .flagstat or .flagstats file, not the flagstat.log file'''
      flagstat_file = self._rg_file_path + os.sep + self._id + '.bam.flagstat'
      d = util.flagstats(open(flagstat_file))
      self._total_reads = d['total_reads'][0]
      self._mapped_reads = d['mapped_reads'][0]
      self._pct_mapped = d['pct_mapped'][0]
    except (IOError, KeyError):
      print >>sys.stderr, "Unable to open/parse .flagstat file: " + flagstat_file

  @rg_file.setter
  def rg_file(self, rgfile):
    '''
    Triggers lots of stats gathering from peripheral files
    '''
    self._rg_file = rgfile
    self._rg_file_path = os.path.dirname(rgfile)
    self._rg_file_name = os.path.basename(rgfile)
    self.process_readgroups()
    self.calc_rmdup_stats()
    self.calc_pct_8x()
    self.calc_depth_of_coverage()
    self.calc_mapability()

