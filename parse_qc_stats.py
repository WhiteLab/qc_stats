import sys
import os
import re
import glob
import argparse
import subprocess
import pysam
import rmdup_calc
import AlignedStats
import util
import HgacRethinkdb
import time
from collections import defaultdict
"""
Generate data accounting, and alignment stats on "
1) HGAC rethinkdb - metadata describing the samples, and what was generated
2) Raw data directory - seq files copied from bionimbus
3) location of unaligned .bam files
4) location of alignment RGfiles.txt files

"""

# globals
TABLE = defaultdict(lambda: defaultdict(int)) # contains gathered metadata
HEADERS = ['BID', 'Sample Type', 'HGAC_files', 'Beagle_files', 'Unaligned_files', 
    'Aligned_files', 'Unaligned_reads', 'Contig_rmdup_pct', 
    'Readgroup_rmdup_pct', 'total_reads', 'mapped', '%mapped', 
    'Sequencing_depth', 'Doc', '%_8x', 'Action']

HGAC_TABLE = {'Sample Type':HEADERS[1], 'Total # of Sequence Files from PE runs (AD)':HEADERS[2]}
TABLE_HGAC = {HEADERS[1]:'Sample Type', HEADERS[2]:'Total # of Sequence Files from PE runs (AD)'}


def count_seq_files(seq_dir):
  """
  Count number of paired-end seq files per library
  """
  for seqfile in glob.glob(seq_dir + '/2*_[1-8]_[12]_sequence.txt.gz'):
    match = re.search(r'(2\d+-\d+)_', os.path.basename(seqfile))
    if match:
      bid = match.group(1)
      TABLE[bid][HEADERS[3]] += 1


def file_counts_by_readgroups(dir_path):
  """
  Count the files by readgroups from the unaligned .bams for each library
  :param dir_path: Path to dir containing unaligned sample directories
  :param seq_dict: dict of library ID's
  :returns: dict of readgroup counts by library
  """
  for bid in TABLE:
    try:
      bam_file = dir_path + "/" + bid + "/" + bid + '.bam'
      bam = pysam.Samfile(bam_file, check_sq=False)
      TABLE[bid][HEADERS[4]] = len(bam.header['RG']) * 2
    except IOError:
      print >>sys.stderr, "[parse_qc_stats] - Unable to count readgroups in bamfile: " + bam_file
   

def parse_unaligned_reads(dir_path):
  for bid in TABLE:
    try:
      flagstats_file = dir_path + os.sep + bid + os.sep + bid + '.bam.flagstats'
      stats = util.flagstats(open(flagstats_file))
      TABLE[bid][HEADERS[6]] = stats['total_reads'][0]
      seq_depth = (stats['total_reads'] * 100.0) / util.GENOMES['hg19']['size']
      TABLE[bid][HEADERS[12]] = "%.2f" % seq_depth
    except ValueError:
      print >>sys.stderr, "[parse_qc_stats] - Looks like flagstats file is empty: " + flagstats_file
    except KeyError:
      print >>sys.stderr, "[parse_qc_stats] - Looks like flagstats file is non-existant: " + flagstats_file


def parse_aligned_stats(rg_file_paths):
  '''
  All aligned datasets will have an RGfiles.txt, provided by the rg_file_paths
  '''
  f = open(rg_file_paths, 'r')
  for file_path in f:
    sys.stderr.write('.')
    file_path = file_path.strip()
    data_path, rg_file = os.path.split(file_path)
    bid = os.path.split(data_path)[1]
    a_stats = AlignedStats.aligned_stats(bid)
    a_stats.rg_file = file_path # triggers loading of many alignment stats

    TABLE[bid][HEADERS[5]] = a_stats.number_of_files # load the TABLE
    TABLE[bid][HEADERS[7]] = "%.4f" % a_stats.contig_rmdup_pct
    TABLE[bid][HEADERS[8]] = "%.4f" % a_stats.readgroup_rmdup_pct
    TABLE[bid][HEADERS[9]] = a_stats.total_reads
    TABLE[bid][HEADERS[10]] = a_stats.mapped_reads
    TABLE[bid][HEADERS[11]] = "%.2f" % a_stats.pct_mapped
    TABLE[bid][HEADERS[13]] = "%.2f" % a_stats.depth_of_coverage
    TABLE[bid][HEADERS[14]] = "%.2f" % a_stats.cov_8x
  sys.stderr.write('\n')

def load_hgac_ids(keyfile):
  '''
  clear the TABLE of metadata
  reload with BIDs and # of HGAC files
  '''
  TABLE.clear()
  hr = HgacRethinkdb.HgacRethinkdb('igsbimg.uchicago.edu', 'hgac', keyfile)
  table = 'FONBC_LibraryTracker'
  cols = ['BID', 'Sample Type', 'Total # of Sequence Files from PE runs (AD)']
  records = hr.get_columns(table, cols) # returns list
  for record in records:
    for i in range(1, len(cols)):
      TABLE[record[cols[0]]][HGAC_TABLE[cols[i]]] = record[cols[i]]

def determine_next_action():
  for bid in TABLE:
    sys.stderr.write('.')
    action = '-'
    if(TABLE[bid][HEADERS[4]] > TABLE[bid][HEADERS[5]]):
      action = 'Align'
    if(TABLE[bid][HEADERS[3]] > TABLE[bid][HEADERS[4]]):
      action = 'Bamify + Align'
    TABLE[bid][HEADERS[15]] = action
  sys.stderr.write('\n')


def regurgitate():
  print "\t".join(HEADERS)
  for bid in sorted(TABLE):
    record = bid
    for H in HEADERS[1:]:
      record += "\t" + str(TABLE[bid][H])
    print record


def main():
  """
  parse command line args
  """
  p = argparse.ArgumentParser()
  p.add_argument('-r', dest='readgroups_file', action='store', required=True, 
      help='List of all RGfiles.txt file paths, used to locate alignment stats per sample')
  p.add_argument('-k', dest='keyfile', action='store', required=True,
      help='RethinkDB keyfile required to query HGAC metadata')
  p.add_argument('-s', dest='seq_dir', action='store', required=True,
      help='Directory containing all raw seq data files')
  p.add_argument('-u', dest='unaligned_dir', action='store', required=True,
      help='Directory containing per-sample directories with unaligned bam files')

  now = time.strftime("%c")

  args = p.parse_args()
  load_hgac_ids(args.keyfile) # into TABLE dict
  print >>sys.stderr, "Counting sequence files"
  count_seq_files(args.seq_dir)
  print >>sys.stderr, "counting unaligned files"
  file_counts_by_readgroups(args.unaligned_dir)
  print >>sys.stderr, "Counting unaligned reads"
  parse_unaligned_reads(args.unaligned_dir)
  print >>sys.stderr, "Gathering aligned stats"
  parse_aligned_stats(args.readgroups_file)
  print >>sys.stderr, "Determining next actions"
  determine_next_action()

  print now
  regurgitate()


if __name__ == '__main__':
  main()
