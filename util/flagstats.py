#!/usr/bin/env python

import re
import logging

REGEX = r'^(?P<pass>\d+) \+ (?P<fail>\d+) (?P<type>.+)$'
REGEX = re.compile(REGEX)

def flagstats(ifs):
  '''
  Parse flagstats from the iterable ifs.
  '''
  # Return value.
  retval = dict()
  fields = {
    0   : 'total_reads',
    1   : 'duplicates',
    2   : 'mapped',
    3   : 'paired',
    4   : 'read1',
    5   : 'read2',
    6   : 'proper',
    7   : 'mate_mapped',
    8   : 'singletons',
    9   : 'mate_chr',
    10  : 'mate_chr_mapq5',
  }

  for i,line in enumerate(ifs):
    match = REGEX.match(line)
    if not match: raise ValueError('line %s malformed' % i)

    try: retval[fields[i]] = (match.group('pass'),match.group('fail'))
    except KeyError: logging.warning('line %s may be extra' % i)

  return retval

if __name__ == '__main__':
  import sys
  import json
  import argparse

  parser = argparse.ArgumentParser(description='flagstats parser')

  parser.add_argument('file',help='Flagstats file for parsing,')
  parser.add_argument('-d','--debug',dest='loglevel',action='store_const',
                      const=logging.DEBUG,default=logging.INFO,
                      help='Set logging level to DEBUG.')

  args = parser.parse_args()

  logging.basicConfig(
    level = args.loglevel,
    format = '%(asctime)s %(name)-6s %(levelname)-4s %(message)s',
  )

  sys.stdout.write('%s\n' % json.dumps(flagstats(open(args.file))))
