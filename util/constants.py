# Pre-calculated sizes of genomes.
GENOMES = {
  # Human Reference Genome
  'hg19' : {
    'size' : 3101804739, # 3,101,804,739 bp
    'chrs' : map(lambda x: 'chr%s' % x, [i for i in range(1,23)]+['X','Y','M']),
  },
}
