import PyHVCF
import argparse

argparser = argparse.ArgumentParser(description = 'Creates HVCF file from provided VCF file.')
argparser.add_argument('--import-gzvcf', metavar = 'file', dest = 'importGZVCFs', nargs = '+', required = True, help = 'Input VCF compressed with gzip.')
argparser.add_argument('--import-populations', metavar = 'file', dest = 'importPopulations', required = False, help = 'Input file with tab-delimited columns sample, pop, super_pop, gender')
argparser.add_argument('--out-hvcf', metavar = 'file', dest = 'outHVCF', required = True, help = 'Output HVCF.')

populations = dict()

def load_populations(file_name):
   with open(file_name) as infile:
      line = infile.readline()
      for line in infile:
         sample, subpop, pop, gender = line.rstrip().split('\t')
         if pop not in populations:
            populations[pop] = PyHVCF.NamesVector()
         populations.get(pop).append(sample)

if __name__ == '__main__':
   args = argparser.parse_args()

   if args.importPopulations:
      load_populations(args.importPopulations)

   hvcf = PyHVCF.HVCF()
   hvcf.create(args.outHVCF)

   for importGZVCF in args.importGZVCFs:
      hvcf.import_vcf(importGZVCF)   
      print 'Imported', importGZVCF

   for pop, samples in populations.iteritems():
      hvcf.create_sample_subset(pop, samples)
      print 'Populations created'

   hvcf.close()
