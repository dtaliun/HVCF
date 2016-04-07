import PyHVCF

def load_populations(subsets_file_name):
   populations = dict()

   with open(subsets_file_name) as infile:
      line = infile.readline()
      for line in infile:
         sample, subpop, pop, gender = line.rstrip().split('\t')
         if pop not in populations:
            populations[pop] = PyHVCF.NamesVector()
         populations.get(pop).append(sample)

   return populations

if __name__ == '__main__':
   vcf_name = '../test/1000G_phase3.ALL.chr20.LD_test.vcf.gz'
   subsets_file_name = '../test/integrated_call_samples_v3.20130502.ALL.panel'
   hvcf_name = 'test.h5'

   populations = load_populations(subsets_file_name)

   hvcf = PyHVCF.HVCF()
   hvcf.create(hvcf_name)
   hvcf.import_vcf(vcf_name)

   for pop, samples in populations.iteritems():
      hvcf.create_sample_subset(pop, samples)

   hvcf.close()

   hvcf.open(hvcf_name)
   pairs = PyHVCF.PairsVector()
   hvcf.compute_ld("20", "EUR", "20:11650214_G/A", 14403183, 55378791, pairs)

   print 'NAME_1 POSITION_1 NAME_2 POSITION_2 R R^2'
   for pair in pairs:
      print pair.name1, pair.position1, pair.name2, pair.position2, pair.r, pair.rsquare   

   hvcf.close()
