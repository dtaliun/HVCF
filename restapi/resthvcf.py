from flask import Flask, request, abort, jsonify
import PyHVCF
import time

app = Flask(__name__)

hvcf_file = 'test.h5'

hvcf = PyHVCF.HVCF()
hvcf.open(hvcf_file)

@app.route('/', methods = ['GET'])
def get_available_datasets():
   pass

@app.route('/populations', methods = ['GET'])
def get_populations():
   names = hvcf.get_sample_subsets()
   result = { 'names': [name for name in names] }
   j = jsonify(result)
   return j

@app.route('/populations/<population>', methods = ['GET'])
def get_samples_in_population(population):
   names = hvcf.get_samples_in_subset(str(population))
   result = { 'names': [name for name in names] }
   j = jsonify(result)
   return j

@app.route('/chromosomes', methods = ['GET'])
def get_chromosomes():
   names = hvcf.get_chromosomes()
   result = { 'names': [name for name in names] }
   j = jsonify(result)
   return j

@app.route('/chromosomes/<chromosome>', methods = ['GET'])
def get_chromosome(chromosome):
   if hvcf.has_chromosome(str(chromosome)):
      n_variants = hvcf.get_n_variants_in_chromosome(str(chromosome))
      start = hvcf.get_chromosome_start(str(chromosome))
      end = hvcf.get_chromosome_end(str(chromosome))
      result = { 'name': chromosome, 'n_variants': n_variants, 'first_variant_bp': start, 'last_variant_bp': end } 
   else:
      result = {}
   j = jsonify(result)
   return j

@app.route('/chromosomes/<chromosome>/variants', methods = ['GET'])
def get_variants_in_region(chromosome):
   start = request.args['start']
   end = request.args['end']

   variants = PyHVCF.VariantsVector()
   
   hvcf.extract_variants(str(chromosome), long(start), long(end), variants)

   result = {
      'chromosome': [str(chromosome)] * len(variants),
      'variant': [None] * len(variants),
      'alt': [None] * len(variants),
      'ref': [None] * len(variants),
      'position': [None] * len(variants)
   }

   for i, variant in enumerate(variants):
      result['variant'][i] = variant.name
      result['position'][i] = variant.position
      result['alt'][i] = variant.alt
      result['ref'][i] = variant.ref

   j = jsonify(result) 

   return j

@app.route('/chromosomes/<chromosome>/af', methods = ['GET'])
def get_af_in_region(chromosome):
   return 'af'

@app.route('/chromosomes/<chromosome>/<population>/haplotypes', methods = ['GET'])
def get_variant_haplotypes(chromosome, population):
   variant = request.args['variant']

   haplotypes = PyHVCF.UCharVector()
   
   hvcf.extract_haplotypes(str(chromosome), str(population), str(variant), haplotypes)

   result = { 'haplotypes' : [allele for allele in haplotypes] }
   j = jsonify(result)

   return j

@app.route('/chromosomes/<chromosome>/haplotypes', methods = ['GET'])
def get_sample_haplotypes(chromosome):
   sample = request.args['sample']
   start = request.args['start']
   end = request.args['end']

   haplotypes = PyHVCF.UCharVector()
   
   hvcf.extract_haplotypes(str(chromosome), str(sample), long(start), long(end), haplotypes)  

   result = { 'haploptypes' : [allele for allele in haplotypes] }
   j = jsonify(result)

   return j

@app.route('/chromosomes/<chromosome>/<population>/ld', methods = ['GET'])
def get_ld_in_region(chromosome, population):
   variant = request.args['variant']
   start = request.args['start']
   end = request.args['end']

   pairs = PyHVCF.PairsVector() 

   start_time = time.time()
   hvcf.compute_ld(str(chromosome), str(population), str(variant), long(start), long(end), pairs)
   elapsed_time = time.time() - start_time
   print 'HVCF request executed in ', elapsed_time, ' sec (', len(pairs), ')'

   start_time = time.time()
   result = { 
      'chromosome': [str(chromosome)] * len(pairs),
      'variant': [None] * len(pairs),
      'position': [None] * len(pairs),
      'r': [None] * len(pairs),
      'rsquare': [None] * len(pairs)
   }
   for i, pair in enumerate(pairs):
      result['variant'][i] = pair.name2
      result['position'][i] = pair.position2
      result['r'][i] = pair.r
      result['rsquare'][i] = pair.rsquare

   j = jsonify(result)       
   elapsed_time = time.time() - start_time
   print 'Response formatting executed in ', elapsed_time, ' sec (', len(pairs), ')' 

   return j

if __name__ == '__main__':
   app.run(host='0.0.0.0', port=5000)
