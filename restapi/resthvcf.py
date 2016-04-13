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

@app.route('/populations/<name>', methods = ['GET'])
def get_samples_in_population(name):
   names = hvcf.get_samples_in_subset(str(name))
   result = { 'names': [name for name in names] }
   j = jsonify(result)
   return j

@app.route('/chromosomes', methods = ['GET'])
def get_chromosomes():
   names = hvcf.get_chromosomes()
   result = { 'names': [name for name in names] }
   j = jsonify(result)
   return j

@app.route('/chromosomes/<name>', methods = ['GET'])
def get_chromosome(name):
   if hvcf.has_chromosome(str(name)):
      n_variants = hvcf.get_n_variants_in_chromosome(str(name))
      start = hvcf.get_chromosome_start(str(name))
      end = hvcf.get_chromosome_end(str(name))
      result = { 'name': name, 'n_variants': n_variants, 'first_variant_bp': start, 'last_variant_bp': end } 
   else:
      result = {}
   j = jsonify(result)
   return j

@app.route('/chromosomes/<name>/variants', methods = ['GET'])
def get_variants_in_region(name):
   start = request.args['start']
   end = request.args['end']

   variants = PyHVCF.VariantsVector()
   
   hvcf.extract_variants(str(name), long(start), long(end), variants)

   result = {
      'chromosome': [str(name)] * len(variants),
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

@app.route('/chromosomes/<name>/af', methods = ['GET'])
def get_af_in_region(name):
   return 'af'

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
