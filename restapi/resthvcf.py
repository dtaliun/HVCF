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
   result = { 'populations': [name for name in names] }
   j = jsonify(result)
   return j

@app.route('/populations/<population>', methods = ['GET'])
def get_samples_in_population(population):
   names = hvcf.get_samples_in_subset(str(population))
   result = { 'samples': [name for name in names] }
   j = jsonify(result)
   return j

@app.route('/chromosomes', methods = ['GET'])
def get_chromosomes():
   names = hvcf.get_chromosomes()
   result = { 'chromosomes': [name for name in names] }
   j = jsonify(result)
   return j

@app.route('/chromosomes/<chromosome>', methods = ['GET'])
def get_chromosome(chromosome):
   if hvcf.has_chromosome(str(chromosome)):
      start_bp = request.args.get('startbp', None)
      end_bp = request.args.get('endbp', None)
      
      if start_bp and end_bp:
         variants = PyHVCF.Variants()

         hvcf.extract_variants(str(chromosome), long(start_bp), long(end_bp), variants)
         result = { 
            'chromosome' : str(chromosome),
            'region_start_bp': long(start_bp),
            'region_end_bp': long(end_bp),
            'number_of_variants': len(variants),
            'variants': [None] * len(variants)
         }

         for i, variant in enumerate(variants):
            result['variants'][i] = {
               'name': variant.name,
               'position': variant.position,
               'reference_allele': variant.ref,
               'alternate_allele': variant.alt
            }
      else:
         n_variants = hvcf.get_n_variants_in_chromosome(str(chromosome))
         start = hvcf.get_chromosome_start(str(chromosome))
         end = hvcf.get_chromosome_end(str(chromosome))
         result = { 'number_of_variants': n_variants, 'first_variant_bp': start, 'last_variant_bp': end } 
   else:
      result = {}

   j = jsonify(result)
   return j

@app.route('/frequency', methods = ['GET'])
def get_frequencies():
   population = request.args['population']
   chromosome = request.args['chromosome']
   start_bp = request.args['startbp']
   end_bp = request.args['endbp']

   frequencies = PyHVCF.Frequencies()

   hvcf.compute_frequencies(str(population), str(chromosome), long(start_bp), long(end_bp), frequencies)

   result = {
      'population': str(population),
      'chromosome': str(chromosome),
      'region_start_bp': long(start_bp),
      'region_end_bp': long(end_bp),
      'number_of_variants': len(frequencies),
      'variants': [None] * len(frequencies)
   }

   for i, variant in enumerate(frequencies):
      result['variants'][i] = {
         'name': variant.name,
         'position': variant.position,
         'reference_allele': variant.ref,
         'alternate_allele': variant.alt,
         'reference_frequency': variant.ref_af,
         'alternate_frequency': variant.alt_af
      }

   j = jsonify(result)

   return j

@app.route('/ld', methods = ['GET'])
def get_ld():
   population = request.args['population']
   chromosome = request.args['chromosome']
   start_bp = request.args['startbp']
   end_bp = request.args['endbp']
   lead_variant = request.args.get('variant', None)
   compact = request.args.get('compact', None)

   pairs = PyHVCF.Pairs()

   start_time = time.time()
   if not lead_variant:
      hvcf.compute_ld(str(population), str(chromosome), long(start_bp), long(end_bp), pairs)
   else:
      hvcf.compute_ld(str(population), str(chromosome), str(lead_variant), long(start_bp), long(end_bp), pairs)
   elapsed_time = time.time() - start_time
   print 'HVCF request executed in ', elapsed_time, ' sec (', len(pairs), ')'

   start_time = time.time()
   if compact is None:
      result = {
         'population': str(population),
         'chromosome': str(chromosome),
         'region_start_bp': long(start_bp),
         'region_end_bp': long(end_bp),
         'number_of_pairs': len(pairs),
         'pairs': [None] * len(pairs)
      }

      for i, pair in enumerate(pairs):
         result['pairs'][i] = { 
            'name1': pair.name1,
            'position1': pair.position1,
            'name2': pair.name2,
            'position2': pair.position2,
            'r': pair.r,
            'rsquare': pair.rsquare
         }
   else:
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

@app.route('/haplotypes/variant', methods = ['GET'])
def get_variant_haplotypes():
   population = request.args['population']
   chromosome = request.args['chromosome']
   name = request.args['name']

   haplotypes = PyHVCF.VariantHaplotypes() 

   hvcf.extract_haplotypes(str(population), str(chromosome), str(name), haplotypes)

   result = {
      'population': str(population),
      'number_of_samples': len(haplotypes),
      'chromosome': str(chromosome),
      'name': str(name),
      'samples': [None] * len(haplotypes)
   }

   for i, haplotype in enumerate(haplotypes):
      result['samples'][i] = {
         'sample': haplotype.sample,
         'allele1': haplotype.allele1,
         'allele2': haplotype.allele2
      }

   j = jsonify(result)
   return j

@app.route('/haplotypes/sample', methods = ['GET'])
def get_sample_haplotypes():
   chromosome = request.args['chromosome']
   start_bp = request.args['startbp']
   end_bp = request.args['endbp']
   name = request.args['name']

   haplotypes = PyHVCF.SampleHaplotypes()

   hvcf.extract_haplotypes(str(name), str(chromosome), long(start_bp), long(end_bp), haplotypes)

   result = {
      'sample': str(name),
      'chromosome': str(chromosome),
      'region_start_bp': long(start_bp),
      'region_end_bp': long(end_bp),
      'number_of_variants': len(haplotypes),
      'variants': [None] * len(haplotypes)
   }

   for i, haplotype in enumerate(haplotypes):
      result['variants'][i] = {
         'name': haplotype.name,
         'position': haplotype.position,
         'allele1': haplotype.allele1,
         'allele2': haplotype.allele2
      }

   j = jsonify(result)
   return j

if __name__ == '__main__':
   app.run(host='0.0.0.0', port=5000)
