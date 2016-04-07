from flask import Flask, request, abort, jsonify
import PyHVCF
import time

app = Flask(__name__)

hvcf_file = 'test.h5'

hvcf = PyHVCF.HVCF()
hvcf.open(hvcf_file)

@app.route('/statistic/pair/LD/', methods = ['GET'])
def get_available_datasets():
   pass

@app.route('/LD/results', methods = ['GET'])
def get_ld():
   population = request.args['population']
   variant = request.args['variant']
   chromosome = request.args['chromosome']
   start = request.args['start']
   end = request.args['end']

   pairs = PyHVCF.PairsVector() 

   start_time = time.time()
   hvcf.compute_ld(str(chromosome), str(population), str(variant), long(start), long(end), pairs)
   elapsed_time = time.time() - start_time
   print 'HVCF request executed in ', elapsed_time, ' sec'

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
   print 'Response formatting executed in ', elapsed_time, ' sec' 

   return j

if __name__ == '__main__':
   app.run(host='0.0.0.0', port=5000)
