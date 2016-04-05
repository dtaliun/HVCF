from flask import Flask, request, abort
import re

app = Flask(__name__)

operator_pattern = r"(eq|le|ge)"
variable_name_pattern = r"([a-zA-Z][a-zA-Z0-9_-]*)"
string_value_pattern = r"(?:'[^']*')"
uint_value_pattern = r"(?:0|[1-9]\d*)"
value_pattern = "".join(["(", string_value_pattern, "|", uint_value_pattern, ")"])
and_pattern = r"\s+and\s+"
operation_pattern = "".join([variable_name_pattern, "\s+", operator_pattern, "\s+", value_pattern]) 

string_value_regex = re.compile(string_value_pattern)
uint_value_regex = re.compile(uint_value_pattern)
operation_regex = re.compile(operation_pattern)
and_regex = re.compile(and_pattern)

def parse_filter_argument(filter):
   start = 0
   end = len(filter)

   filters = []

   while True:
      match = operation_regex.match(filter, start)
      if not match:
         raise Exception()

      filters.append(match.groups())

      start = match.end()
      if start >= end:
         break

      match = and_regex.match(filter, start) 
      if not match:
         raise Exception()

      start = match.end()

   return filters

@app.route('/statistic/pair/LD/', methods = ['GET'])
def get_available_datasets():
   pass

@app.route('/statistic/pair/LD/results', methods = ['GET'])
def get_ld():
   try: 
      for arg in request.args:
         if arg == 'filter':
            filters = parse_filter_argument(request.args[arg])
   except Exception:
      abort(501, "Incorrect arguments provided.")

   reference = None
   variant1 = None
   chromosome2 = None
   position2_start = None
   position2_end = None

   for filter in filters:
      variable, operator, value = filter
      if variable == "reference":
         if operator != "eq":
            abort(501, "Unsupported operator with 'reference' variable in filter parameter")
         if not uint_value_regex.match(value):
            abort(501, "Unsupported data type for 'reference' variable in filter parameter") 
         reference = value
      elif variable == "chromosome1":
         if operator != "eq":
            abort(501, "Unsupported operator with 'chromosome1' variable in filter parameter")     
         if not string_value_regex.match(value):
            abort(501, "Unsupported data type for 'chromosome1' variable in filter parameter")
         chromosome1 = value
      elif variable == "position1":
         if not uint_value_regex.match(value):
            abort(501, "Unsupported data type for 'position1' variable in filter parameter")
         if operator == "ge":
            position1_start = long(value)
         elif operator == "le":
            position1_end = long(value)
         else: 
            abort(501, "Unsupported operator with 'position1' variable in filter parameter")
      elif variable == "variant1":
         if operator != "eq":
            abort(501, "Unsupported operator with 'variant1' variable in filter parameter")
         if not string_value_regex.match(value):
            abort(501, "Unsupported data type for 'variant1' variable in filter parameter") 
         variant1 = value
      elif variable == "chromosome2":
         if operator != "eq":
            abort(501, "Unsupported operator with 'chromosome2' variable in filter parameter")
         if not string_value_regex.match(value):
            abort(501, "Unsupported data type for 'chromosome2' variable in filter parameter") 
         chromosome2 = value
      elif variable == "position2":
         if not uint_value_regex.match(value):
            abort(501, "Unsupported data type for 'position2' variable in filter parameter") 
         if operator == "ge":
            position2_start = long(value)
         elif operator == "le":
            position2_end = long(value)
         else: 
            abort(501, "Unsupported operator with 'position2' variable in filter parameter")
      elif variable == "variant2":
         if operator != "eq":
            abort(501, "Unsupported operator with 'variant2' variable in filter parameter")
         if not string_value_regex.match(value):
            abort(501, "Unsupported data type for 'variant2' variable in filter parameter") 
         variant2 = value
      else:
         abort(501, "Unsupported variable in filter parameter")

   if (reference is None) or (variant1 is None) or (chromosome2 is None) or (position2_start is None) or (position2_end is None):
      abort(501, "This filter configuration is not supported yet") 

   print reference, variant1, chromosome2, position2_start, position2_end

   return 'oops'

if __name__ == '__main__':
   app.run(host='0.0.0.0', port=5000)
