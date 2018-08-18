import sys
import json
import os
import pandas as pd

script, run_id = sys.argv

run_annots = pd.DataFrame()

# First, get the run details
srr_file = run_id + ".ncbi.meta.json"
print("Reading file " + srr_file)
SRR = json.loads(open(srr_file).read())

# Second, get the experiment details. There should only
# be one of these in the directory
srx_id = SRR["EXPERIMENT_REF"]["@accession"]
srx_file = srx_id + ".ncbi.meta.json"
print("Reading file " + srx_file)
SRX = json.loads(open(srx_file).read())
 
# Second, get the sample details. There should only be one
# of these in the directory
sample_id = SRX["SAMPLE"]["@accession"]
sample_file = sample_id + ".ncbi.meta.json"
print("Reading file " + sample_file)
SRS = json.loads(open(sample_file).read())

# Now that we have all the metadata loaded create the non-nested 
# annotation dictionary
annots = {}

# Run info
annots['data:2091'] = SRR['@accession']

# The Experiment
annots['local:SRX_id'] = SRX['EXPERIMENT']['@accession']

# Biological Sample
annots['sep:00195'] = {}
annots['sep:00195']['data:2091'] = SRS['@accession']
annots['sep:00195']['schema:title'] = SRS['TITLE']
annots['sep:00195']['schema:name'] = SRS['@alias']
annots['sep:00195']['obi:organism'] =  {}
annots['sep:00195']['obi:organism']['rdfs:label'] =  SRS['SAMPLE_NAME']['SCIENTIFIC_NAME']
annots['sep:00195']['obi:organism']['NCIT:C43459'] =  SRS['SAMPLE_NAME']['SCIENTIFIC_NAME']
annots['sep:00195']['obi:organism']['data:1179'] =  SRS['SAMPLE_NAME']['TAXON_ID']

# Set defaults
annots['sep:00195']['NCIT:C25150'] = ''
annots['sep:00195']['NCIT:C16631'] = ''
annots['sep:00195']['NCIT:C12801'] = ''
annots['sep:00195']['NCIT:C43531'] = ''
annots['NCIT:C25206'] = ''
annots['EFO:0000721'] = ''
annots['EFO:0000727'] = ''

# Itearte through the sample attributes
attrs = SRS['SAMPLE_ATTRIBUTES']['SAMPLE_ATTRIBUTE']
for attr in attrs:

  # Add the cultivar
  if attr['TAG'] == 'cultivar':
    if attr['VALUE'] != 'missing':
      annots['sep:00195']['obi:organism']['local:infraspecific_type'] = 'cultivar'
      annots['sep:00195']['obi:organism']['TAXRANK:0000045'] = attr['VALUE']
    continue

  # Add the age 
  if attr['TAG'] == 'age':
    if attr['VALUE'] != 'missing':
      annots['sep:00195']['NCIT:C25150'] = attr['VALUE']
    continue
 
  # Add the genotype
  if attr['TAG'] == 'Genotype' or attr['TAG'] == 'genotype':
    if attr['VALUE'] != 'missing':
      annots['sep:00195']['NCIT:C16631'] = attr['VALUE']
    continue
  
  # Add the tissue
  if attr['TAG'] == 'tissue':
    if attr['VALUE'] != 'missing':
      annots['sep:00195']['NCIT:C12801'] = attr['VALUE']
    continue

  # Add the developmental stage
  if attr['TAG'] == 'dev_stage':
    if attr['VALUE'] != 'missing':
      annots['sep:00195']['NCIT:C43531'] = attr['VALUE']
    continue

  # Add the temperature
  if attr['TAG'] == 'temp':
    if attr['VALUE'] != 'missing':
      annots['NCIT:C25206'] = attr['VALUE']
    continue

  # Add the time
  if attr['TAG'] == 'time':
    if attr['VALUE'] != 'missing':
      annots['EFO:0000721'] = attr['VALUE']
    continue

  # Add the treatment
  if attr['TAG'] == 'treatment':
    if attr['VALUE'] != 'missing':
      annots['EFO:0000727'] = attr['VALUE']
    continue

  print("Unhandled sample attribute: '" + attr['TAG'] + "', value: " + attr['VALUE'])

# Save the heirarchical JSON metadata from GEMmaker
jsonfilename = annots['local:SRX_id'] + '.GEMmaker.meta.json'
print("Writing file: " + jsonfilename)
with open(jsonfilename, 'w') as jsonfile:
  json.dump(annots, jsonfile)

# Save a flattened CSV file
tabfilename = annots['local:SRX_id'] + '.GEMmaker.meta.tab'
print("Writing file: " + tabfilename)
with open(tabfilename, 'w') as tabfile:
  tabfile.write(
    annots["data:2091"]
    + "\t" + annots["local:SRX_id"]
    + "\t" + annots["sep:00195"]["data:2091"]
    + "\t" + annots["sep:00195"]["schema:title"]
    + "\t" + annots["sep:00195"]["schema:name"]
    + "\t" + annots["sep:00195"]["obi:organism"]["NCIT:C43459"]
    + "\t" + annots["sep:00195"]["NCIT:C25150"]
    + "\t" + annots["sep:00195"]["NCIT:C43531"]
    + "\t" + annots["sep:00195"]["NCIT:C16631"]
    + "\t" + annots["sep:00195"]["NCIT:C12801"]
    + "\t" + annots["NCIT:C25206"]
    + "\t" + annots["EFO:0000721"]
    + "\t" + annots['EFO:0000727']
    + "\n"
  )

