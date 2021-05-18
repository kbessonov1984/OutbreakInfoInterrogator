# OutbreakInfoInterrogator for SARS-COV-2
Get report for each query mutation from `outbreak.info` resource and find out in which lineages this mutation is present and at what abudance. The output is a tab-delimited file for each query.

The script queries `outbreak.info` and stores cached data. The maximum age of cached file is 24hrs after which a fresh data snapshot is downloaded

## Requirements
- python 3
- requests >= 2.0

## Install
Minimum requirements are needed to run this script. To quickly install all necessary libraries run `pip install -r requirements.txt`


## Usage
Simply run as any Python script by specifying mutation in the following formats: a) for substitution <Gene>:<Substituion> (`s:n501y`); b) deletions <Gene>:<Deletion>(`s:del69/70` or `orf1a:del17`)

```
usage: outbreakinfointerrogator.py uses https://api.outbreak.info/ 
to get information on query mutations supported by GISAID database data

optional arguments:
  -h, --help            show this help message and exit
  -m MUTATION_NAME, --mutation_name MUTATION_NAME
                        Mutation name to query (e.g. s:del69/70, s:d614g,orf1a:del17)
```


## Output

|Mutation|#PangoLineages|PangoNames(prevalence in lineage)|
|--------|--------------|---------------------------------|
|s:d614g|1196	|A.18(1.0),A.19(1.0),B.1.1.3(1.0),B.1.1.4(1.0)|