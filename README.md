# hw_4 Graph de Bruijn (vizualization)

## Description
This program draws graph de Bruijn with full/short labels.  

## Usage and arguments  
usage: graph_de_brujn_klass.py [-h] -i Str -o Str [-k Int] [-f] [-s]  

Graph de Brujin  

optional arguments:  
  -h, --help            show this help message and exit  
  -i Str, --infile Str  Input file  
  -o Str, --outfile Str Input file  
  -k Int, --kmersize Int
                        Kmer size  
  -f, --full            Output = full graph  
  -s, --short           Output = short graph  

## Example  
python3 graph_de_brujn_klass.py -i dataset.fasta -o graph_short.gv -s  
output - see supplementary  

