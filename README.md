# svmerge2vcf
# Tool for Transforming SVMerge to VCF Format

### Overview 
The SVMerge to VCF converter tool is a standalone python script that converts insertion, deletion, and translocation breakpoint calls in an SVMerge file to VCF v4.2 format. There is one required input to the script: the SVMerge file to convert. The QUAL score is calculated as -10 time the log base 10 of (1 minus the confidence)) where confidence is the SVMerge confidence score for the given call. The QUAL ceiling is set at 20. This output VCF file can be used for further downstream analysis using any tools that take VCF files as input.

###Usage

usage: svmerge_to_vcf.py [-h] [-s SVMERGE_PATH] [-n SAMPLE] [-o OUTPUT_PREFIX]
                         [-a REF_ACCESSION] [-b HUMAN_BOOL]

Stand-alone script to convert Bionano SVMerge file format to VCF.

optional arguments:

    -h, --help        show this help message and exit
  
    -s SVMERGE_PATH   Path to SVMerge file to convert (required)
  
    -n SAMPLE         Sample ID name for genotype data (optional, default
                      "Sample1")
                    
    -o OUTPUT_PREFIX  Prefix for output vcf (optional, default to be same as
                      input svmerge)
                    
    -a REF_ACCESSION  RefSeq assembly accession version (optional, default
                      "GCA_000001405.1")
                    
    -b HUMAN_BOOL     Whether sample is human (optional, default "True")

Note:  `python smap_to_vcf.py -h` to see usage on command line

### Requirements
This tool was designed to run with Python 2.7.  

### License
We offer this tool for open source use under the [MIT Software License](https://opensource.org/licenses/MIT). 
