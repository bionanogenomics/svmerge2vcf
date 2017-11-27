#SVMerge-to-toVCF converter

#Version 1 supports conversion of insertions, deletions, and translocation breakpoints. It does not report breakpoint uncertainty.
#Version 2 adds support for duplications and inversions

description='''Stand-alone script to convert Bionano SVMerge file format to VCF.'''

vcfheader='''##fileformat=VCFv4.2
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=TRA,Description="Translocation breakpoint">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=DUP,Description="Duplication">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'''

import os
import argparse
import datetime
import math

#define accepted SV types and their corresponding converted types in VCF
class SV(object):
    def __init__(self, sv_type, vcf_type, svlen=None):
        self.sv_type = sv_type
        self.vcf_type = vcf_type
        self.svlen = svlen #default SV length

sv_types = [SV("insertion", "INS"), SV("deletion","DEL"),
            SV("translocation_intrachr","TRA",0),
            SV("translocation_interchr","TRA",0),
            SV("inversion", "INV"), SV("inversion_partial", "INV"), #SV("inversion_paired", "INV"), #paired are not in merged output
            SV("duplication", "DUP"), SV("duplication_split", "DUP"), SV("duplication_inverted", "INVDUP")]
accepted_svtypes = dict([(sv.sv_type,sv.vcf_type) for sv in sv_types])
default_svlen = dict([(sv.sv_type,sv.svlen) for sv in sv_types])

#check file; return bool (copied from utilities.py from assembly pipeline)
def checkFile(filepath, filesuff="", checkReadable=True) :
    try :
        valid = os.path.isfile(filepath)
        if filesuff != "" :
            valid = (valid and filepath.endswith(filesuff))
        if checkReadable :
            valid = (valid and os.access(filepath, os.R_OK))
    except :
        valid = False
    #print "checkFile:", filepath, valid #debug
    return valid

#main parser function
def svmerge_to_vcf(svmerge_path, sample, output_prefix, ref_accession, human_bool, vcfh) :

    doconf = True #put Confidence in QUAL field
    defaultconf = "." #missing value in vcf format: used if !doconf
    maxconf = "20" #if doconf, this is maximum reported: equivalent to ppv of 0.99

    si = ";" #separator for 'INFO' field
    colhead = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t" #this is complete if no genotype
    colform = "FORMAT\t%s\n" % sample #if genotype, add this

    outpath = svmerge_path.replace(".txt", ".vcf") #create default .vcf at svmerge_path (replacing suffix '.txt' with '.vcf')
    if output_prefix is not None:
        outpath = output_prefix+".vcf"

    today = datetime.date.today() #prepare fileDate for header

    fileDate = "##fileDate=%s" % today

    fout = open(outpath,"w")
    fout.write(fileDate+"\n")
    fout.write(vcfh+"\n")

    reference_header = "##reference=%s\n" % ref_accession

    if human_bool:
        reference_contigids = range(1,25)
        reference_contigids = [str(i) for i in reference_contigids] #make into string
        reference_contigids[22] = 'X'
        reference_contigids[23] = 'Y'

        contig_header = ["##contig=<ID=chr%s>\n" % elem for elem in reference_contigids]

    f1 = open(svmerge_path)
    dogeno = False #SVMerge has genotype info: add to vcf
    nent = 0 #number of processed entries
    conf = defaultconf #if !doconf, this is used
    invdata = {}
    for line in f1 :
        if line[0] == "#" :
            continue
        tokens = line.split("\t")
#         if len(tokens) < 64:
#             print "ERROR: line incomplete, terminating:\n%s\n" % line
#             break

        svtype = tokens[1]
        if svtype not in accepted_svtypes: #only process certain SV types
            continue
        vcftype = accepted_svtypes[svtype]

        if nent == 0 :
            dogeno = True
            fout.write(reference_header)
            fout.write("".join(contig_header))
            fout.write(colhead+colform)

        svmergeid = tokens[0] #SVIndex
        ref = tokens[2] #RefcontigID1
        ref1 = tokens[3] #RefcontigID2
        pos = int(float(tokens[4])) #RefStartPos
        end = int(float(tokens[5])) #RefEndPos
        
        if svtype!="translocation_intrachr" and svtype!="translocation_interchr":
            if end > 0 and end < pos : #make sure they are sorted (partials have end == -1, do not change this case)
                a = end
                end = pos
                pos = a

        if doconf :
            conf = tokens[6] #Confidence
            try :
                conf = float(conf)
                if conf==-1.0 :
                    conf = defaultconf
                elif conf >= 1.0 :
                    conf = maxconf
                else :
                    conf = str(abs(round(-10*math.log10(1-conf),2))) #abs to prevent '-0.0'
            except :
                conf = defaultconf

        if dogeno :
            zygosity = tokens[9]

            if zygosity== "homozygous" :
                gt = "1/1"
            elif zygosity== "heterozygous" :
                gt = "0/1"
            else: #unknown
                gt = "./."

        svlen = default_svlen[svtype]
        if svlen is None:
            svlen = round(float(tokens[8]))

        if human_bool:
            if ref=='23':
                ref='X'
            elif ref=='24':
                ref='Y'

            if ref1=='23':
                ref1='X'
            elif ref1=='24':
                ref1='Y'

        if svtype == "inversion" :
            invdata = {"refstart"   :pos,
                       "refstop"    :end}
            continue #need partial line (next line)
        
        elif svtype == "inversion_partial" :
            #for partial, end (RefEndPos) is always -1 (so ignore it)
            end = max(invdata["refstart"], invdata["refstop"], pos)
            pos = min(invdata["refstart"], invdata["refstop"], pos)
            svlen = round(end-pos)

        elif svtype.startswith("duplication") : #all duplications should be same
            svlen = round(end-pos)
        
        if svtype=="translocation_intrachr" or svtype=="translocation_interchr":
            fout.write(("chr%s\t%i\tSVMerge%s\tN\t<%s>\t%s\tPASS\tSVTYPE=%s"+si+"CHR2=chr%s"+si+"END=%i"+si+"SVLEN=%i") % (ref, pos, svmergeid,vcftype, conf, vcftype, ref1, end, svlen))
        else : #if svtype=="insertion" or svtype=="deletion": #should be same for all remaining types
            fout.write(("chr%s\t%i\tSVMerge%s\tN\t<%s>\t%s\tPASS\tSVTYPE=%s"+si+"END=%i"+si+"SVLEN=%i") % (ref, pos, svmergeid, vcftype, conf, vcftype, end, svlen))
        
        nent += 1
        
        if dogeno:
            fout.write("\tGT\t%s" % gt)
        fout.write("\n")

    #end loop on input svmerge
    f1.close()
    fout.close()

    print "Created vcf with %i entries: %s" % (nent, outpath)
#end svmerge_to_vcf

def getArgs() :
    parser = argparse.ArgumentParser(description=description)

    #required input
    parser.add_argument('-s', dest='svmerge_path', help='Path to SVMerge file to convert (required)', type=str)

    #optional input
    defsamp = "Sample1"
    parser.add_argument('-n', dest='sample', help='Sample ID name for genotype data (optional, default "%s")'%defsamp, type=str, default=defsamp)
    output_prefix = None
    parser.add_argument('-o', dest='output_prefix', help='Prefix for output vcf (optional, default to be same as input svmerge)', type=str, default=output_prefix)
    ref_accession = "GCA_000001405.1"
    parser.add_argument('-a', dest='ref_accession', help='RefSeq assembly accession version (optional, default "%s")'%ref_accession, type=str, default=ref_accession)
    human_bool = True
    parser.add_argument('-b', dest='human_bool', help='Whether sample is human (optional, default "%s")'%human_bool, type=str, default=human_bool)

    result = parser.parse_args()

    #unpack
    svmerge_path = result.svmerge_path #required

    sample = result.sample #optional
    output_prefix = result.output_prefix
    ref_accession = result.ref_accession
    human_bool = result.human_bool

    #checking files
    if not checkFile(svmerge_path,"_mergedSV.txt") :
        print "ERROR: SVMerge file does not exist, is not readable, or does not end with '_mergedSV.txt':", svmerge_path
        return None

    return svmerge_path, sample, output_prefix, ref_accession, human_bool

def run_svmerge_to_vcf():
    getargs = getArgs()
    if getargs != None :
        svmerge_to_vcf(getargs[0], getargs[1], getargs[2], getargs[3], getargs[4], vcfh=vcfheader)

if __name__ == "__main__" :
    run_svmerge_to_vcf()
