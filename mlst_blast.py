#!/usr/bin/env python

"""calculates MLST types from assemblies using BLAST.
If the gene is truncated, it will report a "T" and if
the gene has no blast hit, it will report a "X".
Your reference allele names must all end in "fasta"
and must contain a "_" between gene name and number.
The only external dependency is blast+ - tested version
is 2.2.31 - and igs multithreading tools"""

import optparse
import glob
import subprocess
import os
import re
import sys
try:
    from igs.utils import functional as func
    from igs.utils import logging
    from igs.threading import functional as p_func
except:
    print "igs utilities need to be in your directory"
    sys.exit()

def get_gene_order(reference):
    firstLine = open(reference).readline()
    fields = firstLine.split()
    ordered = fields[1:]
    return ordered

def get_seq_name(in_fasta):
    """used for renaming the sequences"""
    return os.path.basename(in_fasta)

def run_blast(directory,genome):
    name = get_seq_name(genome)
    """Directory is the reference alleles"""
    for infile in glob.glob(os.path.join(directory, '*.fasta')):
        """infile is now the individual allele"""
        allele_name = get_seq_name(infile)
        allele_redux_name = allele_name.replace('.fasta','')
        """This will do the blast analysis"""
        outfile = open("%s_%s_blast.out" % (allele_redux_name,name), "w")
        subprocess.check_call("makeblastdb -in %s -dbtype nucl > /dev/null 2>&1" % genome, shell=True)
        cmd = ["blastn",
               "-task", "blastn",
               "-query", infile,
               "-db", genome,
               "-out", "%s_blast_%s.out" % (allele_redux_name,name),
               "-evalue", "0.0001",
               "-dust", "no",
               "-outfmt", "6"]
        subprocess.check_call(cmd)
        new_infile = open("%s_blast_%s.out" % (allele_redux_name,name))
        new_outfile = open("%s_blast_%s_fixed.out" % (allele_redux_name, name), "w")
        for line in new_infile:
            fields=line.split()
            print >> new_outfile, "\t".join(fields)
        new_infile.close()
        new_outfile.close()

def process_results(genome_names,ordered,trunc,new,reference):
    outfile = open("profile.txt", "w")
    pad_order = list(ordered)
    pad_order.insert(0,"genome")
    pad_order.insert(8,"ST")
    print >> outfile, "\t".join(pad_order)
    for genome in genome_names:
        alleles = [ ]
        for marker in ordered:
            marker_redux = marker.replace(".fasta","")
            new_list = [line.strip().split("\t") for line in open("%s_blast_%s_fixed.out" % (marker_redux, genome), "rU")]
            new_list=sorted(new_list, key=lambda v: float(v[-1]), reverse=True)
            try:
                fields=new_list[0]
                if int(fields[3])/int(fields[7])>=trunc and float(fields[2])>=float(new):
                    alleles.append(fields[0].split("_")[1])
                elif int(fields[3])/int(fields[7])>=trunc and float(fields[2])<float(new):
                    alleles.append("N")
                else:
                    alleles.append("T")
            except:
                alleles.append("X")
        alleles.insert(0, genome)
        in_ref = open(reference, "rU")
        for line in in_ref:
            ref_fields=line.split()
            sum = [ ]
            for i in range(1,8):
                if ref_fields[i]==alleles[i]:
                    sum.append(i)
            if int(len(sum))==7:
                alleles.insert(8, ref_fields[0])
        print >> outfile, "\t".join(alleles)
    outfile.close()

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print '%s cannot be opened' % option
        sys.exit()

def test_dir(option, opt_str, value, parser):
    if os.path.exists(value):
        setattr(parser.values, option.dest, value)
    else:
        print "%s cannot be found" % option
        sys.exit()

def main(directory,reference,fastas,trunc,new,processors):
    curr_dir = os.getcwd()
    ordered = get_gene_order(reference)
    file_dir = glob.glob(os.path.join(fastas,"*fasta"))
    files_and_temp_names = [(str(idx), os.path.join(curr_dir, f))
                            for idx, f in enumerate(file_dir)]
    def _perform_workflow(data):
        tn, f = data
        run_blast(directory,f)
    results = set(p_func.pmap(_perform_workflow,
                              files_and_temp_names,
                              num_workers=processors))
    print "blast finished!"
    genome_names = []
    for infile in glob.glob(os.path.join(fastas, '*.fasta')):
        name = get_seq_name(infile)
        genome_names.append(name)
    process_results(genome_names,ordered,trunc,new,reference)
    subprocess.check_call("rm *.out", shell=True)

if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-d", "--seq_dir", dest="directory",
                      help="/path/to/reference_seq_directory [REQUIRED]",
                      action="callback", callback=test_dir, type="string")
    parser.add_option("-r", "--ref file", dest="reference",
                      help="/path/to/ST_reference_file [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-f", "--fasta dir", dest="fastas",
                      help="/path/to/your_genomes [REQUIRED]",
                      action="callback", callback=test_dir, type="string")
    parser.add_option("-t", "--truncated value", dest="trunc",
                      help="gene is considered truncated if below this value [default = 0.95]",
                      default="0.95", type="float")
    parser.add_option("-n", "--new", dest="new",
                      help="under this value is considered to be a new allele",
                      default="1.0", type="float")
    parser.add_option("-p", "--processors", dest="processors",
                      help="Number of processors to use, defaults to 2",
                      default="2", type="int", action="store")
    options, args = parser.parse_args()

    mandatories = ["directory", "reference", "fastas"]
    for m in mandatories:
        if not getattr(options, m, None):
            print "\nMust provide %s.\n" %m
            parser.print_help()
            exit(-1)
    main(options.directory,options.reference,options.fastas,options.trunc,options.new,options.processors)
