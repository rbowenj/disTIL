# Script to get tumour sample name from a VCF
# Tumour sample identified from #CHROM header line, as the sample NOT ending with "G[0-9]" i.e. not germline

import sys
import os
from argparse import ArgumentParser
import re
import subprocess

def createArgParser():
    # Create argument parser
    description = "Get the tumour and normal sample names from a VCF."
    argparser = ArgumentParser(description = description)
    argparser.add_argument('vcf', type=str, 
                           help='Input VCF containing tumour and normal sample.')
    return argparser

def checkInputs(args):
    if not os.path.exists(args.vcf):
        raise Exception("VCF does not exist")

def execute_shell_command(command):
    po = subprocess.Popen(command,
                          shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)

    stdout, stderr = po.communicate()
    po.wait()

    returnCode = po.returncode

    if returnCode != 0:
        raise Exception("Command {} got return code {}.\nSTDOUT: {}\nSTDERR: {}".format(command, returnCode, stdout, stderr))

    return stdout

def get_sample_names(vcf_path):
    vcf = args.vcf.replace(" ", "\ ")
    sample_ids = execute_shell_command(f"bcftools query -l {vcf}").splitlines()

    tumour_sample = ''
    normal_sample = ''
    for sample in sample_ids:
        # Strip quotes and "b" at the start
        sample = str(sample).strip("b, '")

        # Look for the tumour sample (doesn't have G as second last letter)
        if not re.search(".*G.$", sample):
            tumour_sample = sample
        else:
            normal_sample = sample

    return tumour_sample, normal_sample

if __name__=='__main__':
    argparser = createArgParser()
    args = argparser.parse_args()
    
    tumour_sample, normal_sample = get_sample_names(args.vcf)

    # Write sample names to separate text files

    with open('tumour.txt', 'w+') as tf:
        tf.write(tumour_sample)

    with open('normal.txt', 'w+') as nf:
        nf.write(normal_sample)
