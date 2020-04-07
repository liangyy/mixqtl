from cyvcf2 import VCF
import argparse
import os
import gzip

def _write(f, comps):
    f.write("{}\n".format("\t".join(comps)).encode())

########################################################################################################################
#Functions to de-multiplex variants into one chromosome per file
current_chromosome = None
def activate_output(output, chromosome, samples):
    if chromosome is None or samples is None:
        return None, None, None
    print(f"Activating {chromosome}")
    variant_annotation_ = f'{output}_{chromosome}_variant_annotation.txt.gz'
    variant_annotation = gzip.open(variant_annotation_, "w")
    _write(variant_annotation, ["variant", "chromosome", "position", "non_effect_allele", "effect_allele"])

    hap_1_ = f'{output}_{chromosome}_hap1.txt.gz'
    hap_1 = gzip.open(hap_1_, "w")
    _write(hap_1, ["variant"] + samples)

    hap_2_ = f'{output}_{chromosome}_hap2.txt.gz'
    hap_2 = gzip.open(hap_2_, "w")
    _write(hap_2, ["variant"] + samples)
    return variant_annotation, hap_1, hap_2

def deactivate_output(variant_annotation, hap_1, hap_2):
    if variant_annotation is not None: variant_annotation.close()
    if hap_1 is not None: hap_1.close()
    if hap_2 is not None: hap_2.close()

########################################################################################################################
# Process each vcf variant
def convert_allele(a):
    if a == 0:
        return "0"
    elif a == 1:
        return "1"
    else:
        return "NA"

def convert_sample(s):
    if len(s) == 3 and s[2] == True:
                                                                                                                                e = (convert_allele(s[0]), convert_allele(s[1]))
    else:
        e = ("NA", "NA")

    return e


########################################################################################################################
# Run
parser = argparse.ArgumentParser(description="Convert a vcf to mixqtl formats")
parser.add_argument("-vcf")
parser.add_argument("-output_prefix")
args = parser.parse_args()

folder = os.path.split(args.output_prefix)[0]
if not os.path.exists(folder):
    os.makedirs(folder)

variant_annotation, hap_1, hap_2 = activate_output(args.output_prefix, current_chromosome, None)
vcf = VCF(args.vcf)

print("processing")
for i,variant in enumerate(vcf):

    chromosome = variant.CHROM
    if chromosome != current_chromosome:
        deactivate_output(variant_annotation, hap_1, hap_2)
        current_chromosome =  chromosome
        variant_annotation, hap_1, hap_2 = activate_output(args.output_prefix, current_chromosome, vcf.samples)

    _write(variant_annotation, [variant.ID, variant.CHROM, str(variant.POS), variant.REF, variant.ALT[0]])
    samples = [(variant.ID, variant.ID)] + [convert_sample(s) for s in variant.genotypes]
    hap_1_line, hap_2_line = zip(*samples)

    _write(hap_1, hap_1_line)
    _write(hap_2, hap_2_line)

print("Closing")
deactivate_output(variant_annotation, hap_1, hap_2)

print("done")

