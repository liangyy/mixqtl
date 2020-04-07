import gzip
import argparse

def _write(f, comps):
    f.write("{}\n".format("\t".join(comps)).encode())

parser = argparse.ArgumentParser(description="Convert a gencode file")
parser.add_argument("-gencode")
parser.add_argument("-output")
args = parser.parse_args()

with gzip.open(args.gencode) as gencode:
    for i in range(0,5):
        gencode.readline()

    with gzip.open(args.output, "w") as o:
        _write(o, ["gene_id", "gene_name", "gene_type", "chromosome", "start", "end"])
        for line in gencode:
            comps = line.decode().strip().split()

            if comps[2] != "gene":
                continue

            gene_id=None
            gene_type = None
            gene_name=None
            key_values = {comps[i]:comps[i+1].replace('"', '').replace(';', '') for i in range(8, len(comps), 2)}

            l = [key_values["gene_id"], key_values["gene_name"], key_values["gene_type"], comps[0].replace("chr", ""), comps[3], comps[4]]
            _write(o, l)