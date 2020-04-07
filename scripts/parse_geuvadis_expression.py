import gzip
import argparse

def _write(f, comps):
    f.write("{}\n".format("\t".join(comps)).encode())

parser = argparse.ArgumentParser(description="Convert a geuvados expression file")
parser.add_argument("-expression")
parser.add_argument("-output")
args = parser.parse_args()

with gzip.open(args.expression) as e:
    header = e.readline().decode().strip().split()
    with gzip.open(args.output, "w") as o:
        _write(o, ["gene_list"]+header[4:])
        for line in e:
            comps = line.decode().strip().split()
            _write(o, comps[3:])
