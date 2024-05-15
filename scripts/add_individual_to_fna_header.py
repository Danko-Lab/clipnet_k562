import argparse
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input",
        type=str,
        help="a fasta file. Will auto-detect and uncompress .gz files",
    )
    args = parser.parse_args()

    fn = os.path.split(args.input)[-1]
    id = fn.replace(".gz", "")
    id = fn.replace(".fna", "")
    id = id.replace(".fa", "")
    id = id.replace(".fasta", "")

    if fn.split(".")[-1] == "gz":
        cmd = 'zcat %s | awk \'/>/{sub(">","&""%s""_");sub(/\.fna/,x)}1\'' % (
            args.input,
            id,
        )
    else:
        cmd = 'awk \'/>/{sub(">","&""%s""_");sub(/\.fna/,x)}1\' %s' % (id, args.input)
    os.system(cmd)
