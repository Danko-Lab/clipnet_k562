import argparse
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input", type=str, help="a csv file. Will auto-detect and uncompress .gz files"
    )
    args = parser.parse_args()

    fn = os.path.split(args.input)[-1]
    pref = fn.replace(".gz", "")
    pref = pref.replace(".csv", "")
    pref = pref.replace(".pl", "")
    pref = pref.replace(".mn", "")

    if fn.split(".")[-1] == "gz":
        cmd = "zcat %s | awk -F',' '$1=\"%s_\"$1' | sed -e 's/ /,/g'" % (
            args.input,
            pref,
        )
    else:
        cmd = "awk -F',' '$1=\"%s_\"$1' %s | sed -e 's/ /,/g'" % (pref, args.input)
    os.system(cmd)
