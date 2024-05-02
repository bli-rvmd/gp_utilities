# script to generate random candidate barcodes that satisfy criteria below
# author: Biao Li
# criteria:
    ## 50% or higher GC
    ## no > 3 poly
    ## hamming distance > 10 between any pair of barcodes
# usage:
    ## output candidate barcodes in FASTA format

import random
import numpy as np
import time
from datetime import timedelta
import logging

def gen_barcodes(length = 40,
                 num = 50,
                 out_file = "barcodes.fasta",
                 min_hamming_dist = 10,
                 min_gc_content = 0.5,
                 max_poly = 3,
                 debug = False):
    """
    :param length: desired length of barcodes, default at 40
    :param num: desired number of barcodes, default at 1000
    :param out_file: output file name suffix in .fasta format
    :param min_hamming_dist: minimum hamming distance between barcodes, default at 10
    :param min_gc_content: minimum gc content of barcodes, default at 50%
    :param max_poly: maximum polygon length of any tandem repeats in a barcode, default at 3
    :param debug: if True turn on debugging mode, default at False
    :return: candidate barcodes in fasta format
    """

    res = []
    nucleotides = ['A', 'C', 'G', 'T']

    if debug:
        logging.basicConfig(filename = 'info.log', level=logging.DEBUG)

    while len(res) < num:
        ## initiate barcode string
        bc = ""
        gc_count = 0

        while len(bc) < length:
            nuc = random.choice(nucleotides)
            bc += nuc
            ## check for tandem repeats of poly
            if len(bc) >= max_poly:
                if len(set(bc[-max_poly:])) == 1:
                    # reset and re-initiate
                    bc = ""
                    gc_count = 0
                    continue
            ## track generated GC
            if nuc in ['G', 'C']:
                gc_count += 1

        ## check if GC content satisfies criterion
        if float(gc_count) / length < min_gc_content:
            continue

        else:
            ## edge case for first generated barcode
            if len(res) == 0:
                res.append(bc)
                continue

            ## check if generated barcode having hamming distance > 10 paired with any existing one
            bc_rep = [bc] * len(res)
            hamming_dists = [(np.array(list(x)) != np.array(list(y))).sum() for x, y in zip(res, bc_rep)]

            if debug:
                logging.info("{}\n{}".format(bc, hamming_dists))

            if all([x >= min_hamming_dist for x in hamming_dists]):
                ## found a candidate
                res.append(bc)

            else:
                continue

        # track progress
        if len(res) % 5 == 0:
            print("generated {} candidate barcodes".format(len(res)))

    # output barcode candidates in fasta format
    with open(out_file, "w") as f:
        for bc_idx, bc in enumerate(res):
            f.write(">BC{}_{}\n".format(length, bc_idx + 1))
            f.write("{}\n".format(bc))

    return (res)


if __name__ == '__main__':

    debug = False

    length = 40
    num = 10000
    out_file = "barcodes.fasta"
    min_hamming_dist = 10
    min_gc_content = 0.5
    max_poly = 3

    # track total run time
    start = time.time()
    barcodes = gen_barcodes(length = length,
                       num = num,
                       out_file = out_file,
                       min_hamming_dist = min_hamming_dist,
                       min_gc_content = min_gc_content,
                       max_poly = max_poly,
                       debug = debug)
    end = time.time()
    elapsed = end - start
    print("total run time: {}".format(timedelta(seconds=elapsed)))


