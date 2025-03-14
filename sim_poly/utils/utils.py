import time
import numpy as np
from enum import Enum


class CDS_Status(Enum):
    MISSING_START = 1
    MISSING_STOP = 2
    NOT_MULTIPLE_3 = 3
    EARLY_STOP = 4
    VALID_CDS = 5


def time_print(info):
    print(
        "\033[32m%s\033[0m %s"
        % (time.strftime("[%H:%M:%S]", time.localtime(time.time())), info)
    )


def check_cds(in_seq):
    start_codon = {"ATG"}
    stop_codon = {"TAG", "TAA", "TGA"}
    if in_seq[:3] not in start_codon:
        return CDS_Status.MISSING_START, 0
    for _ in range(3, len(in_seq) - 3, 3):
        if in_seq[_: _ + 3] in stop_codon:
            return CDS_Status.EARLY_STOP, _ + 3
    if in_seq[-3:] not in stop_codon:
        return CDS_Status.MISSING_STOP, 0
    if len(in_seq) % 3 != 0:
        return CDS_Status.NOT_MULTIPLE_3, len(in_seq) % 3
    return CDS_Status.VALID_CDS, 0


def gen_seq(frag_len):
    base = "ATGC"
    return ''.join([base[np.random.randint(4)] for _ in range(frag_len)])


def gen_snp(ori_base):
    np.random.seed()
    base = list("ATGC")
    np.random.shuffle(base)
    for _ in base:
        if _ != ori_base:
            return _


def rev_seq(seq):
    base_db = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return ''.join([base_db[_] if _ in base_db else _ for _ in seq[::-1]])
