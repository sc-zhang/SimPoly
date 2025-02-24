import time
import numpy as np


def time_print(info):
    print(
        "\033[32m%s\033[0m %s"
        % (time.strftime("[%H:%M:%S]", time.localtime(time.time())), info)
    )


def check_cds(in_seq):
    start_codon = {"ATG"}
    stop_codon = {"TAG", "TAA", "TGA"}
    if len(in_seq) % 3 != 0:
        return False
    if in_seq[:3] not in start_codon or in_seq[-3:] not in stop_codon:
        return False
    for _ in range(3, len(in_seq) - 3, 3):
        if in_seq[_: _ + 3] in stop_codon:
            return False
    return True


def gen_seq(frag_len):
    base = "ATGC"
    return ''.join([base[np.random.randint(4)] for _ in range(frag_len)])


def gen_snp(ori_base):
    base = list("ATGC")
    np.random.shuffle(base)
    for _ in base:
        if _ != ori_base:
            return _


def rev_seq(seq):
    base_db = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return ''.join([base_db[_] if _ in base_db else _ for _ in seq[::-1]])
