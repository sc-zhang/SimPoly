import argparse


def get_opts():
    group = argparse.ArgumentParser()
    group.add_argument(
        "-g", "--genome", help="input genome file, gz supported", required=True
    )
    group.add_argument(
        "-a", "--annotation", help="input annotation file, gz supported", required=True
    )
    group.add_argument(
        "-s",
        "--snp",
        type=float,
        help="snp ratio of whole genome, percentage, default: 0.01",
        default=0.01,
    )
    group.add_argument(
        "-i",
        "--insertion",
        type=float,
        help="insertion ratio of whole genome, percentage, default: 0.01",
        default=0.01,
    )
    group.add_argument(
        "--insert_length",
        type=int,
        help="max length of insertion, default: 10",
        default=10,
    )
    group.add_argument(
        "-d",
        "--deletion",
        type=float,
        help="deletion ratio of whole genome, percentage, default: 0.01",
        default=0.01,
    )
    group.add_argument(
        "--delete_length",
        type=int,
        help="max length of deletion, default: 10",
        default=10,
    )
    group.add_argument(
        "--random_length",
        action="store_true",
        help="use this argument for generate random length of indels",
        default=False,
    )
    group.add_argument("-o", "--out", help="output prefix", required=True)

    return group.parse_args()


def main():
    opts = get_opts()
    ref_genome = opts.genome
    gff3_file = opts.annotation
    snp_ratio = opts.snp
    ins_ratio = opts.insertion
    del_ratio = opts.deletion
    ins_len = opts.insert_length
    del_len = opts.delete_length
    random_length = opts.random_length
