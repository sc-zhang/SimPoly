import argparse
from sim_poly.utils.simulator import simulate
from sim_poly.utils.io import *
from sim_poly.utils.utils import time_print
import os
import time


def get_opts():
    group = argparse.ArgumentParser()
    group.add_argument(
        "-g", "--genome", help="input genome file, gz supported", required=True
    )
    group.add_argument(
        "-f", "--gff3", help="input gff3 file, gz supported"
    )
    group.add_argument("-c", "--config", help="input config file", required=True)
    group.add_argument("-o", "--out", help="output directory", required=True)

    return group.parse_args()


def main():
    st = time.time()
    opts = get_opts()
    ref_genome = opts.genome
    ref_gff3 = opts.gff3
    conf = opts.config
    out_dir = opts.out
    params = load_conf(conf)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    time_print("Loading reference")
    ref_fa_db = load_fasta(ref_genome)
    ref_gff3_db = load_gff3_keep_primary(ref_gff3) if ref_gff3 else {}

    time_print("Simulating")
    hap_genome_list = [{} for _ in range(params.ploidy)]
    hap_gff3_list = [{} for _ in range(params.ploidy)]
    hap_cds_list = [{} for _ in range(params.ploidy)]
    sim_cnt = 0
    for idx in range(params.ploidy):
        if params.structure[idx] == 0:
            hap_genome_list[idx], hap_gff3_list[idx], hap_cds_list[idx] = simulate(out_dir, idx + 1, ref_fa_db,
                                                                                   ref_gff3_db,
                                                                                   params.snp[idx],
                                                                                   params.insertion[idx],
                                                                                   params.deletion[idx],
                                                                                   params.insertion_length[idx],
                                                                                   params.deletion_length[idx])
            sim_cnt += 1

    while sim_cnt < params.ploidy:
        add_cnt = 0
        for idx in range(params.ploidy):
            if params.structure[idx] == 0:
                continue
            if hap_genome_list[params.structure[idx] - 1]:
                hap_genome_list[idx], hap_gff3_list[idx], hap_cds_list[idx] = simulate(out_dir, idx + 1,
                                                                                       hap_genome_list[
                                                                                           params.structure[idx] - 1],
                                                                                       hap_gff3_list[
                                                                                           params.structure[idx] - 1],
                                                                                       params.snp[idx],
                                                                                       params.insertion[idx],
                                                                                       params.deletion[idx],
                                                                                       params.insertion_length[idx],
                                                                                       params.deletion_length[idx])
                sim_cnt += 1
                add_cnt += 1
        if add_cnt == 0:
            raise "Structure error!"

    time_print("Saving")
    save_fasta(hap_genome_list, os.path.join(out_dir, "sim.fasta"))
    if ref_gff3:
        save_gff3(hap_gff3_list, os.path.join(out_dir, "sim.gff3"))
        save_fasta(hap_cds_list, os.path.join(out_dir, "sim.cds"))
    et = time.time()

    time_print("Total cost: %ds" % (et - st))
    time_print("Done")
