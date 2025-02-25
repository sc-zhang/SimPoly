import copy
import os
from sim_poly.utils.utils import *


def simulate(out_dir, hap_idx, ref_fa_db, ref_gff3_db, snp_ratio, ins_ratio, del_ratio, max_ins_len, max_del_len):
    np.random.seed()
    new_fa_db = {}
    new_gff3_db = {}
    new_cds_db = {}

    tmp_gff3_db = {}
    snp_db = {}
    ins_db = {}
    del_db = {}
    time_print("Simulating for haplotype %d" % hap_idx)
    for gid in ref_fa_db:
        time_print("\tSimulating %s" % gid)
        new_fa_db[gid] = list(ref_fa_db[gid])
        snp_db[gid] = []
        ins_db[gid] = []
        del_db[gid] = []
        seq_len = len(ref_fa_db[gid])
        time_print("\tReference length: %d" % seq_len)
        pos_list = [1 for _ in range(seq_len)]
        snp_cnt = int(seq_len * snp_ratio)
        time_print("\tSimulating SNP: %d" % snp_cnt)
        while len(snp_db[gid]) < snp_cnt:
            pos = np.random.randint(seq_len)
            if pos_list[pos] == 1:
                snp_base = gen_snp(ref_fa_db[gid][pos])
                snp_db[gid].append([pos, snp_base])
                pos_list[pos] = 0
                new_fa_db[gid][pos] = snp_base

        total_del_len = int(seq_len * del_ratio)
        time_print("\tSimulating deletion length: %d" % total_del_len)
        cur_del_len = 0
        while cur_del_len < total_del_len:
            del_len = np.random.randint(1, max_del_len + 1)
            pos = np.random.randint(seq_len)
            if pos_list[pos] == 1:
                pos_list[pos] = 0
                right_pos = pos
                for _ in range(1, del_len + 1):
                    if pos + _ >= seq_len or pos_list[pos + _] == 0:
                        break
                    else:
                        right_pos = pos + _
                        pos_list[right_pos] = 0
                for _ in range(pos, right_pos + 1):
                    new_fa_db[gid][_] = ''
                del_len = right_pos - pos + 1
                del_db[gid].append([pos, del_len])
                cur_del_len += del_len

        total_ins_len = int(seq_len * ins_ratio)
        time_print("\tSimulating insertion length: %d" % total_ins_len)
        cur_ins_len = 0
        while cur_ins_len < total_ins_len:
            ins_len = np.random.randint(1, max_ins_len + 1)
            pos = np.random.randint(seq_len)
            if pos_list[pos] == 1:
                pos_list[pos] = 0
                ins_seq = gen_seq(ins_len)
                ins_db[gid].append([pos, ins_len, ins_seq])
                new_fa_db[gid][pos] = ''.join([new_fa_db[gid][pos], ins_seq])
            cur_ins_len += ins_len

    for gid in new_fa_db:
        new_fa_db[gid] = ''.join(new_fa_db[gid])

    for gid in ref_fa_db:
        snp_db[gid] = sorted(snp_db[gid])
        ins_db[gid] = sorted(ins_db[gid])
        del_db[gid] = sorted(del_db[gid])

    time_print("\tAdjusting gff3")
    for gid in ref_gff3_db:
        tmp_gff3_db[gid] = []
        ins_idx = 0
        ins_offset = 0
        del_idx = 0
        del_offset = 0
        for sp, ep, record_type, gene_id, data in ref_gff3_db[gid]:
            while ins_idx < len(ins_db[gid]) and sp >= ins_db[gid][ins_idx][0] + 1:
                ins_offset += ins_db[gid][ins_idx][1]
                ins_idx += 1
                if ins_idx >= len(ins_db[gid]):
                    break

            while del_idx < len(del_db[gid]) and sp >= del_db[gid][del_idx][0] + 1:
                del_offset += del_db[gid][del_idx][1]
                del_idx += 1
                if del_idx >= len(del_db[gid]):
                    break
            sp = sp + ins_offset - del_offset

            while ins_idx < len(ins_db[gid]) and ep >= ins_db[gid][ins_idx][0] + 1:
                ins_offset += ins_db[gid][ins_idx][1]
                ins_idx += 1
                if ins_idx >= len(ins_db[gid]):
                    break
            while del_idx < len(del_db[gid]) and ep >= del_db[gid][del_idx][0] + 1:
                del_offset += del_db[gid][del_idx][1]
                del_idx += 1
                if del_idx >= len(del_db[gid]):
                    break
            ep = ep + ins_offset - del_offset
            tmp_gff3_db[gid].append([gene_id, sp, ep, record_type, copy.deepcopy(data)])

    time_print("\tFiltering gff3 with CDS checking")
    retain_genes = set()
    gene_cnt = 0
    for gid in tmp_gff3_db:
        cds_db = {}
        for gene_id, sp, ep, record_type, data in tmp_gff3_db[gid]:
            if gene_id not in cds_db:
                cds_db[gene_id] = []
                gene_cnt += 1
            if record_type == 'CDS':
                seq = new_fa_db[gid][sp - 1: ep]
                if data[6] == '-':
                    seq = rev_seq(seq)
                cds_db[gene_id].append([sp, ep, data[6], seq])

        for gene_id in cds_db:
            seq = ""
            if cds_db[gene_id]:
                if cds_db[gene_id][0][2] == '+':
                    seq_list = [_[3] for _ in sorted(cds_db[gene_id])]
                else:
                    seq_list = [_[3] for _ in sorted(cds_db[gene_id], reverse=True)]
                seq = ''.join(seq_list)
            if check_cds(seq):
                retain_genes.add(gene_id)
                new_cds_db[gene_id] = seq

    time_print("\tTotal: %d, Valid: %d" % (gene_cnt, len(retain_genes)))

    for gid in tmp_gff3_db:
        new_gff3_db[gid] = []
        for gene_id, sp, ep, record_type, data in tmp_gff3_db[gid]:
            if gene_id not in retain_genes:
                continue
            data[3] = str(sp)
            data[4] = str(ep)
            new_gff3_db[gid].append([sp, ep, record_type, gene_id, data])

    time_print("\tWriting variants")
    with open(os.path.join(out_dir, "snp%d.txt" % hap_idx), 'w') as fout:
        fout.write("#CHR\tPOS\tREF\tALT\n")
        for gid in sorted(snp_db):
            for pos, snp_base in snp_db[gid]:
                fout.write("%s\t%d\t%s\t%s\n" % (gid, pos, ref_fa_db[gid][pos], snp_base))

    with open(os.path.join(out_dir, "ins%d.txt" % hap_idx), 'w') as fout:
        fout.write("#CHR\tPOS\tLEN\tSEQ\n")
        for gid in sorted(ins_db):
            for pos, ins_len, ins_seq in ins_db[gid]:
                fout.write("%s\t%d\t%d\t%s\n" % (gid, pos, ins_len, ins_seq))
    with open(os.path.join(out_dir, "del%d.txt" % hap_idx), 'w') as fout:
        fout.write("#CHR\tPOS\tLEN\tSEQ\n")
        for gid in sorted(del_db):
            for pos, del_len in del_db[gid]:
                fout.write("%s\t%d\t%d\t%s\n" % (gid, pos, del_len, ref_fa_db[gid][pos: pos + del_len]))

    return new_fa_db, new_gff3_db, new_cds_db
