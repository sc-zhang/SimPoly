import gzip
import configparser


def load_fasta(in_fa):
    if in_fa.endswith(".gz"):
        fin = gzip.open(in_fa, "rt")
    else:
        fin = open(in_fa, "r")
    fa_db = {}
    gid = ""
    for line in fin:
        if line[0] == ">":
            gid = line.strip().split()[0][1:]
            fa_db[gid] = []
        else:
            fa_db[gid].append(line.strip().upper())

    for gid in fa_db:
        fa_db[gid] = "".join(fa_db[gid])
    fin.close()
    return fa_db


def save_fasta(hap_list, out_fa_file):
    with open(out_fa_file, 'w') as fout:
        for idx in range(len(hap_list)):
            for gid in hap_list[idx]:
                fout.write(">%s%s\n%s\n" % (gid, chr(ord("A") + idx), hap_list[idx][gid]))


def load_gff3_keep_primary(in_gff3):
    if in_gff3.endswith(".gz"):
        fin = gzip.open(in_gff3, "rt")
    else:
        fin = open(in_gff3, "r")
    gff_db = {}
    gene_id = ""
    mrna_cnt = 0
    for line in fin:
        if not line.strip() or line[0] == "#":
            continue
        data = line.strip().split()
        chrn = data[0]
        record_type = data[2]
        sp = int(data[3])
        ep = int(data[4])
        if record_type == 'gene':
            gene_id = data[8].split(";")[0].split("=")[-1]
            mrna_cnt = 0
        if chrn not in gff_db:
            gff_db[chrn] = []
        if record_type == 'mRNA':
            mrna_cnt += 1
        if mrna_cnt <= 1:
            gff_db[chrn].append([sp, ep, record_type, gene_id, line.strip().split()])
    fin.close()
    return gff_db


def save_gff3(hap_list, out_gff3_file):
    gene_id = ""
    with open(out_gff3_file, 'w') as fout:
        for idx in range(len(hap_list)):
            for gid in hap_list[idx]:
                for _, _, record_type, _, data in hap_list[idx][gid]:
                    if record_type == 'gene':
                        gene_id = data[8].split(";")[0].split("=")[-1]
                        new_gene_id = "%s%s" % (gene_id, chr(ord("A") + idx))
                        data[8] = "ID=%s;Name=%s" % (new_gene_id, new_gene_id)
                        idx_db = {}
                    else:
                        if record_type not in idx_db:
                            idx_db[record_type] = 1
                        new_rec_id = "%s%s-%s%d" % (gene_id, chr(ord("A") + idx), record_type, idx_db[record_type])
                        data[8] = "ID=%s;Name=%s;Parent=%s" % (new_rec_id, new_rec_id, new_gene_id)
                        idx_db[record_type] += 1
                    data[0] = "%s%s" % (data[0], chr(ord("A") + idx))
                    fout.write("%s\n" % ("\t".join(data)))


class _Parameters:
    ploidy = 2
    snp = [0.01, 0.01]
    insertion = [0.01, 0.01]
    deletion = [0.01, 0.01]
    insertion_length = [10, 10]
    deletion_length = [10, 10]
    structure = [0, 0, 0, 0]

    def __str__(self):
        return str(self.__dict__)


def load_conf(conf_file):
    cf = configparser.ConfigParser()
    cf.read(conf_file)

    params = _Parameters()
    params.ploidy = int(cf["parameters"]["ploidy"])
    params.snp = [float(_) for _ in cf["parameters"]["snp"].split(',')]
    params.insertion = [float(_) for _ in cf["parameters"]["insertion"].split(',')]
    params.deletion = [float(_) for _ in cf["parameters"]["deletion"].split(',')]
    params.insertion_length = [int(float(_)) for _ in cf["parameters"]["insertion_length"].split(',')]
    params.deletion_length = [int(float(_)) for _ in cf["parameters"]["deletion_length"].split(',')]
    params.structure = [int(_) for _ in cf["parameters"]["structure"].split(',')]
    return params
