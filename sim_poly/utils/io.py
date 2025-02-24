import gzip


def load_fasta(in_fa):
    change_db = {}
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
            fa_db[gid].append(line.strip())

    for gid in fa_db:
        fa_db[gid] = "".join(fa_db[gid])
        change_db[gid] = [1 for _ in range(len(fa_db[gid]))]
    fin.close()
    return fa_db, change_db


def load_gff3(in_gff3):
    if in_gff3.endswith(".gz"):
        fin = gzip.open(in_gff3, "rt")
    else:
        fin = open(in_gff3, "r")
    gff_db = {}
    for line in fin:
        if not line.strip() or line[0] == "#":
            continue
        data = line.strip().split()
        chrn = data[0]
        record_type = data[2]
        sp = int(data[3])
        ep = int(data[4])
        if chrn not in gff_db:
            gff_db[chrn] = []
        else:
            gff_db[chrn].append([sp, ep, record_type, line.strip()])
    fin.close()
    return gff_db
