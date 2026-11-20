#!/usr/bin/env python3
"""Generate the public mock demo data (deterministic, seed 42).

WT VA RNA I / II come from the public GenBank record OP218818.1 (HAdV-5):
  VA I  region 10566-10843 (278 nt), mature RNA = region positions 52-214
  VA II region 10828-11083 (256 nt), mature RNA = region positions 49-208
All "MOCK*" variants are INVENTED for demonstration; they are not real
adenovirus mutants.

Outputs (in demo/):
  VA_RNA_I_mock.fa, VA_RNA_II_mock.fa   annotated FASTA ({XY} = replacement, '-' = deletion)
  fastq/DEMO_<genotype>_<rep>_S1_L001_R2_001.fastq.gz   small mock 10x-style R2 files

Run from the repo root:  python demo/make_demo_data.py
"""

import gzip
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
TSO = "AAGCAGTGGTATCAACGCAGAGTACATGGG"
READ_LEN = 90
COMP = str.maketrans("ACGT", "TGCA")
rc = lambda s: s.translate(COMP)[::-1]  # noqa: E731

genome = "".join(l.strip() for l in (HERE / "OP218818.1.fasta").read_text().splitlines()[1:]).upper()

# name -> (1-based start, 1-based end, mature start0, mature end0-exclusive) within region
REGIONS = {"VAI": (10566, 10843, 51, 214), "VAII": (10828, 11083, 48, 208)}
WT = {r: genome[a - 1 : b] for r, (a, b, _, _) in REGIONS.items()}

# INVENTED mutations: (pos0 in WT region, length, kind) ; kind "sub" or "del"
MOCK = {
    "VAI": {"MOCKA": [(100, 4, "sub")], "MOCKAB": [(100, 4, "sub")], "MOCKC": [(160, 8, "del")]},
    "VAII": {"MOCKB": [(120, 6, "del")], "MOCKAB": [(120, 6, "del")], "MOCKC": [(90, 3, "sub")]},
}
GENOTYPES = ["WT", "MOCKA", "MOCKB", "MOCKAB", "MOCKC"]


def annotate(rna: str, muts) -> tuple[str, str]:
    """Return (annotated string, plain variant sequence)."""
    wt = WT[rna]
    _, _, ms, me = REGIONS[rna]
    case = lambda i, b: b.upper() if ms <= i < me else b.lower()  # noqa: E731
    ann, plain, i = [], [], 0
    for pos, ln, kind in sorted(muts):
        ann += [case(j, wt[j]) for j in range(i, pos)]
        plain.append(wt[i:pos])
        if kind == "sub":
            new = wt[pos : pos + ln].translate(COMP)  # complement: guaranteed different bases
            ann.append("{" + new + "}")
            plain.append(new)
        else:
            ann.append("-" * ln)
        i = pos + ln
    ann += [case(j, wt[j]) for j in range(i, len(wt))]
    plain.append(wt[i:])
    return "".join(ann), "".join(plain)


def wrap(s, w=80):
    return "\n".join(s[i : i + w] for i in range(0, len(s), w))


seqs = {}
for rna, fname in (("VAI", "VA_RNA_I_mock.fa"), ("VAII", "VA_RNA_II_mock.fa")):
    out = [f"; MOCK DEMO DATA. WT = OP218818.1 {REGIONS[rna][0]}-{REGIONS[rna][1]}; MOCK* variants are invented."]
    for g in GENOTYPES:
        ann, plain = annotate(rna, MOCK[rna].get(g, []))
        seqs[(rna, g)] = plain
        out.append(f">{rna}_{g}\n{wrap(ann)}")
    # FASTA readers skip lines before the first '>' — the comment line is informational
    (HERE / fname).write_text("\n".join(out) + "\n")

rng = random.Random(42)


def va_read(seq: str) -> str:
    """90-nt read from a VA transcript: half start with TSO + 5' end (as seen in real
    10x data), the rest are random windows. Random strand, 20% carry 1 sequencing error."""
    if rng.random() < 0.5:
        r = (TSO + seq)[:READ_LEN]
    else:
        s = rng.randrange(0, len(seq) - READ_LEN + 1)
        r = seq[s : s + READ_LEN]
    if rng.random() < 0.5:
        r = rc(r)
    if rng.random() < 0.2:
        p = rng.randrange(READ_LEN)
        r = r[:p] + rng.choice([b for b in "ACGT" if b != r[p]]) + r[p + 1 :]
    return r


def locus_read(seq: str, centre: int) -> str:
    """Read whose window covers position `centre` (used for the contamination sample)."""
    s = max(0, min(len(seq) - READ_LEN, centre - rng.randrange(20, READ_LEN - 20)))
    return seq[s : s + READ_LEN]


def background(n_random=4000, n_viral=300):
    reads = ["".join(rng.choice("ACGT") for _ in range(READ_LEN)) for _ in range(n_random)]
    for _ in range(n_viral):  # viral reads from outside the VA locus
        s = rng.randrange(15000, len(genome) - READ_LEN)
        r = genome[s : s + READ_LEN]
        reads.append(r if rng.random() < 0.5 else rc(r))
    return reads


# (sample, genotype, n VA reads, contaminating WT VA II reads at the MOCKB deletion)
SAMPLES = [
    ("DEMO_WT_1", "WT", 40, 0),
    ("DEMO_MOCKA_1", "MOCKA", 40, 0),
    ("DEMO_MOCKB_1", "MOCKB", 40, 0),
    ("DEMO_MOCKAB_1", "MOCKAB", 40, 0),
    ("DEMO_MOCKC_1", "MOCKC", 40, 0),
    ("DEMO_MOCKB_2", "MOCKB", 30, 6),  # contaminated with WT reads -> should be flagged
    ("DEMO_MOCKA_3", "MOCKA", 0, 0),   # no VA reads captured -> uninformative
]

fq_dir = HERE / "fastq"
fq_dir.mkdir(exist_ok=True)
for sample, g, n_va, n_wt in SAMPLES:
    reads = background()
    for _ in range(n_va):
        rna = rng.choice(["VAI", "VAII"])
        reads.append(va_read(seqs[(rna, g)]))
    for _ in range(n_wt):
        reads.append(locus_read(WT["VAII"], 123))
    rng.shuffle(reads)
    with gzip.GzipFile(fq_dir / f"{sample}_S1_L001_R2_001.fastq.gz", "wb", mtime=0) as fh:
        for i, r in enumerate(reads, 1):
            fh.write(f"@{sample}:{i}\n{r}\n+\n{'I' * len(r)}\n".encode())
    print(f"{sample:15s} {g:7s} reads={len(reads)} VA={n_va} WT-contam={n_wt}")
