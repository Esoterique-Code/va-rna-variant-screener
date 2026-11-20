"""Core logic for the VA-RNA Variant Screener.

Pure functions (parsing, k-mer panel design, running the screen, interpreting
hits) kept separate from the Streamlit UI in streamlit_app.py.
"""

from __future__ import annotations

import datetime as dt
import gzip
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from importlib import metadata
from pathlib import Path

import pandas as pd

# 10x Genomics template-switch oligo (reads starting right at a transcript 5' end)
TSO = "AAGCAGTGGTATCAACGCAGAGTACATGGG"

# {...} groups (mutations, may contain dashes) or bare runs of dashes (deletions)
MUT = re.compile(r"\{([^}]*)\}|(-+)")
_COMP = str.maketrans("ACGT", "TGCA")


# ====================================================================== parsing
def read_fasta_text(text: str):
    """Yield (name, seq) from FASTA text; name = first word after '>'.
    Lines starting with ';' are comments and are ignored."""
    name, buf = None, []
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(";"):
            continue
        if line.startswith(">"):
            if name:
                yield name, "".join(buf)
            name, buf = line[1:].split()[0], []
        elif line and name:
            buf.append(re.sub(r"\s", "", line))
    if name:
        yield name, "".join(buf)


def parse_annotated(text: str):
    """Parse annotated VA FASTA.

    Format: lowercase = genomic flank, UPPERCASE = mature VA RNA,
    {XY} = mutated bases (replacement), '-' = deleted base.
    Returns (full: dict, roi: dict, mutations: DataFrame, warnings: list).
    """
    full, roi, rows, warnings = {}, {}, [], []
    for name, s in read_fasta_text(text):
        pos, clean, last = 0, [], 0
        for m in MUT.finditer(s):
            seg = s[last : m.start()]
            clean.append(seg)
            pos += len(seg)
            body = m.group(1) if m.group(1) is not None else m.group(2)
            ins = body.replace("-", "")
            rows.append((name, pos, ins.upper(), body.count("-")))
            clean.append(ins)
            pos += len(ins)
            last = m.end()
        clean.append(s[last:])
        seq = "".join(clean)
        bad = set(re.sub(r"[acgtACGT]", "", seq))
        if bad:
            warnings.append(f"{name}: unexpected characters {sorted(bad)}")
        full[name] = seq.upper()
        unwrapped = MUT.sub(lambda m: (m.group(1) or "").replace("-", ""), s)
        up = re.search(r"[A-Z][A-Za-z]*[A-Z]", unwrapped)
        roi[name] = up.group(0).upper() if up else ""
    muts = pd.DataFrame(rows, columns=["variant", "pos0", "inserted", "deleted_len"])
    return full, roi, muts, warnings


def rna_of(name: str) -> str:
    return "VAII" if name.startswith("VAII") else "VAI"


def genotype_of(variant: str) -> str:
    """'VAI_MOCKA(alias1)' -> 'MOCKA'; 'VAII_WT' -> 'WT'."""
    return re.sub(r"\(.*\)", "", variant.split("_", 1)[1])


def genotypes(full: dict) -> list[str]:
    """Genotypes in the order they appear, WT first."""
    seen = []
    for n in full:
        g = genotype_of(n)
        if g not in seen:
            seen.append(g)
    return sorted(seen, key=lambda g: g != "WT")


def aliases(full: dict) -> dict[str, str]:
    """Parenthesised aliases: 'VAI_MOCKA(alias1)' -> {'ALIAS1': 'MOCKA'}."""
    out = {}
    for n in full:
        m = re.search(r"\((.+)\)", n)
        if m:
            out[m.group(1).upper()] = genotype_of(n)
    return out


def length_checks(full: dict, muts: pd.DataFrame) -> pd.DataFrame:
    wt_len = {core_rna: len(full[f"{core_rna}_WT"]) for core_rna in ("VAI", "VAII") if f"{core_rna}_WT" in full}
    rows = []
    for n, s in full.items():
        r = rna_of(n)
        exp = wt_len.get(r)
        m = muts[muts["variant"] == n]
        dels, ins = int(m["deleted_len"].sum()), int(m["inserted"].str.len().sum())
        if exp is None:
            check = "no WT to compare"
        elif len(s) == exp:
            check = "ok"
        elif len(s) == exp - dels + ins:
            check = f"{'shorter' if len(s) < exp else 'longer'} as annotated (−{dels} / +{ins})"
        else:
            check = "⚠ length does not match annotation"
        rows.append(dict(name=n, rna=r, length=len(s), wt_length=exp, check=check))
    return pd.DataFrame(rows)


# ====================================================================== k-mers
def rc(s: str) -> str:
    return s.translate(_COMP)[::-1]


def _canon(s: str, k: int) -> set[str]:
    return {min(x, rc(x)) for x in (s[i : i + k] for i in range(len(s) - k + 1))}


def design_panel(full: dict, k: int = 25) -> pd.DataFrame:
    """mutant = k-mers only in the variant (in no WT sequence);
    WT = WT k-mers of the same RNA that the variant lacks. Canonical k-mers."""
    wt = {"VAI": full["VAI_WT"], "VAII": full["VAII_WT"]}
    wt_all = _canon(wt["VAI"], k) | _canon(wt["VAII"], k)
    rows = []
    for name, s in full.items():
        if name.endswith("_WT"):
            continue
        rna = rna_of(name)
        own = _canon(s, k)
        rows += [(rna, name, "mutant", x) for x in sorted(own - wt_all)]
        rows += [(rna, name, "WT", x) for x in sorted(_canon(wt[rna], k) - own)]
    return pd.DataFrame(rows, columns=["rna", "variant", "class", "kmer"])


def panel_summary(panel: pd.DataFrame, full: dict) -> pd.DataFrame:
    variants = [n for n in full if not n.endswith("_WT")]
    counts = (
        panel.groupby(["variant", "class"]).size().unstack(fill_value=0)
        .reindex(index=variants, columns=["mutant", "WT"], fill_value=0)
    )
    counts.columns = ["variant-specific", "WT-specific"]
    return counts.reset_index()


def write_panel(panel: pd.DataFrame, tsv: Path, fa: Path) -> None:
    """TSV + FASTA. FASTA header = rna|variant|class|NR, NR = line number in the TSV."""
    Path(fa).parent.mkdir(parents=True, exist_ok=True)
    panel.to_csv(tsv, sep="\t", index=False)
    with open(fa, "w") as fh:
        rows = zip(panel["rna"], panel["variant"], panel["class"], panel["kmer"])
        for nr, (rna, var, cls, kmer) in enumerate(rows, start=2):
            fh.write(f">{rna}|{var}|{cls}|{nr}\n{kmer}\n")


# ====================================================================== samples
SAMPLE_RE = re.compile(r"^(?P<sample>.+?)_R2_001\.f(ast)?q\.gz$")


def guess_genotype(sample: str, known: list[str], alias: dict[str, str] | None = None) -> str:
    tokens = re.split(r"[_\-.]", sample.upper())
    for g in sorted(known, key=len, reverse=True):  # longest first: MOCKAB before MOCKA
        G = g.upper()
        if any(t == G or (t.startswith(G) and (t[len(G):].isdigit() or t[len(G):].startswith("DL")))
               for t in tokens):
            return g
    for a, g in (alias or {}).items():
        if a in tokens:
            return g
    return "unknown"


def discover_samples(raw_dir: Path, known: list[str], alias=None) -> pd.DataFrame:
    rows = []
    for p in sorted(Path(raw_dir).expanduser().glob("*_R2_001.f*q.gz")):
        m = SAMPLE_RE.match(p.name)
        if m:
            s = m.group("sample")
            rows.append(dict(run=False, sample=s, genotype=guess_genotype(s, known, alias), r2=str(p)))
    return pd.DataFrame(rows, columns=["run", "sample", "genotype", "r2"])


# ====================================================================== engines
def bbduk_exe() -> str | None:
    """bbduk.sh (bioconda) or bbduk (some distro packages), if installed."""
    return shutil.which("bbduk.sh") or shutil.which("bbduk")


def default_mem_gb() -> int:
    try:
        total = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / 2**30
        return max(1, min(16, int(total // 2)))
    except (ValueError, OSError, AttributeError):
        return 4


def bbduk_cmd(r2, panel_fa, outdir, threads=8, mem_gb=16, k=25, max_reads=None):
    outdir = Path(outdir)
    cmd = [
        bbduk_exe() or "bbduk.sh", f"-Xmx{mem_gb}g", f"in={r2}", f"ref={panel_fa}",
        f"k={k}", "hdist=0", "rcomp=t", f"threads={threads}",
        f"stats={outdir / 'bbduk_stats.txt'}",
        f"outm={outdir / 'R2_panel_hits.fq.gz'}", "overwrite=t",
    ]
    if max_reads:
        cmd.append(f"reads={int(max_reads)}")
    return cmd


def run_bbduk(cmd, log_cb=None) -> int:
    """Run bbduk, streaming its output lines to log_cb. Returns exit code."""
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    for line in proc.stdout:
        if log_cb:
            log_cb(line.rstrip())
    return proc.wait()


def run_python_matcher(r2, panel_fa, outdir, k=25, max_reads=None, log_cb=None) -> int:
    """Pure-Python stand-in for bbduk (hdist=0, rcomp=t), for small files only.

    A read matches if any of its k-mers (either strand) is in the panel. Like
    bbduk's stats, each matched read is credited to ONE panel entry: the one
    with the most k-mer hits in that read (ties and duplicate k-mers -> lowest NR).
    Writes bbduk_stats.txt and R2_panel_hits.fq.gz in bbduk's format.
    """
    outdir = Path(outdir)
    names, kmer_to = {}, {}
    name = None
    for line in Path(panel_fa).read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]
        elif name:
            nr = int(name.rsplit("|", 1)[1])
            names[nr] = name
            c = min(line.strip(), rc(line.strip()))
            kmer_to[c] = min(nr, kmer_to.get(c, nr))
    total, counts = 0, {}
    with gzip.open(r2, "rt") as fin, gzip.open(outdir / "R2_panel_hits.fq.gz", "wt") as fout:
        while True:
            rec = [fin.readline() for _ in range(4)]
            if not rec[0]:
                break
            total += 1
            seq = rec[1].strip().upper()
            hits = {}
            for i in range(len(seq) - k + 1):
                x = seq[i : i + k]
                nr = kmer_to.get(min(x, rc(x)))
                if nr is not None:
                    hits[nr] = hits.get(nr, 0) + 1
            if hits:
                best = min(hits, key=lambda n: (-hits[n], n))
                counts[best] = counts.get(best, 0) + 1
                fout.writelines(rec)
            if max_reads and total >= max_reads:
                break
            if log_cb and total % 100_000 == 0:
                log_cb(f"{total:,} reads processed")
    matched = sum(counts.values())
    pct = lambda n: f"{100 * n / total:.5f}%" if total else "0.00000%"  # noqa: E731
    with open(outdir / "bbduk_stats.txt", "w") as fh:
        fh.write(f"#File\t{r2}\n#Total\t{total}\n#Matched\t{matched}\t{pct(matched)}\n#Name\tReads\tReadsPct\n")
        for nr, n in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])):
            fh.write(f"{names[nr]}\t{n}\t{pct(n)}\n")
    if log_cb:
        log_cb(f"Done: {matched:,} / {total:,} reads matched")
    return 0


def tool_versions(engine: str) -> dict:
    v = {"python": sys.version.split()[0]}
    if engine == "bbduk":
        try:
            out = subprocess.run([bbduk_exe(), "--version"], capture_output=True, text=True, timeout=60)
            m = re.search(r"BBTools version ([\d.]+)", out.stdout + out.stderr)
            v["bbduk"] = m.group(1) if m else "unknown"
        except Exception as e:  # noqa: BLE001
            v["bbduk"] = f"error: {e}"
    for pkg in ("streamlit", "pandas", "plotly"):
        try:
            v[pkg] = metadata.version(pkg)
        except metadata.PackageNotFoundError:
            v[pkg] = "not installed"
    return v


def write_run_record(outdir, *, engine, cmd, sample, genotype, r2, panel_fa, k, max_reads,
                     exit_code, seconds, sequences_text):
    rec = dict(
        sample=sample, expected_genotype=genotype, r2=str(r2), panel_fasta=str(panel_fa),
        k=k, max_reads=max_reads, engine=engine,
        command=shlex.join(map(str, cmd)) if cmd else "python matcher (core.run_python_matcher)",
        exit_code=exit_code, runtime_s=round(seconds, 1),
        timestamp=dt.datetime.now().isoformat(timespec="seconds"),
        versions=tool_versions(engine), sequences=sequences_text,
    )
    Path(outdir, "run.json").write_text(json.dumps(rec, indent=2))
    return rec


# ====================================================================== results
def parse_stats(path):
    """Parse bbduk stats=. Returns (total, matched, hits DataFrame)."""
    total = matched = 0
    rows = []
    for line in Path(path).read_text().splitlines():
        if line.startswith("#Total"):
            total = int(line.split("\t")[1])
        elif line.startswith("#Matched"):
            matched = int(line.split("\t")[1])
        elif line and not line.startswith("#"):
            name, reads = line.split("\t")[:2]
            rna, variant, cls, nr = name.split("|")
            rows.append((rna, variant, cls, nr, int(reads)))
    return total, matched, pd.DataFrame(rows, columns=["rna", "variant", "class", "nr", "reads"])


def genotype_kmer_sets(full: dict, k: int) -> dict[str, set[str]]:
    """Canonical k-mers carried by each genotype (its VA I + VA II)."""
    sets: dict[str, set[str]] = {}
    for name, seq in full.items():
        sets.setdefault(genotype_of(name), set()).update(_canon(seq, k))
    return sets


def panel_kmers(fa) -> dict[str, str]:
    """Canonical panel k-mer -> RNA (VAI / VAII), from the panel FASTA."""
    out, rna = {}, None
    for line in Path(fa).read_text().splitlines():
        if line.startswith(">"):
            rna = line[1:].split("|", 1)[0]
        elif rna:
            x = line.strip()
            out[min(x, rc(x))] = rna
    return out


def read_level(hits_fq, panel_fa, full: dict, k: int) -> pd.DataFrame:
    """One row per matched read, using EVERY panel k-mer the read contains.

    compatible_genotypes = genotypes that carry all of the read's panel k-mers.
    This does not depend on which single panel entry bbduk credited the read to
    (that choice is arbitrary when a read holds several k-mers, and identical
    k-mers sit under several panel labels), so bbduk and the python matcher give
    the same answer.
    """
    pk = panel_kmers(panel_fa)
    sets = genotype_kmer_sets(full, k)
    order = genotypes(full)
    rows = []
    with gzip.open(hits_fq, "rt") as fh:
        while True:
            rec = [fh.readline() for _ in range(4)]
            if not rec[0]:
                break
            seq = rec[1].strip().upper()
            found = {c for c in (min(x, rc(x)) for x in (seq[i : i + k] for i in range(len(seq) - k + 1)))
                     if c in pk}
            comp = [g for g in order if found <= sets[g]]
            rows.append(dict(
                read_id=rec[0][1:].split()[0], rna="+".join(sorted({pk[c] for c in found})),
                n_panel_kmers=len(found), compatible_genotypes=",".join(comp) or "none",
                starts_with_TSO=seq.startswith(TSO[:20]), read=seq,
            ))
    return pd.DataFrame(rows, columns=["read_id", "rna", "n_panel_kmers", "compatible_genotypes",
                                       "starts_with_TSO", "read"])


def collect_results(results_dir: Path, default_full: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per-sample summary + per-read table for every results_dir/<sample>/.

    Uses the sequences stored in each run.json when present, so a result is
    always interpreted against the sequences it was screened with.
    """
    summ, reads_all = [], []
    for stats in sorted(Path(results_dir).glob("*/bbduk_stats.txt")):
        d = stats.parent
        rec = json.loads((d / "run.json").read_text()) if (d / "run.json").exists() else {}
        full = parse_annotated(rec["sequences"])[0] if rec.get("sequences") else default_full
        panel_fa = rec.get("panel_fasta")
        if not full or not panel_fa or not Path(panel_fa).exists():
            continue
        all_g = genotypes(full)
        total, matched, _ = parse_stats(stats)
        k = int(rec.get("k", 25))
        hits_fq = d / "R2_panel_hits.fq.gz"
        if not hits_fq.exists():
            continue
        reads = read_level(hits_fq, panel_fa, full, k)
        g = rec.get("expected_genotype", "unknown")
        comp = reads["compatible_genotypes"].map(lambda s: s.split(","))
        compatible = [x for x in all_g if comp.map(lambda c, x=x: x in c).all()]
        fit = int(comp.map(lambda c: g in c).sum())
        summ.append(dict(
            sample=d.name, expected_genotype=g, engine=rec.get("engine", "bbduk"),
            total_reads=total, matched=matched,
            per_million=round(matched / total * 1e6, 3) if total else 0.0,
            reads_fit_expected=fit,
            reads_contradict_expected=(len(reads) - fit) if g != "unknown" else 0,
            compatible_with_all_hits=",".join(compatible) or "none",
            subsampled=rec.get("max_reads") or "", timestamp=rec.get("timestamp", ""),
        ))
        if len(reads):
            reads_all.append(reads.assign(sample=d.name))
    cols = ["sample", "read_id", "rna", "n_panel_kmers", "compatible_genotypes", "starts_with_TSO", "read"]
    reads_df = pd.concat(reads_all)[cols] if reads_all else pd.DataFrame(columns=cols)
    return pd.DataFrame(summ), reads_df


def interpret(row) -> str:
    """Sample-level, cautious wording only."""
    g, comp = row.expected_genotype, row.compatible_with_all_hits.split(",")
    if row.matched == 0:
        return "no panel hits: uninformative (VA RNAs are non-poly(A))"
    if g == "unknown":
        return f"hits compatible with: {row.compatible_with_all_hits}"
    if row.reads_contradict_expected:
        return (f"{row.reads_contradict_expected} read(s) carry k-mers absent from {g}: "
                "check sample identity / contamination / sequencing error")
    others = [x for x in comp if x != g]
    if not others:
        return f"consistent with {g}; no other genotype fits all hits"
    return f"consistent with {g}; not distinguished from {', '.join(others)}"
