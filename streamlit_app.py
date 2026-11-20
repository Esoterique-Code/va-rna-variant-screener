"""VA-RNA Variant Screener — Streamlit GUI.

Demo (mock data):   streamlit run streamlit_app.py
Own data: set any of these environment variables before launching:
  VA_SEQ_DIR       folder with the annotated VA FASTA files (*.fa / *.fasta)
  VA_FASTQ_DIR     folder searched for *_R2_001.fastq.gz
  VA_RESULTS_DIR   where results are written (default: a per-session folder in results/)
  VA_EXTRA_PANELS  extra panel FASTA files offered in tab 2, separated by ':'
"""

from __future__ import annotations

import hashlib
import os
import sys
import time
import uuid
from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import core  # noqa: E402

env = lambda k: os.environ.get(k)  # noqa: E731
DEMO = not env("VA_SEQ_DIR")
SEQ_DIR = Path(env("VA_SEQ_DIR") or ROOT / "demo").expanduser()
FASTQ_DIR = Path(env("VA_FASTQ_DIR") or ROOT / "demo" / "fastq").expanduser()
EXTRA_PANELS = [Path(p).expanduser() for p in (env("VA_EXTRA_PANELS") or "").split(":") if p]
if "session" not in st.session_state:
    st.session_state.session = uuid.uuid4().hex[:8]
# Public deployments share one filesystem, so each browser session gets its own folder.
RESULTS = Path(env("VA_RESULTS_DIR") or ROOT / "results" / f"session_{st.session_state.session}").expanduser()
PANELS = RESULTS / "panels"
RESULTS.mkdir(parents=True, exist_ok=True)

st.set_page_config(page_title="VA-RNA Variant Screener", layout="wide")
st.title("VA-RNA Variant Screener")
st.caption("k-mer panel screen of 10x 3′ scRNA-seq R2 reads for adenovirus VA RNA I/II "
           "WT vs mutant sequence. Sample-level evidence only.")
if DEMO:
    st.info("**Demo mode.** The sequences are the public HAdV-5 WT VA RNAs (GenBank OP218818.1) with "
            "**invented** `MOCK*` mutations, and the FASTQs are small simulated files. See the About tab.",
            icon="🧪")

ABOUT_MD = r"""
### What this tool does
Some adenovirus mutants differ from wild type **only inside the VA RNA I / II genes**. This tool
asks whether the reads of a 10x 3′ single-cell RNA-seq sample contain sequence evidence for one
genotype or another. It works in four steps:
1. It builds a panel of short sequences (**k-mers**, default 25 nt) that tell the variants apart.
2. It searches each sample's cDNA reads (**R2**) for those k-mers, with exact matches on either strand.
3. For every matched read, it asks **which genotypes carry all of the panel k-mers in that read**.
4. It reports, per sample, which genotypes are compatible with all of the matched reads.

It was developed as a proof of concept for an MSc bioinformatics project. **This public version ships
mock data only.**

> **Read results cautiously.** VA RNAs are **not polyadenylated**, and 3′ chemistry captures them
> poorly. Zero hits are **uninformative**, not evidence of absence. Evidence is **sample-level
> only**: there are no per-cell genotype calls, and wording is "consistent with …".

### The demo data
- **Sequences:** `demo/VA_RNA_I_mock.fa`, `demo/VA_RNA_II_mock.fa`. WT = OP218818.1 10566–10843
  (VA I) and 10828–11083 (VA II). The invented variants are `MOCKA` (4-nt substitution in VA I),
  `MOCKB` (6-nt deletion in VA II), `MOCKAB` (both), and `MOCKC` (8-nt deletion in VA I + 3-nt
  substitution in VA II).
- **Reads:** `demo/fastq/`, 7 simulated samples of about 4,300 × 90-nt reads each. They contain random
  background, viral reads from outside the VA locus, and a few VA reads. About half the VA reads start
  with the 10x TSO, and 20 % carry one sequencing error.
- **Expected outcome:** each `DEMO_<G>_1` sample comes out *consistent with G*.
  `DEMO_MOCKB_2` is spiked with WT VA II reads at the MOCKB deletion, so it should be **flagged**.
  `DEMO_MOCKA_3` has no VA reads, so it should be **uninformative**.
- Regenerate everything with `python demo/make_demo_data.py` (seed 42, deterministic).

---
### Tab 1 · Panel design
**Inputs:** annotated FASTA for VA I and VA II. By default these come from the sequence folder;
you can also upload or paste them. Names are `VAI_<genotype>` / `VAII_<genotype>`, and
`VAI_WT` + `VAII_WT` are required.
Format: lowercase = genomic flank, UPPERCASE = mature VA RNA, `{XY}` = replacement bases,
`-` = one deleted base. Lines starting with `;` are comments.
Then choose the **k-mer length**.

| table | column | meaning |
|---|---|---|
| lengths | `name`, `rna` | sequence name; VAI or VAII |
| | `length` / `wt_length` | length after the mutations / length of that RNA's WT |
| | `check` | `ok`, `shorter/longer as annotated`, or ⚠ if the length doesn't match the annotation |
| mutations | `pos0` | 0-based position in the variant sequence |
| | `inserted` / `deleted_len` | replacement bases written in `{…}` / number of `-` |
| panel summary | `variant-specific` | k-mers found only in this variant, in no WT sequence (panel class `mutant`) |
| | `WT-specific` | WT k-mers this variant lacks (panel class `WT`) |

A variant with 0 + 0 k-mers is unmutated in that RNA, which is expected (e.g. `VAI_MOCKB`).
k-mers are canonical, meaning the two strands are collapsed. The panel is saved automatically to the
`panels/` folder of your results and can also be downloaded.

---
### Tab 2 · Screen samples
**Inputs:**
- **Sample sheet:** built from `*_R2_001.fastq.gz` in the FASTQ folder. R2 holds the cDNA; R1 is
  only the barcode + UMI. Tick `run`, and **check `expected genotype`**: it is guessed from the file
  name, and the interpretation depends on it.
- **Panel FASTA:** the panel from tab 1 (default), or any extra panel configured locally. *k* must
  match the panel.
- **Engine:**
  - **bbduk** (BBTools) is used if it is installed: `hdist=0 rcomp=t`, fast, and suitable for real
    10x files.
  - **python** is a pure-Python matcher with the same matching rule. It's for small files such as the
    demo, and it writes output in bbduk's format. Both engines match the **same reads**; on the demo
    their matched-read files are identical.
- **Quick test** screens only the first N reads. **Re-run** overwrites existing results.

**Outputs:** one folder per sample, containing:
- `bbduk_stats.txt`: hits per panel entry.
- `R2_panel_hits.fq.gz`: the matched reads.
- `run.json`: the engine, exact command, inputs, k, versions, runtime, and the sequences used.

---
### Tab 3 · Results
| column | meaning |
|---|---|
| `sample`, `expected_genotype` | sample name; what it *should* be |
| `engine` | bbduk or python |
| `total_reads`, `matched`, `per_million` | reads screened; reads with ≥1 panel k-mer; matched per million |
| `reads_fit_expected` | matched reads compatible with the expected genotype |
| `reads_contradict_expected` | matched reads holding a panel k-mer **absent** from the expected genotype. These should be 0 |
| `compatible_with_all_hits` | genotypes compatible with **every** matched read. The data can't rule these out |
| `subsampled` | N for a quick test, blank for a full run |
| `interpretation` | cautious, sample-level wording |

**Interpretation wording**
- *no panel hits: uninformative*: nothing can be concluded.
- *consistent with G; not distinguished from X*: X carries the same k-mers that were seen.
- *consistent with G; no other genotype fits all hits*: the strongest statement the tool makes.
- *n read(s) carry k-mers absent from G*: check sample identity, contamination or errors.

**Heatmap:** rows = samples; columns = `RNA: genotypes compatible with the read`; cell = number of reads.

**Per-read table** (shown under *Inspect sample*, and available as a download)
| column | meaning |
|---|---|
| `read_id` | read name from the FASTQ |
| `rna` | which RNA(s) the read's panel k-mers belong to |
| `n_panel_kmers` | how many distinct panel k-mers the read contains |
| `compatible_genotypes` | genotypes that carry **all** of them |
| `starts_with_TSO` | the read begins with the 10x template-switch oligo, i.e. right at a transcript's 5′ end |
| `read` | the sequence |

**Why not use the panel labels or bbduk's stats directly?**
1. A `X|WT` label means "a WT k-mer that X lacks". Every genotype that keeps WT sequence there
   carries it too, so WT-class hits are *expected* in many mutant samples.
2. bbduk's stats credit each read to **one** panel entry, even when the read holds several panel
   k-mers. Which entry it picks is engine-specific, and identical k-mers sit under several labels.

The per-read analysis uses every panel k-mer in every matched read, so it is engine-independent. The
raw stats are still shown under *Inspect sample*, for reference only.
"""

tab_about, tab_panel, tab_screen, tab_results = st.tabs(
    ["About", "1 · Panel design", "2 · Screen samples", "3 · Results"])


def seq_files(d: Path):
    return sorted([*d.glob("*.fa"), *d.glob("*.fasta")]) if d.exists() else []


# ======================================================================== About
with tab_about:
    st.markdown(ABOUT_MD)

# ======================================================================== 1
full = None
seq_text = ""
current_panel = None
with tab_panel:
    st.subheader("Annotated VA sequences")
    files = [f for f in seq_files(SEQ_DIR) if f.name != "OP218818.1.fasta"]
    src = st.radio("Source", [f"Sequence folder ({', '.join(f.name for f in files) or 'empty'})",
                              "Upload", "Paste"], horizontal=True)
    if src.startswith("Sequence folder"):
        seq_text = "\n".join(f.read_text() for f in files)
    elif src == "Upload":
        ups = st.file_uploader("FASTA file(s)", accept_multiple_files=True)
        seq_text = "\n".join(u.getvalue().decode() for u in ups or [])
    else:
        seq_text = st.text_area("Paste annotated FASTA", height=200)

    if seq_text.strip():
        full, roi, muts, warns = core.parse_annotated(seq_text)
        for w in warns:
            st.warning(w)
        missing = {"VAI_WT", "VAII_WT"} - full.keys()
        if missing:
            st.error(f"Missing WT reference(s): {', '.join(sorted(missing))}")
            full = None
    if full:
        c1, c2 = st.columns(2)
        c1.markdown("**Sequence lengths**")
        c1.dataframe(core.length_checks(full, muts), hide_index=True, width="stretch")
        c2.markdown("**Mutations**")
        c2.dataframe(muts, hide_index=True, width="stretch")

        k = st.slider("k-mer length", 15, 41, 25, help="25 is the default used so far.")
        panel = core.design_panel(full, k)
        summ = core.panel_summary(panel, full)
        st.markdown(f"**Panel: {len(panel)} k-mers (k={k})**")
        st.dataframe(summ, hide_index=True, width="stretch")
        zero = summ.loc[(summ["variant-specific"] == 0) & (summ["WT-specific"] == 0), "variant"].tolist()
        if zero:
            st.info("0 k-mers (this RNA is unmutated in the variant, which is expected): " + ", ".join(zero))
        with st.expander("Show k-mers"):
            st.dataframe(panel, hide_index=True, width="stretch")

        h = hashlib.sha1(seq_text.encode()).hexdigest()[:8]
        current_panel = PANELS / f"kmer_panel_k{k}_{h}.fa"
        if not current_panel.exists():
            core.write_panel(panel, current_panel.with_suffix(".tsv"), current_panel)
        st.caption(f"Panel saved as `{current_panel.name}` (+ `.tsv`). It is used by default in tab 2.")
        d1, d2 = st.columns(2)
        d1.download_button("Download panel TSV", current_panel.with_suffix(".tsv").read_bytes(),
                           current_panel.with_suffix(".tsv").name)
        d2.download_button("Download panel FASTA", current_panel.read_bytes(), current_panel.name)
    elif not seq_text.strip():
        st.info("Load or paste annotated sequences to build a panel.")

# ======================================================================== 2
with tab_screen:
    if not full:
        st.warning("Load valid sequences in tab 1 first.")
    else:
        st.subheader("Sample sheet")
        known = core.genotypes(full)
        raw_dir = Path(st.text_input("FASTQ folder (searched for *_R2_001.fastq.gz)", str(FASTQ_DIR)))
        key = (str(raw_dir), tuple(known))
        if st.session_state.get("sheet_key") != key:
            st.session_state.sheet = core.discover_samples(raw_dir, known, core.aliases(full)).assign(run=DEMO)
            st.session_state.sheet_key = key
        sheet = st.data_editor(
            st.session_state.sheet, hide_index=True, width="stretch", key="sheet_editor",
            column_config={
                "run": st.column_config.CheckboxColumn("run"),
                "genotype": st.column_config.SelectboxColumn("expected genotype", options=known + ["unknown"]),
                "sample": st.column_config.TextColumn(disabled=True),
                "r2": st.column_config.TextColumn("R2 path", disabled=True),
            },
        )
        if sheet.empty:
            st.warning("No R2 FASTQs found in that folder.")

        panels = [current_panel] + [p for p in EXTRA_PANELS if p.exists()]
        c1, c2, c3, c4 = st.columns(4)
        panel_fa = c1.selectbox("Panel FASTA", panels, format_func=lambda p: p.name)
        k_run = c2.number_input("k", 15, 41, int(panel_fa.stem.split("_k")[1].split("_")[0])
                                if "_k" in panel_fa.stem else 25)
        engines = (["bbduk"] if core.bbduk_exe() else []) + ["python"]
        engine = c3.selectbox("Engine", engines,
                              help="bbduk for real data; python matcher for small files (same matching rule).")
        threads = c4.number_input("threads (bbduk)", 1, 64, min(8, os.cpu_count() or 1))
        c5, c6, c7 = st.columns(3)
        mem = c5.number_input("memory GB (bbduk)", 1, 128, core.default_mem_gb())
        quick = c6.checkbox("Quick test: first N reads only")
        max_reads = c6.number_input("N reads", 1_000, 500_000_000, 5_000_000, step=1_000_000) if quick else None
        rerun = c7.checkbox("Re-run samples that already have results")
        if engine == "python":
            st.caption("The python matcher processes roughly 10⁵ reads per second. Use bbduk for full 10x files.")

        todo = sheet[sheet["run"]]
        if st.button(f"Screen {len(todo)} sample(s)", type="primary", disabled=todo.empty):
            for row in todo.itertuples(index=False):
                out = RESULTS / row.sample
                if (out / "bbduk_stats.txt").exists() and not rerun:
                    st.info(f"{row.sample}: results exist, skipped (tick re-run to overwrite).")
                    continue
                out.mkdir(parents=True, exist_ok=True)
                with st.status(f"{row.sample} ({row.genotype})", expanded=True) as status:
                    box, lines = st.empty(), []

                    def log(line, _box=box, _lines=lines):
                        _lines.append(line)
                        _box.code("\n".join(_lines[-25:]))

                    t0 = time.time()
                    if engine == "bbduk":
                        cmd = core.bbduk_cmd(row.r2, panel_fa, out, threads, mem, k_run, max_reads)
                        st.code(" ".join(map(str, cmd)), language="bash")
                        rc = core.run_bbduk(cmd, log)
                    else:
                        cmd = None
                        rc = core.run_python_matcher(row.r2, panel_fa, out, k_run, max_reads, log)
                    core.write_run_record(out, engine=engine, cmd=cmd, sample=row.sample,
                                          genotype=row.genotype, r2=row.r2, panel_fa=panel_fa, k=k_run,
                                          max_reads=max_reads, exit_code=rc, seconds=time.time() - t0,
                                          sequences_text=seq_text)
                    if rc == 0:
                        total, matched, _ = core.parse_stats(out / "bbduk_stats.txt")
                        status.update(label=f"{row.sample}: {matched:,} / {total:,} reads matched",
                                      state="complete", expanded=False)
                    else:
                        status.update(label=f"{row.sample}: failed (exit {rc})", state="error")
            st.success("Done. See the Results tab.")

# ======================================================================== 3
with tab_results:
    st.subheader("Panel hits across samples")
    st.warning(
        "VA RNAs are non-polyadenylated, and 3′ chemistry captures them poorly: **absence of hits "
        "does not mean absence of the variant.** Hits are sample-level evidence only; no per-cell "
        "genotype calls are made."
    )
    summ, reads = core.collect_results(RESULTS, full or {})
    if summ.empty:
        st.info("No results yet. Screen samples in tab 2.")
    else:
        summ["interpretation"] = summ.apply(core.interpret, axis=1)
        st.dataframe(summ, hide_index=True, width="stretch")
        st.download_button("Download summary TSV", summ.to_csv(sep="\t", index=False),
                           "va_screen_summary.tsv", "text/tab-separated-values")

        if not reads.empty:
            reads["group"] = reads["rna"] + ": " + reads["compatible_genotypes"]
            mat = reads.pivot_table(index="sample", columns="group", values="read_id",
                                    aggfunc="count", fill_value=0)
            fig = px.imshow(mat, text_auto=True, aspect="auto", color_continuous_scale="Blues",
                            labels=dict(color="reads", x="RNA: genotypes compatible with the read"))
            fig.update_layout(height=150 + 40 * len(mat))
            st.plotly_chart(fig, width="stretch")
            st.caption("Each matched read is placed in the column of the genotypes that carry **all** of its "
                       "panel k-mers. A read rules out every genotype not listed.")
            st.download_button("Download per-read TSV", reads.drop(columns="group").to_csv(sep="\t", index=False),
                               "va_screen_reads.tsv", "text/tab-separated-values")

        pick = st.selectbox("Inspect sample", summ["sample"])
        d = RESULTS / pick
        if (d / "run.json").exists():
            with st.expander("Run record (run.json)"):
                st.json((d / "run.json").read_text())
        with st.expander("Raw stats: reads credited per panel entry (bbduk_stats.txt)"):
            st.caption("Each read is credited to ONE panel entry and the choice is engine-specific, so this table "
                       "is kept for reference only. The interpretation uses the per-read analysis below.")
            st.dataframe(core.parse_stats(d / "bbduk_stats.txt")[2], hide_index=True, width="stretch")
        sel = reads[reads["sample"] == pick].drop(columns=["sample", "group"], errors="ignore")
        with st.expander(f"Matched reads ({len(sel)})"):
            st.dataframe(sel, hide_index=True, width="stretch")
