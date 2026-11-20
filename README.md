# VA-RNA Variant Screener

A small Streamlit app that screens **10x 3′ single-cell RNA-seq reads** for sequence evidence that
distinguishes **wild-type adenovirus from VA RNA I / II mutants**. It builds a k-mer panel from
annotated VA sequences, searches the cDNA reads (R2) for it, and reports, per sample, which
genotypes are compatible with every matched read.

It's a proof of concept from an MSc bioinformatics project. **This repository ships mock data only**:
the WT sequences are public (GenBank OP218818.1, HAdV-5), and all `MOCK*` variants are invented.

> **Read results cautiously.** VA RNAs are not polyadenylated, and 3′ chemistry captures them poorly.
> Zero hits are uninformative, not evidence of absence. Evidence is sample-level only: there are no
> per-cell genotype calls.

## Try it
- **Online:** *(add the Streamlit Community Cloud link here once deployed)*
- **Locally:**
  ```bash
  conda env create -f environment.yml      # or: python -m pip install -r requirements.txt
  conda activate va-rna-variant-screener
  streamlit run streamlit_app.py
  ```
  Then click **2 · Screen samples → Screen 7 sample(s) → 3 · Results**.

## How it works
| Tab | Job |
|---|---|
| **About** | Purpose, inputs, and a glossary of every column |
| **1 · Panel design** | Parses the annotated FASTA (lowercase = flank, UPPERCASE = mature RNA, `{XY}` = replacement, `-` = deletion), checks lengths, and builds canonical k-mers. `mutant` = k-mers only in the variant; `WT` = WT k-mers the variant lacks |
| **2 · Screen samples** | Sample sheet from `*_R2_001.fastq.gz`. Runs either **bbduk** (`hdist=0 rcomp=t`, for real-size files) or a **pure-Python matcher** with the same rule (for small files). Each sample gets `bbduk_stats.txt`, `R2_panel_hits.fq.gz` and `run.json`, which records the provenance: command, versions and the sequences used |
| **3 · Results** | For every matched read: which genotypes carry **all** of its panel k-mers. Per sample: reads that fit or contradict the expected genotype, genotypes compatible with all hits, and cautious wording. Includes a heatmap and TSV downloads |

**Why per read, and not bbduk's stats?** bbduk credits each read to one panel entry, even when the
read contains several panel k-mers, and identical k-mers sit under several variant labels. The
per-read analysis uses every k-mer in every matched read, so the result doesn't depend on the engine.

## Demo data (`demo/`)
- **Sequences:** `VA_RNA_I_mock.fa` and `VA_RNA_II_mock.fa`. WT is OP218818.1 10566–10843 (VA I) and
  10828–11083 (VA II).
  - `MOCKA`: 4-nt substitution in VA I
  - `MOCKB`: 6-nt deletion in VA II
  - `MOCKAB`: both
  - `MOCKC`: 8-nt deletion in VA I + 3-nt substitution in VA II
- **Reads:** `fastq/` holds 7 simulated R2 files of about 4,300 × 90-nt reads each, about 1 MB in total.
  They contain random background, off-locus viral reads, and a few VA reads (half starting with the
  10x TSO, 20 % with one error).
- **Expected results:**
  - `DEMO_<G>_1` → *consistent with G; no other genotype fits all hits*
  - `DEMO_MOCKB_2` (spiked with 6 WT reads) → flagged: 6 reads contradict MOCKB
  - `DEMO_MOCKA_3` (no VA reads) → uninformative
- **Regenerate:** `python demo/make_demo_data.py` (deterministic, seed 42).
- **Self-test:** `python tests/test_demo.py` checks every expected outcome. It also checks that bbduk
  and the Python matcher agree, when bbduk is installed.

## Using your own data
Set environment variables before launching:

| variable | meaning |
|---|---|
| `VA_SEQ_DIR` | folder with your annotated VA FASTA files (`VAI_<genotype>`, `VAII_<genotype>`, including `_WT`) |
| `VA_FASTQ_DIR` | folder searched for `*_R2_001.fastq.gz` |
| `VA_RESULTS_DIR` | where results go. Default: a per-session folder under `results/` (gitignored) |
| `VA_EXTRA_PANELS` | extra panel FASTA files to offer in tab 2, separated by `:` |

Install bbduk (`conda install -c bioconda bbmap`) for full-size 10x files. The Python matcher handles
roughly 10⁵ reads per second.

## Deploying on Streamlit Community Cloud
Go to share.streamlit.io and choose **Create app**. Pick this repository, branch `main`, and main file
`streamlit_app.py`; Python 3.11 or newer. bbduk isn't installed there, so the app uses the Python
matcher, which is fine for the demo. Each browser session gets its own results folder, which is
cleared when the app restarts.
