#!/usr/bin/env python3
"""Self-test on the mock demo: screens all demo samples with the python matcher
(and bbduk, if installed), then checks every expected outcome.

Run from the repo root:  python tests/test_demo.py
"""

import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import core  # noqa: E402

EXPECTED = {  # sample prefix -> (compatible_with_all_hits, contradicting reads)
    "DEMO_WT_1": ("WT", 0),
    "DEMO_MOCKA_1": ("MOCKA", 0),
    "DEMO_MOCKB_1": ("MOCKB", 0),
    "DEMO_MOCKAB_1": ("MOCKAB", 0),
    "DEMO_MOCKC_1": ("MOCKC", 0),
    "DEMO_MOCKB_2": ("none", 6),  # 6 spiked WT reads
    "DEMO_MOCKA_3": ("WT,MOCKA,MOCKB,MOCKAB,MOCKC", 0),  # no hits
}

text = "\n".join(f.read_text() for f in sorted((ROOT / "demo").glob("VA_RNA_*_mock.fa")))
full = core.parse_annotated(text)[0]
assert core.genotypes(full) == ["WT", "MOCKA", "MOCKB", "MOCKAB", "MOCKC"]

engines = ["python"] + (["bbduk"] if core.bbduk_exe() else [])
summaries = {}
with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    fa = tmp / "panel.fa"
    core.write_panel(core.design_panel(full, 25), tmp / "panel.tsv", fa)
    sheet = core.discover_samples(ROOT / "demo" / "fastq", core.genotypes(full), core.aliases(full))
    assert len(sheet) == len(EXPECTED), sheet
    for eng in engines:
        for r in sheet.itertuples():
            out = tmp / eng / r.sample
            out.mkdir(parents=True)
            if eng == "bbduk":
                cmd = core.bbduk_cmd(r.r2, fa, out, 2, 1)
                rc = core.run_bbduk(cmd)
            else:
                cmd, rc = None, core.run_python_matcher(r.r2, fa, out, 25)
            assert rc == 0, (eng, r.sample)
            core.write_run_record(out, engine=eng, cmd=cmd, sample=r.sample, genotype=r.genotype, r2=r.r2,
                                  panel_fa=fa, k=25, max_reads=None, exit_code=rc, seconds=0,
                                  sequences_text=text)
        summ, _ = core.collect_results(tmp / eng, {})
        for row in summ.itertuples():
            prefix = row.sample.split("_S1_")[0]
            comp, contra = EXPECTED[prefix]
            assert row.compatible_with_all_hits == comp, (eng, prefix, row.compatible_with_all_hits)
            assert row.reads_contradict_expected == contra, (eng, prefix, row.reads_contradict_expected)
        summaries[eng] = summ.drop(columns=["engine", "timestamp"])
        print(f"{eng}: all {len(summ)} demo samples as expected")

if len(summaries) == 2:
    assert summaries["python"].equals(summaries["bbduk"]), "engines disagree"
    print("python and bbduk engines give identical summaries")
print("OK")
