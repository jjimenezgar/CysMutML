"""Audit enzyme coverage before an ESM pilot. Does not train or deploy a model."""
import argparse
import io
import json
import re
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(args.artifact) as archive:
        df = pd.read_csv(
            archive.open("data/processed/fireprotdb_mutations_aggregated.csv"),
            low_memory=False,
        )
    ids = sorted({
        x for x in df.uniprot_id.dropna().astype(str)
        if re.fullmatch(r"[A-Z0-9]+", x)
    })
    frames = []
    for start in range(0, len(ids), 30):
        query = " OR ".join("accession:" + x for x in ids[start:start + 30])
        params = urllib.parse.urlencode({
            "query": "(" + query + ")", "format": "tsv",
            "fields": "accession,reviewed,ec,sequence", "size": 500,
        })
        with urllib.request.urlopen(
            "https://rest.uniprot.org/uniprotkb/search?" + params, timeout=45
        ) as response:
            frames.append(pd.read_csv(io.BytesIO(response.read()), sep="\t"))
    annotations = pd.concat(frames).drop_duplicates("Entry")
    annotations.to_csv(out / "uniprot_annotations.csv", index=False)
    annotated = df.merge(
        annotations, left_on="uniprot_id", right_on="Entry",
        how="left", validate="many_to_one",
    )
    enzymes = annotated[
        annotated["EC number"].notna() & annotated["Reviewed"].eq("reviewed")
    ].copy()

    def matches(row):
        seq = row["canonical_sequence"]
        if not isinstance(seq, str) or seq != row["Sequence"]:
            return False
        position = int(row["position"])
        return 1 <= position <= len(seq) and seq[position - 1] == row["wt_aa"]

    eligible = enzymes.loc[enzymes.apply(matches, axis=1)].copy() if len(enzymes) else enzymes
    eligible.to_csv(out / "eligible_enzyme_mutations.csv", index=False)
    report = {
        "training_rows": len(df), "reviewed_ec_rows": len(enzymes),
        "eligible_rows": len(eligible),
        "eligible_cys_rows": int(eligible.mut_aa.eq("C").sum()),
        "eligible_sequences": eligible.canonical_sequence.nunique(),
        "excluded_mapping_rows": len(enzymes) - len(eligible),
        "policy": "Reviewed UniProt + EC; exact sequence and WT match; no offset guessing",
        "esm_measured": False,
    }
    (out / "feasibility.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
