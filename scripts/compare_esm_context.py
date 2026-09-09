"""Paired, homology-grouped experiment; does not replace the production model."""
import argparse
import hashlib
import json
import platform
import subprocess
import time
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import sklearn
import torch
import transformers
from sklearn.decomposition import PCA
from sklearn.dummy import DummyRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from transformers import AutoModel, AutoTokenizer

from cysmutml.amino_acids import physicochemical_features
from cysmutml.evaluation.metrics import regression_metrics
from cysmutml.models.pipeline import preprocessor

MODEL = "facebook/esm2_t6_8M_UR50D"
REVISION = "c731040fcd8d73dceaa04b0a8e6329b345b0f5df"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", required=True)
    parser.add_argument("--annotations", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    raw = pd.read_csv(args.features, low_memory=False)
    annotation = pd.read_csv(args.annotations).drop_duplicates("Entry")
    df = raw.merge(annotation, left_on="uniprot_id", right_on="Entry",
                   how="left", validate="many_to_one")
    reasons = []
    for row in df.itertuples():
        seq = row.canonical_sequence
        if not isinstance(seq, str) or not seq:
            reasons.append("no_sequence")
        elif not isinstance(row.Sequence, str) or seq != row.Sequence:
            reasons.append("unverified_reference")
        elif len(seq) > 1000 or set(seq) - set("ACDEFGHIKLMNPQRSTVWY"):
            reasons.append("unsupported_length_or_alphabet")
        elif not (1 <= row.position <= len(seq)) or seq[int(row.position)-1] != row.wt_aa:
            reasons.append("position_or_wt_mismatch")
        elif not np.isfinite(row.destabilization_ddg_kcal_mol):
            reasons.append("nonfinite_target")
        else:
            reasons.append("accepted")
    df["eligibility"] = reasons
    coverage = df.eligibility.value_counts().to_dict()
    df = df[df.eligibility.eq("accepted")].copy().reset_index(drop=True)
    if len(df) < 30:
        raise ValueError("Insufficient verified rows")
    sequences = sorted(df.canonical_sequence.unique())
    seq_ids = {s: "s" + hashlib.sha256(s.encode()).hexdigest() for s in sequences}
    df["seq_id"] = df.canonical_sequence.map(seq_ids)
    df["is_enzyme"] = df["Reviewed"].eq("reviewed") & df["EC number"].notna()
    # Duplicate exact-sequence mutations must not receive conflicting row identities.
    fasta = out / "sequences.fasta"
    fasta.write_text("".join(">" + seq_ids[s] + "\n" + s + "\n" for s in sequences))
    prefix = out / "mmseqs"
    command = ["mmseqs", "easy-cluster", str(fasta), str(prefix), str(out / "mmseqs_tmp"),
               "--min-seq-id", "0.30", "-c", "0.80", "--cov-mode", "0", "--threads", "2"]
    subprocess.run(command, check=True)
    clusters = pd.read_csv(str(prefix) + "_cluster.tsv", sep="\t",
                           names=["representative", "member"])
    mapping = dict(zip(clusters.member, clusters.representative, strict=True))
    if set(seq_ids.values()) != set(mapping):
        raise ValueError("Incomplete cluster mapping")
    # Join clusters that share a protein ID or accession, keeping sequence variants together.
    parent = {x: x for x in mapping.values()}

    def root(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    df["group"] = df.seq_id.map(mapping)
    for col in ["protein_id", "uniprot_id"]:
        for name, frame in df.dropna(subset=[col]).groupby(col, sort=True):
            if str(name).strip().lower() in {"", "unknown", "nan", "none"}:
                continue
            groups = sorted(frame.group.unique())
            for group in groups[1:]:
                parent[root(group)] = root(groups[0])
    df["group"] = df.group.map(root)
    if df.group.nunique() < 3:
        raise ValueError("Need at least three independent groups")
    df["row_id"] = np.arange(len(df))
    df[["row_id", "protein_id", "uniprot_id", "mutation", "seq_id", "group",
        "is_enzyme"]].to_csv(out / "manifest.csv", index=False)

    torch.set_num_threads(2)
    tokenizer = AutoTokenizer.from_pretrained(MODEL, revision=REVISION)
    model = AutoModel.from_pretrained(MODEL, revision=REVISION).eval()
    cache = out / "embeddings"
    cache.mkdir(exist_ok=True)
    embeddings = {}
    start = time.perf_counter()
    for seq in sequences:
        file = cache / (seq_ids[seq] + ".npy")
        if file.exists():
            emb = np.load(file, allow_pickle=False)
        else:
            with torch.inference_mode():
                emb = model(**tokenizer(seq, return_tensors="pt")).last_hidden_state[
                    0, 1:len(seq)+1
                ].numpy()
            np.save(file, emb)
        if emb.shape != (len(seq), 320) or not np.isfinite(emb).all():
            raise ValueError("Invalid embedding")
        embeddings[seq] = emb
    embedding_seconds = time.perf_counter() - start
    e = np.stack([embeddings[r.canonical_sequence][int(r.position)-1]
                  for r in df.itertuples()])
    # Explicit allowlist: no target statistics, IDs, EC annotations or fold labels.
    allowed = list(physicochemical_features("A", "C"))
    numeric = [c for c in allowed if pd.api.types.is_numeric_dtype(df[c])]
    categorical = [c for c in allowed if c not in numeric]
    x = df[allowed]
    y = df.destabilization_ddg_kcal_mol.to_numpy()
    oof = []
    folds = []
    for fold, (train, test) in enumerate(GroupKFold(3).split(x, y, df.group), 1):
        assert not set(df.group.iloc[train]) & set(df.group.iloc[test])
        assert not set(df.seq_id.iloc[train]) & set(df.seq_id.iloc[test])
        prep = preprocessor(numeric, categorical, True)
        xp = prep.fit_transform(x.iloc[train])
        xt = prep.transform(x.iloc[test])
        dim = min(16, len(train)-1, e.shape[1])
        reduction = make_pipeline(
            SimpleImputer(strategy="median"), StandardScaler(),
            PCA(n_components=dim, svd_solver="full"), StandardScaler(),
        )
        ep = reduction.fit_transform(e[train])
        et = reduction.transform(e[test])
        variants = {
            "dummy": (DummyRegressor(), xp, xt),
            "ridge_physchem": (Ridge(alpha=1.0), xp, xt),
            "ridge_esm": (Ridge(alpha=1.0), ep, et),
            "ridge_physchem_esm": (Ridge(alpha=1.0),
                                  np.hstack([xp, ep]), np.hstack([xt, et])),
        }
        folds.append(dict(fold=fold, train_rows=len(train), test_rows=len(test),
                          train_groups=df.group.iloc[train].nunique(),
                          test_groups=df.group.iloc[test].nunique()))
        for name, (estimator, a, b) in variants.items():
            estimator.fit(a, y[train])
            pred = estimator.predict(b)
            frame = df.iloc[test][["row_id", "mutation", "mut_aa", "group", "is_enzyme"]].copy()
            frame["fold"], frame["model"] = fold, name
            frame["observed"], frame["predicted"] = y[test], pred
            oof.append(frame)
            joblib.dump(dict(preprocessor=prep, embedding_reduction=reduction,
                             estimator=estimator, features=allowed, esm_model=MODEL,
                             esm_revision=REVISION, model_variant=name),
                        out / (name + "_fold" + str(fold) + ".joblib"))
    pred = pd.concat(oof, ignore_index=True)
    pred.to_csv(out / "out_of_fold_predictions.csv", index=False)
    metrics = []
    for name, frame in pred.groupby("model"):
        for population, mask in {
            "all": np.ones(len(frame), dtype=bool),
            "cys": frame.mut_aa.eq("C"),
            "enzymes": frame.is_enzyme,
            "enzymes_cys": frame.is_enzyme & frame.mut_aa.eq("C"),
        }.items():
            subset = frame.loc[mask]
            if len(subset) >= 2:
                metrics.append(dict(model=name, population=population, n=len(subset),
                                    groups=subset.group.nunique(),
                                    **regression_metrics(subset.observed, subset.predicted)))
    metrics = pd.DataFrame(metrics)
    metrics.to_csv(out / "metrics.csv", index=False)
    pd.DataFrame(folds).to_csv(out / "folds.csv", index=False)
    report = dict(
        source_features_sha256=digest(args.features),
        annotations_sha256=digest(args.annotations), eligibility=coverage,
        rows=len(df), sequences=len(sequences), groups=df.group.nunique(),
        enzyme_rows=int(df.is_enzyme.sum()),
        enzyme_cys_rows=int((df.is_enzyme & df.mut_aa.eq("C")).sum()),
        esm_model=MODEL, esm_revision=REVISION, embedding_seconds=embedding_seconds,
        pca_components=16, ridge_alpha=1.0, folds=3, mmseqs_command=command,
        mmseqs_version=subprocess.check_output(["mmseqs", "version"], text=True).strip(),
        python=platform.python_version(), sklearn=sklearn.__version__,
        torch=torch.__version__, transformers=transformers.__version__,
        git_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        production_model_changed=False,
        caveats=["Exploratory comparison; not an independent model-selection test.",
                 "Homology groups are operational MMseqs clusters, not proof of no remote homology.",
                 "Sequence identity and WT checks cannot resolve all construct-numbering ambiguity.",
                 "Missing enzyme annotation does not imply non-enzyme.",
                 "ESM pretraining overlap is not ruled out."],
    )
    (out / "run.json").write_text(json.dumps(report, indent=2, default=int))
    print(json.dumps(report, indent=2, default=int), flush=True)
    print(metrics.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
