# -*- coding: utf-8 -*-
"""Control-size check for Fig. 2 (Methods S3, Table S3), called from Paper_figures.ipynb 3.1.2.
Runs pooled_check.R (all weeds pooled vs wheat x wheat; each Fig. 2 call vs the other weeds pooled) and vs_alone.R (each weed vs wheat
alone), and adds three sheets
to the enriched-ASV workbook written by da_core_check.py, saved as Table_S3_enriched_ASVs.xlsx."""
import os, subprocess, shutil
import pandas as pd


def run(matrix_dir, fig2_values_csv, out_dir, rscript=None):
    rscript = rscript or shutil.which("Rscript") or "Rscript"
    r = subprocess.run([rscript, "pooled_check.R", matrix_dir, fig2_values_csv, out_dir], capture_output=True, text=True)
    if r.returncode != 0: raise RuntimeError(r.stderr[-2000:])
    r = subprocess.run([rscript, "vs_alone.R", matrix_dir, out_dir], capture_output=True, text=True)
    if r.returncode != 0: raise RuntimeError(r.stderr[-2000:])
    return pd.read_csv(os.path.join(out_dir, "pooled_weeds_vs_wxw.csv")), pd.read_csv(os.path.join(out_dir, "fig2_neighbour_vs_other_weeds.csv"))


def classify(fig2, pooled, spec, alone=None):
    k = ["kingdom", "soil", "ASV"]
    f = fig2.merge(pooled[k + ["n_methods"]], on=k, how="left")
    f["called_all_weeds_vs_WxW"] = f.n_methods.fillna(0) >= 2
    f = f.merge(spec[k + ["neighbour", "DESeq2_l", "DESeq2_q", "edgeR_l", "edgeR_q"]], on=k + ["neighbour"], how="left")
    f["specific_both_methods"] = (f.DESeq2_q < 0.05) & (f.edgeR_q < 0.05)
    f["specific_one_method"] = (f.DESeq2_q < 0.05) | (f.edgeR_q < 0.05)
    f = f.drop(columns="n_methods")
    if alone is not None:
        a = alone[k + ["neighbour", "n_methods", "DESeq2_l"]].rename(columns={"n_methods": "methods_vs_alone", "DESeq2_l": "DESeq2_l_vs_alone"})
        f = f.merge(a, on=k + ["neighbour"], how="left"); f["enriched_vs_wheat_alone"] = f.methods_vs_alone.fillna(0) >= 2
    return f


def summary(calls):
    a = calls.drop_duplicates(["kingdom", "soil", "ASV"])
    s = dict(calls=len(calls), ASVs=len(a), ASVs_called_all_weeds=int(a.called_all_weeds_vs_WxW.sum()),
             calls_specific_both=int(calls.specific_both_methods.sum()), calls_specific_one=int(calls.specific_one_method.sum()))
    if "enriched_vs_wheat_alone" in calls:
        v = calls[calls.enriched_vs_wheat_alone]
        s.update(calls_vs_alone=len(v), ASVs_vs_alone=v.drop_duplicates(["kingdom", "soil", "ASV"]).shape[0],
                 share_positive_vs_alone=round(float((calls.DESeq2_l_vs_alone > 0).mean()), 2))
    return s


def table_s3(enriched_xlsx, calls, pooled, out_xlsx):
    x = pd.read_excel(enriched_xlsx, sheet_name=None)
    readme = x.pop("README"); readme.columns = ["Table S3"]
    readme.loc[len(readme)] = ["Control-size check (Methods S3): control_size_calls, each Fig. 2 call against the other 27 weeds pooled "
                               "(DESeq2, edgeR; one-sided; Benjamini-Hochberg within kingdom x soil) and whether its ASV is called when all "
                               "weeds are pooled against wheat x wheat, and (enriched_vs_wheat_alone) whether it is also enriched relative to wheat alone with the "
                               "Fig. 2 rule; pooled_weeds_vs_WxW, ASVs called by at least one method in that contrast. metagenomeSeq P values are "
                               "Bonferroni-adjusted, DESeq2 and edgeR Benjamini-Hochberg."]
    c = calls.copy(); p = pooled[pooled.n_methods >= 1].copy()
    for d in (c, p): d["ASV"] = d.ASV.str.replace("Cluster_", "ASV_")
    with pd.ExcelWriter(out_xlsx) as w:
        readme.to_excel(w, sheet_name="README", index=False)
        for k, v in x.items(): v.to_excel(w, sheet_name=k, index=False)
        c.round(4).to_excel(w, sheet_name="control_size_calls", index=False)
        p.round(4).to_excel(w, sheet_name="pooled_weeds_vs_WxW", index=False)
    return out_xlsx
