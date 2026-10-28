"""Publication-style figures for every level (matplotlib, headless)."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

TYPE_COLORS = {"A": "#1a9850", "B": "#91cf60", "C": "#fee08b", "D": "#fc8d59", "E": "#4575b4", "F": "#d73027",
               "model_inadmissible": "#999999", "not_tested": "#555555",
               "cfa_misfit": "#bdbdbd", "cfa_fit_not_testable": "#d9d9d9", "model_not_estimable": "#999999"}
S_MARKERS = {"scalar_invariance": "o", "metric_invariance": "s", "configural_only": "^", "congruence_supported": "P",
             "structurally_discrepant": "D", "model_inadmissible": "x", "not_tested": "+",
             "configural_fit": "^", "cfa_misfit": "v", "cfa_fit_not_testable": "h", "model_not_estimable": "x"}
S_LABELS = {"scalar_invariance": "Scalar", "metric_invariance": "Metric only", "configural_only": "Configural (CFA)",
            "congruence_supported": "Congruence (gated ESEM)", "structurally_discrepant": "Discrepant (gated ESEM)",
            "model_inadmissible": "Model inadmissible", "not_tested": "Not tested",
            "configural_fit": "Configural fit (CFA)", "cfa_misfit": "Baseline CFA misfit",
            "cfa_fit_not_testable": "Baseline fit not testable", "model_not_estimable": "Model not estimable"}


def _save(fig, path: Path, dpi: int):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".png"), dpi=dpi, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def typology_map(typ: pd.DataFrame, labels: dict, out: Path, dpi: int = 200):
    """Scales in the (AUC, CCC) and (r, CCC) planes; colour = type, marker = structural status."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 6), sharey=True)
    for ax, x, xl in zip(axes, ["AUC", "r"], ["Scale-level AUC", "Pearson r"]):
        for _, row in typ.iterrows():
            ax.scatter(row[x], row["CCC"], color=TYPE_COLORS.get(row["final_type"], "grey"),
                       marker=S_MARKERS.get(row["axis_S"], "o"), s=70, alpha=0.85,
                       edgecolor="black" if S_MARKERS.get(row["axis_S"]) not in ("x", "+") else None, linewidth=0.4)
        ax.set_xlabel(xl)
    axes[0].set_ylabel("CCC")
    axes[0].set_title("(a) AUC x CCC")
    axes[1].set_title("(b) r x CCC")
    present = [k for k in TYPE_COLORS if k in set(typ["final_type"])]
    h1 = [Patch(facecolor=TYPE_COLORS[k], edgecolor="black", label=(f"{k}. " if len(k) == 1 else "") + labels.get(k, k)) for k in present]
    used = [k for k in S_MARKERS if k in set(typ["axis_S"])]
    h2 = [Line2D([0], [0], marker=S_MARKERS[k], color="grey", linestyle="", markersize=8, label=S_LABELS[k]) for k in used]
    fig.legend(handles=h1, loc="lower center", ncol=2, fontsize=8, bbox_to_anchor=(0.5, -0.16), title="Type")
    fig.legend(handles=h2, loc="lower center", ncol=4, fontsize=7.5, bbox_to_anchor=(0.5, -0.30), title="Structural status (marker)")
    _save(fig, out, dpi)


def association_heatmap(assoc: pd.DataFrame, out: Path, dpi: int = 200):
    """Spearman rho between item-level / structural indicators (rows) and scale-level accuracy (columns)."""
    if assoc is None or assoc.empty:
        return
    cols = ["r", "AUC", "CCC", "Cb"]
    a = assoc[assoc["y"].isin(cols)]
    rows = list(dict.fromkeys(a["x"]))
    if not rows:
        return
    mat = pd.DataFrame(np.nan, index=rows, columns=cols)
    sig = pd.DataFrame(False, index=rows, columns=cols)
    for _, r in a.iterrows():
        mat.loc[r["x"], r["y"]] = r["rho"]
        sig.loc[r["x"], r["y"]] = bool(pd.notna(r.get("p_BH")) and r["p_BH"] < 0.05)
    fig, ax = plt.subplots(figsize=(7.5, 0.6 * len(rows) + 2))
    cmap = plt.cm.RdBu_r.copy()
    cmap.set_bad("#d9d9d9")
    im = ax.imshow(np.ma.masked_invalid(mat.values.astype(float)), cmap=cmap, vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(len(cols)), cols)
    ax.set_yticks(range(len(rows)), rows, fontsize=8)
    for i in range(len(rows)):
        for j in range(len(cols)):
            v = mat.values[i, j]
            ax.text(j, i, "n/a" if np.isnan(v) else f"{v:.2f}{'*' if sig.values[i, j] else ''}", ha="center", va="center",
                    fontsize=8, color="white" if not np.isnan(v) and abs(v) > 0.55 else "black")
    fig.colorbar(im, ax=ax, fraction=0.04).set_label("Spearman rho")
    ax.set_title("Associations between levels (* = BH q < .05 within family)", fontsize=9)
    _save(fig, out, dpi)


def profile_ccc_hist(per_person: pd.DataFrame, out: Path, dpi: int = 200):
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.hist(per_person["cross_ccc_full_mean"], bins=40, color="#B0B7C6", alpha=0.85, label="With another person's twin")
    ax.hist(per_person["own_ccc_full"], bins=40, color="#4C8C4A", alpha=0.85, label="With own twin")
    ax.axvline(per_person["own_ccc_full"].mean(), color="#2E5C2C", ls="--")
    ax.axvline(per_person["cross_ccc_full_mean"].mean(), color="#5A6270", ls="--")
    ax.set_xlabel("Profile CCC (human-based z, equal instrument weights)")
    ax.set_ylabel("Respondents")
    ax.legend()
    _save(fig, out, dpi)


def scale_fidelity_scatter(master: pd.DataFrame, out: Path, dpi: int = 200):
    """r versus C_b for every scale (the descriptive two-component classification)."""
    fig, ax = plt.subplots(figsize=(7, 5.5))
    names = master["instrument_name"].unique()
    cmap = plt.cm.tab20
    for i, n in enumerate(names):
        g = master[master.instrument_name == n]
        ax.scatter(g["pearson_r"], g["cb"], s=45, color=cmap(i % 20), label=n, edgecolor="white", linewidth=0.5)
    ax.axvline(0.5, color="grey", ls=":")
    ax.axhline(0.8, color="grey", ls=":")
    ax.set_xlabel("Pearson r (individuation)")
    ax.set_ylabel("C_b (accuracy component of CCC)")
    ax.legend(fontsize=7, loc="lower right")
    _save(fig, out, dpi)


def item_fidelity_map(item_master: pd.DataFrame, out: Path, dpi: int = 200):
    fig, ax = plt.subplots(figsize=(7, 5.5))
    names = item_master["instrument_name"].unique()
    for i, n in enumerate(names):
        g = item_master[item_master.instrument_name == n]
        ax.scatter(g["spearman_rho"], g["sd_ratio"], s=20 + 300 * g["mae_range"], alpha=0.7, color=plt.cm.tab20(i % 20), label=n)
    ax.axhline(1, color="grey", ls=":")
    ax.set_xlabel("Item-level Spearman rho")
    ax.set_ylabel("SD ratio (twin / human)")
    ax.legend(fontsize=7)
    ax.set_title("Point size = normalised MAE", fontsize=9)
    _save(fig, out, dpi)


def variance_decomposition_bars(decomp: pd.DataFrame, out: Path, dpi: int = 200):
    d = decomp.copy()
    d["label"] = d["instrument_name"].str.slice(0, 28) + " | " + d["side"]
    d = d.sort_values(["instrument_name", "side"])
    fig, ax = plt.subplots(figsize=(8, 0.35 * len(d) + 1.5))
    left = np.zeros(len(d))
    for col, name, colr in (("pct_person", "Person", "#4C72B0"), ("pct_item", "Item", "#DD8452"),
                            ("pct_person_item_residual", "Person x item / residual", "#8C8C8C")):
        ax.barh(d["label"], d[col], left=left, label=name, color=colr)
        left += d[col].to_numpy()
    ax.set_xlabel("% of variance")
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3)
    _save(fig, out, dpi)


def structural_phi_forest(table: pd.DataFrame, phi_congruence: float, out: Path, dpi: int = 200):
    """Tucker's phi (gated point estimate) with its bootstrap CI, one row per scale with an admissible paired ESEM.

    Colour marks the structural route (congruence supported / structurally discrepant / other); a filled circle is a
    bootstrap-stable estimate (>= threshold share of valid replicates), a faint cross is not.
    """
    d = table[table["tucker_point_available"] == True].sort_values(["instrument_id", "scale_id"])
    if not len(d):
        return
    colour = {"congruence_supported": "#1a9850", "structurally_discrepant": "#fc8d59"}
    fig, ax = plt.subplots(figsize=(7, 0.4 * len(d) + 1.5))
    y = np.arange(len(d))[::-1]
    for i, row in zip(y, d.itertuples()):
        c = colour.get(row.structural_evidence_level, "#4575b4")
        stable = bool(row.tucker_bootstrap_stable)
        if pd.notna(row.tucker_phi_ci_lo):
            ax.plot([row.tucker_phi_ci_lo, row.tucker_phi_ci_hi], [i, i], color=c, lw=2 if stable else 1, alpha=1.0 if stable else 0.35)
        ax.scatter([row.tucker_phi_primary], [i], color=c, marker="o" if stable else "x", s=50, zorder=3,
                  alpha=1.0 if stable else 0.6)
    ax.axvline(phi_congruence, color="grey", ls=":", label=f"congruence threshold ({phi_congruence:g})")
    ax.set_yticks(y, d["scale_id"], fontsize=8)
    ax.set_xlim(0, 1.02)
    ax.set_xlabel("Tucker's phi (gated: both Human and Twin ESEM admissible)")
    handles = [Patch(facecolor=colour["congruence_supported"], label="Congruence supported"),
              Patch(facecolor=colour["structurally_discrepant"], label="Structurally discrepant"),
              Patch(facecolor="#4575b4", label="Other structural route"),
              Line2D([0], [0], marker="x", color="grey", linestyle="", label="Bootstrap not stable (< threshold valid replicates)")]
    ax.legend(handles=handles, fontsize=7, loc="lower right")
    _save(fig, out, dpi)
