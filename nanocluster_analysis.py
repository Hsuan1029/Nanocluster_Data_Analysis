"""Analyse a Quantum Cluster Database CSV and export standalone figures."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


FEATURES = [
    "energy_n_minus_one", "energy_n_plus_one", "homo_lumo_gap",
    "n_val_electrons", "magnetic_moment", "energy_dft", "energy_relative",
]
FAMILIES = {
    "alkali-metals": "Li Na K Rb Cs Fr",
    "alkaline-earth-metals": "Be Mg Ca Sr Ba Ra",
    "post-transition-metals": "Al Ga In Tl Sn Pb Fl Bi Po Lv",
    "metalloids": "B Si Ge As Sb Te Po At",
    "nonmetals": "H C N O F P S Cl Se Br I He Ne Ar Kr Xe Rn",
    "transition-metals": "Sc Ti V Cr Mn Fe Co Ni Cu Zn Y Zr Nb Mo Tc Ru Rh Pd Ag Cd "
                         "Hf Ta W Re Os Ir Pt Au Hg Rf Db Sg Bh Hs Mt Ds Rg Cn",
}


def save(fig, directory, name):
    fig.savefig(directory / f"{name}.png", dpi=160, bbox_inches="tight")
    plt.close(fig)


def draw_biplot(frame, loadings, variance, name, output, by_size=False):
    fig, ax = plt.subplots(figsize=(10, 7))
    palette = plt.get_cmap("tab20")
    for index, (element, rows) in enumerate(frame.groupby("element_symbol")):
        area = 1.8 * rows["n_atoms"] if by_size else 6
        ax.scatter(rows.PC1, rows.PC2, s=area, alpha=0.52, linewidths=0,
                   color=palette(index % 20), label=element)
    ax.axhline(0, color="0.83", lw=0.8)
    ax.axvline(0, color="0.83", lw=0.8)
    offsets = [(5, -11), (5, 9), (-20, -13), (-45, 10), (-42, -5), (4, -11), (5, 0)]
    for feature, point, offset in zip(FEATURES, loadings, offsets):
        ax.annotate("", xy=point, xytext=(0, 0), arrowprops={"arrowstyle": "->", "lw": 1})
        ax.annotate(feature, xy=point, xytext=offset, textcoords="offset points", fontsize=8)
    if name == "transition-metals":
        ax.set(xlim=(-1.55, 1.4), ylim=(-4.5, 5.3))
    elif name == "all-elements":
        ax.set(xlim=(-3, 5), ylim=(-6, 7.5))
    ax.set(xlabel=f"PC1 ({variance[0]:.1f}%)", ylabel=f"PC2 ({variance[1]:.1f}%)",
           title=name.replace("-", " ").capitalize())
    ax.legend(title="Element", loc="upper left", bbox_to_anchor=(1.01, 1),
              ncol=max(1, (frame.element_symbol.nunique() + 19) // 20), fontsize=8)
    save(fig, output, name)


def exploratory_plots(frame, output):
    import seaborn as sns

    numeric = frame[FEATURES + ["n_atoms"]]
    for column in numeric:
        fig, ax = plt.subplots()
        ax.hist(numeric[column], bins=20)
        ax.set(xlabel=column, ylabel="Clusters")
        save(fig, output, f"distribution-{column}")

    fig, ax = plt.subplots(figsize=(13, 4))
    frame.element_symbol.value_counts().plot.bar(ax=ax)
    ax.set(xlabel="Element", ylabel="Clusters")
    save(fig, output, "element-frequency")

    for column in ["n_atoms", "element_symbol"]:
        fig, ax = plt.subplots(figsize=(15, 6))
        sns.violinplot(data=frame, x=column, y="energy_dft", cut=0,
                       order=sorted(frame[column].unique()), ax=ax)
        ax.tick_params(axis="x", rotation=90)
        save(fig, output, f"energy-dft-by-{column}")

    fig, ax = plt.subplots(figsize=(9, 7))
    sns.heatmap(numeric.corr(), annot=True, fmt=".2f", center=0, cmap="vlag", ax=ax)
    save(fig, output, "property-correlations")
    grid = sns.pairplot(numeric, diag_kind="hist", plot_kws={"s": 2, "alpha": 0.25})
    save(grid.fig, output, "property-scatter-matrix")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path)
    parser.add_argument("--output", type=Path, default=Path(__file__).parent / "images")
    parser.add_argument("--eda", action="store_true", help="Also generate descriptive plots")
    args = parser.parse_args()
    columns = FEATURES + ["n_atoms", "element_symbol"]
    frame = pd.read_csv(args.csv, usecols=columns)
    total = len(frame)
    for column in FEATURES + ["n_atoms"]:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame["element_symbol"] = frame.element_symbol.astype("string").str.strip().replace("", pd.NA)
    frame = frame.replace([np.inf, -np.inf], np.nan).dropna().copy()
    if len(frame) < 6 or frame[FEATURES].nunique().max() < 2:
        parser.error("PCA needs at least six complete rows and variation in the features.")
    if (frame.n_atoms <= 0).any():
        parser.error("Atom counts must be positive.")
    scaled = StandardScaler().fit_transform(frame[FEATURES])
    model = PCA(n_components=6, svd_solver="full")
    scores = model.fit_transform(scaled)
    loadings = model.components_[:2].T * np.sqrt(model.explained_variance_[:2])
    # Component signs are arbitrary. Orient the axes to match the report.
    signs = np.sign([loadings[FEATURES.index("energy_relative"), 0],
                     loadings[FEATURES.index("magnetic_moment"), 1]])
    signs[signs == 0] = 1
    frame[["PC1", "PC2"]] = scores[:, :2] * signs
    loadings *= signs
    variance = 100 * model.explained_variance_ratio_
    print(f"Analysed {len(frame):,} of {total:,} rows. PC1 + PC2: {variance[:2].sum():.1f}%")
    args.output.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots()
    bars = ax.bar([f"PC{i}" for i in range(1, 7)], variance, color="skyblue")
    ax.bar_label(bars, fmt="%.1f%%", padding=3)
    ax.set(xlabel="Principal component", ylabel="Variance explained (%)", ylim=(0, max(variance)*1.15))
    save(fig, args.output, "scree")

    draw_biplot(frame, loadings, variance, "all-elements", args.output, by_size=True)
    for name, symbols in FAMILIES.items():
        group = frame[frame.element_symbol.isin(symbols.split())]
        if not group.empty:
            draw_biplot(group, loadings, variance, name, args.output)
    group = frame[frame.element_symbol.isin(["Y", "Sc", "Ni"])]
    if not group.empty:
        draw_biplot(group, loadings, variance, "transition-metals-by-size", args.output, by_size=True)
    if args.eda:
        exploratory_plots(frame, args.output)


if __name__ == "__main__":
    main()
