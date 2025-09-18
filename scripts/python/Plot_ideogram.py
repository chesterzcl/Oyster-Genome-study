import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path
import matplotlib.patheffects as pe 


# === File paths ===
dir="/Users/lizicheng/Desktop/Data/dog_genome/oyster/"
coords_file = dir+"cv_vs_cv_all.filtered.coords"
fai_old = dir+"cv_old.fa.fai"
fai_new = dir+"cv_new.fa.fai"
out_png=dir+"cv_cv_ideogram.png"

# Parameters
min_len_bp   = 10_000          # drop blocks <10 kb
min_identity = 90.0            # drop blocks <90 %
bar_width    = 1.5e6           # width of chromosome bars
track_gap    = 6e6             # horizontal gap between assemblies
row_gap      = 5e7             # vertical gap between chromosomes
curve_alpha  = 0.35            # transparency of chords
curve_lw     = 0.4             # line width

# ---------- 1. load chromosome sizes ----------
def load_fai(path, top=10):
    df = pd.read_csv(
        path, sep="\t",
        names=["chr", "len", "off", "lineB", "lineW"],
        usecols=[0, 1]
    ).sort_values("len", ascending=False).head(top)
    return df.reset_index(drop=True)

ref_df  = load_fai(fai_old)
qry_df  = load_fai(fai_new)

ref_chrs = ref_df["chr"].tolist()
qry_chrs = qry_df["chr"].tolist()

ref_len  = dict(zip(ref_df["chr"], ref_df["len"]))
qry_len  = dict(zip(qry_df["chr"], qry_df["len"]))

# mapping chr → baseline y-offset
ref_y = {c: i * row_gap for i, c in enumerate(ref_chrs)}
qry_y = {c: i * row_gap for i, c in enumerate(qry_chrs)}

# ---------- 2. load .coords ----------
cols = [
    "ref_start", "ref_end", "qry_start", "qry_end",
    "ref_block", "qry_block", "identity",
    "ref_id", "qry_id", "ref_chr", "qry_chr"
]
coords = pd.read_csv(coords_file, sep=r"\s+", header=None, names=cols)

# ---------- 3. filter alignments ----------
coords = coords[
    (coords["ref_chr"].isin(ref_chrs)) &
    (coords["qry_chr"].isin(qry_chrs)) &
    (coords["identity"] >= min_identity) &
    ((coords["ref_end"] - coords["ref_start"]) >= min_len_bp) &
    ((coords["qry_end"] - coords["qry_start"]) >= min_len_bp)
]

print(f"Retaining {len(coords):,} alignments after filtering.")

# ---------- 4. plotting ----------
fig, ax = plt.subplots(figsize=(12, 10))
ax.set_axis_off()

# --- 4a. draw chromosome bars with patheffects for subtle outline
for chr_, length in ref_len.items():
    y = ref_y[chr_]
    bar = patches.Rectangle(
        (0, y), bar_width, length,
        color="#8e8e8e", ec="black", lw=0.3,
        path_effects=[pe.Stroke(linewidth=1, foreground='black'), pe.Normal()]
    )
    ax.add_patch(bar)
    ax.text(-1e6, y + length / 2, chr_, va="center", ha="right", fontsize=9, fontweight="bold")

for chr_, length in qry_len.items():
    y = qry_y[chr_]
    bar = patches.Rectangle(
        (track_gap, y), bar_width, length,
        color="#8e8e8e", ec="black", lw=0.3,
        path_effects=[pe.Stroke(linewidth=1, foreground='black'), pe.Normal()]
    )
    ax.add_patch(bar)
    ax.text(track_gap + bar_width + 1e6, y + length / 2, chr_,
            va="center", ha="left", fontsize=9, fontweight="bold")

# --- 4b. draw chords
for _, row in coords.iterrows():
    # midpoints for smoothness, or use start for precise locus
    y0 = ref_y[row.ref_chr] + (row.ref_start + row.ref_end) / 2
    y1 = qry_y[row.qry_chr] + (row.qry_start + row.qry_end) / 2
    path = Path(
        [(bar_width, y0),
         ((bar_width + track_gap) / 2, (y0 + y1) / 2),
         (track_gap, y1)],
        [Path.MOVETO, Path.CURVE3, Path.CURVE3]
    )
    col = "#2b83ba" if row.ref_start < row.ref_end else "#d7191c"
    patch = patches.PathPatch(
        path, lw=curve_lw, edgecolor=col, alpha=curve_alpha, facecolor="none"
    )
    ax.add_patch(patch)

# ---------- 5. finalize ----------
ax.set_xlim(-2e6, track_gap + bar_width + 2e6)
ax.set_ylim(-row_gap, max(list(ref_y.values()) + list(qry_y.values())) + max(ref_len.values()) + row_gap)
ax.set_title("Linear Ideogram of Synteny Between Old and New Assembly", pad=20, fontsize=14)
plt.tight_layout()
plt.savefig(out_png, dpi=600)
print(f"[✓] Saved ideogram to {out_png}")
plt.show()