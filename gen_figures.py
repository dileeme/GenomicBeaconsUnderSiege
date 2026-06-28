import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os

OUT = "figures"

# Palette
MUTED   = "#7B9EB0"   # steel blue — default bars
ACCENT  = "#2A6496"   # dark blue — Config B highlight
TEXT    = "#2B2B2B"
BG      = "white"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.titlesize": 10,
    "axes.titleweight": "bold",
    "figure.facecolor": BG,
    "axes.facecolor": BG,
})

# ── 1. fhe_performance_overhead.png ──────────────────────────────────────────
fig, ax = plt.subplots(figsize=(3.4, 2.8))

labels  = ["Plaintext\nQuery", "Encryption", "FHE\nComputation"]
values  = [0.003, 11.3, 76.5]
colors  = [MUTED, MUTED, MUTED]

bars = ax.bar(labels, values, color=colors, width=0.5, edgecolor="none")

ax.set_yscale("log")
ax.set_ylabel("Time (seconds)", color=TEXT)
ax.set_title("FHE Performance Overhead")
ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(
    lambda x, _: f"{x:g}s"
))

# Value labels
for bar, val in zip(bars, values):
    label = f"{val:.3f}s" if val < 0.01 else f"{val:.1f}s"
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() * 1.4,
        label,
        ha="center", va="bottom", fontsize=8, color=TEXT, fontweight="bold"
    )

# Annotation — placed in lower-right, below the bars
ax.annotate(
    "29,268× end-to-end\n(~87.8s total)",
    xy=(1.0, 0.0), xycoords="axes fraction",
    xytext=(-6, 6), textcoords="offset points",
    ha="right", va="bottom", fontsize=7.5,
    color="#C0392B", fontweight="bold",
    bbox=dict(boxstyle="round,pad=0.3", fc="#FDECEA", ec="#C0392B", lw=0.8)
)

ax.set_ylim(bottom=0.001)
fig.tight_layout()
path1 = os.path.join(OUT, "fhe_performance_overhead.png")
fig.savefig(path1, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {path1}")

# ── 2. latency_optimization.png ──────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(3.4, 2.8))

configs  = ["Config A", "Config B\n(Recommended)", "Config C", "Config D"]
latency  = [75.6, 38.6, 75.0, 39.3]
bar_cols = [MUTED, ACCENT, MUTED, MUTED]

bars = ax.bar(configs, latency, color=bar_cols, width=0.55, edgecolor="none")

ax.set_ylabel("Latency (seconds)", color=TEXT)
ax.set_title("FHE Latency by Parameter Config")
ax.set_ylim(0, max(latency) * 1.25)

for bar, val in zip(bars, latency):
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() + 0.8,
        f"{val:.1f}s",
        ha="center", va="bottom", fontsize=8, color=TEXT, fontweight="bold"
    )

legend_patch = mpatches.Patch(color=ACCENT, label="Config B (Recommended)")
ax.legend(handles=[legend_patch], fontsize=7.5, loc="upper right",
          frameon=False)

fig.tight_layout()
path2 = os.path.join(OUT, "latency_optimization.png")
fig.savefig(path2, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {path2}")

# ── 3. storage_optimization.png ──────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(3.4, 2.8))

storage  = [326.6, 229.8, 238.0, 155.6]
bar_cols = [MUTED, ACCENT, MUTED, MUTED]

bars = ax.bar(configs, storage, color=bar_cols, width=0.55, edgecolor="none")

ax.set_ylabel("Storage (KB)", color=TEXT)
ax.set_title("FHE Storage by Parameter Config")
ax.set_ylim(0, max(storage) * 1.25)

for bar, val in zip(bars, storage):
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() + 3,
        f"{val:.1f} KB",
        ha="center", va="bottom", fontsize=8, color=TEXT, fontweight="bold"
    )

legend_patch = mpatches.Patch(color=ACCENT, label="Config B (Recommended)")
ax.legend(handles=[legend_patch], fontsize=7.5, loc="upper right",
          frameon=False)

fig.tight_layout()
path3 = os.path.join(OUT, "storage_optimization.png")
fig.savefig(path3, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {path3}")
