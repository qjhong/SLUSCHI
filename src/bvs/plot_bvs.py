# mirrored_bvs_plot.py
import re
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict, OrderedDict

# ---- User: set input filename here (paste your long site listing into this file) ----
INPUT_FILE = "bvs.txt"

# ---- Parsing regex ----
# expects lines like: --- site 0 Sm    frac: [0.876086 0.881428 0.253157] ---
# and a line containing: BVS (fractional) = 2.5887
site_line_re = re.compile(r'---\s*site\s*(\d+)\s+(\w+)', re.IGNORECASE)
bvs_re = re.compile(r'BVS \(fractional\)\s*=\s*([0-9.+-eE]+)')

# ---- Read and parse file ----
records = []  # list of tuples (site_index:int, element:str, bvs:float)
with open(INPUT_FILE, 'r', encoding='utf8') as f:
    content = f.read().splitlines()

current_site = None
current_elem = None
for ln in content:
    m = site_line_re.search(ln)
    if m:
        current_site = int(m.group(1))
        current_elem = m.group(2)
        # continue to next lines to find BVS on a later line (some inputs have it on same line)
        # but we'll continue the loop to locate bvs_re
        continue
    m2 = bvs_re.search(ln)
    if m2 and current_site is not None:
        bvs_val = float(m2.group(1))
        records.append((current_site, current_elem, bvs_val))
        current_site = None
        current_elem = None

if len(records) == 0:
    raise SystemExit("No site/BVS records found. Make sure bvs_sm.txt contains the block you pasted.")

# ---- Collapse symmetry pairs (site 0+1 -> pair 0, 2+3 -> pair 1, etc.) ----
pairs = OrderedDict()  # key: pair_index, value: dict { 'sites':[i,j], 'elem':elem, 'bvs':[v_i,v_j] }
for site_idx, elem, bvs in records:
    pair_idx = site_idx // 2
    if pair_idx not in pairs:
        pairs[pair_idx] = {'sites': [], 'elem': elem, 'bvs': []}
    pairs[pair_idx]['sites'].append(site_idx)
    pairs[pair_idx]['bvs'].append(bvs)
    # store element name from the first site in the pair (assumes both in a pair share element family)
    # (If pair members can have different element tags, you can change this behavior.)
    # Keep the element as the one from the lowest site index encountered:
    if site_idx < min(pairs[pair_idx]['sites']):
        pairs[pair_idx]['elem'] = elem

# average the BVS values per pair
pair_indices = []
pair_bvs = []
pair_elem = []
for pair_idx, info in pairs.items():
    avg_bvs = float(np.mean(info['bvs']))
    pair_indices.append(pair_idx)
    pair_bvs.append(avg_bvs)
    pair_elem.append(info['elem'])

pair_indices = np.array(pair_indices)
pair_bvs = np.array(pair_bvs)
pair_elem = np.array(pair_elem)

# ---- Prepare plotting order & element block locations (to label Sm / P / O) ----
# We'll keep the pair order as-is (0 .. n-1).
n = len(pair_indices)
y_pos = np.arange(n)  # y positions (will plot from top->bottom after inversion)

# Find contiguous blocks of same element for big right-side labels
blocks = []
start = 0
for i in range(1, n+1):
    if i == n or pair_elem[i] != pair_elem[i-1]:
        blocks.append((pair_elem[i-1], start, i-1))
        start = i

# ---- Plot ----
plt.style.use('default')
fig, ax = plt.subplots(figsize=(6.5, 12))  # tall figure like your reference

# mirrored bars:
# left = negative (mirror), right = positive
left_vals = -pair_bvs
right_vals = pair_bvs

bar_height = 0.8

ax.barh(y_pos, left_vals, height=bar_height, align='center', color='#1f77b4', edgecolor='k')   # blue
ax.barh(y_pos, right_vals, height=bar_height, align='center', color='#ff7f0e', edgecolor='k')  # orange

# y-axis labels: use the pair index but you can replace with site labels if you prefer
ax.set_yticks(y_pos)
ax.set_yticklabels([str(int(i)) for i in pair_indices], fontsize=8)

# Invert y to have site 0 at top (matching your reference)
ax.invert_yaxis()

# Vertical dashed lines at nominal valences (both directions)
nominal = {'Sm': 3.0, 'P': 5.0, 'O': 2.0}
# Draw lines symmetric about 0
max_bvs = max(np.max(pair_bvs), max(nominal.values())*1.2)
xlim = max_bvs * 1.15
ax.set_xlim(-xlim, xlim)

for name, val in nominal.items():
    ax.axvline(val, color='red', linestyle='--', linewidth=1)
    ax.axvline(-val, color='red', linestyle='--', linewidth=1)

# central spine / labels
ax.set_xlabel("BVS (mirrored)")
ax.set_title("Bond-Valence Sums — unique sites (collapsed symmetry pairs)", fontsize=12)

# Add big element labels on the right side at the center of each block
for elem_name, start_idx, end_idx in blocks:
    y_center = (start_idx + end_idx) / 2.0
    ax.text(0.99, y_center, elem_name, transform=ax.get_yaxis_transform(), 
            ha='left', va='center', fontsize=16, fontweight='bold',
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.0),
            clip_on=False)

# Tweak layout
ax.grid(False)
ax.set_axisbelow(True)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# Add small ticks at nominal valence positions to echo reference's red ticks
ax.set_xticks(list(ax.get_xticks()) + [nominal['Sm'], nominal['P'], nominal['O'], -nominal['Sm'], -nominal['P'], -nominal['O']])
ax.tick_params(axis='x', labelsize=10)

plt.tight_layout()
plt.show()

# ---- Save figure (optional) ----
fig.savefig("bvs.png", dpi=300, bbox_inches='tight')
#fig.savefig("bvs_mirrored_from_text.svg", format='svg')
print("Saved bvs_mirrored_from_text.png")
