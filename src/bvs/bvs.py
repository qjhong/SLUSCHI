#!/usr/bin/env python3
"""
bvs_pymatgen_fixed_unique.py

Compute fractional bond-valence sums (BVS) using pymatgen.
Optional --unique flag: only report one representative per symmetry-equivalent group.

This version prints out the R0/B parameters actually used during the main calculation
so you can trace which parameter entries were selected for each pair/bond.
"""
import sys, math, argparse
try:
    import pymatgen
    from pymatgen.core import Structure
    from pymatgen.analysis import bond_valence as bvmod
    # symmetry analyzer (may be missing in very old pymatgen)
    try:
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    except Exception:
        SpacegroupAnalyzer = None
except Exception as e:
    print("ERROR: pymatgen required. Install via conda install -c conda-forge pymatgen")
    raise

import load_bvs_param as lb

import os

def get_sluschi_src():
    """
    Read sluschipath from ~/.sluschi.rc and return absolute SLUSCHI src path.
    """
    rc_path = os.path.expanduser("~/.sluschi.rc")
    if not os.path.exists(rc_path):
        raise FileNotFoundError("~/.sluschi.rc not found")

    sluschipath = None
    with open(rc_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("set sluschipath"):
                # supports sluschipath=/path or sluschipath = /path
                sluschipath = line.split("=", 1)[1].strip()
                break

    if not sluschipath:
        raise RuntimeError("sluschipath not defined in ~/.sluschi.rc")

    # canonical SLUSCHI src directory
    src_path = sluschipath
    if not os.path.isdir(src_path):
        raise RuntimeError(f"SLUSCHI src directory not found: {src_path}")

    return src_path

with open(get_sluschi_src() + '/bvs/bvparam2020.cif', 'r', encoding='utf-8') as fh:
    cif_text = fh.read()

# Load builtin parameters (may detect requested pairs from the CIF text)
BUILTIN_PARAMS = lb.get_bvs_params(cif_text, pairs=None, format='dict')

# Print the loaded parameter entries (these are the base parameters available)
if BUILTIN_PARAMS:
    print("\nLoaded BUILTIN_PARAMS entries from bvparam file:")
    #for (A, B), (Ro_val, Bval) in sorted(BUILTIN_PARAMS.items(), key=lambda x: (x[0][0], x[0][1])):
    #    print(f"  Pair {A}-{B}: Ro={Ro_val}, B={Bval}")
else:
    print("  (no parameters loaded)")

def get_primary_symbol(site):
    """Return the primary element symbol for a (Periodic)Site in a robust way."""
    if hasattr(site, "specie") and site.specie is not None:
        try:
            return site.specie.symbol
        except Exception:
            pass
    if hasattr(site, "species"):
        try:
            sp = list(site.species)
            if len(sp) > 0:
                try:
                    return sp[0].symbol
                except Exception:
                    return str(sp[0])
        except Exception:
            pass
    if hasattr(site, "species_string"):
        try:
            s = site.species_string
            return s.split()[0].split(":")[0].split("0")[0] or s
        except Exception:
            pass
    return str(site).split()[0]

def choose_R0_b(elem1, elem2):
    """
    Choose R0 and b for the pair (elem1, elem2). Also print out the decision so the
    main calculation can show which parameters were actually used.
    """
    # prefer exact ordering (elem1, elem2) then unordered
    if (elem1, elem2) in BUILTIN_PARAMS:
        R0, b = BUILTIN_PARAMS[(elem1, elem2)]
        print(f"PARAM SELECT: ({elem1},{elem2}) exact -> R0={R0}, b={b}")
        return (R0, b)
    if (elem2, elem1) in BUILTIN_PARAMS:
        R0, b = BUILTIN_PARAMS[(elem2, elem1)]
        print(f"PARAM SELECT: ({elem1},{elem2}) via ({elem2},{elem1}) -> R0={R0}, b={b}")
        return (R0, b)
    # fallback default
    print(f"PARAM SELECT: ({elem1},{elem2}) not found in BUILTIN_PARAMS -> using fallback R0=1.8, b=0.37")
    return (1.8, 0.37)

def manual_bv_sum_for_site(structure, site_index, cutoff=5.0):
    """
    Compute BVS manually and return total and per-bond details.

    per_bonds entries are tuples:
      (neighbor_index, neighbor_element, numeric_dist, s_bond, R0_used, b_used)
    """
    site = structure[site_index]
    raw_nn = structure.get_neighbors(site, cutoff)
    total = 0.0
    per_bonds = []
    for nsite in raw_nn:
        try:
            # pymatgen's Neighbor distance is an attribute, sometimes a method depending on version
            if hasattr(nsite, "distance") and callable(getattr(nsite, "distance")):
                n_dist = float(nsite.distance(site))
            else:
                n_dist = float(getattr(nsite, "distance", None))
        except Exception:
            n_dist = None

        best_j = None
        best_diff = 1e9
        for j, s2 in enumerate(structure):
            if j == site_index:
                continue
            try:
                d_ij = structure.get_distance(site_index, j)
            except Exception:
                d_ij = math.dist(site.coords, s2.coords)
            if n_dist is None:
                # prefer same element
                try:
                    s2_ss = s2.species_string
                except Exception:
                    s2_ss = get_primary_symbol(s2)
                nsite_ss = getattr(nsite, "species_string", None) or get_primary_symbol(nsite)
                if nsite_ss and s2_ss == nsite_ss:
                    if d_ij < best_diff:
                        best_diff = d_ij; best_j = j
                else:
                    if d_ij < best_diff:
                        best_diff = d_ij; best_j = j
            else:
                diff = abs(d_ij - n_dist)
                if diff < best_diff and s2.species_string == getattr(nsite, "species_string", s2.species_string):
                    best_diff = diff; best_j = j
        if best_j is None:
            continue
        neighbor_idx = best_j
        try:
            numeric_dist = structure.get_distance(site_index, neighbor_idx)
        except Exception:
            numeric_dist = n_dist if n_dist is not None else 0.0
        elem_i = get_primary_symbol(site)
        elem_j = get_primary_symbol(structure[neighbor_idx])
        R0, b = choose_R0_b(elem_i, elem_j)
        if b == 0:
            b = 0.37
        s = math.exp((R0 - numeric_dist) / b)
        total += s
        per_bonds.append((neighbor_idx, elem_j, numeric_dist, s, R0, b))
    return total, per_bonds

def average_std_bvs_by_element(bvs_results):
    """
    Compute average and standard deviation of BVS grouped by element.

    Parameters
    ----------
    bvs_results : list of tuples
        Each entry is (element_symbol: str, bvs_value: float)

    Returns
    -------
    dict
        { element_symbol: (mean_bvs, std_bvs) }
    """
    from collections import defaultdict
    import math

    acc = defaultdict(list)
    for elem, val in bvs_results:
        if val is not None and not math.isnan(val):
            acc[elem].append(val)

    stats = {}
    for elem, vals in acc.items():
        n = len(vals)
        mean = sum(vals) / n
        var = sum((v - mean) ** 2 for v in vals) / n
        std = math.sqrt(var)
        stats[elem] = (mean, std)

    return stats

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("cif")
    parser.add_argument("--cutoff", type=float, default=3.2)
    parser.add_argument("--unique", action="store_true",
                        help="Only compute/report one representative per symmetry-equivalent group")
    args = parser.parse_args()

    print("pymatgen version:", getattr(pymatgen, "__version__", "unknown"))
    struct = Structure.from_file(args.cif)
    print("structure formula:", struct.formula)
    print("lattice lengths (Å):", struct.lattice.abc)
    for i, site in enumerate(struct[:6]):
        print(f"  site {i}: {get_primary_symbol(site):6s}  frac: {site.frac_coords}")

    # determine which indices to process
    indices_to_process = list(range(len(struct)))
    if args.unique:
        if SpacegroupAnalyzer is None:
            print("Warning: pymatgen.symmetry.analyzer.SpacegroupAnalyzer not available; cannot compute unique sites.")
        else:
            try:
                sga = SpacegroupAnalyzer(struct, symprec=1e-3)  # adjust symprec if needed
                symm_struct = sga.get_symmetrized_structure()
                eq_groups = symm_struct.equivalent_indices  # list of lists of indices
                # pick first index of each group as representative
                indices_to_process = [group[0] for group in eq_groups]
                indices_to_process.sort()
                print(f"Detected {len(eq_groups)} unique symmetry groups; processing indices: {indices_to_process}")
            except Exception as e:
                print("Warning: symmetry analysis failed:", e)
                indices_to_process = list(range(len(struct)))

    bvs_results = []
    # Try using calculate_bv_sum(site) if available
    calc_sum = getattr(bvmod, "calculate_bv_sum", None)
    print("\n--- computing BVS ---")
    for i in indices_to_process:
        site = struct[i]
        print(f"\n--- site {i} {get_primary_symbol(site):3s}   frac: {site.frac_coords} ---")
        bvs_val = None
        bonds = None
        if calc_sum is not None:
            try:
                # try site-only API
                bvs_val = calc_sum(site)
            except TypeError:
                try:
                    # try structure,index API
                    bvs_val = calc_sum(struct, i)
                except Exception:
                    bvs_val = None
            except Exception:
                bvs_val = None

        if bvs_val is None:
            try:
                bvs_val, bonds = manual_bv_sum_for_site(struct, i, cutoff=args.cutoff)
            except Exception as e:
                print("manual calculation failed:", e)
                bvs_val = float("nan")

        # Print the BVS and the per-bond parameter choices if available
        print(f"BVS (fractional) = {bvs_val:.4f}")
        if bonds:
            print("  Per-bond details (neighbor_idx, elem, dist(Å), bond_valence, R0_used, b_used):")
            for nb_idx, nb_elem, dist, s_bond, R0_used, b_used in bonds:
                print(f"    -> neighbor {nb_idx:4d} {nb_elem:6s} dist={dist:.4f}  bv={s_bond:.4f}  R0={R0_used:.4f} b={b_used:.4f}")

        bvs_results.append((get_primary_symbol(site), bvs_val))

    stats = average_std_bvs_by_element(bvs_results)
    print("\nBVS statistics by element (mean, std):")
    print(stats)

if __name__ == "__main__":
    main()
