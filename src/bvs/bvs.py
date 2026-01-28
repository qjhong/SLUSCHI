#!/usr/bin/env python3
"""
bvs_pymatgen_fixed_unique.py

Compute fractional bond-valence sums (BVS) using pymatgen.
Optional --unique flag: only report one representative per symmetry-equivalent group.
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

# small builtin R0/b table (extend as needed)
BUILTIN_PARAMS = {
    ("Y", "O"): (2.019, 0.37),
    ("P", "O"): (1.617, 0.37),
    ("V", "O"): (1.803, 0.37),
    #("O", "O"): (0.019, 0.37),
    ("O", "O"): (1.406, 0.37),
    ("Gd", "O"): (2.031, 0.37),
    ("Sm", "O"): (2.055, 0.37),
    # add others you need...
    # add others you need...
}

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
    if (elem1, elem2) in BUILTIN_PARAMS:
        return BUILTIN_PARAMS[(elem1, elem2)]
    if (elem2, elem1) in BUILTIN_PARAMS:
        return BUILTIN_PARAMS[(elem2, elem1)]
    return (1.8, 0.37)

def manual_bv_sum_for_site(structure, site_index, cutoff=5.0):
    site = structure[site_index]
    raw_nn = structure.get_neighbors(site, cutoff)
    total = 0.0
    per_bonds = []
    for nsite in raw_nn:
        try:
            if hasattr(nsite, "distance") and callable(getattr(nsite, "distance")):
                n_dist = float(nsite.distance(site))
            else:
                n_dist = float(getattr(nsite, "distance", None))
        except Exception:
            n_dist = None

        best_j = None
        best_diff = 1e9
        for j, s2 in enumerate(structure):
            if j == site_index: continue
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
        if b == 0: b = 0.37
        s = math.exp((R0 - numeric_dist) / b)
        total += s
        per_bonds.append((neighbor_idx, elem_j, numeric_dist, s))
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
                bvs_val = calc_sum(site)
            except TypeError:
                try:
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
        print(f"BVS (fractional) = {bvs_val:.4f}")
        bvs_results.append((get_primary_symbol(site), bvs_val))
        #if bonds:
        #    for nb_idx, nb_elem, dist, s in bonds:
        #        print(f"    -> neighbor {nb_idx:4d} {nb_elem:6s} dist={dist:.4f}  bv={s:.4f}")
    stats = average_std_bvs_by_element(bvs_results)
    print(stats)

if __name__ == "__main__":
    main()
