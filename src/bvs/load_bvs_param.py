#!/usr/bin/env python3
"""
load_bvs_param.py

Load BVS parameters from a CIF-like text file.

Features / behavior:
 - More permissive candidate-line detection (accepts '?' '+' etc).
 - Robust parsing for space-separated "bvparam" files like the example:
     1 Ac 3   O  -2    2.24     0.37     b   ?
 - Attempts three parsing strategies for each candidate line:
     A) indexed rows: idx elA oxA elB oxB Ro B ...
     B) unindexed but same column order: elA oxA elB oxB Ro B ...
     C) fallback regex search for el1 el2 ... Ro ... B
 - If pairs is None:
     1) try to detect explicitly requested pairs inside CIF (various heuristics)
     2) otherwise build Cartesian product from element symbols found in CIF
 - Returns either a dict {(A,B):(Ro,B)} or a builtin_list [(A,B,Ro,B), ...]
 - CLI has --debug to print candidate lines and parsed pool diagnostics.

Usage:
    python3 load_bvs_param.py bvparam2020.cif --debug
"""

from __future__ import annotations
import re
import statistics
import argparse
from typing import List, Tuple, Dict, Optional

Pair = Tuple[str, str]

# --- Helper regexes ---
# regex to capture 'element' tokens possibly with trailing digits (e.g. O, O1)
_el_token_re = re.compile(r'^[A-Z][a-z]?\d*$')
# plain element symbol without digits
_el_strip_digits_re = re.compile(r'^([A-Z][a-z]?)')
# fallback regex to find el1, el2, Ro, B in arbitrary order
_ro_b_re = re.compile(r"""
    (?P<el1>[A-Z][a-z]?)      # first element symbol
    \s*[,\-/_]?\s*
    (?P<el2>[A-Z][a-z]?)      # second element symbol
    .*?
    (?P<Ro>-?\d+(?:\.\d+)?)   # Ro value (int or float, possibly negative)
    .*?
    (?P<B>\d+(?:\.\d+)?)      # B value (positive float)
    """, re.VERBOSE)

def _strip_trailing_digits(tok: str) -> str:
    """Strip trailing digits from tokens like 'O1' -> 'O'"""
    m = _el_strip_digits_re.match(tok)
    return m.group(1) if m else tok

def _find_bvparm_lines(text: str) -> List[str]:
    """Return lines likely containing BVPARM-style entries.

    More permissive than older version: accepts '?', '+', '*' and looks for at least
    two element-like tokens anywhere in the line (not necessarily adjacent).
    """
    lines = []
    # allow ?, + and some other punctuation found in parameter files
    _line_candidate_re = re.compile(r'^[\s\w\-\.,:;()\/\?\+\*\%]+$')
    # accept two element-like tokens anywhere in the line (not necessarily adjacent)
    el_pair_anywhere_re = re.compile(r'([A-Z][a-z]?)')
    for raw in text.splitlines():
        s = raw.strip()
        if not s:
            continue
        # quick sanity: line must be composed of common table characters
        if not _line_candidate_re.match(s):
            continue
        # require at least two element-like tokens somewhere in the line
        els = re.findall(el_pair_anywhere_re, s)
        if len(els) < 2:
            continue
        # require at least two numeric-looking tokens (covers Ro and B or charges)
        if not (re.search(r'-?\d+(?:\.\d+)?', s) and re.search(r'\d+(?:\.\d+)?', s)):
            continue
        lines.append(s)
    return lines

def _extract_pairs_from_text(text: str) -> List[str]:
    """Conservative extraction of element symbols present in the CIF text."""
    tokens = re.findall(r'\b([A-Z][a-z]?)\b', text)
    seen = []
    for t in tokens:
        if t not in seen:
            seen.append(t)
    return seen

# --- Requested-pairs detection heuristics ---
def _split_elements_from_token(token: str) -> List[str]:
    """Given a token like 'Y,O' or 'Y-O' or 'Y O' return element tokens ['Y','O'] if valid."""
    sep_token = re.sub(r'[\s;/]+', ',', token.strip())
    sep_token = sep_token.replace('-', ',').replace('/', ',')
    parts = [p.strip() for p in sep_token.split(',') if p.strip()]
    el_re = re.compile(r'^[A-Z][a-z]?$')
    parts = [p for p in parts if el_re.match(p)]
    return parts

def _parse_loop_block(lines: List[str], start_idx: int) -> (List[str], List[str], int):
    """
    Parse a CIF 'loop_' block starting at start_idx.
    Returns (labels, data_rows, next_index).
    """
    labels = []
    data_rows = []
    i = start_idx + 1
    # collect labels
    while i < len(lines):
        ln = lines[i].strip()
        if ln.startswith('_'):
            labels.append(ln)
            i += 1
            continue
        break
    # collect data rows
    while i < len(lines):
        ln = lines[i].strip()
        if not ln:
            break
        if ln.lower().startswith('loop_') or ln.startswith('data_') or ln.startswith('_'):
            break
        data_rows.append(ln)
        i += 1
    return (labels, data_rows, i)

def _find_requested_pairs(text: str) -> List[Pair]:
    """
    Heuristics to detect explicitly requested pairs in CIF text.
    Returns ordered list of (A,B) pairs if found, else [].
    """
    lines = text.splitlines()
    found_pairs: List[Pair] = []

    # 1) tags like _bvs_requested_pairs: Y,O; P,O  or _requested_pairs 'Y O','P O'
    tag_re = re.compile(r'^_(?P<tag>[A-Za-z0-9_]*requested[_A-Za-z0-9]*)\s+(?P<val>.+)$', re.IGNORECASE)
    for raw in lines:
        m = tag_re.match(raw.strip())
        if m:
            val = m.group('val').strip().strip("'\"")
            for chunk in re.split(r'[;]+', val):
                elems = _split_elements_from_token(chunk)
                if len(elems) >= 2:
                    found_pairs.append((elems[0], elems[1]))
            if found_pairs:
                return found_pairs

    # 2) loop_ blocks: look for loops where column labels contain 'pair', 'element', 'elem', 'species', 'type'
    label_match_re = re.compile(r'pair|elem|element|species|type', re.IGNORECASE)
    i = 0
    while i < len(lines):
        ln = lines[i].strip()
        if ln.lower().startswith('loop_'):
            labels, data_rows, next_i = _parse_loop_block(lines, i)
            if any(label_match_re.search(lbl) for lbl in labels):
                for row in data_rows:
                    tokens = re.split(r'[\s,;]+', row.strip())
                    el_re = re.compile(r'^[A-Z][a-z]?$')
                    els = [t for t in tokens if el_re.match(t)]
                    if len(els) >= 2:
                        found_pairs.append((els[0], els[1]))
                    else:
                        quoted = re.findall(r"'([^']+)'|\"([^\"]+)\"", row)
                        for q in quoted:
                            qv = q[0] or q[1]
                            elems = _split_elements_from_token(qv)
                            if len(elems) >= 2:
                                found_pairs.append((elems[0], elems[1]))
                if found_pairs:
                    return found_pairs
            i = next_i
        else:
            i += 1

    # 3) free-form lines like "requested pairs: Y O; P O" or "pairs: Y,O P,O"
    free_re = re.compile(r'(requested\s+pairs|pairs)\s*[:=]\s*(.+)', re.IGNORECASE)
    for raw in lines:
        m = free_re.search(raw)
        if m:
            val = m.group(2)
            for chunk in re.split(r'[;,]+', val):
                elems = _split_elements_from_token(chunk)
                if len(elems) >= 2:
                    found_pairs.append((elems[0], elems[1]))
            if found_pairs:
                return found_pairs

    return []

# --- Main loader function ---
def get_bvs_params(cif_text: str,
                   pairs: Optional[List[Pair]] = None,
                   format: str = 'dict'
                   ) -> Dict[Pair, Tuple[float, float]] or List[Tuple[str,str,float,float]]:
    """
    Load BVS parameters from text.

    If pairs is None:
      - try to detect explicit requested pairs in CIF
      - otherwise form full Cartesian pairs from elements found in CIF
    """
    candidate_lines = _find_bvparm_lines(cif_text)

    # allow a global DEBUG knob (set by CLI below)
    debug = globals().get('DEBUG', False)
    if debug:
        print(f"DEBUG: {len(candidate_lines)} candidate bvparm lines found")
        for li in candidate_lines[:200]:
            print("  ", li)

    # parse candidate lines into a pool mapping unordered key -> list of (Ro,B)
    pool: Dict[Tuple[str,str], List[Tuple[float,float]]] = {}

    for ln in candidate_lines:
        tokens = ln.split()
        parsed = False

        # Case A: indexed rows like '1 Ac 3 O -2 2.24 0.37 ...'
        if len(tokens) >= 7 and tokens[0].isdigit():
            cand_elA = tokens[1]
            cand_elB = tokens[3]
            if _el_token_re.match(cand_elA) and _el_token_re.match(cand_elB):
                try:
                    Ro = float(tokens[5])
                    Bval = float(tokens[6])
                    e1 = _strip_trailing_digits(cand_elA)
                    e2 = _strip_trailing_digits(cand_elB)
                    key_u = tuple(sorted((e1, e2)))
                    pool.setdefault(key_u, []).append((Ro, Bval))
                    parsed = True
                except Exception:
                    parsed = False

        # Case B: unindexed rows 'Ac 3 O -2 2.24 0.37 ...' or similar
        if not parsed and len(tokens) >= 6:
            for i in range(len(tokens)-5):
                # tokens[i] = elA, tokens[i+2] = elB, tokens[i+3] = charge, tokens[i+4] = Ro, tokens[i+5] = B
                if _el_token_re.match(tokens[i]) and re.match(r'^[\d\-]+$', tokens[i+1]) and _el_token_re.match(tokens[i+2]):
                    try:
                        Ro = float(tokens[i+4])
                        Bval = float(tokens[i+5])
                        e1 = _strip_trailing_digits(tokens[i])
                        e2 = _strip_trailing_digits(tokens[i+2])
                        key_u = tuple(sorted((e1, e2)))
                        pool.setdefault(key_u, []).append((Ro, Bval))
                        parsed = True
                        break
                    except Exception:
                        continue

        # Case C: fallback to regex-based search
        if not parsed:
            m = _ro_b_re.search(ln)
            if m:
                try:
                    e1 = m.group('el1')
                    e2 = m.group('el2')
                    Ro = float(m.group('Ro'))
                    Bval = float(m.group('B'))
                    key_u = tuple(sorted((_strip_trailing_digits(e1), _strip_trailing_digits(e2))))
                    pool.setdefault(key_u, []).append((Ro, Bval))
                    parsed = True
                except Exception:
                    parsed = False

        if not parsed and debug:
            print("DEBUG: failed to parse candidate line:", ln)

    if debug:
        print("DEBUG: parsed pool keys and counts:")
        for k,v in pool.items():
            print("  ", k, "->", len(v), "entries")

    # If pairs not provided, try detecting requested pairs in CIF else form cartesian product
    if pairs is None:
        requested = _find_requested_pairs(cif_text)
        if requested:
            pairs = requested
            if debug:
                print("DEBUG: using requested pairs found in CIF:", pairs)
        else:
            elements = _extract_pairs_from_text(cif_text)
            # ensure we strip trailing digits from any element-like tokens discovered
            elements = [_strip_trailing_digits(e) for e in elements]
            # unique while preserving order
            uniq = []
            for e in elements:
                if e not in uniq:
                    uniq.append(e)
            pairs = []
            for a in uniq:
                for b in uniq:
                    pairs.append((a, b))
            if debug:
                print("DEBUG: falling back to Cartesian pairs for elements:", uniq)

    results: Dict[Pair, Tuple[float, float]] = {}
    for (A, B) in pairs:
        A_s = _strip_trailing_digits(A)
        B_s = _strip_trailing_digits(B)
        key_u = tuple(sorted((A_s, B_s)))
        values = pool.get(key_u, [])
        if not values:
            print(f"Warning: no BVPARM entry found for pair ({A},{B}) in CIF.")
            continue
        Ro_list = [v[0] for v in values]
        B_list = [v[1] for v in values]
        Ro_med = float(statistics.median(Ro_list))
        B_med = float(statistics.median(B_list))
        results[(A, B)] = (Ro_med, B_med)

    if format == 'dict':
        return results
    elif format == 'builtin_list':
        return [ (A,B,results[(A,B)][0], results[(A,B)][1]) for (A,B) in results ]
    else:
        raise ValueError("format must be 'dict' or 'builtin_list'")

# Convenience file loader
def load_bvs_params_from_file(path: str, pairs: Optional[List[Pair]] = None, format: str = 'dict'):
    with open(path, 'r', encoding='utf-8') as fh:
        text = fh.read()
    return get_bvs_params(text, pairs=pairs, format=format)

# ----- CLI -----
if __name__ == '__main__':
    p = argparse.ArgumentParser(description="Load BV parameters from a CIF-like text file")
    p.add_argument('cif', help='path to CIF/text file')
    p.add_argument('--format', choices=['dict','builtin_list'], default='dict')
    p.add_argument('--pairs', nargs='*', help='optional explicit pairs A,B A,B ... (e.g. Y,O P,O)')
    p.add_argument('--debug', action='store_true', help='print debug info about candidate lines and parsed pool')
    args = p.parse_args()

    # set module-global DEBUG for functions to read
    globals()['DEBUG'] = args.debug

    pairs = None
    if args.pairs:
        pairs = []
        for item in args.pairs:
            if ',' in item:
                a,b = item.split(',',1)
            elif '/' in item:
                a,b = item.split('/',1)
            else:
                parts = item.split()
                if len(parts) >= 2:
                    a,b = parts[0], parts[1]
                else:
                    raise ValueError("Bad pair format in --pairs argument.")
            pairs.append((a.strip(), b.strip()))

    out = load_bvs_params_from_file(args.cif, pairs=pairs, format=args.format)
    if args.format == 'dict':
        print("Loaded BUILTIN_PARAMS entries:")
        for k,v in out.items():
            print(f"  {k}: Ro={v[0]}, B={v[1]}")
    else:
        print(out)
