import json
import os

import pandas as pd
import yaml
from outbreak_data import authenticate_user, outbreak_data as od
from outbreak_tools import crumbs

OUTPUT_DIR = 'validated'
CACHE_FILE = '.api_cache.json'

#authenticate_user.authenticate_new_user()
lineage_key = crumbs.get_alias_key()

# ── barcode + lineage hierarchy (local, no API) ───────────────────────────────

_barcodes = pd.read_feather('usher_barcodes.feather').set_index('index')

with open('lineages.yml') as f:
    _lineage_list = yaml.safe_load(f)

_parent_of: dict = {}
_children_of: dict = {}
for _e in _lineage_list:
    _n = _e['name']
    _p = _e.get('parent')
    _parent_of[_n] = _p
    _children_of.setdefault(_n, [])
    if _p:
        _children_of.setdefault(_p, []).append(_n)

# ── per-run caches ────────────────────────────────────────────────────────────
_ldm_cache: dict = {}
_lineage_prev_cache: dict = {}
_pm_cache: dict = {}
_pml_cache: dict = {}


def _load_caches() -> None:
    if not os.path.isfile(CACHE_FILE):
        return
    with open(CACHE_FILE) as f:
        data = json.load(f)
    for k, v in data.get('lineage_prev', {}).items():
        _lineage_prev_cache[tuple(json.loads(k))] = tuple(v)
    for k, v in data.get('pm', {}).items():
        parts = json.loads(k)
        _pm_cache[(frozenset(parts[0]), parts[1], parts[2])] = v
    for k, v in data.get('pml', {}).items():
        parts = json.loads(k)
        _pml_cache[(frozenset(parts[0]), parts[1], parts[2], parts[3])] = v
    print(f'Loaded cache from {CACHE_FILE}')


def _save_caches() -> None:
    data = {
        'lineage_prev': {json.dumps(list(k)): list(v) for k, v in _lineage_prev_cache.items()},
        'pm': {json.dumps([sorted(k[0]), k[1], k[2]]): v for k, v in _pm_cache.items()},
        'pml': {json.dumps([sorted(k[0]), k[1], k[2], k[3]]): v for k, v in _pml_cache.items()},
    }
    with open(CACHE_FILE, 'w') as f:
        json.dump(data, f)
    print(f'Saved cache to {CACHE_FILE}')


_load_caches()


# ── LDM helpers (barcode-based) ───────────────────────────────────────────────

def to_pango(lineage: str) -> str:
    return lineage.removesuffix('.X')


def _get_all_descendants(pango_lin: str) -> set:
    result, queue = set(), [pango_lin]
    while queue:
        node = queue.pop()
        result.add(node)
        queue.extend(_children_of.get(node, []))
    return result


def _barcode_muts(lineage: str) -> set:
    if lineage not in _barcodes.index:
        return set()
    row = _barcodes.loc[lineage]
    return set(row.index[row == 1.0])


def get_ldms(pango_lin: str) -> set:
    """Union of barcode mutations across pango_lin and all descendants, minus parent mutations."""
    if pango_lin in _ldm_cache:
        return _ldm_cache[pango_lin]

    descendants = _get_all_descendants(pango_lin)
    muts = set().union(*(_barcode_muts(d) for d in descendants))

    parent = _parent_of.get(pango_lin)
    if parent:
        muts -= _barcode_muts(parent)

    _ldm_cache[pango_lin] = muts
    return muts


# ── API helpers ───────────────────────────────────────────────────────────────

def get_lineage_prevalence(pango_lin: str, datemin: str, datemax: str) -> tuple:
    """Returns (p_lineage, lineage_count_sum) for Bayes computation."""
    key = (pango_lin, datemin, datemax)
    if key not in _lineage_prev_cache:
        df = od.lineage_cl_prevalence(
            pango_lin, descendants=True, location='USA',
            datemin=datemin, datemax=datemax, lineage_key=lineage_key,
        )
        lc = float(df['lineage_count'].sum())
        tc = float(df['total_count'].sum())
        _lineage_prev_cache[key] = (lc / tc if tc > 0 else 0.0, lc)
    return _lineage_prev_cache[key]


def get_p_mutations(muts: list, datemin: str, datemax: str) -> float:
    """P(M) — overall prevalence of this mutation set across all lineages."""
    key = (frozenset(muts), datemin, datemax)
    if key not in _pm_cache:
        return 0.0 
        # df = od.lineage_cl_prevalence(
        #     '.', descendants=True, mutations=muts, location='USA',
        #     datemin=datemin, datemax=datemax, lineage_key=lineage_key,
        # )
        # tc = float(df['total_count'].sum())
        # _pm_cache[key] = float(df['lineage_count'].sum()) / tc if tc > 0 else 0.0
    return _pm_cache[key]


def _multiquery_to_df(data: dict) -> pd.DataFrame:
    return pd.concat(
        [pd.DataFrame(v).assign(query=k) for k, v in data['results'].items()], axis=0
    )


def get_p_mutations_given_lineage(
    muts: list, pango_lin: str, lineage_count: float, datemin: str, datemax: str,
) -> float:
    """P(M|L) — fraction of lineage-L sequences that carry this mutation set."""
    key = (frozenset(muts), pango_lin, datemin, datemax)
    if key not in _pml_cache:
        fmt = ' AND '.join(f"mutations:{m.replace(':', '?')}" for m in muts)
        url = (
            'https://api.outbreak.info/genomics/prevalence-by-location'
            f'?lineages=None&q=pangolin_lineage_crumbs:*;{pango_lin};* AND {fmt}'
            f'&cumulative=false&min_date={datemin}&max_date={datemax}&location_id=USA'
        )
        resp = od.requests.get(url, headers=od._get_user_authentication())
        df = _multiquery_to_df(resp.json())
        _pml_cache[key] = float(df['lineage_count'].sum()) / lineage_count if lineage_count > 0 else 0.0
    return _pml_cache[key]


# ── covar cluster utilities ───────────────────────────────────────────────────

def parse_nt_muts(nt_mutations_val) -> list:
    """Parse nt_mutations cell from covar TSV into a deduplicated list."""
    return list(dict.fromkeys(
        m for m in str(nt_mutations_val).split() if m and m != 'nan'
    ))


def parse_aa_muts(aa_mutations_val) -> list:
    """Parse aa_mutations cell from covar TSV into a deduplicated list."""
    return list(dict.fromkeys(m.lower() for m in str(aa_mutations_val).split() if ':' in m))


def count_lineage_clusters(covariants: pd.DataFrame, ldm_set: set) -> dict:
    counts = {
        'n_clusters_1_ldm': 0,
        'cluster_depth_1_ldm': 0,
        'n_clusters_2plus_ldm': 0,
        'cluster_depth_2plus_ldm': 0,
    }
    for _, row in covariants.iterrows():
        nt_muts = parse_nt_muts(row.get('nt_mutations', ''))
        n_ldm = sum(1 for m in nt_muts if m in ldm_set)
        if n_ldm == 0:
            continue
        depth = int(row['cluster_depth'])
        if n_ldm == 1:
            counts['n_clusters_1_ldm'] += 1
            counts['cluster_depth_1_ldm'] += depth
        else:
            counts['n_clusters_2plus_ldm'] += 1
            counts['cluster_depth_2plus_ldm'] += depth
    return counts


def bayes_lineage_probability(
    covariants: pd.DataFrame,
    ldm_set: set,
    pango_lin: str,
    p_lineage: float,
    lineage_count: float,
    datemin: str,
    datemax: str,
) -> float:
    """
    Depth-weighted mean of P(L|M) across all covar clusters containing ≥1 LDM.

    P(L|M) = P(M|L) * P(L) / P(M)   [Bayes, mirrors validate_detects.ipynb]

    Returns nan when no LDM-containing clusters exist or all API calls fail.
    """
    if covariants.empty or p_lineage <= 0:
        return float('nan')

    weighted_sum = 0.0
    total_depth = 0

    for _, row in covariants.iterrows():
        nt_muts = parse_nt_muts(row.get('nt_mutations', ''))
        if not nt_muts or not any(m in ldm_set for m in nt_muts):
            continue

        aa_muts = parse_aa_muts(row.get('aa_mutations', ''))
        if not aa_muts:
            continue

        depth = int(row['cluster_depth'])
        try:
            p_m = get_p_mutations(aa_muts, datemin, datemax)
            if p_m <= 0:
                continue
            p_m_given_l = get_p_mutations_given_lineage(
                aa_muts, pango_lin, lineage_count, datemin, datemax,
            )
            p_l_given_m = min(1.0, p_m_given_l * p_lineage / p_m)
        except Exception:
            continue

        weighted_sum += depth * p_l_given_m
        total_depth += depth

    if total_depth == 0:
        return float('nan')
    return weighted_sum / total_depth


# ── main ──────────────────────────────────────────────────────────────────────

os.makedirs(OUTPUT_DIR, exist_ok=True)

for file in sorted(os.listdir('samples')):
    if not file.endswith('.csv'):
        continue

    lineage = file.split('detections_')[1].split('.csv')[0]
    pango_lin = to_pango(lineage)
    print(f'Processing {lineage} ({pango_lin})')

    samples = pd.read_csv(f'samples/{file}')

    dates = pd.to_datetime(samples['collection_date_ww'], errors='coerce').dropna()
    datemin = dates.min().strftime('%Y-%m-%d')
    datemax = dates.max().strftime('%Y-%m-%d')

    ldm_set = get_ldms(pango_lin)
    p_lineage, lineage_count = get_lineage_prevalence(pango_lin, datemin, datemax)
    print(f'  LDMs: {len(ldm_set)}  P(L): {p_lineage:.4f}  dates: {datemin} – {datemax}')

    rows = []
    for _, sample in samples.iterrows():
        acc = sample['accession']
        covar_path = f'covar/{acc}.covar.tsv'
        if os.path.isfile(covar_path):
            covariants = pd.read_csv(covar_path, sep='\t')
            counts = count_lineage_clusters(covariants, ldm_set)
            p_present = bayes_lineage_probability(
                covariants, ldm_set, pango_lin, p_lineage, lineage_count, datemin, datemax,
            )
        else:
            counts = {
                'n_clusters_1_ldm': 0,
                'cluster_depth_1_ldm': 0,
                'n_clusters_2plus_ldm': 0,
                'cluster_depth_2plus_ldm': 0,
            }
            p_present = float('nan')
        rows.append({
            'collection_date_ww': sample['collection_date_ww'],
            'state': sample['state'],
            'accession': acc,
            **counts,
            'p_lineage_present': p_present,
        })

    out_name = file.replace('initial_', 'validated_')
    pd.DataFrame(rows).to_csv(f'{OUTPUT_DIR}/{out_name}', index=True)
    print(f'  → {OUTPUT_DIR}/{out_name}')

_save_caches()
