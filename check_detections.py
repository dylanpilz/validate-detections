import os

import pandas as pd
from outbreak_data import authenticate_user, outbreak_data as od
from outbreak_tools import crumbs

OUTPUT_DIR = 'validated'

#authenticate_user.authenticate_new_user()
lineage_key = crumbs.get_alias_key()

# ── per-run caches ────────────────────────────────────────────────────────────
_ldm_cache: dict[str, set[str]] = {}
_lineage_prev_cache: dict[tuple, tuple[float, float]] = {}
_pm_cache: dict[tuple, float] = {}
_pml_cache: dict[tuple, float] = {}


# ── API helpers ───────────────────────────────────────────────────────────────

def to_pango(lineage: str) -> str:
    return lineage.removesuffix('.X')


def get_ldms(pango_lin: str) -> set[str]:
    """Lineage-defining aa_mutations from outbreak API (freq ≥ 0.8 in lineage)."""
    if pango_lin not in _ldm_cache:
        df = od.lineage_mutations(pango_lin=pango_lin, descendants=True, lineage_key=lineage_key)
        _ldm_cache[pango_lin] = set(df.index.to_list())
    return _ldm_cache[pango_lin]


def get_lineage_prevalence(pango_lin: str, datemin: str, datemax: str) -> tuple[float, float]:
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


def get_p_mutations(muts: list[str], datemin: str, datemax: str) -> float:
    """P(M) — overall prevalence of this mutation set across all lineages."""
    key = (frozenset(muts), datemin, datemax)
    if key not in _pm_cache:
        df = od.lineage_cl_prevalence(
            '.', descendants=True, mutations=muts, location='USA',
            datemin=datemin, datemax=datemax, lineage_key=lineage_key,
        )
        tc = float(df['total_count'].sum())
        _pm_cache[key] = float(df['lineage_count'].sum()) / tc if tc > 0 else 0.0
    return _pm_cache[key]


def _multiquery_to_df(data: dict) -> pd.DataFrame:
    return pd.concat(
        [pd.DataFrame(v).assign(query=k) for k, v in data['results'].items()], axis=0
    )


def get_p_mutations_given_lineage(
    muts: list[str], pango_lin: str, lineage_count: float, datemin: str, datemax: str,
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

def parse_aa_muts(aa_mutations_val) -> list[str]:
    """
    Parse aa_mutations cell from covar TSV into a deduplicated list.
    Keeps only tokens that contain ':' (valid GENE:CHANGE format).
    """
    return list(dict.fromkeys(m.lower() for m in str(aa_mutations_val).split() if ':' in m))


def count_lineage_clusters(covariants: pd.DataFrame, ldm_set: set[str]) -> dict[str, int]:
    counts = {
        'n_clusters_1_ldm': 0,
        'cluster_depth_1_ldm': 0,
        'n_clusters_2plus_ldm': 0,
        'cluster_depth_2plus_ldm': 0,
    }
    for _, row in covariants.iterrows():
        aa_muts = parse_aa_muts(row.get('aa_mutations', ''))
        n_ldm = sum(1 for m in aa_muts if m in ldm_set)
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
    ldm_set: set[str],
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
        aa_muts = parse_aa_muts(row.get('aa_mutations', ''))
        if not aa_muts or not any(m in ldm_set for m in aa_muts):
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
