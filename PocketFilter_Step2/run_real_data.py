"""Repackage the real-data workflow used in the initial Windows experiment."""
import ast
import csv
import json
from pathlib import Path
import numpy as np
from scipy import stats
from ont_fast5_api.fast5_interface import get_fast5_file
from golden_model import score_preprocessed

LENGTH, MAX_READS = 1000, 10
ROOT = Path(__file__).resolve().parent

def original_helpers():
    nb = json.loads((ROOT / 'upstream_sdtw_analysis.ipynb').read_text(encoding='utf-8'))
    wanted = {'get_fasta', 'rev_comp', 'load_model', 'discrete_normalize',
              'ref_signal', 'get_stall_end', 'trim'}
    definitions = {}
    for cell in nb['cells']:
        if cell['cell_type'] != 'code':
            continue
        source = ''.join(cell['source'])
        try:
            tree = ast.parse(source)
        except SyntaxError:
            continue
        for node in tree.body:
            if isinstance(node, ast.FunctionDef) and node.name in wanted:
                definitions[node.name] = ast.get_source_segment(source, node)
    if set(definitions) != wanted:
        raise RuntimeError('Missing original preprocessing helpers')
    scope = {'np': np, 'stats': stats, 'k': 6, 'prefix_lengths': np.array([LENGTH])}
    for code in definitions.values():
        exec(code, scope)
    return scope

def main():
    out = ROOT / 'results'
    out.mkdir(exist_ok=True)
    h = original_helpers()
    fasta = h['get_fasta'](ROOT / 'data/covid/reference.fasta')
    model = h['load_model'](ROOT / 'data/dna_kmer_model.txt')
    reference = np.concatenate([h['ref_signal'](fasta, model),
        h['ref_signal'](h['rev_comp'](fasta), model)]).astype(np.int32)
    np.save(out / 'reference.npy', reference)
    rows = []
    seen = set()
    for label in ('covid', 'human'):
        accepted = scanned = skipped = 0
        for filename in sorted((ROOT / f'data/{label}/fast5').rglob('*.fast5')):
            with get_fast5_file(str(filename), mode='r') as handle:
                for read in handle.get_reads():
                    if accepted >= MAX_READS:
                        break
                    scanned += 1
                    read_id = str(read.read_id)
                    if read_id in seen:
                        skipped += 1
                        continue
                    raw = np.asarray(read.get_raw_data(), dtype=np.float64)
                    if len(raw) < LENGTH or not np.isfinite(raw).all() or np.std(raw) == 0:
                        skipped += 1
                        continue
                    trimmed, stall_end = h['trim'](raw)
                    if len(trimmed) < LENGTH:
                        skipped += 1
                        continue
                    query = np.empty(LENGTH, dtype=np.float64)
                    valid = True
                    for start in range(0, LENGTH, 500):
                        end = min(start + 500, LENGTH)
                        prefix = trimmed[:end]
                        mean = int(np.mean(prefix))
                        if int(np.mean(np.abs(prefix - mean))) == 0:
                            valid = False
                            break
                        query[start:end] = h['discrete_normalize'](prefix)[start:end]
                    if not valid:
                        skipped += 1
                        continue
                    print(f'{label}: scoring {read_id}', flush=True)
                    cost = float(score_preprocessed(reference, query)[0])
                    accepted += 1
                    seen.add(read_id)
                    vector = out / f'{label}_{accepted:03d}.npz'
                    np.savez(vector, query=query.astype(np.int32), cost=np.array(cost),
                        read_id=np.array(read_id), label=np.array(label), length=np.array(LENGTH),
                        stall_end=np.array(stall_end), source_file=np.array(str(filename)))
                    rows.append(dict(label=label, read_id=read_id, length=LENGTH,
                        cost=cost, vector_file=str(vector), source_file=str(filename)))
                    print(f'cost={cost}', flush=True)
            if accepted >= MAX_READS:
                break
        print(f'{label}: accepted={accepted}, scanned={scanned}, skipped={skipped}')
        if accepted == 0:
            raise RuntimeError(f'No valid {label} reads; check data paths')
    with (out / 'initial_scores.csv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for label in ('covid', 'human'):
        v = [r['cost'] for r in rows if r['label'] == label]
        print(f'{label}: n={len(v)}, min={min(v)}, median={np.median(v)}, max={max(v)}')

if __name__ == '__main__':
    main()
