from pathlib import Path
import argparse
import csv
import hashlib
import json

from analyze_cycle_conservation import cycle_balances


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    report = {'schema': 'primsf_active_smoke_six_balance_normalizations_v1',
              'scope': 'Whole observed interval from first to last complete ledger row; these are not oscillation cycles.',
              'normalizations': 'max over ledger of abs(raw balance), divided separately by initial inventory and end-of-window absolute inlet+outlet processed mass; both limits1e-6.',
              'runs': []}
    for case in ('D_SS1', 'D_SS3'):
        path = args.root / case / 'output/conservation.csv'
        if not path.exists():
            continue
        with path.open() as stream:
            complete = (line for line in stream if line.endswith('\n'))
            data = [{key: float(value) for key, value in row.items()}
                    for row in csv.DictReader(complete) if None not in row and None not in row.values()]
        windows = cycle_balances(data, [data[0]['time_s'], data[-1]['time_s']])
        result = windows[0]
        result['case'] = case
        result['ledger_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
        result['accepted_rows'] = len(data)
        result['first_step'] = int(data[0]['step'])
        result['last_step'] = int(data[-1]['step'])
        result['all_six_pass_both_1e6'] = all(row['passes_inventory_1e6'] and row['passes_processed_1e6']
                                             for row in result['balances'].values())
        report['runs'].append(result)
        print(case, result['all_six_pass_both_1e6'])
        for name, row in result['balances'].items():
            print(name, row['relative_inventory'], row['relative_processed'],
                  row['passes_inventory_1e6'], row['passes_processed_1e6'])
    (args.root / 'balance_normalizations.json').write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
