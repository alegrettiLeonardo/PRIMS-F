from pathlib import Path
import argparse
import csv
import json


def read_rows(path):
    with path.open() as stream:
        return [row for row in csv.DictReader(line for line in stream if line.endswith('\n'))
                if None not in row and None not in row.values()]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    report = {'scope': 'Accepted steps after78.5s with zero additional linear solves; indicates initial residual stopping before Newton update, not proof of a globally stationary flow.', 'runs': []}
    names = ('Ml_kg', 'Mg_kg', 'Ml_pipe_kg', 'Mg_pipe_kg', 'Ml_riser_kg', 'Mg_riser_kg')
    for case in ('D_SS1', 'D_SS3'):
        output = args.root / case / 'output'
        if not output.exists():
            continue
        steps = read_rows(output / 'steps.csv')
        ledger = {int(row['step']): row for row in read_rows(output / 'conservation.csv')}
        observed, zero = [], []
        for previous, row in zip(steps, steps[1:]):
            if float(row['time_s']) <= 78.5 or int(row['step']) != int(previous['step']) + 1:
                continue
            observed.append(row)
            if int(row['linear_solves_total']) != int(previous['linear_solves_total']):
                continue
            before = ledger.get(int(previous['step']))
            after = ledger.get(int(row['step']))
            unchanged = None
            if before is not None and after is not None:
                unchanged = all(float(before[key]) == float(after[key]) for key in names)
            zero.append({'step': int(row['step']), 'time_s': float(row['time_s']),
                         'dt_s': float(row['dt_s']), 'residual_norm': float(row['residual_norm']),
                         'six_inventories_bitwise_equal_as_serialized': unchanged})
        record = {'case': case, 'observed_steps': len(observed), 'zero_solve_accepted_steps': len(zero),
                  'observed_time_increment_s': sum(float(row['dt_s']) for row in observed),
                  'zero_solve_time_increment_s': sum(row['dt_s'] for row in zero),
                  'unchanged_six_inventories': sum(row['six_inventories_bitwise_equal_as_serialized'] is True for row in zero),
                  'first_examples': zero[:3], 'last_examples': zero[-3:]}
        report['runs'].append(record)
        print(json.dumps(record))
    (args.root / 'accepted_zero_updates.json').write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
