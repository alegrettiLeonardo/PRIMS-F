from pathlib import Path
import json
import numpy as np

report = Path('/home/leonardoab/Documentos/PRIMS-F-numerics-20261119/validation/diagnostics/numerics_p0p1_20261119')
provenance = json.loads((report / 'provenance.json').read_text())


def read_operators(path):
    lines = path.read_text().splitlines()[1:]
    records = []
    cursor = 0
    while cursor < len(lines):
        name = lines[cursor]
        profile, cell, closure, status = map(int, lines[cursor + 1].split())
        angle = float(lines[cursor + 2])
        arrays = [np.array(list(map(float, line.split()))) for line in lines[cursor + 3:cursor + 6]]
        record = dict(name=name, profile=profile, cell=cell, closure=closure, status=status,
                      angle=angle, U=arrays[0], v=arrays[1], F=arrays[2])
        cursor += 6
        if status == 0:
            record['A'] = np.array([list(map(float, line.split())) for line in lines[cursor:cursor + 3]])
            record['wr'] = np.array(list(map(float, lines[cursor + 3].split())))
            record['wi'] = np.array(list(map(float, lines[cursor + 4].split())))
            record['R'] = np.array([list(map(float, line.split())) for line in lines[cursor + 5:cursor + 8]])
            cursor += 8
        records.append(record)
    return records


def pressure_gradient(record):
    alpha, rl, rg, p, ul, ug = record['v']
    cl2 = 1498.0**2 * 997.05 / 101325.0
    cg2 = 343.0**2 * 997.05 / 101325.0
    denominator = record['U'][0] / (cl2 * rl**2) + record['U'][1] / (cg2 * rg**2)
    return np.array([1 / rl, 1 / rg, 0]) / denominator


summary = {'kind': 'offline frozen 9f2ed21 operator audit', 'runs': {}}
for case, run in provenance['runs'].items():
    records = read_operators(report / case / 'state_operators.dat')
    speed_scale = run['scales'][-1]
    item = {'states': [], 'linear_characteristic_interface_diagnostics': []}
    for record in records:
        assert record['status'] == 0, record
        assert np.all(record['wi'] == 0), record
        wr = record['wr']
        residual = record['A'] @ record['R'] - record['R'] * wr[None, :]
        entry = {key: record[key] for key in ('name', 'profile', 'cell', 'closure', 'status', 'angle')}
        entry['primitive_SI'] = dict(zip(
            ['alpha', 'rho_l_kg_m3', 'rho_g_kg_m3', 'p_Pa', 'ul_m_s', 'ug_m_s'],
            (record['v'] * np.array([1, 997.05, 997.05, 101325, speed_scale, speed_scale])).tolist()))
        entry['eigenvalues_m_s'] = sorted((wr * speed_scale).tolist())
        entry['eigen_imaginary_max_m_s'] = float(np.max(abs(record['wi'])) * speed_scale)
        entry['eigen_residual_relative'] = float(np.max(abs(residual)) / max(1, np.max(abs(record['A']))))
        entry['positive_count'] = int(sum(wr > 0))
        entry['negative_count'] = int(sum(wr < 0))
        entry['zero_count'] = int(sum(wr == 0))
        if record['name'].startswith('inlet'):
            entry['incoming_count'] = entry['positive_count']
            incoming = record['R'][:, wr > 0]
            constraint = record['A'][:2, :]
            normalized = constraint / np.linalg.norm(constraint, axis=1)[:, None]
            singular = np.linalg.svd(normalized @ incoming, compute_uv=False)
            entry['phase_flux_constraints_projected_singular_values'] = singular.tolist()
            entry['jl_relative_to_prescribed_volume_change'] = float(
                (1 - record['v'][0]) * record['v'][4] * speed_scale / run['controls']['jl'] - 1)
        if record['name'].startswith('outlet'):
            entry['incoming_count'] = entry['negative_count']
            incoming = record['R'][:, wr < 0]
            gradient = pressure_gradient(record)
            projected = gradient @ incoming / np.linalg.norm(gradient)
            entry['pressure_constraint_projected_abs'] = abs(projected).tolist()
        item['states'].append(entry)
    for profile in (1, 2):
        left = next(r for r in records if r['name'] == 'junction_left_cell' and r['profile'] == profile)
        right = next(r for r in records if r['name'] == 'junction_right_cell' and r['profile'] == profile)
        left_in = left['R'][:, left['wr'] < 0]
        right_in = right['R'][:, right['wr'] > 0]
        cl = np.vstack((left['A'][:2, :], pressure_gradient(left)))
        cr = np.vstack((right['A'][:2, :], pressure_gradient(right)))
        matrix = np.column_stack((cl @ left_in, -cr @ right_in))
        mismatch = np.array([*(left['F'][:2] - right['F'][:2]), left['v'][3] - right['v'][3]])
        amplitudes = np.linalg.solve(matrix, -mismatch)
        du_left = left_in @ amplitudes[:left_in.shape[1]]
        du_right = right_in @ amplitudes[left_in.shape[1]:]
        flux_left = left['F'][:2] + left['A'][:2] @ du_left
        flux_right = right['F'][:2] + right['A'][:2] @ du_right
        p_left = left['v'][3] + pressure_gradient(left) @ du_left
        p_right = right['v'][3] + pressure_gradient(right) @ du_right
        normalized = matrix / np.linalg.norm(matrix, axis=1)[:, None]
        item['linear_characteristic_interface_diagnostics'].append({
            'profile': profile, 'incoming_to_left': left_in.shape[1], 'incoming_to_right': right_in.shape[1],
            'constraint_matrix_shape': list(matrix.shape),
            'row_normalized_singular_values': np.linalg.svd(normalized, compute_uv=False).tolist(),
            'linear_phase_flux_mismatch_nd': (flux_left - flux_right).tolist(),
            'linear_pressure_mismatch_nd': float(p_left - p_right),
            'linear_common_pressure_Pa': float(p_left * 101325),
            'linear_common_phase_mass_flux_kg_m2_s': (flux_left * 997.05 * speed_scale).tolist(),
            'linear_left_U': (left['U'] + du_left).tolist(),
            'linear_right_U': (right['U'] + du_right).tolist(),
            'linear_left_U_relative_increment': (du_left / left['U']).tolist(),
            'linear_right_U_relative_increment': (du_right / right['U']).tolist(),
            'qualification': 'one linear diagnostic only; cell angles differ and source integration/nonlinear wave curves are omitted; not a replacement Riemann solver',
        })
    summary['runs'][case] = item
(report / 'characteristics_summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
for case, run in summary['runs'].items():
    print(case)
    for state in run['states']:
        if state['profile'] == 2:
            print(state['name'], state['eigenvalues_m_s'], 'incoming', state.get('incoming_count'))
    print('linear diagnostic', run['linear_characteristic_interface_diagnostics'][1])
