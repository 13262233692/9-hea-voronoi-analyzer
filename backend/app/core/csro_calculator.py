import numpy as np
from typing import Dict
from .cell_binning import CellBinning
from .md_unpacker import TriclinicBox


class CSROCalculator:
    def __init__(self, cutoff: float = 3.5):
        self.cutoff = cutoff

    def calculate_csro(self, coords: np.ndarray, box: TriclinicBox, types: np.ndarray,
                       type1: int, type2: int, use_perfect_crystal: bool = True) -> Dict:
        n_atoms = len(types)
        cb = CellBinning(coords, box, self.cutoff)

        c1 = np.sum(types == type1) / n_atoms
        c2 = np.sum(types == type2) / n_atoms

        N_AA = 0
        N_BB = 0
        N_AB = 0
        N_BA = 0
        Z_total = 0

        for i in range(n_atoms):
            ti = types[i]
            neighbors = cb.get_neighbors(i)
            z_i = len(neighbors)
            Z_total += z_i

            for j in neighbors:
                tj = types[j]
                if ti == type1 and tj == type1:
                    N_AA += 1
                elif ti == type2 and tj == type2:
                    N_BB += 1
                elif ti == type1 and tj == type2:
                    N_AB += 1
                elif ti == type2 and tj == type1:
                    N_BA += 1

        N_AA /= 2
        N_BB /= 2
        N_AB /= 2
        N_BA /= 2

        Z_avg = Z_total / n_atoms

        if use_perfect_crystal:
            N_AB_random = n_atoms * c1 * c2 * Z_avg
            if N_AB_random > 0:
                alpha = 1 - N_AB / N_AB_random
            else:
                alpha = 0.0
        else:
            p_AB = N_AB / (N_AA + N_AB + 1e-10)
            p_BA = N_BA / (N_BB + N_BA + 1e-10)
            alpha = 1 - (p_AB + p_BA) / (c1 + c2)

        S_ij = -alpha

        bond_counts = {}
        for i in range(n_atoms):
            neighbors = cb.get_neighbors(i)
            for j in neighbors:
                if i < j:
                    pair = f'{min(types[i], types[j])}-{max(types[i], types[j])}'
                    if pair not in bond_counts:
                        bond_counts[pair] = 0
                    bond_counts[pair] += 1

        type_pairs = {
            f'{type1}-{type1}': N_AA,
            f'{type2}-{type2}': N_BB,
            f'{type1}-{type2}': N_AB,
            f'{type2}-{type1}': N_BA
        }

        return {
            'S_ij': float(S_ij),
            'alpha': float(alpha),
            'c1': float(c1),
            'c2': float(c2),
            'Z_avg': float(Z_avg),
            'bond_counts': bond_counts,
            'type_pairs': type_pairs,
            'N_AA': float(N_AA),
            'N_BB': float(N_BB),
            'N_AB': float(N_AB),
            'volume': float(box.volume),
            'cell_type': 'triclinic' if not box.orthorhombic else 'orthorhombic'
        }

    def calculate_warren_cowley(self, coords: np.ndarray, box: TriclinicBox,
                         types: np.ndarray, type1: int, type2: int) -> float:
        result = self.calculate_csro(coords, box, types, type1, type2)
        return result['S_ij']

    def analyze_all_pairs(self, coords: np.ndarray, box: TriclinicBox,
                        types: np.ndarray) -> Dict:
        unique_types = sorted(np.unique(types).tolist())
        results = {}
        for i, t1 in enumerate(unique_types):
            for t2 in unique_types[i:]:
                key = f'{t1}-{t2}'
                results[key] = self.calculate_csro(coords, box, types, t1, t2)
        return results

    def calculate_csro_batch(self, frames: list, type1: int, type2: int) -> Dict:
        s_ij_values = []
        timesteps = []
        for frame in frames:
            result = self.calculate_csro(frame.coords, frame.box, frame.types, type1, type2)
            s_ij_values.append(result['S_ij'])
            timesteps.append(frame.timestep)
        return {
            'timesteps': timesteps,
            'S_ij': s_ij_values,
            'S_ij_avg': float(np.mean(s_ij_values)),
            'S_ij_std': float(np.std(s_ij_values))
        }
