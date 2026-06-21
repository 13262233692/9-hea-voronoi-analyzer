import numpy as np
from typing import Dict, Tuple
from .cell_binning import CellBinning


class RDFCalculator:
    def __init__(self, r_min: float = 0.0, r_max: float = 6.0, n_bins: int = 100):
        self.r_min = r_min
        self.r_max = r_max
        self.n_bins = n_bins
        self.bin_edges = np.linspace(r_min, r_max, n_bins + 1)
        self.bin_centers = (self.bin_edges[:-1] + self.bin_edges[1:]) / 2.0
        self.dr = (r_max - r_min) / n_bins

    def calculate(self, coords: np.ndarray, box: np.ndarray, types: np.ndarray = None) -> Dict:
        n_atoms = coords.shape[0]
        box_lengths = box[:, 1] - box[:, 0]
        volume = np.prod(box_lengths)
        density = n_atoms / volume
        cb = CellBinning(coords, box, self.r_max)
        gr = np.zeros(self.n_bins, dtype=np.float64)
        if types is not None:
            unique_types = np.unique(types)
            partial_gr = {}
            for t1 in unique_types:
                for t2 in unique_types:
                    partial_gr[(int(t1), int(t2))] = np.zeros(self.n_bins, dtype=np.float64)
            type_counts = {t: np.sum(types == t) for t in unique_types}
        for i in range(n_atoms):
            neighbors, distances = cb.get_neighbors_with_distances(i)
            valid = (distances >= self.r_min) & (distances < self.r_max)
            d_valid = distances[valid]
            bin_idx = np.floor((d_valid - self.r_min) / self.dr).astype(int)
            for idx in bin_idx:
                gr[idx] += 1.0
            if types is not None:
                t_i = int(types[i])
                for j, d in zip(neighbors, distances):
                    if d >= self.r_min and d < self.r_max:
                        t_j = int(types[j])
                        bidx = int(np.floor((d - self.r_min) / self.dr))
                        partial_gr[(t_i, t_j)][bidx] += 0.5
        shell_volumes = 4.0 * np.pi * (self.bin_centers ** 2) * self.dr
        norm = density * n_atoms
        gr_normalized = gr / (norm * shell_volumes)
        result = {
            'r': self.bin_centers.tolist(),
            'gr': gr_normalized.tolist(),
            'gr_raw': gr.tolist()
        }
        if types is not None:
            partial_normalized = {}
            for (t1, t2), pg in partial_gr.items():
                n1 = type_counts[t1]
                n2 = type_counts[t2]
                if t1 == t2:
                    norm_p = (n2 / volume) * n1
                else:
                    norm_p = (n2 / volume) * n1
                pg_norm = pg / (norm_p * shell_volumes)
                partial_normalized[f'{t1}-{t2}'] = pg_norm.tolist()
            result['partial_gr'] = partial_normalized
        return result

    def calculate_batch(self, frames: list) -> Dict:
        avg_gr = np.zeros(self.n_bins, dtype=np.float64)
        n_frames = len(frames)
        for frame in frames:
            res = self.calculate(frame.coords, frame.box, frame.types)
            avg_gr += np.array(res['gr'])
        avg_gr /= n_frames
        return {
            'r': self.bin_centers.tolist(),
            'gr_avg': avg_gr.tolist()
        }
