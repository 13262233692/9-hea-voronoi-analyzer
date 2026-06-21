import numpy as np
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from .md_unpacker import TriclinicBox
from .cell_binning import CellBinning


@dataclass
class NonAffineResult:
    d2_min: np.ndarray
    shear_band_atoms: np.ndarray
    threshold: float
    n_yielding: int
    total_atoms: int
    dissipated_energy: float
    strain_tensor: Optional[np.ndarray] = None
    atom_strain_magnitudes: Optional[np.ndarray] = None


class NonAffineAnalyzer:
    def __init__(self, cutoff: float = 4.0, yield_threshold: float = 0.08):
        self.cutoff = cutoff
        self.yield_threshold = yield_threshold

    def analyze_deformation(self,
                            coords_ref: np.ndarray,
                            coords_def: np.ndarray,
                            box: TriclinicBox,
                            types: np.ndarray = None) -> NonAffineResult:
        coords_ref = np.asarray(coords_ref, dtype=np.float64)
        coords_def = np.asarray(coords_def, dtype=np.float64)
        n_atoms = coords_ref.shape[0]

        cb = CellBinning(coords_ref, box, self.cutoff)

        d2_min = np.zeros(n_atoms, dtype=np.float64)
        strain_mags = np.zeros(n_atoms, dtype=np.float64)
        F_tensors = np.zeros((n_atoms, 3, 3), dtype=np.float64)

        for i in range(n_atoms):
            neighbors = cb.get_neighbors(i)
            if len(neighbors) < 4:
                d2_min[i] = 0.0
                strain_mags[i] = 0.0
                continue

            r_ref_neighbors = coords_ref[neighbors]
            r_def_neighbors = coords_def[neighbors]

            dr_ref = box.minimum_image_vector(
                r_ref_neighbors - coords_ref[i]
            )
            dr_def = box.minimum_image_vector(
                r_def_neighbors - coords_def[i]
            )

            F, residual = self._solve_affine_F(dr_ref, dr_def)
            d2_min[i] = residual / len(neighbors)
            F_tensors[i] = F

            C = F.T @ F
            trace_C = max(np.trace(C), 3.0)
            strain_mags[i] = np.sqrt(trace_C / 3.0 - 1.0)

        yielding_mask = d2_min > self.yield_threshold
        yielding_indices = np.where(yielding_mask)[0]

        dissipated_energy = np.sum(d2_min[yielding_mask]) * 0.5

        return NonAffineResult(
            d2_min=d2_min,
            shear_band_atoms=yielding_indices,
            threshold=self.yield_threshold,
            n_yielding=len(yielding_indices),
            total_atoms=n_atoms,
            dissipated_energy=float(dissipated_energy),
            atom_strain_magnitudes=strain_mags
        )

    def _solve_affine_F(self, R: np.ndarray, R_prime: np.ndarray) -> Tuple[np.ndarray, float]:
        A = R.T @ R
        b = R.T @ R_prime
        F_T = np.linalg.lstsq(A, b, rcond=None)[0]
        F = F_T.T

        transformed = R @ F_T
        residuals = np.sum((R_prime - transformed) ** 2, axis=1)
        total_residual = np.sum(residuals)

        return F, total_residual

    def detect_shear_band_clusters(self,
                                   coords: np.ndarray,
                                   yielding_indices: np.ndarray,
                                   box: TriclinicBox,
                                   cluster_cutoff: float = 5.0) -> List[Dict]:
        if len(yielding_indices) == 0:
            return []

        yielding_coords = coords[yielding_indices]
        n_yield = len(yielding_indices)

        visited = np.zeros(n_yield, dtype=bool)
        clusters = []

        cb_yield = CellBinning(yielding_coords, box, cluster_cutoff)

        for i in range(n_yield):
            if visited[i]:
                continue

            stack = [i]
            visited[i] = True
            cluster_local = []

            while stack:
                curr = stack.pop()
                cluster_local.append(curr)
                neighbors = cb_yield.get_neighbors(curr)
                for nb in neighbors:
                    if not visited[nb]:
                        visited[nb] = True
                        stack.append(nb)

            cluster_global = yielding_indices[cluster_local]
            cluster_coords = coords[cluster_global]
            center = np.mean(cluster_coords, axis=0)

            clusters.append({
                'atom_indices': cluster_global.tolist(),
                'n_atoms': len(cluster_global),
                'center': center.tolist(),
                'coords': cluster_coords.tolist()
            })

        clusters.sort(key=lambda c: c['n_atoms'], reverse=True)
        return clusters

    def compute_dissipation_curve(self,
                                  coords_ref: np.ndarray,
                                  coords_def: np.ndarray,
                                  box: TriclinicBox,
                                  n_bins: int = 50) -> Dict:
        result = self.analyze_deformation(coords_ref, coords_def, box)
        d2 = result.d2_min

        hist_edges = np.logspace(-4, 1, n_bins + 1)
        hist_counts, _ = np.histogram(d2, bins=hist_edges)

        cumulative = np.cumsum(hist_counts * hist_edges[:-1] * 0.5)

        return {
            'bins': hist_edges.tolist(),
            'counts': hist_counts.tolist(),
            'cumulative_dissipation': cumulative.tolist(),
            'total_dissipated': float(np.sum(d2) * 0.5),
            'yielding_atoms': int(result.n_yielding),
            'threshold': self.yield_threshold
        }
