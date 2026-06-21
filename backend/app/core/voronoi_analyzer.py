import numpy as np
from scipy.spatial import ConvexHull
from typing import List, Dict
from collections import Counter
from .cell_binning import CellBinning
from .md_unpacker import TriclinicBox


class VoronoiAnalyzer:
    def __init__(self, cutoff: float = 5.0, max_neighbors: int = 24):
        self.cutoff = cutoff
        self.max_neighbors = max_neighbors

    def _compute_voronoi_index(self, central_coord: np.ndarray, neighbor_coords: np.ndarray,
                                box: TriclinicBox) -> np.ndarray:
        n_neighbors = len(neighbor_coords)
        if n_neighbors < 4:
            return np.zeros(10, dtype=np.int32)
        rel_coords = np.zeros_like(neighbor_coords)
        for i in range(n_neighbors):
            rel_coords[i] = box.minimum_image_vector(neighbor_coords[i] - central_coord)
        try:
            hull = ConvexHull(rel_coords)
        except:
            return np.zeros(10, dtype=np.int32)
        face_orders = []
        for simplex in hull.simplices:
            face_orders.append(len(simplex))
        vor_index = np.zeros(10, dtype=np.int32)
        for order in face_orders:
            idx = order - 3
            if 0 <= idx < len(vor_index):
                vor_index[idx] += 1
        return vor_index

    def analyze_atom(self, atom_idx: int, coords: np.ndarray, box: TriclinicBox,
                     cb: CellBinning = None) -> np.ndarray:
        if cb is None:
            cb = CellBinning(coords, box, self.cutoff)
        neighbors, distances = cb.get_neighbor_list_sorted(atom_idx, self.max_neighbors)
        if len(neighbors) == 0:
            return np.zeros(10, dtype=np.int32)
        neighbor_coords = np.array([coords[j] for j in neighbors])
        return self._compute_voronoi_index(coords[atom_idx], neighbor_coords, box)

    def analyze_frame(self, coords: np.ndarray, box: TriclinicBox, types: np.ndarray = None,
                      progress_callback=None) -> Dict:
        n_atoms = coords.shape[0]
        cb = CellBinning(coords, box, self.cutoff)
        vor_indices = np.zeros((n_atoms, 10), dtype=np.int32)
        for i in range(n_atoms):
            neighbors, distances = cb.get_neighbor_list_sorted(i, self.max_neighbors)
            if len(neighbors) > 0:
                neighbor_coords = np.array([coords[j] for j in neighbors])
                vor_indices[i] = self._compute_voronoi_index(coords[i], neighbor_coords, box)
            if progress_callback and (i + 1) % 10000 == 0:
                progress_callback(i + 1, n_atoms)
        vor_index_strings = [tuple(idx) for idx in vor_indices]
        index_counts = Counter(vor_index_strings)
        result = {
            'voronoi_indices': vor_indices.tolist(),
            'counts': {str(k): v for k, v in index_counts.items()},
            'n_atoms': n_atoms
        }
        if types is not None:
            unique_types = np.unique(types)
            type_index_counts = {}
            for t in unique_types:
                mask = types == t
                type_indices = [tuple(vor_indices[i]) for i in range(n_atoms) if mask[i]]
                type_index_counts[int(t)] = {str(k): v for k, v in Counter(type_indices).items()}
            result['by_type'] = type_index_counts
        return result

    def classify_polyhedron(self, vor_index: np.ndarray) -> str:
        idx_tuple = tuple(vor_index[:4])
        if idx_tuple == (0, 0, 12, 0):
            return 'Icosahedron'
        elif idx_tuple == (0, 0, 12, 2):
            return 'Distorted Icosahedron'
        elif idx_tuple == (0, 2, 8, 2):
            return 'BCC-like'
        elif idx_tuple == (0, 3, 6, 4):
            return 'FCC-like'
        elif idx_tuple == (0, 4, 4, 6):
            return 'HCP-like'
        else:
            return f'Other_{vor_index[0]}{vor_index[1]}{vor_index[2]}{vor_index[3]}'

    def get_polyhedron_type(self, vor_index: np.ndarray) -> str:
        return f'<{vor_index[0]},{vor_index[1]},{vor_index[2]},{vor_index[3]}>'

    def get_evolution_matrix(self, frames: list, target_types: List[str] = None) -> Dict:
        timesteps = []
        evolution_data = {}
        for fi, frame in enumerate(frames):
            result = self.analyze_frame(frame.coords, frame.box, frame.types)
            timesteps.append(frame.timestep)
            type_strings = [self.get_polyhedron_type(np.array(idx))
                            for idx in result['voronoi_indices']]
            counts = Counter(type_strings)
            for ptype, count in counts.items():
                if ptype not in evolution_data:
                    evolution_data[ptype] = [0] * fi
                evolution_data[ptype].append(count)
            for ptype in evolution_data:
                if len(evolution_data[ptype]) <= fi:
                    evolution_data[ptype].append(0)
        if target_types:
            filtered = {}
            for t in target_types:
                if t in evolution_data:
                    filtered[t] = evolution_data[t]
                else:
                    filtered[t] = [0] * len(timesteps)
            evolution_data = filtered
        return {
            'timesteps': timesteps,
            'evolution': evolution_data,
            'n_frames': len(frames)
        }
