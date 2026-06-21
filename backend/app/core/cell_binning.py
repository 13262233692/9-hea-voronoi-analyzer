import numpy as np
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass
from .md_unpacker import TriclinicBox


class CellBinning:
    def __init__(self, coords: np.ndarray, box: TriclinicBox, cutoff: float,
                 use_fractional: bool = True):
        self.coords = np.asarray(coords, dtype=np.float64)
        self.box = box
        self.cutoff = cutoff
        self.n_atoms = self.coords.shape[0]
        self.use_fractional = use_fractional
        self._build_cells()

    def _build_cells(self):
        if self.use_fractional and not self.box.orthorhombic:
            self._build_fractional_cells()
        else:
            self._build_cartesian_cells()

    def _build_cartesian_cells(self):
        bounds = self.box.get_bounds()
        lengths = bounds[:, 1] - bounds[:, 0]
        self.n_cells = np.maximum(1, (lengths / self.cutoff).astype(int))
        self.cell_size = lengths / self.n_cells
        self.neighbor_layers = 1 if self.box.orthorhombic else 2
        self.cells: Dict[Tuple[int, int, int], List[int]] = {}
        self.atom_to_cell = np.empty((self.n_atoms, 3), dtype=int)
        for i in range(self.n_atoms):
            ci = self._cart_to_cell_index(self.coords[i])
            self.atom_to_cell[i] = ci
            ci_tuple = tuple(ci)
            if ci_tuple not in self.cells:
                self.cells[ci_tuple] = []
            self.cells[ci_tuple].append(i)

    def _cart_to_cell_index(self, coord: np.ndarray) -> Tuple[int, int, int]:
        bounds = self.box.get_bounds()
        lengths = bounds[:, 1] - bounds[:, 0]
        frac = (coord - bounds[:, 0]) / lengths
        frac = np.mod(frac, 1.0)
        ci = (frac * self.n_cells).astype(int)
        ci = np.minimum(ci, self.n_cells - 1)
        return tuple(ci)

    def _build_fractional_cells(self):
        self.frac_coords = self.box.to_fractional(self.coords)
        self.frac_coords = np.mod(self.frac_coords, 1.0)
        max_cell_size_frac = self.cutoff / np.min(self.box.lengths)
        n_cells_per_dim = int(np.ceil(1.0 / max_cell_size_frac))
        self.n_cells = np.array([n_cells_per_dim] * 3, dtype=int)
        self.cell_size_frac = 1.0 / self.n_cells
        self.neighbor_layers = 2 if not self.box.orthorhombic else 1
        self.cells: Dict[Tuple[int, int, int], List[int]] = {}
        self.atom_to_cell = np.empty((self.n_atoms, 3), dtype=int)
        for i in range(self.n_atoms):
            ci = (self.frac_coords[i] * self.n_cells).astype(int)
            ci = np.minimum(ci, self.n_cells - 1)
            ci = np.maximum(ci, 0)
            self.atom_to_cell[i] = ci
            ci_tuple = tuple(ci)
            if ci_tuple not in self.cells:
                self.cells[ci_tuple] = []
            self.cells[ci_tuple].append(i)

    def _get_neighbor_cell_indices(self, ci: Tuple[int, int, int]) -> List[Tuple[int, int, int]]:
        neighbors = []
        cx, cy, cz = ci
        layers = getattr(self, 'neighbor_layers', 1)
        for dx in range(-layers, layers + 1):
            for dy in range(-layers, layers + 1):
                for dz in range(-layers, layers + 1):
                    nci = (
                        (cx + dx) % self.n_cells[0],
                        (cy + dy) % self.n_cells[1],
                        (cz + dz) % self.n_cells[2]
                    )
                    neighbors.append(nci)
        return neighbors

    def _distance_sq(self, r1: np.ndarray, r2: np.ndarray) -> float:
        return self.box.minimum_image_distance_sq(r1, r2)

    def get_neighbors(self, atom_idx: int) -> List[int]:
        if atom_idx < 0 or atom_idx >= self.n_atoms:
            raise IndexError(f"Atom index {atom_idx} out of range")
        ci = tuple(self.atom_to_cell[atom_idx])
        neighbors = []
        cutoff_sq = self.cutoff * self.cutoff
        coord_i = self.coords[atom_idx]
        for nci in self._get_neighbor_cell_indices(ci):
            if nci not in self.cells:
                continue
            for j in self.cells[nci]:
                if j == atom_idx:
                    continue
                dist_sq = self._distance_sq(coord_i, self.coords[j])
                if dist_sq <= cutoff_sq:
                    neighbors.append(j)
        return neighbors

    def get_neighbors_with_distances(self, atom_idx: int) -> Tuple[List[int], np.ndarray]:
        if atom_idx < 0 or atom_idx >= self.n_atoms:
            raise IndexError(f"Atom index {atom_idx} out of range")
        ci = tuple(self.atom_to_cell[atom_idx])
        neighbors = []
        distances = []
        cutoff_sq = self.cutoff * self.cutoff
        coord_i = self.coords[atom_idx]
        for nci in self._get_neighbor_cell_indices(ci):
            if nci not in self.cells:
                continue
            for j in self.cells[nci]:
                if j == atom_idx:
                    continue
                dist_sq = self._distance_sq(coord_i, self.coords[j])
                if dist_sq <= cutoff_sq:
                    neighbors.append(j)
                    distances.append(np.sqrt(dist_sq))
        return neighbors, np.array(distances, dtype=np.float64)

    def get_all_neighbors(self) -> List[List[int]]:
        all_neighbors = [[] for _ in range(self.n_atoms)]
        cutoff_sq = self.cutoff * self.cutoff
        visited_pairs = set()
        for ci, atoms in self.cells.items():
            for nci in self._get_neighbor_cell_indices(ci):
                if nci not in self.cells:
                    continue
                pair_key = tuple(sorted([ci, nci]))
                if pair_key in visited_pairs and ci != nci:
                    continue
                visited_pairs.add(pair_key)
                atoms_j = self.cells[nci]
                is_same_cell = (ci == nci)
                for idx_i, i in enumerate(atoms):
                    coord_i = self.coords[i]
                    start_j = idx_i + 1 if is_same_cell else 0
                    for j in atoms_j[start_j:]:
                        dist_sq = self._distance_sq(coord_i, self.coords[j])
                        if dist_sq <= cutoff_sq:
                            all_neighbors[i].append(j)
                            all_neighbors[j].append(i)
        return all_neighbors

    def get_pair_distances(self, atom_idx: int, atom_indices: List[int]) -> np.ndarray:
        coord_i = self.coords[atom_idx]
        distances = np.empty(len(atom_indices), dtype=np.float64)
        for k, j in enumerate(atom_indices):
            distances[k] = np.sqrt(self._distance_sq(coord_i, self.coords[j]))
        return distances

    def count_neighbors_shell(self, atom_idx: int, r_inner: float, r_outer: float) -> int:
        neighbors = self.get_neighbors(atom_idx)
        count = 0
        coord_i = self.coords[atom_idx]
        r_inner_sq = r_inner * r_inner
        r_outer_sq = r_outer * r_outer
        for j in neighbors:
            dist_sq = self._distance_sq(coord_i, self.coords[j])
            if r_inner_sq <= dist_sq <= r_outer_sq:
                count += 1
        return count

    def get_neighbor_list_sorted(self, atom_idx: int, max_neighbors: int = None) -> Tuple[List[int], np.ndarray]:
        neighbors, distances = self.get_neighbors_with_distances(atom_idx)
        if len(neighbors) == 0:
            return [], np.array([], dtype=np.float64)
        sorted_idx = np.argsort(distances)
        neighbors = [neighbors[i] for i in sorted_idx]
        distances = distances[sorted_idx]
        if max_neighbors is not None and len(neighbors) > max_neighbors:
            neighbors = neighbors[:max_neighbors]
            distances = distances[:max_neighbors]
        return neighbors, distances
