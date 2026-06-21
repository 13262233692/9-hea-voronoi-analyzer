import numpy as np
from typing import List, Tuple, Dict
from dataclasses import dataclass


@dataclass
class Cell:
    cell_index: Tuple[int, int, int]
    atom_indices: List[int]


class CellBinning:
    def __init__(self, coords: np.ndarray, box: np.ndarray, cutoff: float):
        self.coords = coords
        self.box = box
        self.cutoff = cutoff
        self.n_atoms = coords.shape[0]
        self.box_lengths = box[:, 1] - box[:, 0]
        self.n_cells = np.maximum(1, (self.box_lengths / cutoff).astype(int))
        self.cell_size = self.box_lengths / self.n_cells
        self._build_cells()

    def _get_cell_index(self, coord: np.ndarray) -> Tuple[int, int, int]:
        frac = (coord - self.box[:, 0]) / self.box_lengths
        frac = np.mod(frac, 1.0)
        cell_idx = (frac * self.n_cells).astype(int)
        cell_idx = np.minimum(cell_idx, self.n_cells - 1)
        return tuple(cell_idx)

    def _build_cells(self):
        self.cells: Dict[Tuple[int, int, int], List[int]] = {}
        self.atom_to_cell = np.empty((self.n_atoms, 3), dtype=int)
        for i in range(self.n_atoms):
            ci = self._get_cell_index(self.coords[i])
            self.atom_to_cell[i] = ci
            if ci not in self.cells:
                self.cells[ci] = []
            self.cells[ci].append(i)

    def _get_neighbor_cells(self, ci: Tuple[int, int, int]) -> List[Tuple[int, int, int]]:
        neighbors = []
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                for dz in [-1, 0, 1]:
                    nci = (
                        (ci[0] + dx) % self.n_cells[0],
                        (ci[1] + dy) % self.n_cells[1],
                        (ci[2] + dz) % self.n_cells[2]
                    )
                    neighbors.append(nci)
        return neighbors

    def get_neighbors(self, atom_idx: int) -> List[int]:
        ci = tuple(self.atom_to_cell[atom_idx])
        neighbors = []
        coord_i = self.coords[atom_idx]
        for nci in self._get_neighbor_cells(ci):
            if nci not in self.cells:
                continue
            for j in self.cells[nci]:
                if j == atom_idx:
                    continue
                dr = self.coords[j] - coord_i
                dr = dr - self.box_lengths * np.round(dr / self.box_lengths)
                dist_sq = np.sum(dr * dr)
                if dist_sq <= self.cutoff * self.cutoff:
                    neighbors.append(j)
        return neighbors

    def get_all_neighbors(self) -> List[List[int]]:
        all_neighbors = [[] for _ in range(self.n_atoms)]
        for ci, atoms in self.cells.items():
            for nci in self._get_neighbor_cells(ci):
                if nci not in self.cells:
                    continue
                if nci < ci:
                    continue
                atoms_j = self.cells[nci]
                for i in atoms:
                    coord_i = self.coords[i]
                    for j in atoms_j:
                        if i == j:
                            continue
                        dr = self.coords[j] - coord_i
                        dr = dr - self.box_lengths * np.round(dr / self.box_lengths)
                        dist_sq = np.sum(dr * dr)
                        if dist_sq <= self.cutoff * self.cutoff:
                            all_neighbors[i].append(j)
                            all_neighbors[j].append(i)
        return all_neighbors

    def get_neighbors_with_distances(self, atom_idx: int) -> Tuple[List[int], np.ndarray]:
        neighbors = []
        distances = []
        ci = tuple(self.atom_to_cell[atom_idx])
        coord_i = self.coords[atom_idx]
        for nci in self._get_neighbor_cells(ci):
            if nci not in self.cells:
                continue
            for j in self.cells[nci]:
                if j == atom_idx:
                    continue
                dr = self.coords[j] - coord_i
                dr = dr - self.box_lengths * np.round(dr / self.box_lengths)
                dist = np.sqrt(np.sum(dr * dr))
                if dist <= self.cutoff:
                    neighbors.append(j)
                    distances.append(dist)
        return neighbors, np.array(distances, dtype=np.float64)
