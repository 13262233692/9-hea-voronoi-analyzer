import os
import numpy as np
from typing import Dict, List, Optional, Tuple
from ..core.md_unpacker import MmapMDUnpacker, MDFrame
from ..core.cell_binning import CellBinning
from ..core.rdf_calculator import RDFCalculator
from ..core.voronoi_analyzer import VoronoiAnalyzer
from ..core.csro_calculator import CSROCalculator


class AnalysisService:
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._initialized = False
        return cls._instance

    def __init__(self):
        if self._initialized:
            return
        self._initialized = True
        self.unpacker: Optional[MmapMDUnpacker] = None
        self.current_file: Optional[str] = None
        self.rdf_calc = RDFCalculator(r_min=0.0, r_max=6.0, n_bins=100)
        self.voronoi_analyzer = VoronoiAnalyzer(cutoff=5.0)
        self.csro_calc = CSROCalculator(cutoff=3.5)
        self.element_map = {1: 'Al', 2: 'Co', 3: 'Cr', 4: 'Fe', 5: 'Ni'}

    def load_data(self, filepath: str) -> Dict:
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"File not found: {filepath}")
        if self.unpacker:
            self.unpacker.close()
        self.unpacker = MmapMDUnpacker(filepath)
        self.current_file = filepath
        return {
            'n_frames': self.unpacker.n_frames,
            'file': filepath,
            'n_atoms': self.unpacker[0].n_atoms if self.unpacker.n_frames > 0 else 0
        }

    def get_frame(self, idx: int) -> MDFrame:
        if not self.unpacker:
            raise RuntimeError("No data loaded")
        return self.unpacker[idx]

    def get_frame_info(self, idx: int) -> Dict:
        frame = self.get_frame(idx)
        unique_types, counts = np.unique(frame.types, return_counts=True)
        type_info = {}
        for t, c in zip(unique_types, counts):
            elem = self.element_map.get(int(t), f'Type{t}')
            type_info[elem] = int(c)
        return {
            'timestep': frame.timestep,
            'n_atoms': frame.n_atoms,
            'box': frame.box.tolist(),
            'type_counts': type_info
        }

    def calculate_rdf(self, frame_idx: int, r_max: float = 6.0, n_bins: int = 100) -> Dict:
        frame = self.get_frame(frame_idx)
        calc = RDFCalculator(r_max=r_max, n_bins=n_bins)
        return calc.calculate(frame.coords, frame.box, frame.types)

    def calculate_rdf_average(self, start_frame: int, end_frame: int,
                              r_max: float = 6.0, n_bins: int = 100) -> Dict:
        frames = [self.get_frame(i) for i in range(start_frame, end_frame + 1)]
        calc = RDFCalculator(r_max=r_max, n_bins=n_bins)
        return calc.calculate_batch(frames)

    def analyze_voronoi(self, frame_idx: int) -> Dict:
        frame = self.get_frame(frame_idx)
        result = self.voronoi_analyzer.analyze_frame(frame.coords, frame.box, frame.types)
        vor_indices = np.array(result['voronoi_indices'])
        polyhedron_types = [
            self.voronoi_analyzer.get_polyhedron_type(idx)
            for idx in vor_indices
        ]
        result['polyhedron_types'] = polyhedron_types
        result['classifications'] = [
            self.voronoi_analyzer.classify_polyhedron(idx)
            for idx in vor_indices
        ]
        return result

    def get_voronoi_evolution(self, start_frame: int, end_frame: int,
                               polyhedron_types: List[str] = None) -> Dict:
        frames = [self.get_frame(i) for i in range(start_frame, end_frame + 1)]
        return self.voronoi_analyzer.get_evolution_matrix(frames, polyhedron_types)

    def calculate_csro(self, frame_idx: int, type1: int, type2: int) -> Dict:
        frame = self.get_frame(frame_idx)
        result = self.csro_calc.calculate_csro(
            frame.coords, frame.box, frame.types, type1, type2
        )
        result['element1'] = self.element_map.get(type1, f'Type{type1}')
        result['element2'] = self.element_map.get(type2, f'Type{type2}')
        return result

    def calculate_csro_all_pairs(self, frame_idx: int) -> Dict:
        frame = self.get_frame(frame_idx)
        result = self.csro_calc.analyze_all_pairs(
            frame.coords, frame.box, frame.types
        )
        formatted = {}
        for key, val in result.items():
            t1, t2 = key.split('-')
            elem1 = self.element_map.get(int(t1), f'Type{t1}')
            elem2 = self.element_map.get(int(t2), f'Type{t2}')
            formatted[f'{elem1}-{elem2}'] = val
        return formatted

    def find_frame_by_timestep(self, target_timestep: int) -> int:
        if not self.unpacker:
            raise RuntimeError("No data loaded")
        low, high = 0, self.unpacker.n_frames - 1
        while low <= high:
            mid = (low + high) // 2
            mid_ts = self.unpacker[mid].timestep
            if mid_ts == target_timestep:
                return mid
            elif mid_ts < target_timestep:
                low = mid + 1
            else:
                high = mid - 1
        if high < 0:
            return 0
        if low >= self.unpacker.n_frames:
            return self.unpacker.n_frames - 1
        high_ts = self.unpacker[high].timestep
        low_ts = self.unpacker[low].timestep
        if abs(target_timestep - high_ts) < abs(target_timestep - low_ts):
            return high
        return low

    def get_evolution_stream(self, start_frame: int = 0, end_frame: int = None,
                              polyhedron_types: List[str] = None):
        if not self.unpacker:
            raise RuntimeError("No data loaded")
        if end_frame is None:
            end_frame = self.unpacker.n_frames - 1
        for i in range(start_frame, end_frame + 1):
            frame = self.get_frame(i)
            result = self.voronoi_analyzer.analyze_frame(frame.coords, frame.box, frame.types)
            vor_indices = np.array(result['voronoi_indices'])
            type_strings = [
                self.voronoi_analyzer.get_polyhedron_type(idx)
                for idx in vor_indices
            ]
            from collections import Counter
            counts = dict(Counter(type_strings))
            if polyhedron_types:
                filtered = {}
                for pt in polyhedron_types:
                    filtered[pt] = counts.get(pt, 0)
                counts = filtered
            yield {
                'frame_idx': i,
                'timestep': frame.timestep,
                'time_ps': frame.timestep / 1000.0,
                'counts': counts
            }

    def get_atom_neighbors(self, frame_idx: int, atom_idx: int, cutoff: float = 3.5) -> Dict:
        frame = self.get_frame(frame_idx)
        cb = CellBinning(frame.coords, frame.box, cutoff)
        neighbors, distances = cb.get_neighbors_with_distances(atom_idx)
        neighbor_types = [int(frame.types[j]) for j in neighbors]
        neighbor_elements = [self.element_map.get(t, f'Type{t}') for t in neighbor_types]
        return {
            'atom_idx': atom_idx,
            'atom_type': int(frame.types[atom_idx]),
            'atom_element': self.element_map.get(int(frame.types[atom_idx]), f'Type{frame.types[atom_idx]}'),
            'coord': frame.coords[atom_idx].tolist(),
            'neighbors': neighbors,
            'neighbor_types': neighbor_types,
            'neighbor_elements': neighbor_elements,
            'distances': distances.tolist()
        }

    def close(self):
        if self.unpacker:
            self.unpacker.close()
            self.unpacker = None
