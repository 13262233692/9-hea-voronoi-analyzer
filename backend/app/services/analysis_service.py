import os
import gc
import numpy as np
from typing import Dict, List, Optional, Generator
import weakref
from ..core.md_unpacker import MmapMDUnpacker, MDFrame, TriclinicBox, create_test_md_file
from ..core.cell_binning import CellBinning
from ..core.rdf_calculator import RDFCalculator
from ..core.voronoi_analyzer import VoronoiAnalyzer
from ..core.csro_calculator import CSROCalculator
from ..core.nonaffine_analyzer import NonAffineAnalyzer


class AnalysisService:
    _instance = None
    _instance_ref = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._initialized = False
            cls._instance_ref = weakref.ref(cls._instance)
        return cls._instance

    def __init__(self):
        if self._initialized:
            return
        self._initialized = True
        self._unpacker: Optional[MmapMDUnpacker] = None
        self.current_file: Optional[str] = None
        self.rdf_calc = RDFCalculator(r_min=0.0, r_max=6.0, n_bins=100)
        self.voronoi_analyzer = VoronoiAnalyzer(cutoff=5.0)
        self.csro_calc = CSROCalculator(cutoff=3.5)
        self.nonaffine_analyzer = NonAffineAnalyzer(cutoff=4.5, yield_threshold=0.08)
        self.element_map = {1: 'Al', 2: 'Co', 3: 'Cr', 4: 'Fe', 5: 'Ni'}
        self._frame_cache = {}
        self._max_cache_frames = 5
        self._yield_events = []
        self._monitored_pair = None

    @property
    def unpacker(self) -> Optional[MmapMDUnpacker]:
        return self._unpacker

    def load_data(self, filepath: str) -> Dict:
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"File not found: {filepath}")
        self._safe_close_unpacker()
        self._frame_cache.clear()
        gc.collect()
        self._unpacker = MmapMDUnpacker(filepath)
        self.current_file = filepath
        frame0 = self._unpacker[0] if self._unpacker.n_frames > 0 else None
        return {
            'n_frames': self._unpacker.n_frames,
            'file': filepath,
            'n_atoms': frame0.n_atoms if frame0 else 0,
            'is_triclinic': self._unpacker.is_triclinic() if frame0 else False,
            'volume': float(frame0.box.volume) if frame0 else 0.0,
            'cell_matrix': frame0.box.h_matrix.tolist() if frame0 else None
        }

    def _safe_close_unpacker(self):
        if self._unpacker is not None:
            try:
                self._unpacker.close()
            except Exception:
                pass
            self._unpacker = None
            gc.collect()

    def get_frame(self, idx: int) -> MDFrame:
        if self._unpacker is None:
            raise RuntimeError("No data loaded")
        if idx in self._frame_cache:
            return self._frame_cache[idx]
        frame = self._unpacker[idx]
        if len(self._frame_cache) >= self._max_cache_frames:
            oldest_key = next(iter(self._frame_cache))
            del self._frame_cache[oldest_key]
        self._frame_cache[idx] = frame
        return frame

    def get_frame_info(self, idx: int) -> Dict:
        frame = self.get_frame(idx)
        unique_types, counts = np.unique(frame.types, return_counts=True)
        type_info = {}
        for t, c in zip(unique_types, counts):
            elem = self.element_map.get(int(t), f'Type{t}')
            type_info[elem] = int(c)
        bounds = frame.box.get_bounds()
        return {
            'timestep': frame.timestep,
            'n_atoms': frame.n_atoms,
            'box': bounds.tolist(),
            'box_vectors': {
                'a': frame.box.a_vec.tolist(),
                'b': frame.box.b_vec.tolist(),
                'c': frame.box.c_vec.tolist()
            },
            'is_triclinic': not frame.box.orthorhombic,
            'volume': float(frame.box.volume),
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
        if self._unpacker is None:
            raise RuntimeError("No data loaded")
        low, high = 0, self._unpacker.n_frames - 1
        while low <= high:
            mid = (low + high) // 2
            mid_ts = self._unpacker[mid].timestep
            if mid_ts == target_timestep:
                return mid
            elif mid_ts < target_timestep:
                low = mid + 1
            else:
                high = mid - 1
        if high < 0:
            return 0
        if low >= self._unpacker.n_frames:
            return self._unpacker.n_frames - 1
        high_ts = self._unpacker[high].timestep
        low_ts = self._unpacker[low].timestep
        if abs(target_timestep - high_ts) < abs(target_timestep - low_ts):
            return high
        return low

    def get_evolution_stream(self, start_frame: int = 0, end_frame: int = None,
                              polyhedron_types: List[str] = None):
        if self._unpacker is None:
            raise RuntimeError("No data loaded")
        if end_frame is None:
            end_frame = self._unpacker.n_frames - 1
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
                'counts': counts,
                'is_triclinic': not frame.box.orthorhombic
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
            'atom_element': self.element_map.get(
                int(frame.types[atom_idx]),
                f'Type{frame.types[atom_idx]}'
            ),
            'coord': frame.coords[atom_idx].tolist(),
            'neighbors': neighbors,
            'neighbor_types': neighbor_types,
            'neighbor_elements': neighbor_elements,
            'distances': distances.tolist(),
            'n_neighbors': len(neighbors)
        }

    def analyze_nonaffine(self, frame_idx_ref: int, frame_idx_def: int) -> Dict:
        frame_ref = self.get_frame(frame_idx_ref)
        frame_def = self.get_frame(frame_idx_def)

        result = self.nonaffine_analyzer.analyze_deformation(
            frame_ref.coords, frame_def.coords, frame_def.box, frame_def.types
        )

        clusters = self.nonaffine_analyzer.detect_shear_band_clusters(
            frame_def.coords, result.shear_band_atoms, frame_def.box
        )

        dissipation = self.nonaffine_analyzer.compute_dissipation_curve(
            frame_ref.coords, frame_def.coords, frame_def.box
        )

        is_yielding = result.n_yielding > 0
        if is_yielding:
            event = {
                'frame_ref': frame_idx_ref,
                'frame_def': frame_idx_def,
                'timestep_ref': frame_ref.timestep,
                'timestep_def': frame_def.timestep,
                'n_yielding': result.n_yielding,
                'yield_ratio': result.n_yielding / result.total_atoms,
                'dissipated_energy': result.dissipated_energy,
                'clusters': clusters[:5],
                'top_cluster_center': clusters[0]['center'] if clusters else None,
                'top_cluster_size': clusters[0]['n_atoms'] if clusters else 0,
            }
            self._yield_events.append(event)

        return {
            'frame_ref': frame_idx_ref,
            'frame_def': frame_idx_def,
            'timestep_ref': frame_ref.timestep,
            'timestep_def': frame_def.timestep,
            'd2_min': result.d2_min.tolist(),
            'threshold': result.threshold,
            'n_yielding': result.n_yielding,
            'yield_ratio': result.n_yielding / result.total_atoms,
            'total_atoms': result.total_atoms,
            'dissipated_energy': result.dissipated_energy,
            'shear_band_atoms': result.shear_band_atoms.tolist(),
            'strain_magnitudes': result.atom_strain_magnitudes.tolist(),
            'clusters': clusters,
            'dissipation_curve': dissipation,
            'is_yielding': is_yielding,
        }

    def get_nonaffine_stream(self, start_frame: int = 0, end_frame: int = None,
                              yield_threshold: float = None) -> Generator[Dict, None, None]:
        if self._unpacker is None:
            raise RuntimeError("No data loaded")
        if end_frame is None:
            end_frame = self._unpacker.n_frames - 1
        if yield_threshold is not None:
            self.nonaffine_analyzer.yield_threshold = yield_threshold

        self._yield_events = []

        for i in range(start_frame, end_frame):
            result = self.analyze_nonaffine(i, i + 1)
            yield {
                'type': 'nonaffine',
                'frame_idx': i + 1,
                'timestep': result['timestep_def'],
                'n_yielding': result['n_yielding'],
                'yield_ratio': result['yield_ratio'],
                'dissipated_energy': result['dissipated_energy'],
                'is_yielding': result['is_yielding'],
                'clusters': result['clusters'][:3],
                'shear_band_atoms_count': result['n_yielding'],
            }

            if result['is_yielding']:
                yield {
                    'type': 'yield_event',
                    'frame_idx': i + 1,
                    'timestep': result['timestep_def'],
                    'severity': 'critical' if result['yield_ratio'] > 0.05 else 'warning',
                    'n_yielding': result['n_yielding'],
                    'yield_ratio': result['yield_ratio'],
                    'dissipated_energy': result['dissipated_energy'],
                    'top_cluster_center': result['clusters'][0]['center'] if result['clusters'] else None,
                    'top_cluster_size': result['clusters'][0]['n_atoms'] if result['clusters'] else 0,
                    'shear_band_coords': [
                        a['coords'] for a in result['clusters'][:1]
                    ],
                    'dissipation_curve': result['dissipation_curve'],
                    'd2_min_top': sorted(result['d2_min'], reverse=True)[:10],
                }

    def get_yield_events(self) -> List[Dict]:
        return list(self._yield_events)

    def clear_cache(self):
        self._frame_cache.clear()
        gc.collect()

    def close(self):
        self._safe_close_unpacker()
        self._frame_cache.clear()
        gc.collect()

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass
