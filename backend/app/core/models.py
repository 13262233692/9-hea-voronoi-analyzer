from pydantic import BaseModel
from typing import List, Dict, Optional, Union


class LoadDataRequest(BaseModel):
    filepath: str


class FrameInfoResponse(BaseModel):
    timestep: int
    n_atoms: int
    box: List[List[float]]
    type_counts: Dict[str, int]


class RDFResponse(BaseModel):
    r: List[float]
    gr: List[float]
    gr_raw: Optional[List[float]] = None
    partial_gr: Optional[Dict[str, List[float]]] = None


class VoronoiResponse(BaseModel):
    voronoi_indices: List[List[int]]
    counts: Dict[str, int]
    n_atoms: int
    by_type: Optional[Dict[str, Dict[str, int]]] = None
    polyhedron_types: Optional[List[str]] = None
    classifications: Optional[List[str]] = None


class VoronoiEvolutionResponse(BaseModel):
    timesteps: List[int]
    evolution: Dict[str, List[int]]
    n_frames: int


class CSROResponse(BaseModel):
    S_ij: float
    alpha: float
    c1: float
    c2: float
    Z_avg: float
    element1: str
    element2: str
    bond_counts: Dict[str, int]
    type_pairs: Dict[str, float]
    N_AA: float
    N_BB: float
    N_AB: float


class EvolutionStreamData(BaseModel):
    frame_idx: int
    timestep: int
    time_ps: float
    counts: Dict[str, int]


class NeighborsResponse(BaseModel):
    atom_idx: int
    atom_type: int
    atom_element: str
    coord: List[float]
    neighbors: List[int]
    neighbor_types: List[int]
    neighbor_elements: List[str]
    distances: List[float]


class CSROQueryRequest(BaseModel):
    timestep: int
    type1: int = 2
    type2: int = 3
