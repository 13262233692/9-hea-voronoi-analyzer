import mmap
import struct
import os
from typing import Iterator, Tuple, List, Optional
import numpy as np
import polars as pl


class MDFrame:
    def __init__(self, timestep: int, box: np.ndarray, types: np.ndarray, coords: np.ndarray):
        self.timestep = timestep
        self.box = box
        self.types = types
        self.coords = coords
        self.n_atoms = len(types)


class MmapMDUnpacker:
    _MAGIC = b'MDAT'
    _HEADER_FORMAT = '=4sIq3d3d'
    _HEADER_SIZE = struct.calcsize(_HEADER_FORMAT)
    _ATOM_FORMAT = '=I3f'
    _ATOM_SIZE = struct.calcsize(_ATOM_FORMAT)

    def __init__(self, filepath: str):
        self.filepath = filepath
        self.file = None
        self.mm = None
        self.frame_offsets: List[int] = []
        self.n_frames = 0
        self._open()
        self._build_index()

    def _open(self):
        self.file = open(self.filepath, 'rb')
        self.mm = mmap.mmap(self.file.fileno(), 0, access=mmap.ACCESS_READ)

    def _build_index(self):
        pos = 0
        file_size = os.path.getsize(self.filepath)
        while pos < file_size:
            if pos + self._HEADER_SIZE > file_size:
                break
            header = self.mm[pos:pos + self._HEADER_SIZE]
            magic, n_atoms, timestep, bx, by, bz, lx, ly, lz = struct.unpack(self._HEADER_FORMAT, header)
            if magic != self._MAGIC:
                raise ValueError(f"Invalid magic bytes at offset {pos}")
            self.frame_offsets.append(pos)
            pos += self._HEADER_SIZE + n_atoms * self._ATOM_SIZE
        self.n_frames = len(self.frame_offsets)

    def __len__(self) -> int:
        return self.n_frames

    def _read_frame_at(self, offset: int) -> MDFrame:
        header = self.mm[offset:offset + self._HEADER_SIZE]
        magic, n_atoms, timestep, bx, by, bz, lx, ly, lz = struct.unpack(self._HEADER_FORMAT, header)
        box = np.array([[bx, bx + lx], [by, by + ly], [bz, bz + lz]], dtype=np.float64)
        atom_data_offset = offset + self._HEADER_SIZE
        atom_data = self.mm[atom_data_offset:atom_data_offset + n_atoms * self._ATOM_SIZE]
        types = np.empty(n_atoms, dtype=np.int32)
        coords = np.empty((n_atoms, 3), dtype=np.float32)
        for i in range(n_atoms):
            idx = i * self._ATOM_SIZE
            t, x, y, z = struct.unpack_from(self._ATOM_FORMAT, atom_data, idx)
            types[i] = t
            coords[i, 0] = x
            coords[i, 1] = y
            coords[i, 2] = z
        return MDFrame(timestep=timestep, box=box, types=types, coords=coords)

    def __getitem__(self, idx: int) -> MDFrame:
        if idx < 0:
            idx += self.n_frames
        if idx < 0 or idx >= self.n_frames:
            raise IndexError(f"Frame index {idx} out of range")
        return self._read_frame_at(self.frame_offsets[idx])

    def __iter__(self) -> Iterator[MDFrame]:
        for i in range(self.n_frames):
            yield self[i]

    def get_frame_as_polars(self, idx: int) -> pl.DataFrame:
        frame = self[idx]
        return pl.DataFrame({
            'type': frame.types,
            'x': frame.coords[:, 0],
            'y': frame.coords[:, 1],
            'z': frame.coords[:, 2]
        })

    def close(self):
        if self.mm:
            self.mm.close()
        if self.file:
            self.file.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


def create_test_md_file(filepath: str, n_frames: int = 10, n_atoms: int = 100000):
    element_types = [1, 2, 3, 4, 5]
    box_size = 50.0
    with open(filepath, 'wb') as f:
        for frame_idx in range(n_frames):
            timestep = frame_idx * 1000
            header = struct.pack(
                MmapMDUnpacker._HEADER_FORMAT,
                MmapMDUnpacker._MAGIC,
                n_atoms,
                timestep,
                0.0, 0.0, 0.0,
                box_size, box_size, box_size
            )
            f.write(header)
            types = np.random.choice(element_types, size=n_atoms)
            coords = np.random.uniform(0, box_size, size=(n_atoms, 3)).astype(np.float32)
            for i in range(n_atoms):
                atom_data = struct.pack(
                    MmapMDUnpacker._ATOM_FORMAT,
                    int(types[i]),
                    coords[i, 0], coords[i, 1], coords[i, 2]
                )
                f.write(atom_data)
    print(f"Created test file with {n_frames} frames, {n_atoms} atoms each")
