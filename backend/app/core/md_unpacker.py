import mmap
import struct
import os
import threading
from typing import Iterator, Tuple, List, Optional
from contextlib import contextmanager
import numpy as np
import polars as pl


class TriclinicBox:
    def __init__(self, a_vec: np.ndarray, b_vec: np.ndarray, c_vec: np.ndarray):
        self.a_vec = np.asarray(a_vec, dtype=np.float64)
        self.b_vec = np.asarray(b_vec, dtype=np.float64)
        self.c_vec = np.asarray(c_vec, dtype=np.float64)
        self.h_matrix = np.column_stack([self.a_vec, self.b_vec, self.c_vec])
        self.inv_h_matrix = np.linalg.inv(self.h_matrix)
        self.volume = np.dot(self.a_vec, np.cross(self.b_vec, self.c_vec))
        self.lengths = np.array([
            np.linalg.norm(self.a_vec),
            np.linalg.norm(self.b_vec),
            np.linalg.norm(self.c_vec)
        ], dtype=np.float64)
        self.orthorhombic = self._is_orthorhombic()

    def _is_orthorhombic(self, tol: float = 1e-10) -> bool:
        h = self.h_matrix
        return (
            abs(h[0, 1]) < tol and abs(h[0, 2]) < tol and
            abs(h[1, 0]) < tol and abs(h[1, 2]) < tol and
            abs(h[2, 0]) < tol and abs(h[2, 1]) < tol
        )

    def to_fractional(self, coords: np.ndarray) -> np.ndarray:
        return coords @ self.inv_h_matrix.T

    def to_cartesian(self, frac_coords: np.ndarray) -> np.ndarray:
        return frac_coords @ self.h_matrix.T

    def minimum_image_vector(self, dr: np.ndarray) -> np.ndarray:
        if self.orthorhombic:
            return dr - self.lengths * np.round(dr / self.lengths)
        frac_dr = dr @ self.inv_h_matrix.T
        frac_dr -= np.round(frac_dr)
        return frac_dr @ self.h_matrix.T

    def minimum_image_distance_sq(self, r1: np.ndarray, r2: np.ndarray) -> float:
        dr = r2 - r1
        dr_min = self.minimum_image_vector(dr)
        return np.sum(dr_min * dr_min)

    def minimum_image_distance(self, r1: np.ndarray, r2: np.ndarray) -> float:
        return np.sqrt(self.minimum_image_distance_sq(r1, r2))

    def get_bounds(self) -> np.ndarray:
        if self.orthorhombic:
            return np.array([
                [0.0, self.a_vec[0]],
                [0.0, self.b_vec[1]],
                [0.0, self.c_vec[2]]
            ], dtype=np.float64)
        else:
            corners = np.array([
                [0, 0, 0],
                self.a_vec,
                self.b_vec,
                self.c_vec,
                self.a_vec + self.b_vec,
                self.a_vec + self.c_vec,
                self.b_vec + self.c_vec,
                self.a_vec + self.b_vec + self.c_vec
            ])
            return np.array([
                [corners[:, 0].min(), corners[:, 0].max()],
                [corners[:, 1].min(), corners[:, 1].max()],
                [corners[:, 2].min(), corners[:, 2].max()]
            ], dtype=np.float64)

    @classmethod
    def from_orthorhombic(cls, lx: float, ly: float, lz: float):
        return cls(
            a_vec=np.array([lx, 0.0, 0.0]),
            b_vec=np.array([0.0, ly, 0.0]),
            c_vec=np.array([0.0, 0.0, lz])
        )

    @classmethod
    def from_lammps(cls, xlo: float, xhi: float, ylo: float, yhi: float,
                    zlo: float, zhi: float, xy: float = 0.0, xz: float = 0.0,
                    yz: float = 0.0):
        a_vec = np.array([xhi - xlo, 0.0, 0.0])
        b_vec = np.array([xy, yhi - ylo, 0.0])
        c_vec = np.array([xz, yz, zhi - zlo])
        return cls(a_vec, b_vec, c_vec)


class MDFrame:
    def __init__(self, timestep: int, box: TriclinicBox, types: np.ndarray, coords: np.ndarray):
        self.timestep = timestep
        self.box = box
        self.types = types
        self.coords = coords
        self.n_atoms = len(types)


class MmapMDUnpacker:
    _MAGIC = b'MDAT'
    _VERSION = 2
    _FLAG_ORTHO = 0x00
    _FLAG_TRICLINIC = 0x01
    _HEADER_FORMAT = '=4sHBIq'
    _HEADER_SIZE = struct.calcsize(_HEADER_FORMAT)
    _ATOM_FORMAT = '=I3f'
    _ATOM_SIZE = struct.calcsize(_ATOM_FORMAT)
    _ORTHO_BOX_FORMAT = '=6d'
    _ORTHO_BOX_SIZE = struct.calcsize(_ORTHO_BOX_FORMAT)
    _TRICLINIC_BOX_FORMAT = '=9d'
    _TRICLINIC_BOX_SIZE = struct.calcsize(_TRICLINIC_BOX_FORMAT)
    _instance_count = 0
    _instance_lock = threading.Lock()
    _max_instances = 8

    def __init__(self, filepath: str, validate: bool = True):
        self.filepath = os.path.abspath(filepath)
        self._lock = threading.RLock()
        self._file = None
        self._mm = None
        self._closed = False
        self.frame_offsets: List[int] = []
        self.n_frames = 0
        self._box_cache = {}
        with MmapMDUnpacker._instance_lock:
            MmapMDUnpacker._instance_count += 1
            if MmapMDUnpacker._instance_count > MmapMDUnpacker._max_instances:
                raise RuntimeError(
                    f"Too many MmapMDUnpacker instances ({MmapMDUnpacker._instance_count}). "
                    f"Memory leak detected - please ensure proper cleanup."
                )
        try:
            self._open()
            if validate:
                self._validate_file()
            self._build_index()
        except Exception:
            self.close()
            raise

    def _open(self):
        with self._lock:
            if self._file is not None or self._mm is not None:
                self._close_resources()
            self._file = open(self.filepath, 'rb')
            self._mm = mmap.mmap(self._file.fileno(), 0, access=mmap.ACCESS_READ)

    def _validate_file(self):
        if self._HEADER_SIZE > len(self._mm):
            raise ValueError("File too small to contain a valid header")
        magic, version, flags, n_atoms, timestep = struct.unpack_from(
            self._HEADER_FORMAT, self._mm, 0
        )
        if magic != self._MAGIC:
            raise ValueError(f"Invalid magic bytes: expected {self._MAGIC!r}, got {magic!r}")
        if version > self._VERSION:
            raise ValueError(f"Unsupported file version: {version}")

    def _get_frame_header_size(self, flags: int) -> int:
        if flags & self._FLAG_TRICLINIC:
            return self._HEADER_SIZE + self._TRICLINIC_BOX_SIZE
        else:
            return self._HEADER_SIZE + self._ORTHO_BOX_SIZE

    def _read_box(self, offset: int, flags: int) -> Tuple[TriclinicBox, int]:
        hdr_size = self._HEADER_SIZE
        if flags & self._FLAG_TRICLINIC:
            box_data = struct.unpack_from(self._TRICLINIC_BOX_FORMAT, self._mm, offset + hdr_size)
            ax, ay, az, bx, by, bz, cx, cy, cz = box_data
            box = TriclinicBox(
                a_vec=np.array([ax, ay, az]),
                b_vec=np.array([bx, by, bz]),
                c_vec=np.array([cx, cy, cz])
            )
            box_size = self._TRICLINIC_BOX_SIZE
        else:
            box_data = struct.unpack_from(self._ORTHO_BOX_FORMAT, self._mm, offset + hdr_size)
            bx, by, bz, lx, ly, lz = box_data
            box = TriclinicBox.from_orthorhombic(lx, ly, lz)
            box_size = self._ORTHO_BOX_SIZE
        return box, box_size

    def _build_index(self):
        pos = 0
        file_size = len(self._mm)
        offsets = []
        while pos < file_size:
            if pos + self._HEADER_SIZE > file_size:
                break
            magic, version, flags, n_atoms, timestep = struct.unpack_from(
                self._HEADER_FORMAT, self._mm, pos
            )
            if magic != self._MAGIC:
                raise ValueError(f"Invalid magic bytes at offset {pos}: file may be corrupted")
            frame_hdr_size = self._get_frame_header_size(flags)
            if pos + frame_hdr_size > file_size:
                break
            atom_data_size = n_atoms * self._ATOM_SIZE
            next_pos = pos + frame_hdr_size + atom_data_size
            if next_pos > file_size + 1:
                break
            offsets.append(pos)
            pos = next_pos
        self.frame_offsets = offsets
        self.n_frames = len(offsets)

    def __len__(self) -> int:
        return self.n_frames

    def _read_frame_at(self, offset: int) -> MDFrame:
        with self._lock:
            if self._closed:
                raise RuntimeError("Cannot read from closed unpacker")
            magic, version, flags, n_atoms, timestep = struct.unpack_from(
                self._HEADER_FORMAT, self._mm, offset
            )
            box, box_size = self._read_box(offset, flags)
            frame_hdr_size = self._HEADER_SIZE + box_size
            atom_data_offset = offset + frame_hdr_size
            types = np.empty(n_atoms, dtype=np.int32)
            coords = np.empty((n_atoms, 3), dtype=np.float32)
            for i in range(n_atoms):
                idx = atom_data_offset + i * self._ATOM_SIZE
                t, x, y, z = struct.unpack_from(self._ATOM_FORMAT, self._mm, idx)
                types[i] = t
                coords[i, 0] = x
                coords[i, 1] = y
                coords[i, 2] = z
            return MDFrame(timestep=timestep, box=box, types=types, coords=coords)

    def __getitem__(self, idx: int) -> MDFrame:
        if self._closed:
            raise RuntimeError("Cannot access closed unpacker")
        if idx < 0:
            idx += self.n_frames
        if idx < 0 or idx >= self.n_frames:
            raise IndexError(f"Frame index {idx} out of range (0-{self.n_frames - 1})")
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

    def _close_resources(self):
        if self._mm is not None:
            try:
                self._mm.close()
            except Exception:
                pass
            self._mm = None
        if self._file is not None:
            try:
                self._file.close()
            except Exception:
                pass
            self._file = None

    def close(self):
        with self._lock:
            if self._closed:
                return
            self._closed = True
            self._close_resources()
            self.frame_offsets.clear()
            self._box_cache.clear()
            with MmapMDUnpacker._instance_lock:
                MmapMDUnpacker._instance_count -= 1

    def __del__(self):
        if not self._closed:
            try:
                self.close()
            except Exception:
                pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def is_triclinic(self) -> bool:
        if self.n_frames == 0:
            return False
        magic, version, flags, n_atoms, timestep = struct.unpack_from(
            self._HEADER_FORMAT, self._mm, self.frame_offsets[0]
        )
        return bool(flags & self._FLAG_TRICLINIC)

    @property
    def closed(self) -> bool:
        return self._closed


def create_test_md_file(filepath: str, n_frames: int = 10, n_atoms: int = 100000,
                        triclinic: bool = False, tilt_factor: float = 0.3):
    element_types = [1, 2, 3, 4, 5]
    os.makedirs(os.path.dirname(os.path.abspath(filepath)), exist_ok=True)

    if triclinic:
        a = 50.0
        b = 48.0
        c = 52.0
        a_vec = np.array([a, 0.0, 0.0])
        b_vec = np.array([b * tilt_factor, b, 0.0])
        c_vec = np.array([c * tilt_factor * 0.5, c * tilt_factor * 0.3, c])
        flags = MmapMDUnpacker._FLAG_TRICLINIC
        box = TriclinicBox(a_vec, b_vec, c_vec)
    else:
        box_size = 50.0
        a_vec = np.array([box_size, 0.0, 0.0])
        b_vec = np.array([0.0, box_size, 0.0])
        c_vec = np.array([0.0, 0.0, box_size])
        flags = MmapMDUnpacker._FLAG_ORTHO
        box = TriclinicBox(a_vec, b_vec, c_vec)

    with open(filepath, 'wb') as f:
        for frame_idx in range(n_frames):
            timestep = frame_idx * 1000
            header = struct.pack(
                MmapMDUnpacker._HEADER_FORMAT,
                MmapMDUnpacker._MAGIC,
                MmapMDUnpacker._VERSION,
                flags,
                n_atoms,
                timestep
            )
            f.write(header)

            if triclinic:
                box_data = struct.pack(
                    MmapMDUnpacker._TRICLINIC_BOX_FORMAT,
                    a_vec[0], a_vec[1], a_vec[2],
                    b_vec[0], b_vec[1], b_vec[2],
                    c_vec[0], c_vec[1], c_vec[2]
                )
            else:
                box_data = struct.pack(
                    MmapMDUnpacker._ORTHO_BOX_FORMAT,
                    0.0, 0.0, 0.0,
                    a_vec[0], b_vec[1], c_vec[2]
                )
            f.write(box_data)

            types = np.random.choice(element_types, size=n_atoms)
            frac_coords = np.random.random((n_atoms, 3)).astype(np.float64)
            coords = (frac_coords @ box.h_matrix.T).astype(np.float32)

            for i in range(n_atoms):
                atom_data = struct.pack(
                    MmapMDUnpacker._ATOM_FORMAT,
                    int(types[i]),
                    coords[i, 0], coords[i, 1], coords[i, 2]
                )
                f.write(atom_data)

    print(f"Created test file: {filepath}")
    print(f"  Frames: {n_frames}")
    print(f"  Atoms/frame: {n_atoms}")
    print(f"  Cell type: {'triclinic' if triclinic else 'orthorhombic'}")
    if triclinic:
        print(f"  a = {a_vec}")
        print(f"  b = {b_vec}")
        print(f"  c = {c_vec}")
        print(f"  Volume: {box.volume:.2f}")
