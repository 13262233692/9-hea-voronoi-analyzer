import numpy as np
import time
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from app.core.md_unpacker import create_test_md_file, MmapMDUnpacker
from app.core.cell_binning import CellBinning
from app.core.rdf_calculator import RDFCalculator
from app.core.voronoi_analyzer import VoronoiAnalyzer
from app.core.csro_calculator import CSROCalculator
from app.services.analysis_service import AnalysisService


def test_md_unpacker():
    print("=" * 60)
    print("TEST 1: MD Unpacker (mmap-based streaming reader)")
    print("=" * 60)

    test_file = r"d:\SOLO-0621-1\9-hea-voronoi-analyzer\backend\data\test_heau.md"
    os.makedirs(os.path.dirname(test_file), exist_ok=True)

    n_frames = 5
    n_atoms = 10000

    print(f"Generating test data: {n_frames} frames, {n_atoms} atoms each...")
    t0 = time.time()
    create_test_md_file(test_file, n_frames, n_atoms)
    t1 = time.time()
    print(f"  Test data generated in {t1 - t0:.3f}s")

    print("Opening with mmap unpacker...")
    t0 = time.time()
    with MmapMDUnpacker(test_file) as unpacker:
        t1 = time.time()
        print(f"  Index built in {t1 - t0:.3f}s, {unpacker.n_frames} frames found")

        for i in range(min(3, unpacker.n_frames)):
            t0 = time.time()
            frame = unpacker[i]
            t1 = time.time()
            print(f"  Frame {i}: timestep={frame.timestep}, atoms={frame.n_atoms}, "
                  f"read in {t1 - t0:.3f}s")

        df = unpacker.get_frame_as_polars(0)
        print(f"  Polars DataFrame shape: {df.shape}")

    print("  PASSED\n")
    return test_file


def test_cell_binning():
    print("=" * 60)
    print("TEST 2: Cell Binning (Spatial Partitioning)")
    print("=" * 60)

    n_atoms = 10000
    box_size = 50.0
    cutoff = 3.5

    print(f"Generating {n_atoms} random atoms...")
    coords = np.random.uniform(0, box_size, size=(n_atoms, 3)).astype(np.float32)
    box = np.array([[0, box_size], [0, box_size], [0, box_size]], dtype=np.float64)

    t0 = time.time()
    cb = CellBinning(coords, box, cutoff)
    t1 = time.time()
    print(f"  Cell grid built in {t1 - t0:.3f}s")
    print(f"  Grid dimensions: {cb.n_cells}, cell size: {cb.cell_size}")

    t0 = time.time()
    neighbors = cb.get_neighbors(0)
    t1 = time.time()
    print(f"  Atom 0 neighbors: {len(neighbors)}, query time: {(t1 - t0) * 1000:.3f}ms")

    t0 = time.time()
    all_neighbors = cb.get_all_neighbors()
    t1 = time.time()
    avg_neighbors = sum(len(n) for n in all_neighbors) / n_atoms
    print(f"  All neighbors found in {t1 - t0:.3f}s, avg neighbors: {avg_neighbors:.2f}")

    print("  PASSED\n")


def test_rdf():
    print("=" * 60)
    print("TEST 3: Radial Distribution Function g(r)")
    print("=" * 60)

    n_atoms = 5000
    box_size = 30.0

    print(f"Generating {n_atoms} atoms...")
    coords = np.random.uniform(0, box_size, size=(n_atoms, 3)).astype(np.float32)
    types = np.random.choice([1, 2, 3, 4, 5], size=n_atoms)
    box = np.array([[0, box_size], [0, box_size], [0, box_size]], dtype=np.float64)

    calc = RDFCalculator(r_min=0.0, r_max=6.0, n_bins=100)

    t0 = time.time()
    result = calc.calculate(coords, box, types)
    t1 = time.time()
    print(f"  RDF calculated in {t1 - t0:.3f}s")
    print(f"  g(r) max: {max(result['gr']):.3f}")
    if 'partial_gr' in result:
        print(f"  Partial RDFs computed for {len(result['partial_gr'])} pairs")

    print("  PASSED\n")


def test_voronoi():
    print("=" * 60)
    print("TEST 4: Voronoi Polyhedron Analysis")
    print("=" * 60)

    n_atoms = 2000
    box_size = 20.0

    print(f"Generating {n_atoms} atoms...")
    coords = np.random.uniform(0, box_size, size=(n_atoms, 3)).astype(np.float32)
    types = np.random.choice([1, 2, 3, 4, 5], size=n_atoms)
    box = np.array([[0, box_size], [0, box_size], [0, box_size]], dtype=np.float64)

    analyzer = VoronoiAnalyzer(cutoff=5.0)

    t0 = time.time()
    result = analyzer.analyze_frame(coords, box, types)
    t1 = time.time()
    print(f"  Voronoi analysis completed in {t1 - t0:.3f}s")
    print(f"  Avg time per atom: {(t1 - t0) / n_atoms * 1000:.4f}ms")

    top_indices = sorted(result['counts'].items(), key=lambda x: -x[1])[:5]
    print(f"  Top 5 polyhedron types:")
    for idx, count in top_indices:
        print(f"    {idx}: {count} atoms ({count / n_atoms * 100:.1f}%)")

    test_idx = np.array(result['voronoi_indices'][0])
    print(f"  Atom 0 Voronoi index: {analyzer.get_polyhedron_type(test_idx)}")
    print(f"  Classification: {analyzer.classify_polyhedron(test_idx)}")

    print("  PASSED\n")
    return result


def test_csro():
    print("=" * 60)
    print("TEST 5: Chemical Short-Range Order (CSRO)")
    print("=" * 60)

    n_atoms = 5000
    box_size = 30.0

    print(f"Generating {n_atoms} atoms...")
    coords = np.random.uniform(0, box_size, size=(n_atoms, 3)).astype(np.float32)
    types = np.random.choice([1, 2, 3, 4, 5], size=n_atoms)
    box = np.array([[0, box_size], [0, box_size], [0, box_size]], dtype=np.float64)

    calc = CSROCalculator(cutoff=3.5)

    t0 = time.time()
    result = calc.calculate_csro(coords, box, types, type1=2, type2=3)
    t1 = time.time()
    print(f"  Co-Cr CSRO calculated in {t1 - t0:.3f}s")
    print(f"  Warren-Cowley S_ij: {result['S_ij']:.4f}")
    print(f"  Alpha parameter: {result['alpha']:.4f}")
    print(f"  Average coordination Z_avg: {result['Z_avg']:.2f}")
    print(f"  Concentrations: c_Co={result['c1']:.3f}, c_Cr={result['c2']:.3f}")
    print(f"  Bond counts: Co-Co={result['N_AA']:.0f}, Cr-Cr={result['N_BB']:.0f}, Co-Cr={result['N_AB']:.0f}")

    t0 = time.time()
    all_pairs = calc.analyze_all_pairs(coords, box, types)
    t1 = time.time()
    print(f"  All pair CSRO computed in {t1 - t0:.3f}s, {len(all_pairs)} pairs")

    print("  PASSED\n")


def test_analysis_service():
    print("=" * 60)
    print("TEST 6: Analysis Service (Integration)")
    print("=" * 60)

    service = AnalysisService()
    test_file = r"d:\SOLO-0621-1\9-hea-voronoi-analyzer\backend\data\test_heau.md"

    t0 = time.time()
    info = service.load_data(test_file)
    t1 = time.time()
    print(f"  Data loaded in {t1 - t0:.3f}s")
    print(f"  File: {info['file']}")
    print(f"  Frames: {info['n_frames']}, Atoms/frame: {info['n_atoms']}")

    t0 = time.time()
    frame_info = service.get_frame_info(0)
    t1 = time.time()
    print(f"  Frame info retrieved in {(t1 - t0) * 1000:.2f}ms")
    print(f"  Timestep: {frame_info['timestep']}")
    print(f"  Type counts: {frame_info['type_counts']}")

    t0 = time.time()
    csro_result = service.calculate_csro(0, 2, 3)
    t1 = time.time()
    print(f"  CSRO computed in {t1 - t0:.3f}s, S_ij={csro_result['S_ij']:.4f}")

    t0 = time.time()
    target_ts = 5000
    frame_idx = service.find_frame_by_timestep(target_ts)
    t1 = time.time()
    actual_ts = service.get_frame(frame_idx).timestep
    print(f"  Binary search for timestep {target_ts}: found frame {frame_idx} "
          f"(timestep={actual_ts}) in {(t1 - t0) * 1000:.3f}ms")

    print("\nStreaming evolution data...")
    frame_count = 0
    t0 = time.time()
    for data in service.get_evolution_stream(0, 3):
        frame_count += 1
        print(f"  Frame {data['frame_idx']}: timestep={data['timestep']}, "
              f"time={data['time_ps']:.2f}ps, types={len(data['counts'])}")
    t1 = time.time()
    print(f"  Streamed {frame_count} frames in {t1 - t0:.3f}s")

    service.close()
    print("  PASSED\n")


def main():
    print("\n" + "=" * 60)
    print("  AlCoCrFeNi HEA Microstructure Analysis Workstation")
    print("  Core Module Integration Test Suite")
    print("=" * 60 + "\n")

    total_t0 = time.time()

    try:
        test_file = test_md_unpacker()
        test_cell_binning()
        test_rdf()
        test_voronoi()
        test_csro()
        test_analysis_service()

        total_t1 = time.time()
        print("=" * 60)
        print(f"ALL TESTS PASSED! Total time: {total_t1 - total_t0:.3f}s")
        print("=" * 60)
        return 0
    except Exception as e:
        print(f"\nTEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
