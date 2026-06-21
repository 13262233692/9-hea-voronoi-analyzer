import numpy as np
import sys
sys.path.insert(0, '.')
from app.core.md_unpacker import TriclinicBox
from app.core.nonaffine_analyzer import NonAffineAnalyzer

box = TriclinicBox(
    np.array([30.0, 0.0, 0.0]),
    np.array([0.0, 30.0, 0.0]),
    np.array([0.0, 0.0, 30.0])
)

n_atoms = 200
coords_ref = np.random.random((n_atoms, 3)) * 30.0

coords_def = coords_ref.copy()
band_mask = (coords_def[:, 0] > 12) & (coords_def[:, 0] < 18)
coords_def[band_mask, 1] += 2.0 + np.random.randn(band_mask.sum()) * 0.3
coords_def += np.random.randn(*coords_def.shape) * 0.05

analyzer = NonAffineAnalyzer(cutoff=4.5, yield_threshold=0.08)
result = analyzer.analyze_deformation(coords_ref, coords_def, box)

print('Total atoms:', result.total_atoms)
print('Yielding atoms:', result.n_yielding)
print('Dissipated energy: %.4f' % result.dissipated_energy)
print('D2_min range: %.6f - %.6f' % (result.d2_min.min(), result.d2_min.max()))
print('Yield threshold:', result.threshold)

clusters = analyzer.detect_shear_band_clusters(coords_def, result.shear_band_atoms, box)
print('Clusters found:', len(clusters))
for i, c in enumerate(clusters[:3]):
    print('  Cluster %d: %d atoms, center=%s' % (i, c['n_atoms'], c['center']))

diss_curve = analyzer.compute_dissipation_curve(coords_ref, coords_def, box, n_bins=20)
print('Dissipation curve bins:', len(diss_curve['bins']) - 1)
print('Total dissipated: %.4f' % diss_curve['total_dissipated'])
print('OK')
