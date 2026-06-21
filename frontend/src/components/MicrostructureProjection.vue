<template>
  <div class="projection-container">
    <div class="projection-header">
      <h4>微观结构投影图</h4>
      <div class="projection-controls">
        <label>
          投影面
          <select v-model="projectionPlane" @change="redraw">
            <option value="xy">XY 面</option>
            <option value="xz">XZ 面</option>
            <option value="yz">YZ 面</option>
          </select>
        </label>
        <label>
          原子大小
          <input type="range" min="1" max="8" v-model.number="atomRadius" @input="redraw" />
        </label>
      </div>
    </div>
    <div class="canvas-wrapper" ref="canvasWrapper">
      <canvas ref="bgCanvas" class="projection-canvas"></canvas>
      <canvas ref="atomCanvas" class="projection-canvas overlay"></canvas>
      <canvas ref="shearCanvas" class="projection-canvas overlay shear-layer"></canvas>
    </div>
    <div v-if="yieldEvent" class="yield-badge" :class="yieldEvent.severity">
      <span class="badge-icon">⚠</span>
      <span>检测到塑性屈服</span>
      <span class="yield-count">{{ yieldEvent.n_yielding }} 原子</span>
    </div>
  </div>
</template>

<script setup>import { ref, onMounted, onUnmounted, watch } from 'vue';
const props = defineProps({
 atoms: {
 type: Array,
 default: () => []
 },
 boxVectors: {
 type: Object,
 default: null
 },
 shearBandAtoms: {
 type: Array,
 default: () => []
 },
 yieldEvent: {
 type: Object,
 default: null
 },
 atomTypes: {
 type: Array,
 default: () => []
 }
});
const canvasWrapper = ref(null);
const bgCanvas = ref(null);
const atomCanvas = ref(null);
const shearCanvas = ref(null);
const projectionPlane = ref('xy');
const atomRadius = ref(3);
const typeColors = {
 1: '#c0c0c0',
 2: '#4a90d9',
 3: '#cd853f',
 4: '#daa520',
 5: '#a9a9a9'
};
let resizeObserver = null;
const getCanvasSize = () => {
 if (!canvasWrapper.value)
 return { w: 800, h: 400 };
 const rect = canvasWrapper.value.getBoundingClientRect();
 return { w: rect.width, h: rect.height };
};
const setupCanvases = () => {
 const { w, h } = getCanvasSize();
 [bgCanvas.value, atomCanvas.value, shearCanvas.value].forEach(canvas => {
 if (!canvas)
 return;
 canvas.width = w * window.devicePixelRatio;
 canvas.height = h * window.devicePixelRatio;
 canvas.style.width = w + 'px';
 canvas.style.height = h + 'px';
 const ctx = canvas.getContext('2d');
 ctx.scale(window.devicePixelRatio, window.devicePixelRatio);
 });
 drawBackground();
};
const drawBackground = () => {
 const canvas = bgCanvas.value;
 if (!canvas)
 return;
 const ctx = canvas.getContext('2d');
 const { w, h } = getCanvasSize();
 ctx.clearRect(0, 0, w, h);
 const gradient = ctx.createLinearGradient(0, 0, 0, h);
 gradient.addColorStop(0, 'rgba(20, 20, 40, 0.9)');
 gradient.addColorStop(1, 'rgba(30, 30, 50, 0.9)');
 ctx.fillStyle = gradient;
 ctx.fillRect(0, 0, w, h);
 ctx.strokeStyle = 'rgba(100, 100, 150, 0.15)';
 ctx.lineWidth = 1;
 const gridSize = 40;
 for (let x = 0; x < w; x += gridSize) {
 ctx.beginPath();
 ctx.moveTo(x, 0);
 ctx.lineTo(x, h);
 ctx.stroke();
 }
 for (let y = 0; y < h; y += gridSize) {
 ctx.beginPath();
 ctx.moveTo(0, y);
 ctx.lineTo(w, y);
 ctx.stroke();
 }
};
const projectCoords = (coords) => {
 const plane = projectionPlane.value;
 if (plane === 'xy')
 return [coords[0], coords[1]];
 if (plane === 'xz')
 return [coords[0], coords[2]];
 return [coords[1], coords[2]];
};
const getBounds = () => {
 if (!props.boxVectors) {
 return { min: [0, 0, 0], max: [100, 100, 100] };
 }
 const a = props.boxVectors.a || [100, 0, 0];
 const b = props.boxVectors.b || [0, 100, 0];
 const c = props.boxVectors.c || [0, 0, 100];
 const allX = [0, a[0], b[0], c[0], a[0] + b[0], a[0] + c[0], b[0] + c[0], a[0] + b[0] + c[0]];
 const allY = [0, a[1], b[1], c[1], a[1] + b[1], a[1] + c[1], b[1] + c[1], a[1] + b[1] + c[1]];
 const allZ = [0, a[2], b[2], c[2], a[2] + b[2], a[2] + c[2], b[2] + c[2], a[2] + b[2] + c[2]];
 return {
 min: [Math.min(...allX), Math.min(...allY), Math.min(...allZ)],
 max: [Math.max(...allX), Math.max(...allY), Math.max(...allZ)]
 };
};
const worldToScreen = (x, y) => {
 const bounds = getBounds();
 const { w, h } = getCanvasSize();
 const boundsMin = projectCoords(bounds.min);
 const boundsMax = projectCoords(bounds.max);
 const rangeX = boundsMax[0] - boundsMin[0] || 1;
 const rangeY = boundsMax[1] - boundsMin[1] || 1;
 const padding = 20;
 const scaleX = (w - padding * 2) / rangeX;
 const scaleY = (h - padding * 2) / rangeY;
 const scale = Math.min(scaleX, scaleY);
 const offsetX = (w - rangeX * scale) / 2;
 const offsetY = (h - rangeY * scale) / 2;
 return {
 x: offsetX + (x - boundsMin[0]) * scale,
 y: offsetY + (y - boundsMin[1]) * scale,
 scale
 };
};
const drawAtoms = () => {
 const canvas = atomCanvas.value;
 if (!canvas || !props.atoms || props.atoms.length === 0)
 return;
 const ctx = canvas.getContext('2d');
 const { w, h } = getCanvasSize();
 ctx.clearRect(0, 0, w, h);
 const atoms = props.atoms;
 const types = props.atomTypes;
 for (let i = 0; i < atoms.length; i++) {
 const [x, y] = projectCoords(atoms[i]);
 const { x: sx, y: sy, scale } = worldToScreen(x, y);
 const type = types[i] || 1;
 const color = typeColors[type] || '#888';
 const r = atomRadius.value;
 ctx.beginPath();
 ctx.arc(sx, sy, r, 0, Math.PI * 2);
 ctx.fillStyle = color;
 ctx.fill();
 ctx.strokeStyle = 'rgba(255,255,255,0.2)';
 ctx.lineWidth = 0.5;
 ctx.stroke();
 }
};
const drawShearBand = () => {
 const canvas = shearCanvas.value;
 if (!canvas)
 return;
 const ctx = canvas.getContext('2d');
 const { w, h } = getCanvasSize();
 ctx.clearRect(0, 0, w, h);
 const shearAtoms = props.shearBandAtoms;
 if (!shearAtoms || shearAtoms.length === 0)
 return;
 ctx.save();
 ctx.shadowColor = '#ff00ff';
 ctx.shadowBlur = 15;
 for (let i = 0; i < shearAtoms.length; i++) {
 const [x, y] = projectCoords(shearAtoms[i]);
 const { x: sx, y: sy } = worldToScreen(x, y);
 const r = atomRadius.value * 1.5;
 const gradient = ctx.createRadialGradient(sx, sy, 0, sx, sy, r);
 gradient.addColorStop(0, 'rgba(255, 0, 255, 0.9)');
 gradient.addColorStop(0.5, 'rgba(255, 0, 255, 0.6)');
 gradient.addColorStop(1, 'rgba(255, 0, 255, 0)');
 ctx.beginPath();
 ctx.arc(sx, sy, r, 0, Math.PI * 2);
 ctx.fillStyle = gradient;
 ctx.fill();
 }
 ctx.restore();
 if (shearAtoms.length > 10) {
 drawShearLine(ctx, shearAtoms);
 }
};
const drawShearLine = (ctx, atoms) => {
 const coords2d = atoms.map(a => projectCoords(a));
 let minX = Infinity, maxX = -Infinity;
 let minY = Infinity, maxY = -Infinity;
 for (const [x, y] of coords2d) {
 if (x < minX)
 minX = x;
 if (x > maxX)
 maxX = x;
 if (y < minY)
 minY = y;
 if (y > maxY)
 maxY = y;
 }
 const p1 = worldToScreen(minX, (minY + maxY) / 2);
 const p2 = worldToScreen(maxX, (minY + maxY) / 2);
 ctx.save();
 ctx.shadowColor = '#ff00ff';
 ctx.shadowBlur = 20;
 ctx.strokeStyle = 'rgba(255, 0, 255, 0.8)';
 ctx.lineWidth = 4;
 ctx.lineCap = 'round';
 ctx.beginPath();
 ctx.moveTo(p1.x, p1.y);
 ctx.lineTo(p2.x, p2.y);
 ctx.stroke();
 ctx.strokeStyle = 'rgba(255, 200, 255, 0.9)';
 ctx.lineWidth = 2;
 ctx.beginPath();
 ctx.moveTo(p1.x, p1.y);
 ctx.lineTo(p2.x, p2.y);
 ctx.stroke();
 ctx.restore();
};
const redraw = () => {
 drawAtoms();
 drawShearBand();
};
watch(() => props.atoms, redraw, { deep: true });
watch(() => props.shearBandAtoms, redraw, { deep: true });
watch(() => props.yieldEvent, redraw);
onMounted(() => {
 setupCanvases();
 redraw();
 resizeObserver = new ResizeObserver(() => {
 setupCanvases();
 redraw();
 });
 if (canvasWrapper.value) {
 resizeObserver.observe(canvasWrapper.value);
 }
});
onUnmounted(() => {
 if (resizeObserver) {
 resizeObserver.disconnect();
 }
});
defineExpose({ redraw });
</script>

<style scoped>
.projection-container {
  background: rgba(30, 30, 50, 0.8);
  border-radius: 12px;
  border: 1px solid rgba(100, 100, 150, 0.3);
  overflow: hidden;
  position: relative;
}

.projection-header {
  display: flex;
  justify-content: space-between;
  align-items: center;
  padding: 15px 20px;
  border-bottom: 1px solid rgba(100, 100, 150, 0.3);
}

.projection-header h4 {
  color: #8ab4f8;
  font-size: 16px;
  margin: 0;
}

.projection-controls {
  display: flex;
  gap: 20px;
}

.projection-controls label {
  display: flex;
  align-items: center;
  gap: 8px;
  font-size: 12px;
  color: #a0a0c0;
}

.projection-controls select,
.projection-controls input[type="range"] {
  background: rgba(20, 20, 40, 0.8);
  border: 1px solid rgba(100, 100, 150, 0.4);
  color: #e0e0e0;
  border-radius: 4px;
  padding: 4px 8px;
  font-size: 12px;
}

.canvas-wrapper {
  position: relative;
  width: 100%;
  height: 400px;
}

.projection-canvas {
  position: absolute;
  top: 0;
  left: 0;
  width: 100%;
  height: 100%;
}

.projection-canvas.overlay {
  pointer-events: none;
}

.shear-layer {
  mix-blend-mode: screen;
}

.yield-badge {
  position: absolute;
  top: 20px;
  right: 20px;
  background: linear-gradient(135deg, rgba(255, 0, 128, 0.3), rgba(255, 0, 255, 0.3));
  border: 1px solid #ff00ff;
  border-radius: 8px;
  padding: 8px 16px;
  color: #ff80ff;
  font-size: 13px;
  font-weight: 600;
  display: flex;
  align-items: center;
  gap: 8px;
  box-shadow: 0 0 20px rgba(255, 0, 255, 0.4);
  animation: pulse 1.5s ease-in-out infinite;
}

.yield-badge.critical {
  background: linear-gradient(135deg, rgba(255, 50, 50, 0.4), rgba(255, 0, 100, 0.4));
  border-color: #ff5050;
  color: #ff8080;
  box-shadow: 0 0 20px rgba(255, 50, 50, 0.4);
}

.yield-badge .badge-icon {
  font-size: 16px;
}

.yield-count {
  font-size: 11px;
  opacity: 0.8;
}

@keyframes pulse {
  0%, 100% {
    opacity: 1;
  }
  50% {
    opacity: 0.7;
  }
}
</style>
