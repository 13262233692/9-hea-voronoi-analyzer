<template>
  <Teleport to="body">
    <Transition name="modal-fade">
      <div v-if="visible" class="modal-backdrop" @mousedown.self="">
        <Transition name="modal-scale">
          <div v-if="visible" class="failure-modal" :class="severity">
            <div class="modal-header">
              <div class="header-icon">
                <span v-if="severity === 'critical'">💥</span>
                <span v-else>⚠️</span>
              </div>
              <div class="header-text">
                <h3>材料屈服失效警报</h3>
                <p>检测到局部塑性剪切带形成，晶格发生不可逆撕裂</p>
              </div>
              <button class="close-btn" @click="$emit('close')">×</button>
            </div>

            <div class="modal-body">
              <div class="stats-grid">
                <div class="stat-item">
                  <span class="stat-label">失效时间步</span>
                  <span class="stat-value">{{ formatNumber(timestep) }}</span>
                </div>
                <div class="stat-item">
                  <span class="stat-label">屈服原子数</span>
                  <span class="stat-value highlight">{{ formatNumber(nYielding) }}</span>
                </div>
                <div class="stat-item">
                  <span class="stat-label">屈服比例</span>
                  <span class="stat-value">{{ (yieldRatio * 100).toFixed(2) }}%</span>
                </div>
                <div class="stat-item">
                  <span class="stat-label">耗散能</span>
                  <span class="stat-value energy">{{ dissipatedEnergy.toFixed(3) }}</span>
                </div>
              </div>

              <div class="chart-section">
                <h4>瞬时耗散能积分曲线</h4>
                <div class="chart-container">
                  <canvas ref="chartCanvas" class="dissipation-chart"></canvas>
                </div>
                <div class="chart-legend">
                  <div class="legend-item">
                    <span class="legend-dot threshold"></span>
                    <span>屈服阈值 {{ threshold }} Å²</span>
                  </div>
                  <div class="legend-item">
                    <span class="legend-dot cumulative"></span>
                    <span>累积耗散能</span>
                  </div>
                </div>
              </div>

              <div class="warning-text">
                <p>⚡ <strong>物理临界值已突破</strong></p>
                <p class="sub">D²<sub>min</sub> = 0.08 Å² 表征晶格发生不可逆撕裂与位错交滑移</p>
              </div>
            </div>

            <div class="modal-footer">
              <button class="btn-secondary" @click="$emit('acknowledge')">
                确认并继续监测
              </button>
              <button class="btn-primary danger" @click="$emit('halt')">
                立即停机分析
              </button>
            </div>
          </div>
        </Transition>
      </div>
    </Transition>
  </Teleport>
</template>

<script setup>
import { ref, watch, onMounted, nextTick } from 'vue'

const props = defineProps({
  visible: Boolean,
  severity: {
    type: String,
    default: 'warning'
  },
  timestep: {
    type: Number,
    default: 0
  },
  nYielding: {
    type: Number,
    default: 0
  },
  yieldRatio: {
    type: Number,
    default: 0
  },
  dissipatedEnergy: {
    type: Number,
    default: 0
  },
  threshold: {
    type: Number,
    default: 0.08
  },
  dissipationCurve: {
    type: Object,
    default: () => null
  }
})

defineEmits(['close', 'acknowledge', 'halt'])

const chartCanvas = ref(null)

const formatNumber = (n) => {
  return n?.toLocaleString() || '0'
}

const drawChart = () => {
  const canvas = chartCanvas.value
  if (!canvas) return

  const ctx = canvas.getContext('2d')
  const rect = canvas.getBoundingClientRect()
  const dpr = window.devicePixelRatio || 1
  
  canvas.width = rect.width * dpr
  canvas.height = rect.height * dpr
  ctx.scale(dpr, dpr)

  const w = rect.width
  const h = rect.height
  const padding = { top: 20, right: 20, bottom: 40, left: 60 }
  const chartW = w - padding.left - padding.right
  const chartH = h - padding.top - padding.bottom

  ctx.clearRect(0, 0, w, h)

  const curve = props.dissipationCurve
  if (!curve || !curve.bins || curve.bins.length < 2) {
    ctx.fillStyle = '#606080'
    ctx.font = '12px sans-serif'
    ctx.textAlign = 'center'
    ctx.fillText('暂无耗散能数据', w / 2, h / 2)
    return
  }

  const bins = curve.bins
  const counts = curve.counts
  const cumulative = curve.cumulative_dissipation

  const xMin = bins[0]
  const xMax = bins[bins.length - 1]
  const yMax = Math.max(...cumulative) * 1.1 || 1

  const logX = (x) => {
    return padding.left + chartW * (Math.log10(x) - Math.log10(xMin)) / (Math.log10(xMax) - Math.log10(xMin))
  }

  const yScale = (y) => {
    return padding.top + chartH - (y / yMax) * chartH
  }

  ctx.strokeStyle = 'rgba(100, 100, 150, 0.2)'
  ctx.lineWidth = 1
  for (let i = 0; i <= 5; i++) {
    const y = padding.top + (chartH / 5) * i
    ctx.beginPath()
    ctx.moveTo(padding.left, y)
    ctx.lineTo(w - padding.right, y)
    ctx.stroke()
  }

  ctx.strokeStyle = 'rgba(255, 0, 255, 0.6)'
  ctx.setLineDash([6, 4])
  ctx.lineWidth = 2
  const thresholdX = logX(props.threshold)
  ctx.beginPath()
  ctx.moveTo(thresholdX, padding.top)
  ctx.lineTo(thresholdX, padding.top + chartH)
  ctx.stroke()
  ctx.setLineDash([])

  ctx.fillStyle = '#ff80ff'
  ctx.font = '11px sans-serif'
  ctx.textAlign = 'left'
  ctx.fillText(`阈值 ${props.threshold}`, thresholdX + 5, padding.top + 15)

  const grad = ctx.createLinearGradient(0, padding.top, 0, padding.top + chartH)
  grad.addColorStop(0, 'rgba(255, 0, 255, 0.3)')
  grad.addColorStop(1, 'rgba(255, 0, 255, 0.0)')

  ctx.beginPath()
  ctx.moveTo(logX(bins[0]), yScale(0))
  for (let i = 0; i < cumulative.length; i++) {
    const x = logX(bins[i + 1] || bins[i])
    const y = yScale(cumulative[i])
    ctx.lineTo(x, y)
  }
  ctx.lineTo(logX(bins[bins.length - 1]), yScale(0))
  ctx.closePath()
  ctx.fillStyle = grad
  ctx.fill()

  ctx.beginPath()
  ctx.moveTo(logX(bins[0]), yScale(0))
  for (let i = 0; i < cumulative.length; i++) {
    const x = logX(bins[i + 1] || bins[i])
    const y = yScale(cumulative[i])
    ctx.lineTo(x, y)
  }
  ctx.strokeStyle = '#ff00ff'
  ctx.lineWidth = 2
  ctx.stroke()

  ctx.fillStyle = '#8080a0'
  ctx.font = '11px sans-serif'
  ctx.textAlign = 'center'
  ctx.fillText('D²_min (Å²) - 对数刻度', w / 2, h - 10)

  ctx.save()
  ctx.translate(15, padding.top + chartH / 2)
  ctx.rotate(-Math.PI / 2)
  ctx.textAlign = 'center'
  ctx.fillStyle = '#8080a0'
  ctx.fillText('累积耗散能 (eV?)', 0, 0)
  ctx.restore()
}

watch(() => props.visible, (val) => {
  if (val) {
    nextTick(() => {
      setTimeout(drawChart, 50)
    })
  }
})

watch(() => props.dissipationCurve, () => {
  if (props.visible) {
    nextTick(drawChart)
  }
}, { deep: true })

onMounted(() => {
  if (props.visible) {
    nextTick(drawChart)
  }
})
</script>

<style scoped>
.modal-backdrop {
  position: fixed;
  top: 0;
  left: 0;
  right: 0;
  bottom: 0;
  background: rgba(0, 0, 0, 0.7);
  display: flex;
  align-items: center;
  justify-content: center;
  z-index: 9999;
  backdrop-filter: blur(4px);
}

.failure-modal {
  background: linear-gradient(135deg, rgba(30, 20, 50, 0.98), rgba(50, 20, 60, 0.98));
  border: 2px solid #ff00ff;
  border-radius: 16px;
  width: 560px;
  max-width: 90vw;
  max-height: 90vh;
  overflow: hidden;
  box-shadow: 
    0 0 60px rgba(255, 0, 255, 0.4),
    0 20px 60px rgba(0, 0, 0, 0.5);
  animation: glow 2s ease-in-out infinite alternate;
}

.failure-modal.critical {
  border-color: #ff3030;
  box-shadow: 
    0 0 60px rgba(255, 50, 50, 0.4),
    0 20px 60px rgba(0, 0, 0, 0.5);
}

@keyframes glow {
  from {
    box-shadow: 0 0 40px rgba(255, 0, 255, 0.3), 0 20px 60px rgba(0, 0, 0, 0.5);
  }
  to {
    box-shadow: 0 0 80px rgba(255, 0, 255, 0.6), 0 20px 60px rgba(0, 0, 0, 0.5);
  }
}

.failure-modal.critical {
  animation: glow-critical 2s ease-in-out infinite alternate;
}

@keyframes glow-critical {
  from {
    box-shadow: 0 0 40px rgba(255, 50, 50, 0.3), 0 20px 60px rgba(0, 0, 0, 0.5);
  }
  to {
    box-shadow: 0 0 80px rgba(255, 50, 50, 0.6), 0 20px 60px rgba(0, 0, 0, 0.5);
  }
}

.modal-header {
  display: flex;
  align-items: center;
  gap: 15px;
  padding: 20px 25px;
  border-bottom: 1px solid rgba(255, 0, 255, 0.2);
  background: linear-gradient(90deg, rgba(255, 0, 255, 0.1), transparent);
}

.failure-modal.critical .modal-header {
  background: linear-gradient(90deg, rgba(255, 50, 50, 0.15), transparent);
  border-bottom-color: rgba(255, 50, 50, 0.2);
}

.header-icon {
  font-size: 36px;
  animation: pulse 1s ease-in-out infinite;
}

@keyframes pulse {
  0%, 100% { transform: scale(1); }
  50% { transform: scale(1.1); }
}

.header-text {
  flex: 1;
}

.header-text h3 {
  color: #ff80ff;
  font-size: 20px;
  margin: 0 0 4px 0;
  font-weight: 700;
}

.failure-modal.critical .header-text h3 {
  color: #ff6060;
}

.header-text p {
  color: #a080c0;
  font-size: 13px;
  margin: 0;
}

.close-btn {
  background: none;
  border: none;
  color: #8080a0;
  font-size: 28px;
  cursor: pointer;
  padding: 0;
  width: 32px;
  height: 32px;
  line-height: 1;
  border-radius: 6px;
  transition: all 0.2s;
}

.close-btn:hover {
  background: rgba(255, 0, 255, 0.2);
  color: #ff80ff;
}

.modal-body {
  padding: 20px 25px;
  max-height: 60vh;
  overflow-y: auto;
}

.stats-grid {
  display: grid;
  grid-template-columns: repeat(4, 1fr);
  gap: 12px;
  margin-bottom: 20px;
}

.stat-item {
  background: rgba(40, 20, 60, 0.6);
  border: 1px solid rgba(255, 0, 255, 0.2);
  border-radius: 8px;
  padding: 12px;
  text-align: center;
}

.failure-modal.critical .stat-item {
  border-color: rgba(255, 50, 50, 0.2);
  background: rgba(60, 20, 30, 0.6);
}

.stat-label {
  display: block;
  font-size: 11px;
  color: #8080a0;
  margin-bottom: 6px;
}

.stat-value {
  display: block;
  font-size: 18px;
  font-weight: 700;
  color: #e0e0f0;
}

.stat-value.highlight {
  color: #ff80ff;
}

.failure-modal.critical .stat-value.highlight {
  color: #ff6060;
}

.stat-value.energy {
  color: #ffcc00;
}

.chart-section {
  background: rgba(20, 10, 40, 0.6);
  border: 1px solid rgba(255, 0, 255, 0.15);
  border-radius: 10px;
  padding: 15px;
  margin-bottom: 16px;
}

.chart-section h4 {
  color: #c0a0e0;
  font-size: 14px;
  margin: 0 0 10px 0;
}

.chart-container {
  width: 100%;
  height: 180px;
}

.dissipation-chart {
  width: 100%;
  height: 100%;
}

.chart-legend {
  display: flex;
  justify-content: center;
  gap: 20px;
  margin-top: 8px;
}

.legend-item {
  display: flex;
  align-items: center;
  gap: 6px;
  font-size: 11px;
  color: #8080a0;
}

.legend-dot {
  width: 10px;
  height: 10px;
  border-radius: 2px;
}

.legend-dot.threshold {
  background: #ff00ff;
  height: 10px;
}

.legend-dot.cumulative {
  background: linear-gradient(180deg, rgba(255,0,255,0.5), rgba(255,0,255,0.1));
}

.warning-text {
  background: rgba(255, 200, 0, 0.1);
  border-left: 3px solid #ffcc00;
  padding: 12px 15px;
  border-radius: 6px;
}

.warning-text p {
  margin: 0;
  color: #ffdd66;
  font-size: 13px;
}

.warning-text .sub {
  color: #b0a060;
  font-size: 11px;
  margin-top: 4px;
}

.modal-footer {
  display: flex;
  gap: 12px;
  padding: 18px 25px;
  border-top: 1px solid rgba(255, 0, 255, 0.2);
  background: rgba(20, 10, 40, 0.5);
}

.failure-modal.critical .modal-footer {
  border-top-color: rgba(255, 50, 50, 0.2);
}

.btn-secondary,
.btn-primary {
  flex: 1;
  padding: 12px 20px;
  border-radius: 8px;
  font-size: 14px;
  font-weight: 600;
  cursor: pointer;
  transition: all 0.2s;
  border: 1px solid;
}

.btn-secondary {
  background: rgba(60, 60, 90, 0.5);
  border-color: rgba(100, 100, 150, 0.4);
  color: #b0b0d0;
}

.btn-secondary:hover {
  background: rgba(80, 80, 120, 0.6);
  color: #e0e0f0;
}

.btn-primary {
  background: linear-gradient(135deg, #667eea, #764ba2);
  border-color: #667eea;
  color: white;
}

.btn-primary:hover {
  transform: translateY(-1px);
  box-shadow: 0 4px 15px rgba(102, 126, 234, 0.4);
}

.btn-primary.danger {
  background: linear-gradient(135deg, #ff3030, #ff0066);
  border-color: #ff3030;
}

.btn-primary.danger:hover {
  box-shadow: 0 4px 15px rgba(255, 50, 50, 0.4);
}

.modal-fade-enter-active,
.modal-fade-leave-active {
  transition: opacity 0.3s ease;
}

.modal-fade-enter-from,
.modal-fade-leave-to {
  opacity: 0;
}

.modal-scale-enter-active,
.modal-scale-leave-active {
  transition: transform 0.3s cubic-bezier(0.34, 1.56, 0.64, 1);
}

.modal-scale-enter-from,
.modal-scale-leave-to {
  transform: scale(0.8);
}
</style>
