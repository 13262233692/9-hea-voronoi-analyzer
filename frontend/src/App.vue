<template>
  <div class="app-container">
    <header class="app-header">
      <div class="header-content">
        <h1>
          <span class="title-icon">⚛</span>
          AlCoCrFeNi 高熵合金微观结构分析工作站
        </h1>
        <p class="subtitle">Voronoi 多面体剖分 · RDF · CSRO · 分子动力学可视化</p>
      </div>
      <div class="header-status">
        <div class="status-badge" :class="{ online: backendOnline }">
          {{ backendOnline ? '后端在线' : '后端离线' }}
        </div>
      </div>
    </header>

    <main class="main-content">
      <aside class="sidebar">
        <ControlPanel
          :is-connected="isConnected"
          :is-streaming="isStreaming"
          :current-frame="currentFrame"
          :total-frames="totalFrames"
          :data-info="dataInfo"
          :frame-info="frameInfo"
          :polyhedron-types="evolutionData.polyhedronTypes"
          :get-color="getColor"
          :error="error"
          @load-data="handleLoadData"
          @start-stream="handleStartStream"
          @stop-stream="handleStopStream"
          @reset="handleReset"
          @generate-test-data="handleGenerateTestData"
        />

        <div v-if="frameInfo" class="info-card">
          <h4>当前帧信息</h4>
          <div class="info-grid">
            <div class="info-item">
              <span class="info-label">时间步</span>
              <span class="info-value">{{ frameInfo.timestep?.toLocaleString() }}</span>
            </div>
            <div class="info-item">
              <span class="info-label">原子数</span>
              <span class="info-value">{{ frameInfo.n_atoms?.toLocaleString() }}</span>
            </div>
          </div>
          <div class="type-counts">
            <h5>元素分布</h5>
            <div class="type-bars">
              <div v-for="(count, element) in frameInfo.type_counts" :key="element" class="type-bar-item">
                <span class="element-name">{{ element }}</span>
                <div class="type-bar">
                  <div class="type-bar-fill" :style="{ width: getPercent(count) + '%' }" :class="element"></div>
                </div>
                <span class="count-value">{{ count }}</span>
              </div>
            </div>
          </div>
        </div>
      </aside>

      <section class="content-area">
        <StackedAreaChart
          v-if="evolutionData.timePs.length > 0"
          :evolution-data="evolutionData"
          :polyhedron-types="evolutionData.polyhedronTypes"
          :get-color="getColor"
          :binary-search="binarySearchTimestep"
          :query-c-s-r-o="queryCSROByTimestep"
          @hover="handleChartHover"
        />

        <div v-else class="empty-state">
          <div class="empty-icon">📊</div>
          <h3>等待数据</h3>
          <p>加载MD数据文件并开始播放以查看多面体类型演化</p>
          <div class="quick-actions">
            <button @click="handleGenerateTestData" class="btn-primary">
              生成测试数据
            </button>
          </div>
        </div>

        <div v-if="hoveredData" class="hover-info-panel">
          <h4>悬停时间点详情</h4>
          <div class="hover-grid">
            <div class="hover-item">
              <span>模拟时间</span>
              <strong>{{ hoveredData.timePs.toFixed(3) }} ps</strong>
            </div>
            <div class="hover-item">
              <span>时间步</span>
              <strong>{{ hoveredData.timestep?.toLocaleString() }}</strong>
            </div>
            <div class="hover-item">
              <span>帧索引</span>
              <strong>#{{ hoveredData.frameIndex }}</strong>
            </div>
          </div>
        </div>

        <div class="analysis-cards">
          <div class="analysis-card rdf-card">
            <h4>径向分布函数 g(r)</h4>
            <button @click="analyzeRDF" :disabled="!dataInfo || isAnalyzing">
              {{ isAnalyzing ? '分析中...' : '计算RDF' }}
            </button>
            <div v-if="rdfResult" class="rdf-result">
              <span>第一峰位置: {{ findFirstPeak(rdfResult.r, rdfResult.gr).toFixed(2) }} Å</span>
            </div>
          </div>

          <div class="analysis-card voronoi-card">
            <h4>Voronoi 分析</h4>
            <button @click="analyzeVoronoi" :disabled="!dataInfo || isAnalyzing">
              {{ isAnalyzing ? '分析中...' : 'Voronoi剖分' }}
            </button>
            <div v-if="voronoiResult" class="voronoi-result">
              <span>主要构型: {{ getTopPolyhedron(voronoiResult.counts) }}</span>
            </div>
          </div>

          <div class="analysis-card csro-card">
            <h4>化学短程有序度</h4>
            <button @click="analyzeCSRO" :disabled="!dataInfo || isAnalyzing">
              {{ isAnalyzing ? '分析中...' : '计算CSRO' }}
            </button>
            <div v-if="csroResult" class="csro-result">
              <span>Co-Cr S_ij: {{ csroResult.S_ij?.toFixed(4) }}</span>
            </div>
          </div>
        </div>
      </section>
    </main>
  </div>
</template>

<script setup>
import { ref, computed, onMounted } from 'vue'
import ControlPanel from './components/ControlPanel.vue'
import StackedAreaChart from './components/StackedAreaChart.vue'
import { useEvolutionStream } from './composables/useEvolutionStream'
import { useApi } from './composables/useApi'

const {
  isConnected,
  isStreaming,
  currentFrame,
  totalFrames,
  evolutionData,
  error,
  connect,
  disconnect,
  reset,
  binarySearchTimestep,
  getColor
} = useEvolutionStream()

const {
  loadData,
  getDataInfo,
  getFrameInfo,
  getRDF,
  getVoronoi,
  getCSRO,
  queryCSROByTimestep,
  healthCheck
} = useApi()

const backendOnline = ref(false)
const dataInfo = ref(null)
const frameInfo = ref(null)
const hoveredData = ref(null)
const isAnalyzing = ref(false)
const rdfResult = ref(null)
const voronoiResult = ref(null)
const csroResult = ref(null)

const totalAtoms = computed(() => {
  if (!frameInfo.value?.type_counts) return 0
  return Object.values(frameInfo.value.type_counts).reduce((a, b) => a + b, 0)
})

const getPercent = (count) => {
  const total = totalAtoms.value || 1
  return (count / total) * 100
}

const checkBackend = async () => {
  try {
    await healthCheck()
    backendOnline.value = true
  } catch (e) {
    backendOnline.value = false
  }
}

const handleLoadData = async (filepath) => {
  try {
    const result = await loadData(filepath)
    dataInfo.value = result.data
    frameInfo.value = await getFrameInfo(0)
  } catch (e) {
    console.error('Failed to load data:', e)
  }
}

const handleStartStream = (options) => {
  connect(options)
}

const handleStopStream = () => {
  disconnect()
}

const handleReset = () => {
  reset()
  rdfResult.value = null
  voronoiResult.value = null
  csroResult.value = null
  hoveredData.value = null
}

const handleGenerateTestData = async () => {
  try {
    const response = await fetch('http://localhost:8000/api/generate-test-data', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        filepath: 'd:/SOLO-0621-1/9-hea-voronoi-analyzer/backend/data/test_heau.md',
        n_frames: 20,
        n_atoms: 50000
      })
    })
    if (response.ok) {
      await handleLoadData('d:/SOLO-0621-1/9-hea-voronoi-analyzer/backend/data/test_heau.md')
    }
  } catch (e) {
    console.error('Failed to generate test data:', e)
  }
}

const handleChartHover = (data) => {
  hoveredData.value = data
}

const analyzeRDF = async () => {
  if (!dataInfo.value) return
  isAnalyzing.value = true
  try {
    rdfResult.value = await getRDF(0)
  } finally {
    isAnalyzing.value = false
  }
}

const analyzeVoronoi = async () => {
  if (!dataInfo.value) return
  isAnalyzing.value = true
  try {
    voronoiResult.value = await getVoronoi(0)
  } finally {
    isAnalyzing.value = false
  }
}

const analyzeCSRO = async () => {
  if (!dataInfo.value) return
  isAnalyzing.value = true
  try {
    csroResult.value = await getCSRO(0, 2, 3)
  } finally {
    isAnalyzing.value = false
  }
}

const findFirstPeak = (r, gr) => {
  let maxVal = 0
  let maxIdx = 0
  for (let i = 1; i < gr.length - 1; i++) {
    if (gr[i] > gr[i - 1] && gr[i] > gr[i + 1] && gr[i] > maxVal) {
      maxVal = gr[i]
      maxIdx = i
    }
  }
  return r[maxIdx] || 0
}

const getTopPolyhedron = (counts) => {
  if (!counts) return '-'
  const entries = Object.entries(counts)
  if (!entries.length) return '-'
  entries.sort((a, b) => b[1] - a[1])
  return `${entries[0][0]} (${entries[0][1]})`
}

onMounted(() => {
  checkBackend()
  setInterval(checkBackend, 5000)
})
</script>

<style scoped>
.app-container {
  min-height: 100vh;
  display: flex;
  flex-direction: column;
}

.app-header {
  background: linear-gradient(135deg, rgba(102, 126, 234, 0.2) 0%, rgba(118, 75, 162, 0.2) 100%);
  border-bottom: 1px solid rgba(100, 100, 150, 0.3);
  padding: 20px 30px;
  display: flex;
  justify-content: space-between;
  align-items: center;
}

.header-content h1 {
  font-size: 24px;
  color: #ffffff;
  font-weight: 600;
  display: flex;
  align-items: center;
  gap: 12px;
}

.title-icon {
  font-size: 32px;
}

.subtitle {
  color: #a0a0c0;
  font-size: 13px;
  margin-top: 6px;
}

.status-badge {
  padding: 8px 16px;
  border-radius: 20px;
  font-size: 12px;
  font-weight: 500;
  background: rgba(255, 107, 107, 0.2);
  color: #ff6b6b;
  border: 1px solid rgba(255, 107, 107, 0.3);
}

.status-badge.online {
  background: rgba(78, 205, 196, 0.2);
  color: #4ecdc4;
  border-color: rgba(78, 205, 196, 0.3);
}

.main-content {
  flex: 1;
  display: flex;
  gap: 20px;
  padding: 20px 30px;
}

.sidebar {
  width: 380px;
  flex-shrink: 0;
  display: flex;
  flex-direction: column;
  gap: 20px;
}

.content-area {
  flex: 1;
  display: flex;
  flex-direction: column;
  gap: 20px;
}

.info-card {
  background: rgba(30, 30, 50, 0.8);
  border-radius: 12px;
  padding: 20px;
  border: 1px solid rgba(100, 100, 150, 0.3);
}

.info-card h4 {
  color: #8ab4f8;
  font-size: 16px;
  margin-bottom: 15px;
}

.info-grid {
  display: grid;
  grid-template-columns: 1fr 1fr;
  gap: 10px;
  margin-bottom: 15px;
}

.info-item {
  display: flex;
  flex-direction: column;
  gap: 4px;
}

.info-label {
  font-size: 11px;
  color: #8080a0;
}

.info-value {
  font-size: 16px;
  font-weight: 600;
  color: #e0e0e0;
}

.type-counts h5 {
  color: #a0a0c0;
  font-size: 13px;
  margin-bottom: 10px;
}

.type-bars {
  display: flex;
  flex-direction: column;
  gap: 8px;
}

.type-bar-item {
  display: flex;
  align-items: center;
  gap: 10px;
}

.element-name {
  width: 30px;
  font-size: 12px;
  font-weight: 600;
  color: #8ab4f8;
}

.type-bar {
  flex: 1;
  height: 8px;
  background: rgba(40, 40, 60, 0.8);
  border-radius: 4px;
  overflow: hidden;
}

.type-bar-fill {
  height: 100%;
  transition: width 0.5s ease;
}

.type-bar-fill.Al { background: linear-gradient(90deg, #a8a8a8, #c0c0c0); }
.type-bar-fill.Co { background: linear-gradient(90deg, #4a90d9, #6ab0f9); }
.type-bar-fill.Cr { background: linear-gradient(90deg, #8b4513, #cd853f); }
.type-bar-fill.Fe { background: linear-gradient(90deg, #b87333, #daa520); }
.type-bar-fill.Ni { background: linear-gradient(90deg, #708090, #a9a9a9); }

.count-value {
  width: 60px;
  text-align: right;
  font-size: 12px;
  color: #a0a0c0;
}

.empty-state {
  background: rgba(30, 30, 50, 0.8);
  border-radius: 12px;
  padding: 80px 40px;
  text-align: center;
  border: 1px solid rgba(100, 100, 150, 0.3);
}

.empty-icon {
  font-size: 64px;
  margin-bottom: 20px;
  opacity: 0.5;
}

.empty-state h3 {
  color: #a0a0c0;
  font-size: 20px;
  margin-bottom: 10px;
}

.empty-state p {
  color: #8080a0;
  font-size: 14px;
  margin-bottom: 30px;
}

.quick-actions button {
  padding: 12px 30px;
  font-size: 14px;
}

.btn-primary {
  background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
  color: white;
  border: none;
  border-radius: 8px;
  padding: 10px 20px;
  font-weight: 500;
  cursor: pointer;
  transition: all 0.2s ease;
}

.btn-primary:hover {
  transform: translateY(-1px);
  box-shadow: 0 4px 15px rgba(102, 126, 234, 0.4);
}

.btn-primary:disabled {
  opacity: 0.5;
  cursor: not-allowed;
}

.hover-info-panel {
  background: rgba(30, 30, 50, 0.8);
  border-radius: 12px;
  padding: 20px;
  border: 1px solid rgba(100, 100, 150, 0.3);
}

.hover-info-panel h4 {
  color: #8ab4f8;
  font-size: 16px;
  margin-bottom: 15px;
}

.hover-grid {
  display: grid;
  grid-template-columns: repeat(3, 1fr);
  gap: 20px;
}

.hover-item {
  display: flex;
  flex-direction: column;
  gap: 4px;
}

.hover-item span {
  font-size: 12px;
  color: #8080a0;
}

.hover-item strong {
  font-size: 18px;
  color: #4ecdc4;
}

.analysis-cards {
  display: grid;
  grid-template-columns: repeat(3, 1fr);
  gap: 20px;
}

.analysis-card {
  background: rgba(30, 30, 50, 0.8);
  border-radius: 12px;
  padding: 20px;
  border: 1px solid rgba(100, 100, 150, 0.3);
}

.analysis-card h4 {
  color: #8ab4f8;
  font-size: 16px;
  margin-bottom: 15px;
}

.analysis-card button {
  width: 100%;
  margin-bottom: 15px;
}

.rdf-result,
.voronoi-result,
.csro-result {
  padding-top: 15px;
  border-top: 1px solid rgba(100, 100, 150, 0.3);
  font-size: 13px;
  color: #a0a0c0;
}

.rdf-card { border-left: 3px solid #667eea; }
.voronoi-card { border-left: 3px solid #4ecdc4; }
.csro-card { border-left: 3px solid #ff6b6b; }

@media (max-width: 1400px) {
  .analysis-cards {
    grid-template-columns: 1fr;
  }
  .hover-grid {
    grid-template-columns: 1fr;
  }
}
</style>
