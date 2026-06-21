<template>
  <div class="control-panel">
    <h3>系统控制</h3>
    <div class="control-section">
      <label>数据文件路径</label>
      <div class="file-input-group">
        <input type="text" v-model="filepath" placeholder="输入MD数据文件路径..." />
        <button @click="handleLoadData" :disabled="loading">
          {{ loading ? '加载中...' : '加载数据' }}
        </button>
      </div>
      <div v-if="dataInfo" class="data-info">
        <span>帧数: {{ dataInfo.n_frames }}</span>
        <span>每帧原子数: {{ frameInfo?.n_atoms?.toLocaleString() }}</span>
      </div>
    </div>
    <div class="control-section">
      <label>多面体类型筛选</label>
      <div class="type-selector">
        <label v-for="type in polyhedronTypes" :key="type">
          <input type="checkbox" v-model="selectedTypes" :value="type" />
          <span :style="{ color: getColor(type) }">{{ type }}</span>
        </label>
      </div>
    </div>
    <div class="control-section">
      <label>流式参数</label>
      <div class="param-grid">
        <div class="param-item">
          <span>起始帧</span>
          <input type="number" v-model.number="startFrame" min="0" :max="dataInfo?.n_frames - 1" />
        </div>
        <div class="param-item">
          <span>结束帧</span>
          <input type="number" v-model.number="endFrame" min="0" :max="dataInfo?.n_frames - 1" />
        </div>
        <div class="param-item">
          <span>间隔 (ms)</span>
          <input type="number" v-model.number="intervalMs" min="50" max="2000" />
        </div>
      </div>
    </div>
    <div class="control-section">
      <label>回放控制</label>
      <div class="button-group">
        <button @click="handleStartStream" :disabled="isStreaming || !dataInfo" class="btn-primary">
          {{ isStreaming ? '播放中...' : '开始播放' }}
        </button>
        <button @click="handleStopStream" :disabled="!isStreaming" class="btn-secondary">
          停止
        </button>
        <button @click="handleReset" class="btn-danger">
          重置
        </button>
        <button @click="generateTestData" class="btn-info">
          生成测试数据
        </button>
      </div>
    </div>
    <div class="status-section">
      <div class="status-item">
      <span class="status-label">连接状态</span>
      <span class="status-value" :class="{ connected: isConnected }">
        {{ isConnected ? '已连接' : '未连接' }}
      </span>
    </div>
    <div class="status-item">
      <span class="status-label">当前帧</span>
      <span class="status-value">{{ currentFrame }} / {{ totalFrames || dataInfo?.n_frames }}</span>
    </div>
    <div class="progress-bar">
      <div class="progress-fill" :style="{ width: progressPercent + '%' }"></div>
    </div>
  </div>
    <div v-if="error" class="error-message">
      {{ error }}
    </div>
  </div>
</template>

<script setup>
import { ref, computed, watch } from 'vue'

const props = defineProps({
  isConnected: Boolean,
  isStreaming: Boolean,
  currentFrame: Number,
  totalFrames: Number,
  dataInfo: Object,
  frameInfo: Object,
  polyhedronTypes: {
    type: Array,
    default: () => []
  },
  getColor: {
    type: Function,
    required: true
  },
  error: String
})

const emit = defineEmits(['load-data', 'start-stream', 'stop-stream', 'reset', 'generate-test-data'])

const filepath = ref('d:/SOLO-0621-1/9-hea-voronoi-analyzer/backend/data/test_heau.md')
const loading = ref(false)

const startFrame = ref(0)
const endFrame = ref(9)
const intervalMs = ref(300)
const selectedTypes = ref([])

const progressPercent = computed(() => {
  const total = props.totalFrames || props.dataInfo?.n_frames || 1
  return Math.min(100, (props.currentFrame / total) * 100)
})

watch(() => props.dataInfo, (info) => {
  if (info) {
    endFrame.value = Math.min(9, info.n_frames - 1)
  }
}, { immediate: true })

const handleLoadData = async () => {
  loading.value = true
  try {
    await emit('load-data', filepath.value)
  } finally {
    loading.value = false
  }
}

const handleStartStream = () => {
  emit('start-stream', {
    startFrame: startFrame.value,
    endFrame: endFrame.value,
    types: selectedTypes.value.length > 0 ? selectedTypes.value.join(',') : null,
    intervalMs: intervalMs.value
  })
}

const handleStopStream = () => {
  emit('stop-stream')
}

const handleReset = () => {
  emit('reset')
}

const generateTestData = () => {
  emit('generate-test-data')
}
</script>

<style scoped>
.control-panel {
  background: rgba(30, 30, 50, 0.8);
  border-radius: 12px;
  padding: 20px;
  border: 1px solid rgba(100, 100, 150, 0.3);
}

.control-panel h3 {
  color: #8ab4f8;
  font-size: 18px;
  margin-bottom: 20px;
  padding-bottom: 10px;
  border-bottom: 1px solid rgba(100, 100, 150, 0.3);
}

.control-section {
  margin-bottom: 20px;
}

.control-section label {
  display: block;
  color: #a0a0c0;
  font-size: 13px;
  margin-bottom: 8px;
  font-weight: 500;
}

.file-input-group {
  display: flex;
  gap: 10px;
}

.file-input-group input {
  flex: 1;
  padding: 10px 12px;
  background: rgba(20, 20, 40, 0.8);
  border: 1px solid rgba(100, 100, 150, 0.3);
  border-radius: 6px;
  color: #e0e0e0;
  font-size: 13px;
}

.file-input-group input:focus {
  outline: none;
  border-color: #8ab4f8;
}

button {
  padding: 10px 20px;
  border: none;
  border-radius: 6px;
  font-size: 13px;
  font-weight: 500;
  cursor: pointer;
  transition: all 0.2s ease;
}

button:disabled {
  opacity: 0.5;
  cursor: not-allowed;
}

.btn-primary {
  background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
  color: white;
}

.btn-primary:hover:not(:disabled) {
  transform: translateY(-1px);
  box-shadow: 0 4px 15px rgba(102, 126, 234, 0.4);
}

.btn-secondary {
  background: rgba(100, 100, 150, 0.3);
  color: #e0e0e0;
}

.btn-secondary:hover:not(:disabled) {
  background: rgba(100, 100, 150, 0.5);
}

.btn-danger {
  background: rgba(255, 107, 107, 0.2);
  color: #ff6b6b;
  border: 1px solid rgba(255, 107, 107, 0.5);
}

.btn-danger:hover:not(:disabled) {
  background: rgba(255, 107, 107, 0.3);
}

.btn-info {
  background: rgba(78, 205, 196, 0.2);
  color: #4ecdc4;
  border: 1px solid rgba(78, 205, 196, 0.5);
}

.btn-info:hover:not(:disabled) {
  background: rgba(78, 205, 196, 0.3);
}

.data-info {
  display: flex;
  gap: 20px;
  margin-top: 10px;
  padding: 10px;
  background: rgba(20, 20, 40, 0.6);
  border-radius: 6px;
  font-size: 12px;
  color: #a0a0c0;
}

.type-selector {
  display: flex;
  flex-wrap: wrap;
  gap: 15px;
}

.type-selector label {
  display: flex;
  align-items: center;
  gap: 6px;
  font-size: 12px;
  cursor: pointer;
  font-family: 'Courier New', monospace;
}

.type-selector input {
  cursor: pointer;
}

.param-grid {
  display: grid;
  grid-template-columns: repeat(3, 1fr);
  gap: 10px;
}

.param-item {
  display: flex;
  flex-direction: column;
  gap: 5px;
}

.param-item span {
  font-size: 11px;
  color: #8080a0;
}

.param-item input {
  padding: 8px 10px;
  background: rgba(20, 20, 40, 0.8);
  border: 1px solid rgba(100, 100, 150, 0.3);
  border-radius: 6px;
  color: #e0e0e0;
  font-size: 13px;
  width: 100%;
}

.button-group {
  display: flex;
  gap: 10px;
  flex-wrap: wrap;
}

.status-section {
  background: rgba(20, 20, 40, 0.6);
  border-radius: 8px;
  padding: 15px;
}

.status-item {
  display: flex;
  justify-content: space-between;
  margin-bottom: 10px;
  font-size: 13px;
}

.status-label {
  color: #8080a0;
}

.status-value {
  color: #e0e0e0;
  font-weight: 500;
}

.status-value.connected {
  color: #4ecdc4;
}

.progress-bar {
  height: 6px;
  background: rgba(40, 40, 60, 0.8);
  border-radius: 3px;
  overflow: hidden;
  margin-top: 10px;
}

.progress-fill {
  height: 100%;
  background: linear-gradient(90deg, #667eea 0%, #764ba2 100%);
  transition: width 0.3s ease;
}

.error-message {
  margin-top: 15px;
  padding: 12px;
  background: rgba(255, 107, 107, 0.1);
  border: 1px solid rgba(255, 107, 107, 0.3);
  border-radius: 6px;
  color: #ff6b6b;
  font-size: 13px;
}
</style>
