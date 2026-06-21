import { ref, reactive, onUnmounted } from 'vue'

export function useEvolutionStream() {
  const isConnected = ref(false)
  const isStreaming = ref(false)
  const currentFrame = ref(0)
  const totalFrames = ref(0)
  const evolutionData = reactive({
    timesteps: [],
    timePs: [],
    counts: {},
    polyhedronTypes: []
  })
  const error = ref(null)
  let ws = null

  const colorScale = {
    '<0,0,12,0>': '#FF6B6B',
    '<0,0,12,2>': '#4ECDC4',
    '<0,2,8,2>': '#45B7D1',
    '<0,3,6,4>': '#96CEB4',
    '<0,4,4,6>': '#FFEAA7',
    '<0,1,10,2>': '#DDA0DD',
    '<0,1,10,3>': '#98D8C8',
    '<0,2,8,1>': '#F7DC6F',
    '<0,2,8,4>': '#BB8FCE',
    '<0,3,6,3>': '#85C1E9'
  }

  const getColor = (type) => {
    if (colorScale[type]) return colorScale[type]
    const hash = type.split('').reduce((acc, char) => acc + char.charCodeAt(0), 0)
    const hue = hash % 360
    return `hsl(${hue}, 70%, 60%)`
  }

  const connect = (options = {}) => {
    const {
      startFrame = 0,
      endFrame = null,
      types = null,
      intervalMs = 200
    } = options

    if (ws) {
      ws.close()
    }

    const proto = window.location.protocol === 'https:' ? 'wss:' : 'ws:'
    let url = `${proto}//${window.location.host}/ws/evolution?start_frame=${startFrame}&interval_ms=${intervalMs}`
    if (endFrame !== null) {
      url += `&end_frame=${endFrame}`
    }
    if (types) {
      url += `&types=${encodeURIComponent(types)}`
    }

    ws = new WebSocket(url)

    ws.onopen = () => {
      isConnected.value = true
      isStreaming.value = true
      error.value = null
    }

    ws.onmessage = (event) => {
      const data = JSON.parse(event.data)
      if (data.error) {
        error.value = data.error
        isStreaming.value = false
        return
      }
      if (data.done) {
        isStreaming.value = false
        totalFrames.value = data.frames_processed
        return
      }
      currentFrame.value = data.frame_idx + 1
      evolutionData.timesteps.push(data.timestep)
      evolutionData.timePs.push(data.time_ps)
      for (const [type, count] of Object.entries(data.counts)) {
        if (!evolutionData.counts[type]) {
          evolutionData.counts[type] = []
          evolutionData.polyhedronTypes.push(type)
        }
        evolutionData.counts[type].push(count)
      }
      for (const type of evolutionData.polyhedronTypes) {
        if (!data.counts[type]) {
          evolutionData.counts[type].push(0)
        }
      }
    }

    ws.onerror = (event) => {
      error.value = 'WebSocket connection error'
      isConnected.value = false
      isStreaming.value = false
    }

    ws.onclose = () => {
      isConnected.value = false
      isStreaming.value = false
    }
  }

  const disconnect = () => {
    if (ws) {
      ws.close()
      ws = null
    }
  }

  const reset = () => {
    disconnect()
    currentFrame.value = 0
    totalFrames.value = 0
    evolutionData.timesteps = []
    evolutionData.timePs = []
    evolutionData.counts = {}
    evolutionData.polyhedronTypes = []
    error.value = null
  }

  const binarySearchTimestep = (targetPs) => {
    const times = evolutionData.timePs
    let low = 0
    let high = times.length - 1
    while (low <= high) {
      const mid = Math.floor((low + high) / 2)
      if (Math.abs(times[mid] - targetPs) < 0.001) {
        return { index: mid, timestep: evolutionData.timesteps[mid], timePs: times[mid] }
      }
      if (times[mid] < targetPs) {
        low = mid + 1
      } else {
        high = mid - 1
      }
    }
    if (high < 0) return { index: 0, timestep: evolutionData.timesteps[0], timePs: times[0] }
    if (low >= times.length) return { index: times.length - 1, timestep: evolutionData.timesteps[times.length - 1], timePs: times[times.length - 1] }
    if (Math.abs(targetPs - times[high]) < Math.abs(targetPs - times[low])) {
      return { index: high, timestep: evolutionData.timesteps[high], timePs: times[high] }
    }
    return { index: low, timestep: evolutionData.timesteps[low], timePs: times[low] }
  }

  onUnmounted(() => {
    disconnect()
  })

  return {
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
  }
}
