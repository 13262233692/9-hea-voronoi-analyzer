import { ref, onUnmounted } from 'vue'

export function useNonaffineMonitor() {
  const isMonitoring = ref(false)
  const isPaused = ref(false)
  const currentFrame = ref(-1)
  const totalFrames = ref(0)
  const yieldEvents = ref([])
  const latestNonaffine = ref(null)
  const error = ref(null)
  const connectionStatus = ref('disconnected')

  let ws = null
  let isConnecting = false

  const connect = (options = {}) => {
    if (isConnecting || (ws && ws.readyState === WebSocket.OPEN)) return

    isConnecting = true
    connectionStatus.value = 'connecting'

    const {
      startFrame = 0,
      endFrame = null,
      yieldThreshold = 0.08,
      intervalMs = 200
    } = options

    let url = `${window.location.protocol === 'https:' ? 'wss:' : 'ws:'}//${window.location.host}/ws/nonaffine-monitor?`
    url += `start_frame=${startFrame}&yield_threshold=${yieldThreshold}&interval_ms=${intervalMs}`
    if (endFrame !== null) url += `&end_frame=${endFrame}`

    try {
      ws = new WebSocket(url)
    } catch (e) {
      error.value = e.message
      isConnecting = false
      connectionStatus.value = 'error'
      return
    }

    ws.onopen = () => {
      isConnecting = false
      isMonitoring.value = true
      connectionStatus.value = 'connected'
    }

    ws.onmessage = (event) => {
      try {
        const data = JSON.parse(event.data)

        if (data.type === 'nonaffine') {
          latestNonaffine.value = data
          currentFrame.value = data.frame_idx
        } else if (data.type === 'yield_event') {
          yieldEvents.value.push(data)
        } else if (data.done) {
          totalFrames.value = data.frames_processed
        } else if (data.error) {
          error.value = data.error
        }
      } catch (e) {
        console.warn('Failed to parse WS message:', e)
      }
    }

    ws.onclose = () => {
      isMonitoring.value = false
      isConnecting = false
      connectionStatus.value = 'disconnected'
    }

    ws.onerror = (e) => {
      error.value = 'WebSocket error'
      connectionStatus.value = 'error'
    }
  }

  const disconnect = () => {
    if (ws) {
      try {
        ws.close()
      } catch (e) {
        // ignore
      }
      ws = null
    }
    isMonitoring.value = false
    isConnecting = false
    connectionStatus.value = 'disconnected'
  }

  const reset = () => {
    disconnect()
    yieldEvents.value = []
    latestNonaffine.value = null
    currentFrame.value = -1
    error.value = null
  }

  onUnmounted(() => {
    disconnect()
  })

  return {
    isMonitoring,
    isPaused,
    currentFrame,
    totalFrames,
    yieldEvents,
    latestNonaffine,
    error,
    connectionStatus,
    connect,
    disconnect,
    reset
  }
}
