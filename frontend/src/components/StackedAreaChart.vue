<template>
  <div class="stacked-area-chart">
    <div class="chart-header">
      <h3>多面体类型演化 - 堆叠面积图</h3>
      <div class="legend" ref="legendRef"></div>
    </div>
    <div class="chart-container" ref="chartContainer">
      <svg ref="svgRef"></svg>
      <div class="tooltip" ref="tooltipRef"></div>
      <div class="crosshair" ref="crosshairRef" v-show="showCrosshair">
        <div class="crosshair-line"></div>
        <div class="crosshair-info">
          <div>时间: {{ hoverInfo.timePs.toFixed(2) }} ps</div>
          <div>帧: {{ hoverInfo.frameIndex }}</div>
          <div v-if="csroData" class="csro-popup">
            <strong>Co-Cr CSRO:</strong> {{ csroData.S_ij.toFixed(4) }}
            <div class="csro-detail">
              <div>α: {{ csroData.alpha.toFixed(4) }}</div>
              <div>Z_avg: {{ csroData.Z_avg.toFixed(2) }}</div>
              <div>c_Co: {{ csroData.c1.toFixed(4) }}</div>
              <div>c_Cr: {{ csroData.c2.toFixed(4) }}</div>
            </div>
          </div>
          <div v-else-if="csroLoading" class="csro-loading">
            加载CSRO数据...
          </div>
        </div>
      </div>
    </div>
  </div>
</template>

<script setup>
import { ref, watch, onMounted, onUnmounted, nextTick } from 'vue'
import * as d3 from 'd3'

const props = defineProps({
  evolutionData: {
    type: Object,
    required: true
  },
  polyhedronTypes: {
    type: Array,
    default: () => []
  },
  getColor: {
    type: Function,
    required: true
  },
  binarySearch: {
    type: Function,
    required: true
  },
  queryCSRO: {
    type: Function,
    required: true
  }
})

const emit = defineEmits(['hover'])

const chartContainer = ref(null)
const svgRef = ref(null)
const legendRef = ref(null)
const tooltipRef = ref(null)
const crosshairRef = ref(null)

const showCrosshair = ref(false)
const hoverInfo = ref({ timePs: 0, frameIndex: 0 })
const csroData = ref(null)
const csroLoading = ref(false)
let csroRequestId = 0

let svg = null
let xScale = null
let yScale = null
let colorScale = null
let stack = null
let area = null
let xAxisGroup = null
let yAxisGroup = null
let bisector = null

const margin = { top: 20, right: 150, bottom: 60, left: 80 }
let width = 0
let height = 0

const initChart = () => {
  if (!chartContainer.value || !svgRef.value) return
  const containerRect = chartContainer.value.getBoundingClientRect()
  width = containerRect.width - margin.left - margin.right
  height = 500 - margin.top - margin.bottom
  svg = d3.select(svgRef.value)
    .attr('width', containerRect.width)
    .attr('height', 500)
  const g = svg.append('g')
    .attr('transform', `translate(${margin.left},${margin.top})`)
  xScale = d3.scaleLinear().range([0, width])
  yScale = d3.scaleLinear().range([height, 0])
  colorScale = d3.scaleOrdinal()
  stack = d3.stack()
    .offset(d3.stackOffsetNone)
    .order(d3.stackOrderNone)
  area = d3.area()
    .x(d => xScale(d.data.timePs))
    .y0(d => yScale(d[0]))
    .y1(d => yScale(d[1]))
    .curve(d3.curveMonotoneX)
  xAxisGroup = g.append('g')
    .attr('class', 'x-axis')
    .attr('transform', `translate(0,${height})`)
  yAxisGroup = g.append('g')
    .attr('class', 'y-axis')
  g.append('text')
    .attr('class', 'x-label')
    .attr('x', width / 2)
    .attr('y', height + 45)
    .attr('text-anchor', 'middle')
    .attr('fill', '#a0a0c0')
    .text('模拟时间 (ps)')
  g.append('text')
    .attr('class', 'y-label')
    .attr('transform', 'rotate(-90)')
    .attr('x', -height / 2)
    .attr('y', -55)
    .attr('text-anchor', 'middle')
    .attr('fill', '#a0a0c0')
    .text('原子计数')
  const overlay = g.append('rect')
    .attr('class', 'overlay')
    .attr('width', width)
    .attr('height', height)
    .attr('fill', 'none')
    .attr('pointer-events', 'all')
    .on('mousemove', handleMouseMove)
    .on('mouseleave', handleMouseLeave)
  bisector = d3.bisector(d => d.timePs).center
  updateLegend()
}

const updateLegend = () => {
  if (!legendRef.value || !props.polyhedronTypes.length) return
  const legend = d3.select(legendRef.value).html('')
  const items = legend.selectAll('.legend-item')
    .data(props.polyhedronTypes)
    .enter()
    .append('div')
    .attr('class', 'legend-item')
  items.append('span')
    .attr('class', 'legend-color')
    .style('background-color', d => props.getColor(d))
  items.append('span')
    .attr('class', 'legend-text')
    .text(d => d)
}

const updateChart = () => {
  if (!svg || !props.evolutionData.timePs.length) return
  const timePs = props.evolutionData.timePs
  const counts = props.evolutionData.counts
  const types = props.polyhedronTypes.length > 0
    ? props.polyhedronTypes
    : Object.keys(counts)
  if (types.length === 0) return
  const data = timePs.map((t, i) => {
    const row = { timePs: t }
    types.forEach(type => {
      row[type] = counts[type] ? counts[type][i] || 0 : 0
    })
    return row
  })
  colorScale.domain(types).range(types.map(t => props.getColor(t)))
  xScale.domain(d3.extent(timePs))
  stack.keys(types)
  const stacked = stack(data)
  const maxY = d3.max(stacked, layer => d3.max(layer, d => d[1]))
  yScale.domain([0, maxY || 100])
  const g = svg.select('g')
  const paths = g.selectAll('.area-path')
    .data(stacked, d => d.key)
  paths.exit().remove()
  const pathsEnter = paths.enter()
    .append('path')
    .attr('class', 'area-path')
  pathsEnter.merge(paths)
    .transition()
    .duration(300)
    .attr('d', area)
    .attr('fill', d => colorScale(d.key))
    .attr('opacity', 0.85)
    .style('mix-blend-mode', 'normal')
  const xAxis = d3.axisBottom(xScale)
    .ticks(10)
    .tickFormat(d => d.toFixed(1))
  xAxisGroup.transition().duration(300).call(xAxis)
    .selectAll('text')
    .attr('fill', '#a0a0c0')
  xAxisGroup.selectAll('line').attr('stroke', '#4a4a6a')
  xAxisGroup.select('.domain').attr('stroke', '#4a4a6a')
  const yAxis = d3.axisLeft(yScale)
    .ticks(10)
    .tickFormat(d3.format('~s'))
  yAxisGroup.transition().duration(300).call(yAxis)
    .selectAll('text')
    .attr('fill', '#a0a0c0')
  yAxisGroup.selectAll('line').attr('stroke', '#4a4a6a')
  yAxisGroup.select('.domain').attr('stroke', '#4a4a6a')
  g.selectAll('.grid-line').remove()
  const gridLines = g.append('g').attr('class', 'grid-line')
  yScale.ticks(10).forEach(tick => {
    gridLines.append('line')
      .attr('x1', 0)
      .attr('x2', width)
      .attr('y1', yScale(tick))
      .attr('y2', yScale(tick))
      .attr('stroke', '#2a2a4a')
      .attr('stroke-dasharray', '2,2')
  })
  updateLegend()
}

const handleMouseMove = (event) => {
  if (!xScale || !bisector || !props.evolutionData.timePs.length) return
  const [mx] = d3.pointer(event)
  const x0 = xScale.invert(mx)
  const data = props.evolutionData.timePs.map((t, i) => ({ timePs: t, index: i }))
  const idx = bisector(data, x0)
  const d = data[idx]
  if (!d) return
  showCrosshair.value = true
  const crosshairX = xScale(d.timePs) + margin.left
  d3.select(crosshairRef.value)
    .style('left', `${crosshairX}px`)
  hoverInfo.value = {
    timePs: d.timePs,
    frameIndex: d.index
  }
  fetchCSRO(props.evolutionData.timesteps[d.index])
  emit('hover', {
    timePs: d.timePs,
    frameIndex: d.index,
    timestep: props.evolutionData.timesteps[d.index]
  })
}

const fetchCSRO = async (timestep) => {
  const requestId = ++csroRequestId
  csroLoading.value = true
  csroData.value = null
  try {
    const result = await props.queryCSRO(timestep, 2, 3)
    if (requestId === csroRequestId) {
      csroData.value = result
      csroLoading.value = false
    }
  } catch (e) {
    if (requestId === csroRequestId) {
      csroLoading.value = false
    }
  }
}

const handleMouseLeave = () => {
  showCrosshair.value = false
  csroData.value = null
  csroLoading.value = false
}

const handleResize = () => {
  if (svgRef.value) {
    d3.select(svgRef.value).selectAll('*').remove()
    initChart()
    updateChart()
  }
}

watch(() => props.evolutionData.timePs.length, () => {
  nextTick(() => {
    updateChart()
  })
}, { deep: true })

watch(() => props.polyhedronTypes, () => {
  nextTick(() => {
    updateChart()
  })
}, { deep: true })

onMounted(() => {
  initChart()
  updateChart()
  window.addEventListener('resize', handleResize)
})

onUnmounted(() => {
  window.removeEventListener('resize', handleResize)
})
</script>

<style scoped>
.stacked-area-chart {
  background: rgba(30, 30, 50, 0.8);
  border-radius: 12px;
  padding: 20px;
  border: 1px solid rgba(100, 100, 150, 0.3);
}

.chart-header {
  display: flex;
  justify-content: space-between;
  align-items: flex-start;
  margin-bottom: 15px;
  flex-wrap: wrap;
  gap: 15px;
}

.chart-header h3 {
  color: #8ab4f8;
  font-size: 18px;
  font-weight: 600;
}

.legend {
  display: flex;
  flex-wrap: wrap;
  gap: 12px;
  max-width: 600px;
}

.legend-item {
  display: flex;
  align-items: center;
  gap: 6px;
  font-size: 12px;
  color: #b0b0d0;
}

.legend-color {
  width: 14px;
  height: 14px;
  border-radius: 3px;
  display: inline-block;
}

.legend-text {
  font-family: 'Courier New', monospace;
}

.chart-container {
  position: relative;
  width: 100%;
}

.overlay {
  cursor: crosshair;
}

.x-axis text,
.y-axis text {
  font-size: 11px;
}

.crosshair {
  position: absolute;
  top: 0;
  pointer-events: none;
  z-index: 100;
}

.crosshair-line {
  width: 1px;
  height: 500px;
  background: linear-gradient(to bottom, transparent, #ff6b6b, transparent);
  box-shadow: 0 0 10px rgba(255, 107, 107, 0.5);
}

.crosshair-info {
  position: absolute;
  top: 20px;
  left: 15px;
  background: rgba(20, 20, 40, 0.95);
  border: 1px solid #ff6b6b;
  border-radius: 8px;
  padding: 12px 16px;
  font-size: 13px;
  color: #e0e0e0;
  min-width: 200px;
  box-shadow: 0 4px 20px rgba(0, 0, 0, 0.5);
  backdrop-filter: blur(10px);
}

.csro-popup {
  margin-top: 10px;
  padding-top: 10px;
  border-top: 1px solid rgba(255, 107, 107, 0.3);
}

.csro-popup strong {
  color: #4ecdc4;
  font-size: 15px;
}

.csro-detail {
  margin-top: 6px;
  font-size: 11px;
  color: #a0a0c0;
  display: grid;
  grid-template-columns: 1fr 1fr;
  gap: 4px;
}

.csro-loading {
  margin-top: 10px;
  padding-top: 10px;
  border-top: 1px solid rgba(255, 107, 107, 0.3);
  color: #ffd93d;
  font-size: 12px;
  animation: pulse 1.5s infinite;
}

@keyframes pulse {
  0%, 100% { opacity: 0.6; }
  50% { opacity: 1; }
}

.area-path {
  transition: opacity 0.2s ease;
}

.area-path:hover {
  opacity: 1;
}
</style>
