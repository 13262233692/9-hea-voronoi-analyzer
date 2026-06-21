import axios from 'axios'

const apiClient = axios.create({
  baseURL: '/api',
  headers: {
    'Content-Type': 'application/json'
  },
  timeout: 300000
})

export function useApi() {
  const loadData = async (filepath) => {
    const response = await apiClient.post('/data/load', { filepath })
    return response.data
  }

  const getDataInfo = async () => {
    const response = await apiClient.get('/data/info')
    return response.data
  }

  const getFrameInfo = async (frameIdx) => {
    const response = await apiClient.get(`/frame/${frameIdx}`)
    return response.data
  }

  const getRDF = async (frameIdx, rMax = 6.0, nBins = 100) => {
    const response = await apiClient.get(`/rdf/${frameIdx}`, {
      params: { r_max: rMax, n_bins: nBins }
    })
    return response.data
  }

  const getRDFAverage = async (startFrame, endFrame, rMax = 6.0, nBins = 100) => {
    const response = await apiClient.get('/rdf/average', {
      params: { start_frame: startFrame, end_frame: endFrame, r_max: rMax, n_bins: nBins }
    })
    return response.data
  }

  const getVoronoi = async (frameIdx) => {
    const response = await apiClient.get(`/voronoi/${frameIdx}`)
    return response.data
  }

  const getVoronoiEvolution = async (startFrame, endFrame, types = null) => {
    const params = { start_frame: startFrame, end_frame: endFrame }
    if (types) params.types = types
    const response = await apiClient.get('/voronoi/evolution', { params })
    return response.data
  }

  const getCSRO = async (frameIdx, type1, type2) => {
    const response = await apiClient.get(`/csro/${frameIdx}`, {
      params: { type1, type2 }
    })
    return response.data
  }

  const getCSROAllPairs = async (frameIdx) => {
    const response = await apiClient.get(`/csro/all/${frameIdx}`)
    return response.data
  }

  const queryCSROByTimestep = async (timestep, type1 = 2, type2 = 3) => {
    const response = await apiClient.post('/csro/query', { timestep, type1, type2 })
    return response.data
  }

  const getAtomNeighbors = async (frameIdx, atomIdx, cutoff = 3.5) => {
    const response = await apiClient.get(`/neighbors/${frameIdx}/${atomIdx}`, {
      params: { cutoff }
    })
    return response.data
  }

  const getNonaffineAnalysis = async (frameIdxRef, frameIdxDef) => {
    const response = await apiClient.get(`/nonaffine/${frameIdxRef}/${frameIdxDef}`)
    return response.data
  }

  const getYieldEvents = async () => {
    const response = await apiClient.get('/yield-events')
    return response.data
  }

  const healthCheck = async () => {
    const response = await apiClient.get('/health')
    return response.data
  }

  return {
    loadData,
    getDataInfo,
    getFrameInfo,
    getRDF,
    getRDFAverage,
    getVoronoi,
    getVoronoiEvolution,
    getCSRO,
    getCSROAllPairs,
    queryCSROByTimestep,
    getAtomNeighbors,
    getNonaffineAnalysis,
    getYieldEvents,
    healthCheck
  }
}
