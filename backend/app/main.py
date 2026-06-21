from fastapi import FastAPI, WebSocket, WebSocketDisconnect, HTTPException
from fastapi.middleware.cors import CORSMiddleware
import json
import asyncio
from pydantic import BaseModel
from typing import List, Optional
from .services.analysis_service import AnalysisService
from .core.md_unpacker import create_test_md_file
from .core.models import (
    LoadDataRequest, FrameInfoResponse, RDFResponse,
    VoronoiResponse, VoronoiEvolutionResponse, CSROResponse,
    EvolutionStreamData, NeighborsResponse, CSROQueryRequest
)


class GenerateTestDataRequest(BaseModel):
    filepath: str
    n_frames: int = 10
    n_atoms: int = 100000
    triclinic: bool = False
    tilt_factor: float = 0.3

app = FastAPI(title="AlCoCrFeNi HEA Microstructure Analysis Workstation",
              description="High-entropy alloy microstructure analysis backend",
              version="1.0.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

analysis_service = AnalysisService()


@app.on_event("shutdown")
async def shutdown_event():
    analysis_service.close()


@app.post("/api/data/load")
async def load_data(request: LoadDataRequest):
    try:
        result = analysis_service.load_data(request.filepath)
        return {"status": "success", "data": result}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.get("/api/data/info")
async def get_data_info():
    if not analysis_service.unpacker:
        raise HTTPException(status_code=400, detail="No data loaded")
    return {
        "n_frames": analysis_service.unpacker.n_frames,
        "current_file": analysis_service.current_file
    }


@app.get("/api/frame/{frame_idx}", response_model=FrameInfoResponse)
async def get_frame_info(frame_idx: int):
    try:
        return analysis_service.get_frame_info(frame_idx)
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.get("/api/rdf/{frame_idx}", response_model=RDFResponse)
async def get_rdf(frame_idx: int, r_max: float = 6.0, n_bins: int = 100):
    try:
        return analysis_service.calculate_rdf(frame_idx, r_max, n_bins)
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.get("/api/rdf/average")
async def get_rdf_average(start_frame: int, end_frame: int, r_max: float = 6.0, n_bins: int = 100):
    try:
        return analysis_service.calculate_rdf_average(start_frame, end_frame, r_max, n_bins)
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.get("/api/voronoi/{frame_idx}", response_model=VoronoiResponse)
async def get_voronoi(frame_idx: int):
    try:
        return analysis_service.analyze_voronoi(frame_idx)
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.get("/api/voronoi/evolution", response_model=VoronoiEvolutionResponse)
async def get_voronoi_evolution(start_frame: int, end_frame: int, types: str = None):
    try:
        poly_types = types.split(',') if types else None
        return analysis_service.get_voronoi_evolution(start_frame, end_frame, poly_types)
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.get("/api/csro/{frame_idx}", response_model=CSROResponse)
async def get_csro(frame_idx: int, type1: int, type2: int):
    try:
        return analysis_service.calculate_csro(frame_idx, type1, type2)
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.get("/api/csro/all/{frame_idx}")
async def get_csro_all_pairs(frame_idx: int):
    try:
        return analysis_service.calculate_csro_all_pairs(frame_idx)
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.post("/api/csro/query")
async def query_csro_by_timestep(request: CSROQueryRequest):
    try:
        frame_idx = analysis_service.find_frame_by_timestep(request.timestep)
        result = analysis_service.calculate_csro(frame_idx, request.type1, request.type2)
        result['frame_idx'] = frame_idx
        result['actual_timestep'] = analysis_service.get_frame(frame_idx).timestep
        return result
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.get("/api/neighbors/{frame_idx}/{atom_idx}", response_model=NeighborsResponse)
async def get_atom_neighbors(frame_idx: int, atom_idx: int, cutoff: float = 3.5):
    try:
        return analysis_service.get_atom_neighbors(frame_idx, atom_idx, cutoff)
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.get("/api/nonaffine/{frame_idx_ref}/{frame_idx_def}")
async def get_nonaffine_analysis(frame_idx_ref: int, frame_idx_def: int):
    try:
        return analysis_service.analyze_nonaffine(frame_idx_ref, frame_idx_def)
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.get("/api/yield-events")
async def get_yield_events():
    try:
        return {"events": analysis_service.get_yield_events()}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.websocket("/ws/evolution")
async def websocket_evolution(websocket: WebSocket,
                              start_frame: int = 0,
                              end_frame: int = None,
                              types: str = None,
                              interval_ms: int = 200):
    await websocket.accept()
    try:
        if not analysis_service.unpacker:
            await websocket.send_json({"error": "No data loaded"})
            await websocket.close()
            return
        if end_frame is None:
            end_frame = analysis_service.unpacker.n_frames - 1
        poly_types = types.split(',') if types else None
        for data in analysis_service.get_evolution_stream(start_frame, end_frame, poly_types):
            await websocket.send_json(data)
            await asyncio.sleep(interval_ms / 1000.0)
        await websocket.send_json({"done": True, "frames_processed": end_frame - start_frame + 1})
    except WebSocketDisconnect:
        pass
    except Exception as e:
        try:
            await websocket.send_json({"error": str(e)})
        except:
            pass
    finally:
        try:
            await websocket.close()
        except:
            pass


@app.websocket("/ws/nonaffine-monitor")
async def websocket_nonaffine_monitor(websocket: WebSocket,
                                      start_frame: int = 0,
                                      end_frame: int = None,
                                      yield_threshold: float = 0.08,
                                      interval_ms: int = 200):
    await websocket.accept()
    try:
        if not analysis_service.unpacker:
            await websocket.send_json({"error": "No data loaded"})
            await websocket.close()
            return
        if end_frame is None:
            end_frame = analysis_service.unpacker.n_frames - 1
        for data in analysis_service.get_nonaffine_stream(
            start_frame, end_frame, yield_threshold
        ):
            await websocket.send_json(data)
            await asyncio.sleep(interval_ms / 1000.0)
        await websocket.send_json({
            "done": True,
            "frames_processed": end_frame - start_frame,
            "total_yield_events": len(analysis_service.get_yield_events())
        })
    except WebSocketDisconnect:
        pass
    except Exception as e:
        try:
            await websocket.send_json({"error": str(e)})
        except:
            pass
    finally:
        try:
            await websocket.close()
        except:
            pass


@app.post("/api/generate-test-data")
async def generate_test_data(request: GenerateTestDataRequest):
    try:
        import os
        os.makedirs(os.path.dirname(os.path.abspath(request.filepath)), exist_ok=True)
        create_test_md_file(
            request.filepath,
            request.n_frames,
            request.n_atoms,
            triclinic=request.triclinic,
            tilt_factor=request.tilt_factor
        )
        return {
            "status": "success",
            "filepath": request.filepath,
            "n_frames": request.n_frames,
            "n_atoms": request.n_atoms,
            "triclinic": request.triclinic,
            "tilt_factor": request.tilt_factor
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.get("/api/health")
async def health_check():
    return {"status": "ok", "service": "HEA Analysis Backend"}
