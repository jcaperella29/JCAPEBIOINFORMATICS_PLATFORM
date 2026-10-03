from fastapi import FastAPI

from models import CellStateRequest, CellStateResponse
from service import analyze_cell_states


app = FastAPI(
    title="JCAP Cell-State / Lineage API",
    version="0.1.0",
)


@app.get("/health")
def health():
    return {
        "status": "ok",
        "service": "cell_state_api",
    }


@app.post("/analyze", response_model=CellStateResponse)
def analyze(req: CellStateRequest):
    return analyze_cell_states(req)