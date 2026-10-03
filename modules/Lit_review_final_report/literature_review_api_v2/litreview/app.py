from fastapi import FastAPI
from .models import ReviewRequest, ReviewResponse
from .service import run_review

app = FastAPI(title="JCAP Literature Review API", version="2.0.0")


@app.get("/health")
def health():
    return {"status": "ok"}


@app.post("/review", response_model=ReviewResponse)
async def review(req: ReviewRequest):
    return await run_review(req)
