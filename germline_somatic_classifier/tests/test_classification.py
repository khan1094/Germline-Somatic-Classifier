import math
from somatic_variant_classifier import (
    compute_germline_score,
    compute_somatic_score,
    classify_variant
)

def test_germline_high_af_high_an():
    score = compute_germline_score(0.5, 200_000, 100_000)
    assert score == 0.5

def test_germline_high_af_low_an():
    score = compute_germline_score(0.5, 50_000, 100_000)
    assert score == 0.25

def test_germline_missing_af():
    score = compute_germline_score(None, 100_000, 100_000)
    assert score == 0.0

def test_germline_missing_an():
    score = compute_germline_score(0.5, None, 100_000)
    assert score == 0.0

def test_somatic_zero():
    assert compute_somatic_score(0) == 0.0

def test_somatic_one():
    assert compute_somatic_score(1) == math.log10(2)

def test_somatic_large():
    score = compute_somatic_score(100)
    assert score > 1.0

GERM_THR = 0.1
SOM_THR = 1.0

def test_likely_germline():
    c = classify_variant(0.5, 0.2, GERM_THR, SOM_THR)
    assert c == "LIKELY_GERMLINE"

def test_likely_somatic():
    c = classify_variant(0.01, 2.0, GERM_THR, SOM_THR)
    assert c == "LIKELY_SOMATIC"

def test_conflicting():
    c = classify_variant(0.5, 2.0, GERM_THR, SOM_THR)
    assert c == "CONFLICTING"

def test_unknown():
    c = classify_variant(0.01, 0.2, GERM_THR, SOM_THR)
    assert c == "UNKNOWN"
