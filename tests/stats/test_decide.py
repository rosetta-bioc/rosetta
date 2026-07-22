import pytest

try:
    import rpy2.robjects as ro
    HAS_RPY2 = True
except Exception:
    HAS_RPY2 = False

pytestmark = pytest.mark.skipif(not HAS_RPY2, reason="rpy2 not available or incompatible")

if HAS_RPY2:
    from rosetta.stats.decide import run_decide_tests

def test_run_decide_tests_success(limma_fit_object):
    # 1. Execute wrapper
    results = run_decide_tests(limma_fit_object, p_value=0.05)
    
    # 2. Convert to numpy for validation
    results_np = ro.numpy2ri.rpy2py(results)
    
    # 3. Validation
    assert results_np.shape == (3, 1)

