def test_run_treat_success(limma_fit_object):
    from rosetta.stats.treat import run_treat
    treated_fit = run_treat(limma_fit_object, lfc=1.0)
    assert treated_fit is not None