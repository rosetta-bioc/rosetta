"""Mock-based tests for stats/decide.py and stats/treat.py — no live R session."""
import pytest
from unittest.mock import MagicMock, patch


# --- decide.py ---

def test_run_decide_tests_success_mock():
    """Happy path: mocked limma.decideTests returns a sentinel."""
    import rosetta.stats.decide as decide_mod

    sentinel = object()
    mock_limma = MagicMock()
    mock_limma.decideTests.return_value = sentinel

    with patch.object(decide_mod, "importr", return_value=mock_limma):
        from rosetta.stats.decide import run_decide_tests
        result = run_decide_tests("fake_fit", method="separate", adj="BH", p_value=0.05)

    assert result is sentinel
    mock_limma.decideTests.assert_called_once_with(
        "fake_fit", method="separate", adjust_method="BH", p_value=0.05
    )


def test_run_decide_tests_r_error_raises_runtime():
    """R-side error is wrapped in RuntimeError."""
    import rosetta.stats.decide as decide_mod

    mock_limma = MagicMock()
    mock_limma.decideTests.side_effect = RuntimeError("R error")

    with patch.object(decide_mod, "importr", return_value=mock_limma):
        from rosetta.stats.decide import run_decide_tests
        with pytest.raises(RuntimeError, match="decideTests\\(\\) failed"):
            run_decide_tests("fake_fit")


def test_run_decide_tests_custom_params_passed():
    """Non-default method/adj/p_value flow through correctly."""
    import rosetta.stats.decide as decide_mod

    mock_limma = MagicMock()
    mock_limma.decideTests.return_value = "result"

    with patch.object(decide_mod, "importr", return_value=mock_limma):
        from rosetta.stats.decide import run_decide_tests
        run_decide_tests("fit", method="global", adj="bonferroni", p_value=0.01)

    mock_limma.decideTests.assert_called_once_with(
        "fit", method="global", adjust_method="bonferroni", p_value=0.01
    )


# --- treat.py ---

def test_run_treat_success_mock():
    """Happy path: mocked limma.treat returns a sentinel."""
    import rosetta.stats.treat as treat_mod

    sentinel = object()
    mock_limma = MagicMock()
    mock_limma.treat.return_value = sentinel

    with patch.object(treat_mod, "importr", return_value=mock_limma):
        from rosetta.stats.treat import run_treat
        result = run_treat("fake_fit", lfc=1.5, trend=True)

    assert result is sentinel
    mock_limma.treat.assert_called_once_with("fake_fit", lfc=1.5, trend=True)


def test_run_treat_default_params():
    """Default lfc=1.0, trend=False are passed through."""
    import rosetta.stats.treat as treat_mod

    mock_limma = MagicMock()
    mock_limma.treat.return_value = "fit"

    with patch.object(treat_mod, "importr", return_value=mock_limma):
        from rosetta.stats.treat import run_treat
        run_treat("fit")

    mock_limma.treat.assert_called_once_with("fit", lfc=1.0, trend=False)
