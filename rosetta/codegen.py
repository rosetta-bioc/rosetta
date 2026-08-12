"""rosetta.codegen — Show the equivalent R code being run behind the scenes.

Usage:
    import rosetta as rb
    rb.codegen.enable()
    results = rb.deseq2(counts, meta, design="~ condition")
    # prints equivalent R code as it runs

    # Or get it as a string:
    code = rb.codegen.last()
"""

import threading

_state = threading.local()


def _get_log() -> list[str]:
    if not hasattr(_state, "log"):
        _state.log = []
    return _state.log


def _is_enabled() -> bool:
    return getattr(_state, "enabled", False)


def enable():
    """Enable R code printing for all rosetta calls."""
    _state.enabled = True
    _get_log().clear()


def disable():
    """Disable R code printing."""
    _state.enabled = False


def last() -> str:
    """Return the last recorded R code block as a string."""
    return "\n".join(_get_log())


def clear():
    """Clear the code log."""
    _get_log().clear()


def _emit(line: str):
    if not _is_enabled():
        return
    _get_log().append(line)
    print(f"  \033[2mR>\033[0m \033[36m{line}\033[0m")


def _block(lines: list[str]):
    """Record and optionally print a block of R code."""
    for line in lines:
        _emit(line)
