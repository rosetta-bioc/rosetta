import logging
# Safely try importing rpy2 for keyword filtering if available
try:
    import rpy2.robjects as ro
except Exception:
    ro = None

def filter_kwargs(kwargs, allowed_args):
    """
    Filters kwargs based on allowed_args and converts Python types to R types.
    """
    filtered = {}
    for key, value in kwargs.items():
        if key not in allowed_args:
            logging.warning(f"Parameter '{key}' is not supported and will be ignored.")
            continue
            
        if isinstance(value, bool):
            filtered[key] = value # rpy2 handles bool automatically
        elif isinstance(value, list):
            if all(isinstance(x, bool) for x in value):
                filtered[key] = ro.BoolVector(value)
            elif all(isinstance(x, int) for x in value):
                filtered[key] = ro.IntVector(value)
            elif all(isinstance(x, float) for x in value):
                filtered[key] = ro.FloatVector(value)
            elif all(isinstance(x, str) for x in value):
                filtered[key] = ro.StrVector(value)
            else:
                filtered[key] = ro.StrVector([str(x) for x in value])
        elif value is None:
            if ro is not None:
                filtered[key] = ro.r('NULL')
            else:
                filtered[key] = None
        else:
            filtered[key] = value
            
    return filtered