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
        elif isinstance(value, list): #TODO: this logic can be better
            if all(isinstance(x, int) for x in value):
                filtered[key] = ro.IntVector(value)
            else:
                filtered[key] = ro.FloatVector(value)
        elif value is None:
            filtered[key] = ro.r('NULL')
        else:
            filtered[key] = value
            
    return filtered