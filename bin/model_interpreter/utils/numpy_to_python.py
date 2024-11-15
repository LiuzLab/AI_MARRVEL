import numpy as np

def convert_numpy_to_python(obj):
    """
    Convert Numpy to Python object for JSON storage

    Usage: for json.dumps
    json_string = json.dumps(json_data, indent=2, default=convert_numpy_to_python)

    """
    
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, np.float32):
        return float(obj)
    return obj
