import numpy as np
def refactor_to4chan(v):
    return [np.array(v[i]) for i in range(4)]
def data_refactor(sci:dict):
    return {k:v if len(v) != 4 else refactor_to4chan(v) for k,v in sci.items()}

def adict_info(adict):
    for k,v in adict.items():
        if isinstance(v, np.ndarray):
            print(f"{k}:{v.shape},{v.dtype}")
        if isinstance(v, list) and len(v) == 4:
            print(f"{k}:{[i.shape for i in v]},{v[0].dtype}")
        else:
            print(f"{k}:{v}")