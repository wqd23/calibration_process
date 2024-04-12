import numpy as np
def refactor_to4chan(v):
    return [np.array(v[i]) for i in range(4)]
def data_refactor(sci:dict):
    return {k:v if len(v) != 4 else refactor_to4chan(v) for k,v in sci.items()}
