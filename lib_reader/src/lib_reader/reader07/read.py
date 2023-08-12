from . import gridBasicFunctions02 as basic
import numpy as np
def single_read07(path:str):
    sciExtracted, telExtracted = basic.dataReadout(path, isHex=True, isBin=False, newProgramme=True)
    for k,v in sciExtracted.items():
        sciExtracted[k] = [np.array(v[i]) for i in range(4)]
    for k,v in telExtracted.items():
        telExtracted[k] = [np.array(v[i]) for i in range(4)]
    
    return sciExtracted, telExtracted