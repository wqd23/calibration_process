from . import gridBasicFunctions02 as basic
import numpy as np
from ..reader05.version_lib import data_refactor
def single_read04(path:str):
    sciExtracted, telExtracted = basic.dataReadout(path, isHex=True, isBin=False, newProgramme=True)
    sciExtracted, telExtracted = data_refactor(sciExtracted), data_refactor(telExtracted)
    
    return sciExtracted, telExtracted