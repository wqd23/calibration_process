from ..operation import EC_operation_03B
import os
from .. import util_lib as util
from .__init__ import ec_op
energy = 'data/05B/single_process/ec_energy.json'

dir = 'data/05B/single_process/EC_fit_result'
x_files = [f'{os.path.splitext(x_energy)[0]}.pickle' for x_energy in ec_op.x_list]
src_files = [f'{os.path.splitext(file)[0]}.pickle' for file in ec_op.src_list]
ec_energy = util.json_load(energy)
def get_fit_result(file):
    path = os.path.join(dir, file)
    data = util.pickle_load(path)
    key = f"{os.path.splitext(os.path.basename(data['file']))[0]}"
    if "src" in key:
        key = f"{key}.dat"
    else:
        key = key + ".dat"
    energy = ec_energy[key]
    return energy, data['fit_result']

src_result = []
src_energy = []
for src_file in src_files:
    energy, fit_result = get_fit_result(src_file)
    src_result.append(fit_result)
    src_energy.append(energy)
x_result = []
x_energy = []
for x_file in x_files:
    energy, fit_result = get_fit_result(x_file)
    x_result.append(fit_result)
    x_energy.append(energy)
ec_op.ec_fit(src_result, src_energy, x_result, x_energy)