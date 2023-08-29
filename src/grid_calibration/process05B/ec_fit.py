from .. import util_lib as util
from .__init__ import ec_op

x_res, src_res = util.get_fit_dict(ec_op.save_path, ec_op.energy)

src_result = list(src_res.values())
src_energy = list(src_res.keys())
x_result = list(x_res.values())
x_energy = list(x_res.keys())
ec_op.ec_fit(src_result, src_energy, x_result, x_energy)