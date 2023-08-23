from .__init__ import ec_op
import fire
from .. import cmd
class ec05():
    def __init__(self) -> None:
        # sub command for xray data process
        self.x = cmd.file_list_op(ec_op.to_x_op(), fp_method=None)
        # sub command for src data process
        self.src = cmd.file_list_op(ec_op.to_src_op())
if __name__ == "__main__":
    fire.Fire(ec05)