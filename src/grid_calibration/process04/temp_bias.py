import fire

from .. import cmd
from .__init__ import tb_op

if __name__ == "__main__":
    fire.Fire(cmd.file_list_op(tb_op.to_op()))