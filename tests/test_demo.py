from grid_calibration.process import (
    process03B,
    process04,
    process05B,
    process07,
    process10B,
)
from grid_calibration.cmd import VersionProcessOp
import cProfile


def __test_all(p: VersionProcessOp):
    p.tb.list()
    p.tb.run(0)
    p.tb.run(0, 1)
    p.tbfit()
    p.ec["x"].list()
    p.ec["x"].run(0)
    p.ec["x"].run(0, 1)
    p.ec["src"].list()
    p.ec["src"].run(0)
    p.ec["src"].run(0, 1)
    p.ecfit()


def test_process03B():
    __test_all(process03B)


def test_process04():
    __test_all(process04)


def test_process05B():
    __test_all(process05B)


def test_process07():
    __test_all(process07)
if __name__ == "__main__":
    with cProfile.Profile() as pr:
        process03B.tb.run(0)
    pr.print_stats(sort='cumulative')
    pr.dump_stats("cali.prof")