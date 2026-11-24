"""Optional independent numerical cross-check: pip install python-copasi==4.48.311."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import COPASI as copasi
import numpy as np

from chemflow.network import Network
from chemflow.projects import discover, manifest, steps


def verify_case(net, stage_path, h, c, expected):
    dm = copasi.CRootContainer.addDatamodel()
    try:
        path = stage_path / f"{h}__{c}.xml"
        if not dm.importSBML(str(path)):
            raise ValueError("Cannot import " + str(path))
        task = dm.getTask("Time-Course")
        task.setMethodType(copasi.CTaskEnum.Method_deterministic)
        method = task.getMethod()
        method.getParameter("Relative Tolerance").setDblValue(1e-9)
        method.getParameter("Absolute Tolerance").setDblValue(1e-13)
        problem = task.getProblem()
        problem.setDuration(expected[0])
        problem.setStepNumber(24)
        problem.setTimeSeriesRequested(True)
        if not task.process(True):
            raise ValueError("Simulation failed: " + str(path))
        model = dm.getModel()
        actual = {s.getSBMLId(): s.getConcentration() for s in model.getMetabolites()}
        expected = expected[1]
        np.testing.assert_allclose(
            [actual[n] for n in net.names],
            expected,
            rtol=2e-5,
            atol=2e-8,
            err_msg=str(path),
        )
    finally:
        copasi.CRootContainer.removeDatamodel(dm)


def main(root):
    import csv

    count = 0
    for project in discover():
        for stage in steps(manifest(project)):
            if stage["step"] != "simulate":
                continue
            net = Network(project)
            path = root / project / stage["id"]
            with (path / "trajectories_M.csv").open() as f:
                rows = list(csv.DictReader(f))
            endpoint = max(float(row["time_h"]) for row in rows)
            for row in rows:
                if float(row["time_h"]) == endpoint:
                    expected = [float(row[n]) for n in net.names]
                    verify_case(
                        net,
                        path,
                        row["hypothesis"],
                        row["condition"],
                        (endpoint, expected),
                    )
                    count += 1
    if not count:
        raise ValueError("No simulation cases to verify")
    print(f"{count} SBML cases agree with SciPy at the final sampled time")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
