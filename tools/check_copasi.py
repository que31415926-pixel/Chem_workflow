"""Optional independent numerical cross-check: pip install python-copasi==4.48.311."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import COPASI as copasi
import numpy as np

from chemflow.network import PROJECTS, Network


def main(root):
    count = 0
    for project in PROJECTS:
        net = Network(project)
        for h in net.hypotheses:
            for c in net.conditions:
                dm = copasi.CRootContainer.addDatamodel()
                try:
                    path = root / project / "scenarios" / f"{h}__{c}.xml"
                    if not dm.importSBML(str(path)):
                        raise ValueError("Cannot import " + str(path))
                    task = dm.getTask("Time-Course")
                    task.setMethodType(copasi.CTaskEnum.Method_deterministic)
                    method = task.getMethod()
                    method.getParameter("Relative Tolerance").setDblValue(1e-9)
                    method.getParameter("Absolute Tolerance").setDblValue(1e-13)
                    problem = task.getProblem()
                    problem.setDuration(net.data["end_h"])
                    problem.setStepNumber(24)
                    problem.setTimeSeriesRequested(True)
                    if not task.process(True):
                        raise ValueError("Simulation failed: " + str(path))
                    model = dm.getModel()
                    actual = {
                        s.getSBMLId(): s.getConcentration()
                        for s in model.getMetabolites()
                    }
                    expected = net.simulate(h, c, [net.data["end_h"]])[:, 0]
                    np.testing.assert_allclose(
                        [actual[n] for n in net.names],
                        expected,
                        rtol=2e-5,
                        atol=2e-8,
                        err_msg=str(path),
                    )
                    count += 1
                finally:
                    copasi.CRootContainer.removeDatamodel(dm)
    print(f"{count} SBML cases agree with SciPy at project endpoints")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
