"""Rebuild independent geometry guesses from connectivity (RDKit ETKDG + UFF).
No fitted or literature Cartesian coordinates are distributed.
"""

from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).resolve().parents[1]
ORG = {
    "co": {
        "veratrole": "COc1ccccc1OC",
        "hbpin": "[BH]1OC(C)(C)C(C)(C)O1",
        "hydrogen": "[H][H]",
        "aryl_bpin": "COc1ccc(B2OC(C)(C)C(C)(C)O2)cc1OC",
    },
    "fe": {
        "methyl_rad": "[CH3]",
        "ethyl_rad": "C[CH2]",
        "neopentyl_rad": "CC(C)(C)[CH2]",
        "pivalic_acid": "CC(C)(C)C(=O)O",
        "tert_butyl_rad": "C[C](C)C",
        "benzyl_rad": "[CH2]c1ccccc1",
        "sec_benzyl_rad": "C[CH]c1ccccc1",
        "triphenylpyridine": "n1c(-c2ccccc2)cc(-c2ccccc2)cc1-c1ccccc1",
        "product": "CC(c1ccccc1)C(C)(C)C",
        "salt_primary": "c1ccc(C[n+]2c(-c3ccccc3)cc(-c3ccccc3)cc2-c2ccccc2)cc1",
        "salt_secondary": "CC(c1ccccc1)[n+]1c(-c2ccccc2)cc(-c2ccccc2)cc1-c1ccccc1",
        "salt_mod_primary": "c1ccc(C[n+]2c(C(C)C)cc(-c3ccccc3)cc2-c2ccccc2)cc1",
        "salt_mod_secondary": "CC(c1ccccc1)[n+]1c(C(C)C)cc(-c2ccccc2)cc1-c1ccccc1",
    },
}
OEP = "CCC1=C(C2=CC3=NC(=CC4=NC(=CC5=C(C(=C(N5)C=C1N2)CC)CC)C(=C4CC)CC)C(=C3CC)CC)CC"
PINCER = "Cc1cc(C)c(-n2cc[n+](C(C)C)c2)nc1-n1cc[n+](C(C)C)c1"


def unit(v):
    return v / np.linalg.norm(v)


def prepare(smiles, pincer=False):
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    par = AllChem.ETKDGv3()
    par.randomSeed = 61453
    par.numThreads = 1
    if AllChem.EmbedMolecule(m, par) < 0:
        raise ValueError(smiles)
    ff = AllChem.UFFGetMoleculeForceField(m)
    if pincer:
        nit = [a for a in m.GetAtoms() if a.GetSymbol() == "N" and a.IsInRingSize(6)][
            0
        ].GetIdx()
        donors = [
            a.GetIdx()
            for a in m.GetAtoms()
            if a.GetSymbol() == "C"
            and a.IsInRingSize(5)
            and sum(b.GetSymbol() == "N" for b in a.GetNeighbors()) == 2
        ]
        for i in donors:
            ff.AddDistanceConstraint(nit, i, 2.55, 2.8, 400.0)
        ff.AddDistanceConstraint(*donors, 3.65, 3.85, 400.0)
    ff.Minimize(maxIts=3500)
    return m, m.GetConformer().GetPositions()


def xyz(project, name, syms, coords, desc):
    path = ROOT / "projects" / project / "xyz" / (name + ".xyz")
    path.parent.mkdir(exist_ok=True)
    with path.open("w") as f:
        print(len(syms), desc, sep="\n", file=f)
        for atom, point in zip(syms, coords, strict=True):
            print(atom, *(f"{x:.7f}" for x in point), file=f)


def symbols(m):
    return [a.GetSymbol() for a in m.GetAtoms()]


def metal_seeds():
    m, p = prepare(PINCER, True)
    atoms = list(m.GetAtoms())
    donors = [
        a
        for a in atoms
        if a.GetSymbol() == "C"
        and a.IsInRingSize(5)
        and sum(b.GetSymbol() == "N" for b in a.GetNeighbors()) == 2
    ]
    nh = [a.GetIdx() for a in atoms if a.GetSymbol() == "N" and a.IsInRingSize(6)][0]
    stripped = {
        n.GetIdx() for a in donors for n in a.GetNeighbors() if n.GetSymbol() == "H"
    }
    keep = [i for i in range(len(atoms)) if i not in stripped]
    axis = unit(np.mean([p[a.GetIdx()] for a in donors], axis=0) - p[nh])
    co = p[nh] + 1.92 * axis

    def co_xyz(label, tail, positions):
        xyz(
            "co",
            label,
            [atoms[i].GetSymbol() for i in keep] + ["Co"] + tail,
            np.vstack([p[keep], co, *positions]),
            "Independent idealized 3,5-Me2 / N-iPr ACNC seed; inspect donor geometry",
        )

    co_xyz("co_h", ["H"], [co + 1.55 * axis])
    c = co + 2.00 * axis
    v = unit(np.cross(axis, [0.0, 0.0, 1.0]))
    u = np.cross(axis, v)
    ch = [
        c + 1.09 * (axis / 3 + (8 / 9) ** 0.5 * (v * np.cos(angle) + u * np.sin(angle)))
        for angle in (0, 2 * np.pi / 3, 4 * np.pi / 3)
    ]
    co_xyz("co_me", ["C", "H", "H", "H"], [c, *ch])

    bh, q = prepare(ORG["co"]["hbpin"])
    boron = next(a for a in bh.GetAtoms() if a.GetSymbol() == "B")
    h = next(n.GetIdx() for n in boron.GetNeighbors() if n.GetSymbol() == "H")
    b = boron.GetIdx()
    rot, _ = Rotation.align_vectors((-axis)[None, :], unit(q[h] - q[b])[None, :])
    bpin = rot.apply(q - q[b]) + co + 2.03 * axis
    retain = [i for i in range(bh.GetNumAtoms()) if i != h]
    co_xyz(
        "co_bpin",
        [bh.GetAtomWithIdx(i).GetSymbol() for i in retain],
        list(bpin[retain]),
    )

    m, p = prepare(OEP)
    ns = [a for a in m.GetAtoms() if a.GetSymbol() == "N"]
    stripped = {
        n.GetIdx() for a in ns for n in a.GetNeighbors() if n.GetSymbol() == "H"
    }
    assert len(stripped) == 2 and len(ns) == 4
    keep = [i for i in range(m.GetNumAtoms()) if i not in stripped]
    pts = p[[a.GetIdx() for a in ns]]
    c = pts.mean(axis=0)
    axis = np.linalg.svd(pts - c)[2][-1]
    core = [m.GetAtomWithIdx(i).GetSymbol() for i in keep] + ["Fe"]

    def fe_xyz(label, offset, tail, positions):
        xyz(
            "fe",
            label,
            core + tail,
            np.vstack([p[keep], c + offset * axis, *positions]),
            "Idealized Fe-OEP seed from free-base OEP connectivity (PubChem CID 102311)",
        )

    fe_xyz("fe_oep", 0, [], [])
    fe_xyz("fe_oep_cl", 0.35, ["Cl"], [c + 2.60 * axis])
    for label, rad in [
        ("fe_methyl", "methyl_rad"),
        ("fe_ethyl", "ethyl_rad"),
        ("fe_neopentyl", "neopentyl_rad"),
        ("fe_bn", "benzyl_rad"),
        ("fe_sec_bn", "sec_benzyl_rad"),
    ]:
        mol, q = prepare(ORG["fe"][rad])
        j = next(a.GetIdx() for a in mol.GetAtoms() if a.GetNumRadicalElectrons())
        neighbors = [a.GetIdx() for a in mol.GetAtomWithIdx(j).GetNeighbors()]
        # The new Fe-C bond is initially perpendicular to the radical's local plane.
        normal = np.linalg.svd(q[neighbors] - q[neighbors].mean(axis=0))[2][-1]
        rot, _ = Rotation.align_vectors(axis[None, :], normal[None, :])
        other = rot.apply(q - q[j]) + c + 2.26 * axis
        fe_xyz(label, 0.22, symbols(mol), list(other))


def main():
    for proj, molecules in ORG.items():
        for name, smi in molecules.items():
            m, p = prepare(smi)
            xyz(
                proj,
                name,
                symbols(m),
                p,
                "RDKit ETKDGv3/UFF coordinate guess; SMILES " + smi,
            )
    metal_seeds()
    # Separated fragments give atom order/composition for a manually refined encounter or TS.
    for proj, name, first, second in [
        ("co", "co_bpin_ar", "co_bpin", "veratrole"),
        ("fe", "sh2_sec", "fe_sec_bn", "tert_butyl_rad"),
        ("fe", "sh2_pri", "fe_bn", "tert_butyl_rad"),
    ]:
        rows = []
        for index, part in enumerate([first, second]):
            lines = (
                (ROOT / "projects" / proj / "xyz" / (part + ".xyz"))
                .read_text(encoding="utf-8")
                .splitlines()[2:]
            )
            els = [lin.split()[0] for lin in lines]
            pos = np.array([lin.split()[1:] for lin in lines], dtype=float)
            if index:
                pos[:, 0] += max(q[1][0] for q in rows) - pos[:, 0].min() + 4.0
            rows.extend(zip(els, pos))
        xyz(
            proj,
            name,
            *zip(*rows),
            "Separated fragments only; orient and review before TS search",
        )


if __name__ == "__main__":
    main()
