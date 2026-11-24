"""Generate organic starting coordinates from project-defined SMILES."""

from .projects import folder, identifier, read_json


def smiles_molecule(record):
    try:
        from rdkit import Chem
    except ImportError as exc:
        raise ValueError(
            "SMILES generation requires RDKit (conda environment or pip install rdkit)"
        ) from exc
    mol = Chem.MolFromSmiles(record["smiles"])
    if mol is None:
        raise ValueError(f"Invalid SMILES: {record['smiles']}")
    if Chem.GetFormalCharge(mol) != record["charge"]:
        raise ValueError("SMILES formal charge differs from the registered charge")
    return Chem.AddHs(mol)


def build(project, dest, names=None, seed=137):
    from .qm import ORGANIC_ELEMENTS

    records = read_json(folder(project) / "species.json")
    selected = names or [n for n, r in records.items() if r.get("smiles")]
    if not selected:
        raise ValueError(
            "No species with SMILES declared; supply reviewed XYZs or add smiles records"
        )
    if isinstance(selected, str) or len(set(selected)) != len(selected):
        raise ValueError("Specify a unique list of species IDs")
    dest.mkdir(parents=True, exist_ok=True)
    for name in selected:
        identifier(name)
        rec = records[name]
        if not rec.get("smiles"):
            raise ValueError(f"{name}: missing smiles")
        mol = smiles_molecule(rec)
        from rdkit import Chem
        from rdkit.Chem import AllChem

        if (
            any(a.GetSymbol() not in ORGANIC_ELEMENTS for a in mol.GetAtoms())
            or len(Chem.GetMolFrags(mol)) != 1
        ):
            raise ValueError(
                f"{name}: provide a reviewed XYZ for coordination complexes or separate fragments"
            )
        par = AllChem.ETKDGv3()
        par.randomSeed = seed
        par.numThreads = 1
        if AllChem.EmbedMolecule(mol, par) != 0 or not AllChem.UFFHasAllMoleculeParams(
            mol
        ):
            raise ValueError(f"Failed to create supported organic geometry: {name}")
        if AllChem.UFFOptimizeMolecule(mol, maxIters=1500) != 0:
            raise ValueError(f"UFF failed to converge: {name}")
        mol.SetProp("_Name", name + "; RDKit ETKDGv3/UFF; " + rec["smiles"])
        (dest / (name + ".xyz")).write_text(Chem.MolToXYZBlock(mol), encoding="utf-8")
