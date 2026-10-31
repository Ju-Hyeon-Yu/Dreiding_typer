# Dreiding_typer

**English** | [한국어](README_kor.md)

Dreiding_typer builds a single-molecule moltemplate `.lt` file with DREIDING atom types from a SMILES string,
optionally together with a mol2 file that carries RESP charges. Atom types follow the labeling scheme of the
`force_fields/dreiding.lt` file shipped with moltemplate (Matthew Bone's version). The output has the form
`import "dreiding.lt"` + `MOL inherits DREIDING { ... }`, so it can be pulled straight into a `system.lt`.

- Types are derived from the bond orders in the SMILES, so they do not depend on the conformer.
- Hydrogen-bond labels (`H_HB`, `_hd`, `_ha`) are assigned automatically.
- Bond-order flags (`_b1`/`_b2`) are chosen by reading the rules in `dreiding.lt` directly, so that the bond
  order moltemplate assigns matches the real one and every bond has a dihedral rule.

Program messages and the HTML interface are currently in Korean.

## Installation

Installing into a miniconda environment makes the `dreiding-lt` command available only while that environment
is active.

```bash
conda activate myenv
git clone https://github.com/Ju-Hyeon-Yu/Dreiding_typer.git ~/Dreiding_typer
cd ~/Dreiding_typer
pip install -e .                 # also installs rdkit and moltemplate if missing
```

- `-e` (editable install) means a later `git pull` takes effect without reinstalling.
- The command is created only in that environment's `bin`, so it is not visible from other environments or base.
- If RDKit is already installed with conda (`conda install -c conda-forge rdkit`), pip leaves it alone.
- To remove it, run `pip uninstall dreiding-typer` in that environment.
- From Python code you can use `from dreiding_typer import DreidingFF, assign_labels`.

Python 3.9 or later is required. The label list and rules are read from the installed moltemplate's
`dreiding.lt` on every run.

## Usage

```bash
# With a mol2 that has RESP charges — atom order, coordinates and charges come from the mol2,
# bond orders come from the SMILES
dreiding-lt --smiles "Cn1cc[n+](C)c1.[O-]S(=O)(=O)C(F)(F)F" \
    --mol2 ip.mol2 -o single_molecule.lt --name MOL

# SMILES only — RDKit 3D structure, all charges 0 (a warning is printed)
dreiding-lt --smiles "CCO" -o ethanol.lt --name MOL

# Several molecules with a shell loop (list.csv: name,SMILES)
while IFS=, read name smi; do
  dreiding-lt --smiles "$smi" --mol2 "$name.mol2" -o "$name.lt" --name "$name"
done < list.csv
```

Quote the SMILES on the command line. Characters such as `(`, `[` and `#` are interpreted by the shell
before the program sees them.

| Option | Description |
|---|---|
| `--smiles` | Required. Join the ions of an ion pair with `.` (e.g. `A.B`) |
| `--mol2` | mol2 to take coordinates and charges from. Without it, RDKit builds a 3D structure and charges are 0 |
| `-o`, `--out` | Output `.lt` path |
| `--name` | moltemplate object name. Defaults to the output file name |
| `--dreiding-lt` | `dreiding.lt` to read labels and rules from. Defaults to the installed moltemplate's file |
| `--seed` | Random seed for 3D structure generation when there is no mol2 |

The per-atom labels and charges are printed to the screen. Anything DREIDING cannot represent exactly is
reported on standard error as a warning (`경고:`). If a label cannot be determined, no `.lt` is written and the
program stops.

### Example output

```
# dreiding_lt.py 생성 — SMILES: Cn1cc[n+](C)c1.[O-]S(=O)(=O)C(F)(F)F
# 전하: mol2 (ip.mol2)

import "dreiding.lt"

MOL inherits DREIDING {

  write("Data Atoms") {
    $atom:C1     $mol:f1 @atom:C_3              0.079135     3.44930     1.07820     0.99060
    $atom:N2     $mol:f1 @atom:N_R_d2          -0.239526     2.56160    -0.01380     0.61200
    ...
  }

  write("Data Bond List") {
    $bond:b1    $atom:C1     $atom:N2
    ...
  }

}
```

Each fragment of an ion pair gets its own molecule ID (`$mol:f1`, `$mol:f2`).

### Loading it from system.lt

```
import "single_molecule.lt"
m = new MOL
```

Then run `moltemplate.sh system.lt` to generate the LAMMPS input. The default settings in `dreiding.lt` are
`hybrid/overlay hbond/dreiding/lj ... lj/cut/coul/long` + `kspace_style ewald`, and LAMMPS needs the MOLECULE
package (`hbond/dreiding/lj`, `dihedral charmm`, `improper umbrella`).

## Offline HTML tool (`dreiding_typer.html`)

A single HTML file that does the same job in a browser, without Python. Download `dreiding_typer.html` from the
repository and open it in Edge or Chrome; it works with no internet connection (the RDKit.js engine and the
DREIDING rule tables are embedded, about 3.6 MB).

1. Enter a SMILES (required).
2. Optionally drop a coordinate file (mol2, pdb, xyz). Atom order is matched automatically.
   - mol2: coordinates and RESP charges are used.
   - pdb, xyz: only coordinates are used and charges are set to 0. If a pdb has no CONECT records, or for xyz,
     bonds are determined from interatomic distances.
   - Without a coordinate file, a 3D starting structure is built from 2D coordinates by minimizing DREIDING bond
     and angle terms, and charges are set to 0.
3. Press **타입 할당** (assign types).
   - The 3D viewer labels every atom with its type. Drag to rotate, scroll to zoom, Shift+drag to pan,
     double-click to reset. Labels can show main atoms (heavy atoms and `H_HB`), all atoms or none; atoms can be
     colored by type or by element.
   - The **원자 타입** (atom types) cards show the count and atom numbers for each type. Clicking a card
     highlights only the atoms of that type.
   - Hovering over an atom shows its number, name, fragment and charge. Atoms with warnings have an orange ring.
   - The `.lt` text appears below, with **복사** (copy) and **내려받기** (download, `single_molecule.lt` by
     default). The viewer image can be saved as PNG.

The typing code is a direct port of the Python typer and gives the same results. Tests confirm that labels,
warnings and errors agree for 128 SMILES, and that the full `.lt` agrees with the Python CLI for 12 systems given
as mol2 files with shuffled atom order. A 3D structure built without a coordinate file is only a starting point:
minimize it before running MD. Stereochemistry in the SMILES is not used when building 3D coordinates.

### Rebuilding the HTML

Only needed after changing code in `web/src` or when moltemplate's `dreiding.lt` changes.

```bash
cd web && npm install && cd ..        # RDKit.js (@rdkit/rdkit 2026.9.1)
python3 web/build.py                  # → dreiding_typer.html
cd web && npm test                    # compare against the Python typer (Node)
node tests/test_browser.js            # check the built HTML in a browser (needs playwright)
```

## When you only have a mol2 or PDB file, without SMILES

`--smiles` is required. DREIDING labels (`C_R`/`C_2`, `_b1`, ...) are determined by bond orders, and bond orders
and formal charges are taken only from the SMILES (the bond-order column of a mol2 is not used).

Without a SMILES, first infer bond orders from the coordinates to get a SMILES, then pass it in. The script below
uses RDKit `DetermineBonds`; it ignores any bond information in the file and uses only coordinates.

```python
# guess_smiles.py — usage: python3 guess_smiles.py ip.mol2 0   (last argument: total charge of the system)
import sys
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

path, charge = sys.argv[1], int(sys.argv[2])
if path.endswith(".pdb"):
    mol = Chem.MolFromPDBFile(path, removeHs=False, sanitize=False, proximityBonding=False)
else:
    mol = Chem.MolFromMol2File(path, removeHs=False, sanitize=False, cleanupSubstructures=False)
rw = Chem.RWMol(mol)
for b in list(rw.GetBonds()):
    rw.RemoveBond(b.GetBeginAtomIdx(), b.GetEndAtomIdx())
for a in rw.GetAtoms():
    a.SetFormalCharge(0)
    a.SetNoImplicit(True)
mol = rw.GetMol()
rdDetermineBonds.DetermineBonds(mol, charge=charge)
print(Chem.MolToSmiles(Chem.RemoveHs(mol)))
```

```bash
smi=$(python3 guess_smiles.py ip.mol2 0)
echo "$smi"                                   # always check it by eye
dreiding-lt --smiles "$smi" --mol2 ip.mol2 -o single_molecule.lt --name MOL
```

Points to watch:

- The structure must contain **all hydrogens**. The total charge must be given by hand (usually 0 for an ion pair).
- **Ions with several resonance structures may come out in a different form.** For example, TFSI⁻ is inferred
  with the negative charge on O instead of N and an S=N double bond, so the N label becomes `N_2_d1` instead of
  `N_R_d1`. If the printed SMILES is not the structure you intended, correct it by hand.
- **Coordinates and charges are read only from mol2.** With only a PDB, the script above still gives you a SMILES,
  but running without `--mol2` makes RDKit build new coordinates and sets charges to 0. To use PDB coordinates
  and RESP charges, convert the PDB to mol2 first.

## mol2 requirements

- Elements are read from the atom-type column (SYBYL `C.ar`, GAFF `ca`, ...), falling back to the leading letters
  of the atom name.
- Charges are read from the 9th column. If it is missing, the charge is set to 0 with a warning.
- **Atom order does not need to match.** Atoms of the SMILES and the mol2 are matched by element and connectivity
  (graph isomorphism), and the output `.lt` keeps the mol2 atom order. Only the connectivity in the BOND section
  is used, not its bond-order column. If the atom count, bond count or connectivity differ, the program stops.
- A warning is printed if the sum of charges differs from the formal charge of the SMILES.

## Typing rules in brief

| Case | Label |
|---|---|
| Aromatic C/N/O | `C_R`, `N_R_d1` (2 neighbors), `N_R_d2` (3 neighbors), `O_R` |
| Carboxylate and nitro resonance groups | Center `C_R`/`N_R`, O as `O_R` (bond order 1.5) |
| Three-coordinate N bonded to an sp2/aromatic atom (amide, aniline) | `N_R_d2` |
| Two-coordinate N without multiple bonds (imide anion) | `N_R_d1` |
| O in esters, phenols and acid OH | `O_3` (dihedral rule (i) of the original paper applies) |
| S | `S_3` (the only S type) |
| I, Br, Cl, F | `I`, `Br`, `Cl`, `F` |
| H bonded to N/O/F | `H_HB`, and its heavy atom gets `_hd` |
| Hydrogen-bond acceptors | `_ha` on O, F, and N with a free lone pair |
| sp2 atom pairs joined by a single bond | `_b1`/`_b2` assigned by search |

### Relation to the examples in the label manual

Only labels defined in moltemplate's `dreiding.lt` are used. Some labels differ in name from the figures in the
label manual ([DREIDING_Label_Manual.pdf](https://github.com/m-bone/moltemplate-DREIDING)), but tests confirm that
moltemplate assigns identical bond, angle, dihedral and improper coefficients in those cases
(cyclopentadienone, imidazole N–H).

The one deliberate exception is amides. The manual's example (`C_2` + `N_2_b1_d2_hd`) treats C–N as a double bond
(equilibrium 1.275 Å, rotation barrier 45 kcal/mol), while this tool (`C_2_b1` + `N_R_b1_d2_hd`) treats it as a
single bond (1.31 Å, barrier 5 kcal/mol). The latter is closer to a real amide (about 1.33 Å, 15–20 kcal/mol).

## Known limitations (reported as warnings)

- **Three- and four-coordinate S (sulfonate, sulfone, TFSI)** — DREIDING has no such type, so `S_3`
  (based on H₂S, equilibrium angle 92.1°) is used. S=O gets a single-bond force constant, and the dihedral barrier
  through S is divided assuming two neighbors, so it is three times too large.
- **Four-coordinate N (quaternary ammonium, protonated amines)** — `N_3` assumes three neighbors, so the dihedral
  barrier through this N is 1.5 times too large.
- **Cl** — the Cl mass in `dreiding.lt` is 34.45 (the real value is 35.45).

## Files

| File | Contents |
|---|---|
| `dreiding_lt.py` | CLI (`dreiding-lt` command) — SMILES/mol2 input, atom matching, `.lt` writing |
| `dreiding_typer.py` | Typer — label rules, `dreiding.lt` rule parser, `_b1`/`_b2` search |
| `pyproject.toml` | pip install settings (registers the `dreiding-lt` command) |
| `tests/` | Unit and integration tests |
| `dreiding_typer.html` | Offline single-file HTML tool (build output) |
| `web/src/` | HTML tool sources — `typer.js` (typer), `io.js` (files, atom matching, .lt), `geom.js` (3D build), `pipeline.js`, `viewer.js`, `app.js`, `index.html` |
| `web/build.py`, `web/build_tables.py` | Rule-table extraction and single-file HTML build |
| `web/tests/` | Comparison tests against Python (Node) and browser tests |
| `plan.md`, `checklist.md`, `context-notes.md` | Design plan, task list, and a log of decisions and their reasons (Korean) |

## Tests

```bash
pip install -e ".[test]"         # pytest, lammps (LAMMPS tests are skipped if lammps is missing)
python3 -m pytest -q tests       # about a minute, including moltemplate runs
```

The integration tests run moltemplate on the `.lt` files of 14 molecules, including ion pairs, TFSI/BMIM,
amides, carboxylates, nitro groups, biphenyl and conjugated chains, and check that:

- every bond force constant equals 350 × the expected bond order
- every angle and dihedral receives a type (nothing is dropped)
- dihedral coefficients match an independent implementation of rules (a)–(j) of Mayo 1990
  (except for the limitations above)
- each fragment of an ion pair has its own molecule ID, and mol2 atom order and charges are preserved
- the label-manual examples and this tool's labels give identical moltemplate coefficients
- LAMMPS can read the result and run a single-point calculation

## References

- Mayo, Olafson, Goddard, *J. Phys. Chem.* **94**, 8897 (1990) — DREIDING
- [jewettaij/moltemplate](https://github.com/jewettaij/moltemplate) — `force_fields/dreiding.lt`
- [m-bone/moltemplate-DREIDING](https://github.com/m-bone/moltemplate-DREIDING) — label manual
- [RDKit.js](https://github.com/rdkit/rdkit-js) — chemistry engine of the HTML tool (BSD-3-Clause)
