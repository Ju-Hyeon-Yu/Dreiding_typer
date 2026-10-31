# SMILES 결합 차수로 moltemplate dreiding.lt 라벨(원자 타입)을 할당하는 타이퍼
"""RDKit 분자(명시적 H 포함)에 moltemplate DREIDING 라벨을 붙인다.

라벨 체계: m-bone/moltemplate-DREIDING 의 DREIDING_Label_Manual.pdf
  <원소>_<혼성>[_b1|_b2][_dN][_hd|_ha]
매뉴얼 본문 규칙을 따른다: R 은 방향족 고리 원자에만 쓰고, _b1/_b2 는 짝으로만 쓴다.
결합 차수는 라벨 쌍으로 정해지므로 _b1/_b2 플래그를 제약 만족 탐색으로
골라 실제 결합 차수와 일치시키고, 모든 결합에 이면각 규칙이 있는 조합만 고른다.
판정은 dreiding.lt 의 규칙을 직접 읽어서 한다(DreidingFF).
"""
import os
import re
from collections import deque

from rdkit import Chem


class TypingError(Exception):
    pass


def default_lt_path():
    import moltemplate
    return os.path.join(os.path.dirname(moltemplate.__file__), "force_fields", "dreiding.lt")


class DreidingFF:
    """dreiding.lt 에서 라벨, 결합 규칙, 이면각 규칙을 읽어 실제 할당 결과를 조회한다.

    moltemplate 의 "By Type" 규칙은 와일드카드(*)를 쓰고, 여러 규칙이 맞으면
    파일에서 나중 규칙이 이긴다. 결합 차수는 그 규칙의 힘상수 K/350 으로 얻는다.
    """

    def __init__(self, path=None):
        self.path = path or default_lt_path()
        self.labels = {}
        bond_rules, bond_k, self.dih_rules = [], {}, []
        section = None
        with open(self.path) as f:
            for line in f:
                s = line.strip()
                m = re.match(r'write_once\("([^"]+)"\)', s)
                if m:
                    section = m.group(1)
                    continue
                if s.startswith("}"):
                    section = None
                    continue
                t = s.split()
                if section == "Data Masses" and t and t[0].startswith("@atom:"):
                    self.labels[t[0][6:]] = float(t[1])
                elif section == "Data Bonds By Type" and t and t[0].startswith("@bond:"):
                    bond_rules.append((t[0][6:], _pat(t[1]), _pat(t[2])))
                elif section == "Data Dihedrals By Type" and t and t[0].startswith("@dihedral:"):
                    self.dih_rules.append((_pat(t[2]), _pat(t[3])))
                elif t[:1] == ["bond_coeff"]:
                    bond_k[t[1][6:]] = float(t[2])
        if not self.labels:
            raise TypingError(f"{self.path} 에서 Data Masses 라벨을 찾지 못했습니다")
        self.bond_rules = [(a, b, bond_k[n]) for n, a, b in bond_rules]
        self._order, self._dih = {}, {}

    def bond_order(self, a, b):
        """라벨 쌍에 moltemplate 가 붙이는 결합의 차수 (K/350). 규칙이 없으면 None."""
        key = (a, b) if a <= b else (b, a)
        if key not in self._order:
            k = None
            for pa, pb, kk in self.bond_rules:
                if (pa.match(a) and pb.match(b)) or (pa.match(b) and pb.match(a)):
                    k = kk
            self._order[key] = None if k is None else k / 350.0
        return self._order[key]

    def has_dihedral(self, j, k):
        """J–K 를 중심으로 하는 이면각 규칙이 하나라도 있는지."""
        key = (j, k) if j <= k else (k, j)
        if key not in self._dih:
            self._dih[key] = any((pj.match(j) and pk.match(k)) or (pj.match(k) and pk.match(j))
                                 for pj, pk in self.dih_rules)
        return self._dih[key]


def _pat(token):
    """'@atom:C_R*' → 와일드카드 정규식."""
    body = token.split(":", 1)[1]
    return re.compile("^" + ".*".join(re.escape(x) for x in body.split("*")) + "$")


def sp2_class(label):
    """결합 차수 플래그(_b1/_b2)가 의미 있는 sp2/공명/sp 계열 라벨인지."""
    return any(t in label for t in ("_2", "_R", "_1"))


def _aromatic_rings(mol):
    """방향족 결합마다 그 결합을 지나는 가장 작은 방향족 고리 [(원자 집합, 결합 집합)].

    RDKit SSSR 대신 이 함수를 쓰는 것은 JS 이식(typer.js)과 결과를 똑같이 맞추기 위해서다.
    """
    aro = [b for b in mol.GetBonds() if b.GetIsAromatic()]
    nbrs = {}
    for b in aro:
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        nbrs.setdefault(i, []).append((j, b.GetIdx()))
        nbrs.setdefault(j, []).append((i, b.GetIdx()))
    rings, seen = [], set()
    for b in aro:
        s, t, skip = b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetIdx()
        prev, q = {s: None}, deque([s])
        while q and t not in prev:
            u = q.popleft()
            for v, bi in nbrs[u]:
                if bi != skip and v not in prev:
                    prev[v] = (u, bi)
                    q.append(v)
        if t not in prev:
            continue
        atoms, bonds, v = {t}, {skip}, t
        while prev[v] is not None:
            v, bi = prev[v]
            atoms.add(v)
            bonds.add(bi)
        key = frozenset(bonds)
        if key not in seen:
            seen.add(key)
            rings.append((atoms, bonds))
    return rings


def _dearomatize_exocyclic(mol):
    """고리 밖 C/N/O/B 이중결합이 붙은 방향족 고리(피리돈, 우라실 등)를 케쿨레 구조로 바꾼 복사본.

    방향족 원자(C_R)와 이중결합 상대의 라벨 쌍은 차수 2 를 낼 수 없으므로, 그 고리만 매뉴얼의
    사이클로펜타디에논처럼 비방향족 sp2 고리로 본다. 다른 방향족 고리와 공유한 원자·결합은 방향족으로 둔다.
    반환: (분자, 비방향족으로 바꾼 원자 번호 목록)
    """
    pinned = set()
    for b in mol.GetBonds():
        if b.GetIsAromatic() or b.GetBondTypeAsDouble() != 2.0:
            continue
        for a, o in ((b.GetBeginAtom(), b.GetEndAtom()), (b.GetEndAtom(), b.GetBeginAtom())):
            if a.GetIsAromatic() and o.GetSymbol() in ("C", "N", "O", "B"):
                pinned.add(a.GetIdx())
    if not pinned:
        return mol, []
    rings = _aromatic_rings(mol)
    keep_atoms, keep_bonds, drop_atoms, drop_bonds = set(), set(), set(), set()
    for atoms, bonds in rings:
        if atoms & pinned:
            drop_atoms |= atoms
            drop_bonds |= bonds
        else:
            keep_atoms |= atoms
            keep_bonds |= bonds
    drop_atoms -= keep_atoms
    drop_bonds -= keep_bonds
    kek = Chem.Mol(mol)
    try:
        Chem.Kekulize(kek, clearAromaticFlags=True)
    except Exception:
        return mol, []
    rw = Chem.RWMol(mol)
    for bi in drop_bonds:
        bond = rw.GetBondWithIdx(bi)
        bond.SetBondType(kek.GetBondWithIdx(bi).GetBondType())
        bond.SetIsAromatic(False)
    for ai in drop_atoms:
        rw.GetAtomWithIdx(ai).SetIsAromatic(False)
    return rw.GetMol(), sorted(drop_atoms)


def _is_sp2_class(atom):
    """라벨이 sp2/공명/sp 계열이 되는 원자인지: 방향족이거나 C/N/O/B 의 다중결합.

    S·P 는 S=O, P=O 가 있어도 라벨이 S_3, P_3 이므로 넣지 않는다.
    """
    if atom.GetIsAromatic():
        return True
    if atom.GetSymbol() not in ("C", "N", "O", "B"):
        return False
    return any(b.GetBondTypeAsDouble() > 1.0 for b in atom.GetBonds())


def _base_label(atom, warnings):
    """b/h 플래그를 뺀 라벨 몸통과 d 플래그를 정한다. 예: ('N_2', '_d2').

    매뉴얼은 R 을 "방향족 고리의 공명 구조에 참여하는 sp2 원자"로 정의하므로 고리 밖 공액 원자
    (아마이드·아닐린 N, 카복실레이트·니트로)는 _2 로 둔다.
    """
    sym, deg, idx = atom.GetSymbol(), atom.GetDegree(), atom.GetIdx()
    orders = [b.GetBondTypeAsDouble() for b in atom.GetBonds()]
    has2, has3 = 2.0 in orders, 3.0 in orders
    aromatic = atom.GetIsAromatic()

    if sym == "H":
        return "H", ""
    if sym == "C":
        if aromatic:
            return "C_R", ""
        return {4: "C_3", 3: "C_2", 2: "C_1"}.get(deg, None) or _fail(atom), ""
    if sym == "N":
        if aromatic:
            return "N_R", f"_d{deg - 1}"
        if deg == 4:
            return "N_3", ""
        if has3:
            return "N_1", ""
        if has2:
            return "N_2", f"_d{deg - 1}"
        if deg in (2, 3):
            # 매뉴얼 아마이드 그림처럼 π 결합 원자에 붙은 N 은 N_2, 아니면 N_3
            if any(_is_sp2_class(n) for n in atom.GetNeighbors()):
                return "N_2", f"_d{deg - 1}"
            return "N_3", ""
        _fail(atom)
    if sym == "O":
        if aromatic:
            return "O_R", ""
        if deg == 1:
            nb = atom.GetNeighbors()[0]
            if nb.GetSymbol() in ("S", "P") or has2:
                return "O_2", ""
            if has3:
                return "O_1", ""
            if nb.GetSymbol() in ("C", "N") and not nb.GetIsAromatic() and _is_sp2_class(nb):
                return "O_2", ""     # 카복실레이트·니트로의 -O⁻ → 짝 플래그로 O_2_b1 (단일결합)
            return "O_3", ""
        if deg == 3:
            warnings.append(f"원자 {idx + 1} O: 3배위 O 용 라벨이 없어 O_3 로 둡니다")
        return "O_3", ""
    if sym == "B":
        if deg == 4:
            return "B_3", ""
        if deg in (2, 3):
            return "B_2", f"_d{deg - 1}"
        _fail(atom)
    if sym in ("F", "Cl", "Br", "I"):
        if deg > 1:
            raise TypingError(f"원자 {idx + 1} ({sym}, 이웃 {deg}개): dreiding.lt 의 {sym} 라벨은 이웃 1개용뿐이라 "
                              "초원자가 할로젠(ClO₄⁻, 아이오도늄 등)에는 맞는 라벨이 없습니다")
        return sym, ""
    if sym == "S":
        if deg > 2 or has2:
            warnings.append(f"원자 {idx + 1} S: {deg}배위/다중결합 S 용 라벨이 없어 S_3 로 둡니다 "
                            "(S_3 은 이웃 2개, 평형각 92.1° 를 가정한 라벨)")
        return "S_3", ""
    if sym == "P":
        return "P_3", f"_d{deg - 1}"
    if sym == "Si":
        return "Si_3", ""
    if deg == 0:
        return sym, ""
    _fail(atom)


def _fail(atom):
    raise TypingError(f"원자 {atom.GetIdx() + 1} ({atom.GetSymbol()}, 이웃 {atom.GetDegree()}개): "
                      "DREIDING 라벨 규칙이 없습니다")


def _h_flag(atom):
    """'_hd' (수소결합 주개), '_ha' (받개), '' 중 하나."""
    sym = atom.GetSymbol()
    if sym not in ("N", "O", "F"):
        return ""
    if any(n.GetSymbol() == "H" for n in atom.GetNeighbors()):
        return "_hd"
    if sym in ("O", "F"):
        return "_ha"
    # N: 비공유전자쌍이 남는 경우만 받개
    deg, aromatic = atom.GetDegree(), atom.GetIsAromatic()
    if atom.GetFormalCharge() > 0:
        return ""
    if aromatic:
        return "_ha" if deg == 2 else ""
    orders = [b.GetBondTypeAsDouble() for b in atom.GetBonds()]
    if 3.0 in orders or 2.0 in orders or deg == 2:
        return "_ha"
    if deg == 3 and not any(_is_sp2_class(n) for n in atom.GetNeighbors()):
        return "_ha"
    return ""


def _expected_orders(mol):
    exp = {}
    for b in mol.GetBonds():
        key = tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx())))
        exp[key] = 1.5 if b.GetIsAromatic() else b.GetBondTypeAsDouble()
    return exp


def assign_labels(mol, ff, max_steps=1_000_000):
    """mol(명시적 H 포함, sanitize 완료), ff(DreidingFF) → (라벨 리스트, 경고 리스트)."""
    warnings = []
    mol, dearo = _dearomatize_exocyclic(mol)
    if dearo:
        warnings.append(f"원자 {', '.join(str(i + 1) for i in dearo)}: 고리 밖 이중결합이 붙은 방향족 고리라 "
                        "매뉴얼의 사이클로펜타디에논처럼 비방향족 고리(_2 + b 플래그)로 둡니다")

    n = mol.GetNumAtoms()
    cands = []      # 원자별 후보 라벨 (선호 순서)
    wanted_h = []
    for atom in mol.GetAtoms():
        eh, d = _base_label(atom, warnings)
        h = _h_flag(atom)
        wanted_h.append(h)
        flaggable = sp2_class(eh)
        bflags = ("", "_b1", "_b2") if flaggable else ("",)
        hs = (h, "") if h else ("",)
        c = [eh + b + d + hh for hh in hs for b in bflags if eh + b + d + hh in ff.labels]
        if not c:
            raise TypingError(f"원자 {atom.GetIdx() + 1} ({atom.GetSymbol()}): "
                              f"'{eh}{d}' 계열 라벨이 dreiding.lt 에 없습니다")
        cands.append(c)

    # H_HB: N/O/F 에 붙은 H
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "H" and atom.GetDegree() == 1 \
                and atom.GetNeighbors()[0].GetSymbol() in ("N", "O", "F"):
            cands[atom.GetIdx()] = ["H_HB"]

    exp = _expected_orders(mol)
    multi = [i for i in range(n) if len(cands[i]) > 1]
    adj = {i: [] for i in range(n)}
    for (i, j) in exp:
        adj[i].append(j)
        adj[j].append(i)

    # 탐색 순서: 후보가 여러 개인 원자를 BFS 로 (이웃끼리 연달아 정해 가지치기)
    order, seen = [], set()
    for s in multi:
        if s in seen:
            continue
        q = deque([s])
        seen.add(s)
        while q:
            u = q.popleft()
            order.append(u)
            for v in adj[u]:
                if v not in seen and len(cands[v]) > 1:
                    seen.add(v)
                    q.append(v)

    chosen = [c[0] if len(c) == 1 else None for c in cands]

    def cost(i, lab):
        # 플래그가 적을수록, b1 이 b2 보다, 원하는 h 플래그가 붙을수록 선호
        c = 1.0 if "_b1" in lab else 1.5 if "_b2" in lab else 0.0
        if wanted_h[i] and not lab.endswith(wanted_h[i]):
            c += 10.0
        return c

    deg = [a.GetDegree() for a in mol.GetAtoms()]

    def bond_ok(i, j):
        a, b = chosen[i], chosen[j]
        # 이면각 규칙이 없으면 moltemplate 가 그 이면각을 조용히 빼므로 금지
        if deg[i] > 1 and deg[j] > 1 and not ff.has_dihedral(a, b):
            return False
        if sp2_class(a) and sp2_class(b):
            o = ff.bond_order(a, b)
            return o is not None and abs(o - exp[tuple(sorted((i, j)))]) < 1e-6
        return True

    def paired(a):
        # 매뉴얼: "bX 플래그는 반드시 짝으로" — 같은 플래그를 단 이웃이 하나는 있어야 함
        lab = chosen[a]
        if lab is None:
            return True
        for fl in ("_b1", "_b2"):
            if fl in lab:
                nb = [chosen[j] for j in adj[a]]
                if all(x is not None for x in nb) and not any(fl in x for x in nb):
                    return False
        return True

    def ok(i):
        return (all(bond_ok(i, j) for j in adj[i] if chosen[j] is not None)
                and paired(i) and all(paired(j) for j in adj[i]))

    # 분기 한정 탐색: 결합 차수·이면각·짝 규칙을 모두 만족하는 조합 중 비용 최소
    best, best_cost, steps = None, float("inf"), 0

    def solve(k, acc):
        nonlocal best, best_cost, steps
        if acc >= best_cost:
            return
        if k == len(order):
            if all(paired(a) for a in range(n)):
                best, best_cost = list(chosen), acc
            return
        i = order[k]
        for lab in sorted(cands[i], key=lambda lb: cost(i, lb)):
            steps += 1
            if steps > max_steps:
                break
            chosen[i] = lab
            if ok(i):
                solve(k + 1, acc + cost(i, lab))
        chosen[i] = None

    solve(0, 0.0)
    if best is None:
        bad = [f"{i + 1}-{j + 1}" for (i, j) in exp
               if len(cands[i]) == 1 and len(cands[j]) == 1 and not bond_ok(i, j)]
        raise TypingError("결합 차수·이면각·짝 규칙을 모두 만족하는 _b1/_b2 플래그 조합이 없습니다 "
                          "(DREIDING 라벨 체계로 표현할 수 없는 결합 패턴)"
                          + (f" — 고정 라벨끼리 맞지 않는 결합: {', '.join(bad)}" if bad else ""))
    if steps > max_steps:
        warnings.append("플래그 탐색 단계 한도에 걸려 최소 비용 해가 아닐 수 있습니다(결합 차수는 만족)")
    chosen = best

    # sp3/할로젠/S 쪽이 낀 결합(S=O 등)은 라벨 쌍으로 차수를 나타낼 수 없어 경고
    for (i, j), e in exp.items():
        imp = ff.bond_order(chosen[i], chosen[j])
        if imp is None:
            raise TypingError(f"결합 {i + 1}-{j + 1} ({chosen[i]}–{chosen[j]}): dreiding.lt 에 결합 규칙이 없습니다")
        if abs(imp - e) > 1e-6:
            warnings.append(f"결합 {i + 1}-{j + 1} ({chosen[i]}–{chosen[j]}): 실제 차수 {e:g} 를 "
                            f"DREIDING 에서 {imp:g} 로 둡니다")
    # N_3 라벨은 이웃 3개를 가정함 (암모늄 4배위, TFSI 같은 2배위 N 에 맞는 라벨이 없음)
    for i, lab in enumerate(chosen):
        if lab.startswith("N_3") and deg[i] != 3:
            warnings.append(f"원자 {i + 1} N: {deg[i]}배위 N 용 라벨이 없어 {lab} 로 둡니다 "
                            "(N_3 은 이웃 3개를 가정한 라벨)")
    for i, h in enumerate(wanted_h):
        if h and not chosen[i].endswith(h):
            warnings.append(f"원자 {i + 1}: '{chosen[i]}{h}' 라벨이 없어 수소결합 플래그 {h} 를 뺐습니다")
    return chosen, warnings
