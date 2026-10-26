# SMILES 결합 차수로 moltemplate dreiding.lt 라벨(원자 타입)을 할당하는 타이퍼
"""RDKit 분자(명시적 H 포함)에 moltemplate DREIDING 라벨을 붙인다.

라벨 체계: m-bone/moltemplate-DREIDING 의 DREIDING_Label_Manual.pdf
  <원소>_<혼성>[_b1|_b2][_dN][_hd|_ha]
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


def _resonance_groups(mol):
    """C/N 중심에 붙은 말단 O 가 =O 와 -O⁻ 로 섞여 있으면 공명 묶음(카복실레이트, 니트로)."""
    groups = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in ("C", "N") or atom.GetIsAromatic():
            continue
        term = [(b.GetOtherAtom(atom), b.GetBondTypeAsDouble()) for b in atom.GetBonds()
                if b.GetOtherAtom(atom).GetSymbol() == "O" and b.GetOtherAtom(atom).GetDegree() == 1]
        has_double = any(od == 2.0 for _, od in term)
        has_anion = any(od == 1.0 and o.GetFormalCharge() == -1 for o, od in term)
        if has_double and has_anion:
            groups.append((atom.GetIdx(), [o.GetIdx() for o, _ in term]))
    return groups


def _is_sp2_class(atom, res_atoms):
    """방향족이거나 다중결합을 가진(또는 공명 묶음의) 원자인지."""
    if atom.GetIsAromatic() or atom.GetIdx() in res_atoms:
        return True
    return any(b.GetBondTypeAsDouble() > 1.0 for b in atom.GetBonds())


def _base_label(atom, res_atoms, warnings):
    """b/h 플래그를 뺀 라벨 몸통과 d 플래그를 정한다. 예: ('N_R', '_d2')."""
    sym, deg, idx = atom.GetSymbol(), atom.GetDegree(), atom.GetIdx()
    orders = [b.GetBondTypeAsDouble() for b in atom.GetBonds()]
    has2, has3 = 2.0 in orders, 3.0 in orders
    resonant = atom.GetIsAromatic() or idx in res_atoms

    if sym == "H":
        return "H", ""
    if sym == "C":
        if resonant:
            return "C_R", ""
        return {4: "C_3", 3: "C_2", 2: "C_1"}.get(deg, None) or _fail(atom), ""
    if sym == "N":
        if resonant:
            return "N_R", f"_d{deg - 1}"
        if deg == 4:
            return "N_3", ""
        if has3:
            return "N_1", ""
        if has2:
            return "N_2", f"_d{deg - 1}"
        if deg == 3:
            if any(_is_sp2_class(n, res_atoms) for n in atom.GetNeighbors()):
                return "N_R", "_d2"
            return "N_3", ""
        if deg == 2:
            return "N_R", "_d1"
        _fail(atom)
    if sym == "O":
        if resonant:
            return "O_R", ""
        if deg == 1:
            nb = atom.GetNeighbors()[0]
            if nb.GetSymbol() in ("S", "P") or has2:
                return "O_2", ""
            if has3:
                return "O_1", ""
            return "O_3", ""
        if deg == 3:
            warnings.append(f"원자 {idx + 1} O: 3배위 O 는 DREIDING 에 없어 O_3 로 둡니다")
        return "O_3", ""
    if sym in ("F", "Cl", "Br", "I"):
        if sym == "Cl":
            warnings.append("dreiding.lt 의 Cl 질량이 34.45 로 들어가 있습니다(실제 35.45)")
        return sym, ""
    if sym == "S":
        if deg > 2 or has2:
            warnings.append(
                f"원자 {idx + 1} S: {deg}배위/다중결합 S 는 DREIDING 원 논문에 없어 S_3 로 둡니다 "
                "(평형각 92.1°, 이면각 장벽 나눗셈 가정이 실제와 다름)")
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


def _h_flag(atom, res_atoms):
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
    if atom.GetIdx() in res_atoms:
        return ""
    orders = [b.GetBondTypeAsDouble() for b in atom.GetBonds()]
    if 3.0 in orders or 2.0 in orders or deg == 2:
        return "_ha"
    if deg == 3 and not any(_is_sp2_class(n, res_atoms) for n in atom.GetNeighbors()):
        return "_ha"
    return ""


def _expected_orders(mol, res_bonds):
    exp = {}
    for b in mol.GetBonds():
        key = tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx())))
        if b.GetIsAromatic() or key in res_bonds:
            exp[key] = 1.5
        else:
            exp[key] = b.GetBondTypeAsDouble()
    return exp


def assign_labels(mol, ff, max_steps=1_000_000):
    """mol(명시적 H 포함, sanitize 완료), ff(DreidingFF) → (라벨 리스트, 경고 리스트)."""
    warnings = []
    groups = _resonance_groups(mol)
    res_atoms, res_bonds = set(), set()
    for c, os_ in groups:
        res_atoms.add(c)
        res_atoms.update(os_)
        res_bonds.update(tuple(sorted((c, o))) for o in os_)

    n = mol.GetNumAtoms()
    cands = []      # 원자별 후보 라벨 (선호 순서)
    wanted_h = []
    for atom in mol.GetAtoms():
        eh, d = _base_label(atom, res_atoms, warnings)
        h = _h_flag(atom, res_atoms)
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

    exp = _expected_orders(mol, res_bonds)
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

    def ok(i):
        return all(bond_ok(i, j) for j in adj[i] if chosen[j] is not None)

    # 분기 한정 탐색: 결합 차수·이면각 규칙을 모두 만족하는 조합 중 비용 최소
    best, best_cost, steps = None, float("inf"), 0

    def solve(k, acc):
        nonlocal best, best_cost, steps
        if acc >= best_cost:
            return
        if k == len(order):
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
        raise TypingError("결합 차수와 이면각 규칙을 모두 만족하는 _b1/_b2 플래그 조합이 없습니다 "
                          "(DREIDING 라벨 체계로 표현할 수 없는 결합 패턴)"
                          + (f" — 고정 라벨끼리 맞지 않는 결합: {', '.join(bad)}" if bad else ""))
    if steps > max_steps:
        warnings.append("플래그 탐색 단계 한도에 걸려 최소 비용 해가 아닐 수 있습니다(결합 차수는 만족)")
    chosen = best

    # sp3/할로젠/S 쪽이 낀 결합(S=O 등)은 DREIDING 한계라 경고
    for (i, j), e in exp.items():
        imp = ff.bond_order(chosen[i], chosen[j])
        if imp is None:
            raise TypingError(f"결합 {i + 1}-{j + 1} ({chosen[i]}–{chosen[j]}): dreiding.lt 에 결합 규칙이 없습니다")
        if abs(imp - e) > 1e-6:
            warnings.append(f"결합 {i + 1}-{j + 1} ({chosen[i]}–{chosen[j]}): 실제 차수 {e:g} 를 "
                            f"DREIDING 에서 {imp:g} 로 둡니다")
    # N_3 라벨은 이웃 3개를 가정해 이면각 장벽을 나누므로 4배위 N 은 장벽이 1.5배 커짐
    for i, lab in enumerate(chosen):
        if lab.startswith("N_3") and deg[i] == 4:
            warnings.append(f"원자 {i + 1} N: 4배위 N({lab})은 dreiding.lt 가 이웃 3개를 가정해 "
                            "이 N 을 지나는 이면각 장벽이 1.5배로 들어갑니다")
    for i, h in enumerate(wanted_h):
        if h and not chosen[i].endswith(h):
            warnings.append(f"원자 {i + 1}: '{chosen[i]}{h}' 라벨이 없어 수소결합 플래그 {h} 를 뺐습니다")
    return chosen, warnings
