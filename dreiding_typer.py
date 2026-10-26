# SMILES 결합 차수로 moltemplate dreiding.lt 라벨(원자 타입)을 할당하는 타이퍼
"""RDKit 분자(명시적 H 포함)에 moltemplate DREIDING 라벨을 붙인다.

라벨 체계: m-bone/moltemplate-DREIDING 의 DREIDING_Label_Manual.pdf
  <원소>_<혼성>[_b1|_b2][_dN][_hd|_ha]
결합 차수는 라벨 쌍으로 정해지므로(implied_order), _b1/_b2 플래그를
제약 만족 탐색으로 골라 실제 결합 차수와 일치시킨다.
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


def load_labels(path=None):
    """dreiding.lt 의 Data Masses 블록에서 허용 라벨과 질량을 읽는다."""
    path = path or default_lt_path()
    labels, inside = {}, False
    with open(path) as f:
        for line in f:
            if 'write_once("Data Masses")' in line:
                inside = True
                continue
            if inside:
                if line.strip().startswith("}"):
                    break
                m = re.match(r"\s*@atom:(\S+)\s+(\S+)", line)
                if m:
                    labels[m.group(1)] = float(m.group(2))
    if not labels:
        raise TypingError(f"{path} 에서 Data Masses 라벨을 찾지 못했습니다")
    return labels


def sp2_class(label):
    """결합 차수 플래그(_b1/_b2)가 의미 있는 sp2/공명/sp 계열 라벨인지."""
    return any(t in label for t in ("_2", "_R", "_1"))


def implied_order(a, b):
    """라벨 쌍이 dreiding.lt 에서 갖는 결합 차수 (Dreiding_bond.py 의 판정 순서 그대로)."""
    if "_b1" in a and "_b1" in b:
        return 1.0
    if "_b2" in a and "_b2" in b:
        return 2.0
    if ("_b2" in a and "_2" in b) or ("_2" in a and "_b2" in b):
        return 1.0
    if "_1" in a and "_1" in b:
        return 3.0
    if "_2" in a and "_2" in b:
        return 2.0
    if "_R" in a and "_R" in b:
        return 1.5
    return 1.0


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


def assign_labels(mol, labels, max_steps=1_000_000):
    """mol(명시적 H 포함, sanitize 완료) → (라벨 리스트, 경고 리스트)."""
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
        c = [eh + b + d + hh for hh in hs for b in bflags if eh + b + d + hh in labels]
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

    def ok(i):
        for j in adj[i]:
            if chosen[j] is None or not (sp2_class(chosen[i]) and sp2_class(chosen[j])):
                continue
            key = tuple(sorted((i, j)))
            if abs(implied_order(chosen[i], chosen[j]) - exp[key]) > 1e-6:
                return False
        return True

    # 분기 한정 탐색: 결합 차수를 모두 만족하는 조합 중 비용 최소
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
        raise TypingError("모든 결합 차수를 만족하는 _b1/_b2 플래그 조합이 없습니다 "
                          "(DREIDING 라벨 체계로 표현할 수 없는 결합 패턴)")
    if steps > max_steps:
        warnings.append("플래그 탐색 단계 한도에 걸려 최소 비용 해가 아닐 수 있습니다(결합 차수는 만족)")
    chosen = best

    # sp3/할로젠/S 쪽이 낀 결합(S=O 등)은 DREIDING 한계라 경고, sp2 끼리 불일치는 오류
    for (i, j), e in exp.items():
        imp = implied_order(chosen[i], chosen[j])
        if abs(imp - e) > 1e-6:
            if sp2_class(chosen[i]) and sp2_class(chosen[j]):
                raise TypingError(f"결합 {i + 1}-{j + 1} ({chosen[i]}–{chosen[j]}): 실제 차수 {e:g} 와 "
                                  f"라벨 차수 {imp:g} 가 다릅니다")
            warnings.append(f"결합 {i + 1}-{j + 1} ({chosen[i]}–{chosen[j]}): 실제 차수 {e:g} 를 "
                            f"DREIDING 에서 {imp:g} 로 둡니다")
    for i, h in enumerate(wanted_h):
        if h and not chosen[i].endswith(h):
            warnings.append(f"원자 {i + 1}: '{chosen[i]}{h}' 라벨이 없어 수소결합 플래그 {h} 를 뺐습니다")
    return chosen, warnings
