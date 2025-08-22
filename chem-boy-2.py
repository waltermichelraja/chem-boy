from sympy import *
import math
import re

# author: Walter
init_printing(use_unicode=True)

equation = input("Enter the equation (use -> between reactants and products): \n").replace(" ", "")

# --- Utility Functions ---

def subscript(s: str) -> str:
    """Convert digits to Unicode subscript for pretty display."""
    return "".join(chr(int(ch) + 8320) if ch.isdigit() else ch for ch in s)


def format_compound(formula: str) -> str:
    """Insert implicit '1's into molecular formulas for parsing."""
    f_formula = ""
    for i in range(len(formula)):
        if i + 1 < len(formula):
            if formula[i].isalpha() and (formula[i + 1].isupper() or formula[i + 1] in "([)]"):
                f_formula += formula[i] + "1"
            elif formula[i] in ")]" and formula[i + 1].isalpha():
                f_formula += formula[i] + "1"
            else:
                f_formula += formula[i]
        else:
            f_formula += formula[i] + ("1" if formula[i].isalpha() or formula[i] in ")]" else "")
    return f_formula


def open_brackets(formula: str) -> str:
    """Expand formulas with parentheses, e.g., Ca(OH)2 -> Ca1O2H2."""
    while "(" in formula or ")" in formula:
        match = re.search(r"\(([A-Za-z0-9]+)\)(\d+)", formula)
        if not match:
            break
        inside, mult = match.groups()
        expanded = ""
        for elem, count in zip(re.findall(r"[A-Z][a-z]*", inside), re.findall(r"\d+", inside)):
            expanded += elem + str(int(count) * int(mult))
        formula = formula.replace(match.group(0), expanded, 1)
    return formula


def simplify(formula: str) -> str:
    """Fully simplify molecular formula."""
    formula = format_compound(formula)
    formula = open_brackets(formula.replace("[", "(").replace("]", ")"))
    return formula


# --- Data Structures ---

all_elements = []

class Compound:
    """Represents a chemical compound with element counts."""
    def __init__(self, formula: str):
        self.compound_n = formula
        self.compound_f = simplify(formula)
        self.elements = {}
        for elem, val in zip(re.findall(r"[A-Z][a-z]*", self.compound_f), re.findall(r"\d+", self.compound_f)):
            self.elements[elem] = self.elements.get(elem, 0) + int(val)
            if elem not in all_elements:
                all_elements.append(elem)


# --- Parse Input ---
lhs, rhs = equation.split("->")
compounds = [Compound(f) for f in lhs.split("+")] + [Compound(f) for f in rhs.split("+")]

cols, rows = len(compounds), len(all_elements)
m = zeros(rows, cols)

for c in range(cols):
    for r in range(rows):
        m[r, c] = compounds[c].elements.get(all_elements[r], 0) * (1 if c < len(lhs.split("+")) else -1)

# --- Solve System ---
m_rref, _ = m.rref()
coeffs = []
denoms = []

for val in m_rref[:, -1]:
    num, den = val.as_numer_denom()
    coeffs.append(int(num))
    denoms.append(int(den))

coeffs.append(1)
denoms.append(1)

# scale to integers
lcm = denoms[0]
for d in denoms[1:]:
    lcm = lcm * d // math.gcd(lcm, d)

coeffs = [c * (lcm // d) for c, d in zip(coeffs, denoms)]

# --- Construct Balanced Equation ---
balanced = ""
for i, cmpd in enumerate(compounds):
    balanced += (str(coeffs[i]) if coeffs[i] > 1 else "") + cmpd.compound_n
    if i + 1 == len(lhs.split("+")):
        balanced += " → "
    elif i < len(compounds) - 1:
        balanced += " + "

print("\nBalanced equation:\n" + subscript(balanced))
