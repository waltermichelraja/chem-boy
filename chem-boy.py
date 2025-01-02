from sympy import *
import math
init_printing(use_unicode=True)


def simplify(formula_in):
    def open_brackets(formula):
        if formula.find("(") >= 0:
            formula_br = temp = ""
            m = 1
            for i in range(0, len(formula)):
                if i + 1 < len(formula):
                    if formula[i].isalpha() and (formula[i + 1].isupper() or formula[i + 1] == "(" or formula[i + 1] == ")"):
                        temp += formula[i] + "1"
                    elif formula[i] == ")" and formula[i + 1].isalpha():
                        temp += formula[i] + "1"
                    else:
                        temp += formula[i]
                else:
                    if formula[i].isalpha():
                        temp += formula[i] + "1"
                    elif formula[i] == ")":
                        temp += formula[i] + "1"
                    else:
                        temp += formula[i]
            formula = str(temp)
            for i in range(formula.index(")"), len(formula)):
                if str(formula[i]).isalpha() or str(formula[i]).isspace() or i == len(formula) - 1:
                    if (i == len(formula) - 1):
                        i += 1
                    formula_br = formula[formula.index("("):i]
                    break
            m *= (int(formula_br[formula_br.index(")") + 1:]))
            val = ""
            temp = formula_br
            i = 1
            while i < temp.index(")") + 1:
                if temp[i].isdigit():
                    val += str(temp[i])
                else:
                    if len(val) > 0:
                        temp = temp.replace(
                            temp[1:i], (temp[1:i - len(val)] + str(int(val) * m)))
                        val = ""
                i += 1
            formula = formula.replace(formula_br, temp[1: temp.index(")")])
        return formula


    if formula_in.find(".") >= 0:
        formula_in = formula_in[0: formula_in.index(
            ".")] + "." + open_brackets(formula_in[formula_in.index(".") + 1:])
        d = ""
        for i in range(formula_in.index(".") + 1, len(formula_in)):
            if formula_in[i].isdigit():
                d += formula_in[i]
            else:
                break
        formula_in = formula_in[:formula_in.index(
            ".")] + open_brackets("(" + formula_in[formula_in.index(".") + len(d) + 1:] + ")" + d)
    formula_in = open_brackets(formula_in)
    return formula_in


def find_lcm(num1, num2):  # function to find the lcm of a list of numbers
    lcm = int(int(num1 * num2) / int(math.gcd(num1, num2)))
    return lcm


def sub_script(s):
    for i in range(1, len(s)):
        if s[i].isdigit() and (s[i - 1].isalpha() or s[i - 1] == ")"):
            s = s[0:i] + s[i].replace("0", "\u2080").replace("1", "\u2081").replace("2", "\u2082").replace("3", "\u2083").replace("4", "\u2084").replace(
                "5", "\u2085").replace("6", "\u2086").replace("7", "\u2087").replace("8", "\u2088").replace("9", "\u2089") + s[i + 1:]
    return s


class compound(object):  # A class of compounds. It stores all the relevant data for the compound
    def __init__(self, n_compound):
        self.n_compound = str(n_compound)
        self.f_compound = simplify((str(n_compound)))
        temp = ""
        e = ""
        v = "0"
        self.element = []
        self.val = []
        for i in range(0, len(self.f_compound) - 1):
            if self.f_compound[i].isalpha() and self.f_compound[i + 1].isupper():
                temp += self.f_compound[i] + "1"
            else:
                temp += self.f_compound[i]
        temp += self.f_compound[len(self.f_compound) - 1]
        if temp[len(temp) - 1].isalpha():
            temp += "1"
        self.f_compound = temp
        for i in range(0, len(self.f_compound)):
            if self.f_compound[i].isalpha():
                if v != "0":
                    if e in self.element:
                        self.val[self.element.index(e)] = int(
                            self.val[self.element.index(e)]) + int(v)
                    else:
                        self.element.append(str(e))
                        self.val.append(int(v))
                    e = self.f_compound[i]
                    v = "0"
                    i -= 1
                else:
                    e += self.f_compound[i]
            elif self.f_compound[i].isdigit():
                v += self.f_compound[i]

        if e in self.element:
            self.val[self.element.index(e)] = int(
                self.val[self.element.index(e)]) + int(v)
        else:
            self.element.append(str(e))
            self.val.append(int(v))


equation = str(input("Enter a chemical equation: "))

equation = equation.replace("->", "+").replace(' ', '')
compounds = []  # An array of compounds
elements = []  # An array of all the elements present in the equation

# Assigning a compound name to an object of the class compound and storing it in the array
for i in range(0, len(equation.split("+"))):
    compounds.append(compound(equation.split("+")[i]))
# Checking if element has already been added to the 'elements' array. Addidng it if it has not been added
for i in range(0, len(compounds)):
    for x in range(0, len(compounds[i].element)):
        if not(compounds[i].element[x] in elements):
            elements.append(compounds[i].element[x])

# Number of columns in the matrix. There is one column for every compound in the equation
cols = len(compounds)
# Number of rows in the matrix. There is one row for every element in the equation
rows = len(elements)
m = (zeros(rows, cols))

# The number of atoms of an element present in a particular compound is stored in the matrix
for c in range(0, int(cols)):
    for r in range(0, int(rows)):
        try:
            m[r, c] = compounds[c].val[compounds[c].element.index(elements[r])]
        except:
            m[r, c] = 0

print(m)
print("\nRREF Matrix:\n")
print(m.rref())
m = (list(m.rref()))[0]  # The matrix is converted to Row Reduced Echelon Form

coefficients = []  # Array to store the coefficients of the compounds
denominator = []  # Array to store the denominator of the coefficients of the compounds

for r in range(0, rows):
    if str(m[r, cols - 1]) == "0":
        break
    elif str(m[r, cols - 1]).find("/") > 0:
        coefficients.append(int(str(m[r, cols - 1]).split("/")[0]))
        denominator.append(int(str(m[r, cols - 1]).split("/")[1]))
    else:
        coefficients.append(int(str(m[r, cols - 1])))
        denominator.append(1)
coefficients.append(int(m[0, 0]))
denominator.append(1)


# LCM of the denominators is calculated
lcm = find_lcm(int(denominator[0]), int(denominator[1]))
for i in range(2, len(denominator)):
    lcm = find_lcm(lcm, denominator[i])

x = 0
# The coefficients are multiplied by the LCM
for i in range(0, len(coefficients)):
    if coefficients[i] < 0:
        x += 1
    if x > 1:
        coefficients[i] = abs(coefficients[i])
    coefficients[i] = int(coefficients[i] * lcm / denominator[i])

equation = ""

# Concatenating the coefficients and the molecular formula of the compounds
for i in range(0, len(compounds)):
    if coefficients[i] < 0:
        equation = equation.strip()[:-2].strip() + " \u2794 "
    equation += (str(abs(coefficients[i])) if (abs(coefficients[i])) > 1 else "") + \
        "" + sub_script(compounds[i].n_compound)
    if not(i == len(compounds) - 1):
        equation += " + "
print("\nThe balanced equation is : " + equation)
