from sympy import *
init_printing(use_unicode=True)
import math
import re
# H3PO4 + (NH4)2MoO4 + HNO3 + H2O-> (NH4)3PO4.12MoO3 + NH4NO3


def subscript(s):  # Function to convert number to subscript. 8320 is the ascii of subscript 0
    i = 0
    while i < len(s):
        if s[i].isdigit() and (s[i - 1].isalpha() or s[i - 1] == ")" or ord(s[i - 1]) in range(8320, 8340)):
            s = s[0:i] + "[sub]" + s[i] + "[/sub]" + s[i + 1:]
            print(s)
        i += 1
    s = s.replace("[b]", "[b][color=000080]").replace("[/b]", "[/b][/color]")
    return s


def format_compound(formula):
    f_formula = ""
    for part in formula.split("."):
        for i in range(0, len(part)):
            if i + 1 < len(part):
                if part[i].isalpha() and (part[i + 1].isupper() or part[i + 1] in ["(", "[", ")", "]"]):
                    f_formula += part[i] + "1"
                elif part[i] in [")", "]"] and part[i + 1].isalpha():
                    f_formula += part[i] + "1"
                else:
                    f_formula += part[i]
            else:
                if part[i].isalpha() or part[i] in [")", "]"]:
                    f_formula += part[i] + "1"
                else:
                    f_formula += part[i]
        f_formula += ("." if len(formula.split(".")) > 1 else "")
    if f_formula.find(".") > -1:
        part = f_formula.split(".")
        val = str(re.findall('[0-9][0-9]*', part[1])[0])
        val = ("1" if val == "" else val)
        f_formula = part[0] + \
            "(" + str(re.findall('[A-Za-z]\S*', part[1])[0]) + ")" + val
    return f_formula


def open_brackets(formula):  # Function to open the brackets
    formula_br = re.findall('\((.*?)\)', formula)
    m = re.findall('\)([0-9][0-9]*)', formula)
    for (i, indexc) in zip(formula_br, range(0, len(formula_br))):
        f_formula = ""
        val = re.findall('[0-9][0-9]*', i)
        for (n, indexv) in zip(val, range(0, len(val))):
            f_formula += re.findall('([A-Za-z][A-Za-z]*)', i)[indexv]
            f_formula += str(int(n) * int(m[indexc]))
        formula = formula.replace(re.findall(
            '(\(.*?\)[0-9]*)', formula)[0], f_formula)
    return formula


def simplify(formula):  # Function to simplify the molecular formula by using the above functions
    formula = format_compound(formula)
    formula = open_brackets(formula)
    formula = formula.replace("[", "(").replace("]", ")")
    formula = open_brackets(formula)
    return formula


all_elements = list()  # A list of all the elements in the list


class compound(object):  # A class of compounds. I contains all the attributes of each compound
    def __init__(self, compound_n):
        self.compound_n = compound_n  # It stores the original molecular formula
        # It stores the simplified molecular formula
        self.compound_f = simplify(self.compound_n)
        # Dict containing all the elements and number of atoms of each element
        self.elements = dict()
        for (val, element) in zip(re.findall('[0-9][0-9]*', self.compound_f), re.findall('[A-Za-z][A-Za-z]*', self.compound_f)):
            self.elements[element] = str(
                int(self.elements.get(element, "0")) + int(val))
            if element not in all_elements:
                all_elements.append(element)


def balance(equation):
    if not isValid(equation):
        return "Invalid Equation"

    try:
        equation = equation.replace("->", "+").replace(" ", "")

        compounds = []  # A list of objects of the class 'compound'
        for compound_n in equation.split("+"):
            compounds.append(compound(compound_n))

        cols = len(compounds)
        rows = len(all_elements)
        # Matrix which solves the system of linear equations
        m = (zeros(rows, cols))

        for c in range(0, int(cols)):  # Inserting the values into the matrix
            for r in range(0, int(rows)):
                try:
                    m[r, c] = compounds[c].elements[all_elements[r]]
                except:
                    m[r, c] = 0

        # Converting the matrix to RREF form and removing the list of pivots
        m = list(m.rref())[0]

        coefficients = list()
        denominator = list()

        for r in range(0, rows):  # Storing the coefficints in the lists
            if m[r, cols - 1] == 0:
                break
            else:
                coefficients.append(int(str(m[r, cols - 1]).split("/")[0]))
                try:
                    denominator.append(int(str(m[r, cols - 1]).split("/")[1]))
                except:
                    denominator.append(1)
        coefficients.append(1)
        denominator.append(1)

        lcm = denominator[0]
        for num in denominator[1:]:  # Calculating the lcm of the denominators
            lcm = int(int(lcm * num) / int(math.gcd(lcm, num)))

        rhs = False
        equation = ""
        # Concatenating the coefficients and the molecular formula
        for i in range(0, len(compounds)):
            coefficients[i] = int(coefficients[i] / denominator[i] * lcm)
            if coefficients[i] < 1 and rhs == false:
                equation = equation.strip()[:-2].strip() + " --> "
                rhs = True
            equation += "[b]" + (str(abs(coefficients[i])) if abs(
                coefficients[i]) > 1 else "") + "[/b]"
            equation += compounds[i].compound_n + \
                (" + " if i < len(compounds) - 1 else "")

        equation = subscript(equation)
        return equation

    except:
        return "Invalid Equation"


def isValid(equation):
    for i in equation:
        if i.isalpha() == False and i.isdigit() == False and i not in ['(', ')', '[', ']', '-', '+', ' ', '>', '.']:
            return False
    return True


print(balance("Al2(SO3)3 + NaOH -> Na2SO3 + Al(OH)32"))