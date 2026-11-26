from fractions import Fraction
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent

def read(path):
    return json.loads(path.read_bytes())

def check_paper_values(paths):
    table = read(ROOT / "paper_values.json")
    radii = read(paths["15"] / "summary.json")["endpoints"]
    bounds = read(paths["24"] / "summary.json")["sides"]
    for side, row in table.items():
        def prefix(interval, digits, exponent):
            scale = Fraction(10) ** exponent
            floor = Fraction(digits) * scale
            step = scale / 10 ** len(digits.split(".")[1])
            lo, hi = map(Fraction, interval)
            if not floor <= lo <= hi < floor + step:
                raise ValueError("Printed decimal prefix mismatch: " + side)
        radius = radii[side]["root_interval"]
        prefix([radius["lo"], radius["hi"]], row["radius_prefix"], 0)
        coefficient = read(paths["19"] / (side + "_prefactor.json"))["C"]
        digits = row["coefficient_mantissa"]
        scale = Fraction(10) ** row["coefficient_exponent"]
        printed = Fraction(digits) * scale
        half_unit = scale / (2 * 10 ** len(digits.split(".")[1]))
        if not (printed - half_unit < Fraction(coefficient["lower"]) <= Fraction(coefficient["upper"]) < printed + half_unit):
            raise ValueError("Coefficient differs at the paper's displayed precision: " + side)
        for key in ["lower", "upper"]:
            if Fraction(row["constant_" + key]) != Fraction(bounds[side]["D_" + key]):
                raise ValueError("Constant-term endpoint mismatch: " + side)
    print("Matched the printed radius and prefactor digits and the four constant-term bounds.")
