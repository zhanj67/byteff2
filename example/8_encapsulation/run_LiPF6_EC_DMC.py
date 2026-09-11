#!/usr/bin/env python
"""
Minimal test: 1:10 LiPF6 in EC (or a solvent mixture)
Calculate: conductivity, viscosity, dielectric constant
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from byteff2.toolkit.properties_calculator import PropertiesCalculator

calc = PropertiesCalculator(
    solvent=["EC", "DMC"],
    solvent_ratio="3:7",
    anion="PF6",
    li_count=34,
    salt_to_solvent_ratio_str="1:10",
    base_dir="./test_LiPF6_ECDMC",
)

# Calculate conductivity, viscosity, dielectric (skip density)
results = results = calc.calculate(properties=["density"])

#results = calc.calculate(properties=["conductivity", "viscosity", "dielectric"])
print(json.dumps(results, indent=2))
