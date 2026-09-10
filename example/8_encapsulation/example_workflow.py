#!/usr/bin/env python
"""
Example: Using PropertiesCalculator in a workflow
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from byteff2.toolkit.properties_calculator import PropertiesCalculator


def screen_ionic_liquids():
    """Example: Screen different ionic liquids"""

    systems = [
        ("EC", "PF6", 34, "1:10"),
        ("EC", "TFSI", 34, "1:10"),
        ("DMC", "PF6", 34, "1:10"),
    ]

    results_summary = []

    for solvent, anion, li_count, ratio in systems:
        print(f"\nProcessing: {li_count} Li with {solvent}/{anion} ({ratio})")

        calc = PropertiesCalculator(
            solvent=solvent,
            anion=anion,
            li_count=li_count,
            salt_to_solvent_ratio_str=ratio,
            base_dir=f"./results_{solvent}_{anion}"
        )

        print(f"  Components: {calc.components}")
        # results = calc.calculate(properties=["density", "conductivity"])
        # results_summary.append({"system": f"{solvent}/{anion}", "results": results})

    print(f"\n✓ Workflow complete")


def parametric_study_li_concentration():
    """Example: Vary Li concentration"""

    li_counts = [20, 34, 50]

    for li_count in li_counts:
        print(f"\nLi count: {li_count}")

        calc = PropertiesCalculator(
            solvent="EC",
            anion="PF6",
            li_count=li_count,
            salt_to_solvent_ratio_str="1:10",
            base_dir=f"./results_li{li_count}"
        )

        print(f"  Components: {calc.components}")
        # results = calc.calculate(properties=["density", "conductivity"])


if __name__ == "__main__":
    print("\n" + "="*60)
    print("Example: PropertiesCalculator in Workflows")
    print("="*60)

    print("\nExample 1: Screen ionic liquids")
    screen_ionic_liquids()

    print("\n\nExample 2: Parametric study")
    parametric_study_li_concentration()

    print("\n" + "="*60)
    print("✓ Examples complete")
    print("="*60)
