#!/usr/bin/env python
"""
Test: 1:10 LiPF6 in EC
- Density + Conductivity
- Python interface encapsulation
"""

import sys
import json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from byteff2.toolkit.properties_calculator import PropertiesCalculator

def main():
    print("\n" + "="*70)
    print("Encapsulation Test: 1:10 LiPF6 in EC")
    print("="*70)

    try:
        # Create calculator
        calc = PropertiesCalculator(
            solvent="EC",
            anion="PF6",
            li_count=34,
            salt_to_solvent_ratio_str="1:10",
            temperature=298.0,
            base_dir="./test_LiPF6_EC"
        )

        print(f"\n✓ Calculator initialized")
        print(f"  Components: {calc.components}")
        print(f"  Temperature: {calc.temperature} K")

        # Calculate density and conductivity
        print(f"\n Calculating density and conductivity...")
        results = calc.calculate(properties=["density", "conductivity"])

        print(f"\n✓ Calculation complete")
        print(f"\nResults:")
        print(json.dumps(results, indent=2))

        return 0

    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
