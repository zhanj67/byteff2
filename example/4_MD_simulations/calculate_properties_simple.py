#!/usr/bin/env python
"""
Simple usage examples of PropertiesCalculator with the new simplified interface.

Input format:
- solvent: Solvent name (e.g., "DMC", "EC")
- anion: Anion name (e.g., "PF6", "TFSI")
- ion_count: Number of ion pairs
- Default: 5000 atoms, lithium cation (LI)
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from byteff2.toolkit.properties_calculator import PropertiesCalculator


def example_1_simple():
    """Example 1: Simplest usage."""
    print("\n" + "="*60)
    print("EXAMPLE 1: Simplest Usage")
    print("="*60)

    # Input: solvent, anion, ion_count
    # Defaults: 5000 atoms, LI cation, 298K
    calculator = PropertiesCalculator(
        solvent="DMC",
        anion="PF6",
        ion_count=34,
        base_dir="./example_simple_dmc_pf6",
    )

    print(f"System: DMC + PF6 + LI (default)")
    print(f"Ion count: {calculator.ion_count}")
    print(f"Components: {calculator.components}")

    # Calculate only density
    results = calculator.calculate(properties=["density"])
    print(f"\nDensity Result: {results['density']}")


def example_2_multiple_solvents():
    """Example 2: Test with different solvents."""
    print("\n" + "="*60)
    print("EXAMPLE 2: Multiple Solvent Systems")
    print("="*60)

    systems = [
        ("DMC", "PF6", 34),
        ("EC", "TFSI", 20),
        ("H2O", "TFSI", 10),
    ]

    for solvent, anion, ion_count in systems:
        print(f"\nCalculating: {solvent} + {anion}")
        try:
            calc = PropertiesCalculator(
                solvent=solvent,
                anion=anion,
                ion_count=ion_count,
                base_dir=f"./example_{solvent}_{anion}",
            )
            print(f"  ✓ Valid system")
            print(f"  Components: {calc.components}")
        except Exception as e:
            print(f"  ✗ Error: {e}")


def example_3_auto_ion_count():
    """Example 3: Auto-calculate ion count."""
    print("\n" + "="*60)
    print("EXAMPLE 3: Auto-Calculate Ion Count")
    print("="*60)

    # Don't specify ion_count, let it auto-calculate
    calculator = PropertiesCalculator(
        solvent="EC",
        anion="PF6",
        ion_count=None,  # Will be auto-calculated
        natoms=5000,
        base_dir="./example_auto_ion_count",
    )

    print(f"Auto-calculated ion count: {calculator.ion_count}")
    print(f"Components: {calculator.components}")


def example_4_parametric_study():
    """Example 4: Parametric study with varying ion counts."""
    print("\n" + "="*60)
    print("EXAMPLE 4: Parametric Study (Different Ion Counts)")
    print("="*60)

    ion_counts = [10, 20, 30, 40, 50]

    print("\nParametric study: DMC + PF6 with varying ion counts\n")
    for ion_count in ion_counts:
        calc = PropertiesCalculator(
            solvent="DMC",
            anion="PF6",
            ion_count=ion_count,
            base_dir=f"./parametric_dmc_pf6_{ion_count}",
        )
        print(f"Ion count: {ion_count:2d} | Components: {calc.components}")


if __name__ == "__main__":
    print("\n" + "="*60)
    print("PropertiesCalculator - Simplified Interface Examples")
    print("="*60)

    try:
        example_1_simple()
    except Exception as e:
        print(f"Example 1 error: {e}")

    try:
        example_2_multiple_solvents()
    except Exception as e:
        print(f"Example 2 error: {e}")

    try:
        example_3_auto_ion_count()
    except Exception as e:
        print(f"Example 3 error: {e}")

    try:
        example_4_parametric_study()
    except Exception as e:
        print(f"Example 4 error: {e}")

    print("\n" + "="*60)
    print("CLI Usage Examples:")
    print("="*60)
    print("""
# Basic usage (calculate all properties)
python -m byteff2.toolkit.properties_calculator \\
    --solvent DMC --anion PF6 --ion-count 34 \\
    --work-dir ./my_simulation

# Calculate only density
python -m byteff2.toolkit.properties_calculator \\
    --solvent EC --anion TFSI --ion-count 20 \\
    --calculate density \\
    --work-dir ./density_only

# Custom natoms and temperature
python -m byteff2.toolkit.properties_calculator \\
    --solvent H2O --anion TFSI --ion-count 50 \\
    --natoms 8000 --temperature 313.0 \\
    --work-dir ./aqueous_system

# Auto-calculate ion count
python -m byteff2.toolkit.properties_calculator \\
    --solvent DMC --anion PF6 \\
    --natoms 5000 --work-dir ./auto_ions
    """)
