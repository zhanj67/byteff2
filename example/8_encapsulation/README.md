# PropertiesCalculator Encapsulation

Unified interface for ByteFF2 MD property calculations.

## Quick Start

```python
from byteff2.toolkit.properties_calculator import PropertiesCalculator

# 1:10 LiPF6 in EC (34 Li, 340 EC)
calc = PropertiesCalculator(
    solvent="EC",
    anion="PF6",
    li_count=34,
    salt_to_solvent_ratio_str="1:10",
    base_dir="./results"
)

# Calculate density + conductivity
results = calc.calculate(properties=["density", "conductivity"])
```

### Mixed solvents

`salt_to_solvent_ratio_str` fixes the *total* number of solvent molecules;
`solvent_ratio` splits that total between the solvents (largest-remainder
rounding, so the counts always add up exactly).

```python
# 1:10 LiPF6 in EC:DMC = 3:7 (34 Li, 102 EC, 238 DMC)
calc = PropertiesCalculator(
    solvent=["EC", "DMC"],
    solvent_ratio="3:7",          # or [3, 7]; omit for equal parts
    anion="PF6",
    li_count=34,
    salt_to_solvent_ratio_str="1:10",
    base_dir="./results_ec_dmc"
)

# Equivalent dict form: proportions are the dict values
calc = PropertiesCalculator(solvent={"EC": 3, "DMC": 7}, anion="PF6")
```

## CLI Usage

```bash
# 1:10 LiPF6 in EC (default Li=34)
python -m byteff2.toolkit.properties_calculator \
    --solvent EC --anion PF6 --ratio 1:10 \
    --calculate density,conductivity

# Custom: 1:8 with 50 Li
python -m byteff2.toolkit.properties_calculator \
    --solvent EC --anion TFSI --li-count 50 --ratio 1:8

# Mixed solvent: EC:DMC = 3:7
python -m byteff2.toolkit.properties_calculator \
    --solvent EC,DMC --solvent-ratio 3:7 --anion PF6 --ratio 1:10
```

## Parameters

- `solvent`: Solvent name (DMC, EC, H2O, etc.), a list of names for a mixture, or a dict `{name: proportion}`
- `solvent_ratio`: Molar proportions between the solvents, "3:7" or [3, 7] (default: equal parts; ignored for the dict form)
- `anion`: Anion name (PF6, TFSI, FSI, etc.)
- `li_count`: Number of Li ions (default: 34)
- `salt_to_solvent_ratio_str`: Format "1:10" means 1 salt pair : 10 solvent molecules in total (default: "1:10")
- `temperature`: Temperature in K (default: 298.0)
- `base_dir`: Working directory

## Supported Properties

- `density` - Density calculation (~2-4 hours)
- `conductivity` - Ionic conductivity via Onsager formula (~12-24 hours)
- `viscosity` - Viscosity via NEMD (~12-24 hours)
- `dielectric` - Dielectric constant (~12-24 hours)

## Test Status

✓ Interface initialization working  
✓ Configuration generation working  
✓ MD workflow integration ready  

Full simulations require significant compute time. Start with density calculations for quick validation.
