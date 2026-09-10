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

## CLI Usage

```bash
# 1:10 LiPF6 in EC (default Li=34)
python -m byteff2.toolkit.properties_calculator \
    --solvent EC --anion PF6 --ratio 1:10 \
    --calculate density,conductivity

# Custom: 1:8 with 50 Li
python -m byteff2.toolkit.properties_calculator \
    --solvent EC --anion TFSI --li-count 50 --ratio 1:8
```

## Parameters

- `solvent`: Solvent name (DMC, EC, H2O, etc.)
- `anion`: Anion name (PF6, TFSI, FSI, etc.)
- `li_count`: Number of Li ions (default: 34)
- `salt_to_solvent_ratio_str`: Format "1:10" means 1 salt pair : 10 solvent molecules (default: "1:10")
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
