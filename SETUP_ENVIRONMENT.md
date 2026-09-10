# Environment Setup Guide

## Quick Setup

### Option 1: Using requirements.txt (Recommended)

```bash
# Clone the repo
git clone https://github.com/zhanj67/byteff2.git
cd byteff2

# Create virtual environment
python3.10 -m venv venv310

# Activate
source venv310/bin/activate  # Linux/Mac
# or
venv310\Scripts\activate     # Windows

# Install dependencies
pip install -r requirements_venv310.txt
```

### Option 2: From scratch

```bash
python3.10 -m venv venv310
source venv310/bin/activate
pip install --upgrade pip
pip install numpy scipy pandas openmm ase rdkit bytemol
```

## Test Installation

```bash
source venv310/bin/activate
cd example/8_encapsulation
python example_workflow.py
```

## Using PropertiesCalculator

```bash
source venv310/bin/activate
python -m byteff2.toolkit.properties_calculator --help
```
