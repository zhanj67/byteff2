# Copyright (c) 2025 Bytedance Ltd. and/or its affiliates

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Reducing PropertiesCalculator output to a flat results file."""

import json
from pathlib import Path
from typing import Union

PROPERTIES = ["density", "dielectric", "viscosity", "conductivity"]

UNITS = {
    "density": "g/mL",
    "dielectric": "dimensionless",
    "viscosity": "cP",
    "conductivity": "mS/cm",
}


def flatten(results: dict) -> dict:
    """Lift the four property values out of the protocols' nested output.

    A property whose protocol failed comes out as None rather than going
    missing, so the mapping always has the same shape.
    """
    return {
        "density": (results.get("density") or {}).get("density"),
        "dielectric": (results.get("dielectric") or {}).get("dielectric"),
        "viscosity": results.get("viscosity"),
        "conductivity": results.get("conductivity"),
        "units": UNITS,
    }


def save(results: dict, base_dir: Union[str, Path]) -> dict:
    """Write the flattened results to <base_dir>/results.json."""
    out = flatten(results)
    path = Path(base_dir) / "results.json"
    path.write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))
    print(f"wrote {path}")
    return out
