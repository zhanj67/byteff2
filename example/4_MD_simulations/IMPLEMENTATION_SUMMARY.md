# PropertiesCalculator 实现总结

## 系统概览

已为 ByteFF2 工作流创建了统一的、可复用的属性计算封装，支持：
- ✅ 密度 (Density)
- ✅ 电导率 (Conductivity)  
- ✅ 黏度 (Viscosity)
- ✅ 介电常数 (Dielectric Constant)
- ✅ 压缩率 (Compressibility)

## 关键特性

### 1. 简化的输入格式

```
输入：溶剂名 + 阴离子名 + 离子对数量

默认：
- 原子数: 5000
- 阳离子: LI（锂，固定）
- 温度: 298K
```

### 2. 双重接口支持

**命令行接口 (CLI)**
```bash
python -m byteff2.toolkit.properties_calculator \
    --solvent DMC --anion PF6 --ion-count 34 \
    --work-dir ./my_sim
```

**Python 函数接口**
```python
from byteff2.toolkit.properties_calculator import PropertiesCalculator

calc = PropertiesCalculator(
    solvent="DMC",
    anion="PF6", 
    ion_count=34,
    base_dir="./my_sim"
)
results = calc.calculate()
```

### 3. 灵活的协议集成

- `DensityProtocol` → 计算密度
- `TransportProtocol` → 计算电导率和黏度
- `DielectricProtocol` → 计算介电常数
- `CompressibilityProtocol` → 计算压缩率

## 文件结构

```
/root/byteff2/byteff2/toolkit/
├── properties_calculator.py          ← 主要实现
└── protocol.py                       ← 协议定义（现有）

/root/byteff2/example/4_MD_simulations/
├── QUICK_START.md                    ← 快速入门
├── calculate_properties_simple.py    ← Python 示例
├── SMILES_DATABASE_GUIDE.md          ← 数据库扩展指南
└── IMPLEMENTATION_SUMMARY.md         ← 本文件
```

## 数据库管理

### 当前支持

**溶剂** (6 种)
- DMC, EC, EMC, DEC, PC, H2O

**阴离子** (5 种)
- PF6, BF4, ClO4, TFSI, OTf

**阳离子**（固定）
- LI（锂）

### 后续扩展

SMILES 数据可从以下来源扩展：
1. `/root/byteff2/example/7_paper_source_data/Figure_3_d.csv`
2. PubChem
3. ChemSpider
4. 手动输入

## 使用示例

### 例 1: 基础计算
```python
calc = PropertiesCalculator(
    solvent="EC",
    anion="TFSI",
    ion_count=20
)
results = calc.calculate()  # 计算所有属性
```

### 例 2: 选择性计算
```python
results = calc.calculate(
    properties=["density", "viscosity"]
)
```

### 例 3: 参数化研究
```python
for ion_count in [10, 20, 30, 40, 50]:
    calc = PropertiesCalculator(
        solvent="DMC",
        anion="PF6",
        ion_count=ion_count,
        base_dir=f"./parametric_li{ion_count}"
    )
    results = calc.calculate()
```

### 例 4: 自定义参数
```python
calc = PropertiesCalculator(
    solvent="H2O",
    anion="PF6",
    ion_count=50,
    natoms=8000,           # 自定义原子数
    temperature=313.0,     # 自定义温度
    base_dir="./aqueous"
)
```

## 输出结构

```
base_dir/
├── summary.json                      # 计算总结
├── density_config.json               # 配置文件
├── transport_config.json
├── dielectric_config.json
│
├── density_results/
│   ├── npt_state.csv                # MD 轨迹数据
│   ├── npt.dcd
│   └── density_results.json          # 计算结果
│
├── transport_results/
│   ├── npt_state.csv
│   ├── nvt_state.csv
│   ├── nonequ.dcd
│   └── results.json                  # 电导率、黏度结果
│
└── dielectric_results/
    ├── npt_state.csv
    ├── dipole.csv
    └── dielectric_results.json       # 介电常数结果
```

## 工作流集成

### 在你的代码中使用

```python
from byteff2.toolkit.properties_calculator import PropertiesCalculator

def my_material_screening_workflow():
    """材料筛选工作流"""
    
    candidates = [
        ("DMC", "PF6", 34),
        ("EC", "TFSI", 20),
        ("H2O", "PF6", 50),
    ]
    
    results_table = []
    
    for solvent, anion, ion_count in candidates:
        print(f"Processing {solvent}/{anion}...")
        
        calc = PropertiesCalculator(
            solvent=solvent,
            anion=anion,
            ion_count=ion_count,
            base_dir=f"./results/{solvent}_{anion}"
        )
        
        results = calc.calculate()
        
        # 处理结果
        summary = {
            "solvent": solvent,
            "anion": anion,
            "density": results["density"]["density"],
            "conductivity": results.get("conductivity"),
            "viscosity": results.get("viscosity"),
        }
        results_table.append(summary)
    
    return results_table
```

## 命令行快速参考

```bash
# 基础用法
python -m byteff2.toolkit.properties_calculator \
    --solvent DMC --anion PF6 --ion-count 34

# 只计算密度
python -m byteff2.toolkit.properties_calculator \
    --solvent EC --anion TFSI --ion-count 20 \
    --calculate density

# 自定义参数
python -m byteff2.toolkit.properties_calculator \
    --solvent H2O --anion PF6 --ion-count 50 \
    --natoms 8000 --temperature 313.0

# 多属性计算
python -m byteff2.toolkit.properties_calculator \
    --solvent DMC --anion BF4 --ion-count 30 \
    --calculate density,conductivity,viscosity,dielectric

# 查看帮助
python -m byteff2.toolkit.properties_calculator --help
```

## 配置文件示例

`config.json`:
```json
{
  "solvent": "DMC",
  "anion": "PF6",
  "ion_count": 34,
  "natoms": 5000,
  "temperature": 298.0,
  "base_dir": "./my_simulation"
}
```

## 主要设计决策

| 决策 | 原因 |
|------|------|
| 阳离子固定为 LI | 简化输入，Li 是常用的 |
| 默认 5000 原子 | 平衡计算成本和系统大小 |
| 简化输入为三个参数 | 易用性，减少错误 |
| 支持自动离子数计算 | 灵活性，支持多种工作流 |
| CLI + Python 接口 | 满足不同使用场景 |

## 待办事项

- [ ] 从 Figure_3_d.csv 批量导入溶剂 SMILES
- [ ] 实现配置文件支持
- [ ] 添加结果可视化功能
- [ ] 创建单元测试
- [ ] 添加更多溶剂和阴离子
- [ ] 性能优化

## 测试

运行示例代码：
```bash
python /root/byteff2/example/4_MD_simulations/calculate_properties_simple.py
```

测试 CLI：
```bash
python -m byteff2.toolkit.properties_calculator \
    --solvent DMC --anion PF6 --ion-count 10 \
    --work-dir ./test_sim
```

## 许可

Apache License 2.0

---

## 快速链接

- **快速开始**: [QUICK_START.md](./QUICK_START.md)
- **数据库扩展**: [SMILES_DATABASE_GUIDE.md](./SMILES_DATABASE_GUIDE.md)
- **Python 示例**: [calculate_properties_simple.py](./calculate_properties_simple.py)
- **核心代码**: `/root/byteff2/byteff2/toolkit/properties_calculator.py`

## 联系与反馈

如有问题或建议，欢迎提交 Issue 或 Pull Request。
