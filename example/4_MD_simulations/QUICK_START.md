# PropertiesCalculator - 快速开始指南

## 超简单的输入格式

只需指定三个参数：
1. **溶剂名** (solvent)
2. **阴离子名** (anion)  
3. **离子对数量** (ion_count)

默认：5000 原子，锂离子（LI），298K

## 命令行方式（最简单）

```bash
# 基础用法
python -m byteff2.toolkit.properties_calculator \
    --solvent DMC \
    --anion PF6 \
    --ion-count 34 \
    --work-dir ./my_sim

# 指定计算属性
python -m byteff2.toolkit.properties_calculator \
    --solvent EC \
    --anion TFSI \
    --ion-count 20 \
    --calculate density,viscosity \
    --work-dir ./density_viscosity

# 自动计算离子数（基于 natoms）
python -m byteff2.toolkit.properties_calculator \
    --solvent H2O \
    --anion PF6 \
    --natoms 8000 \
    --work-dir ./water_system
```

## Python 方式（工作流集成）

```python
from byteff2.toolkit.properties_calculator import PropertiesCalculator

# 最简单的用法
calc = PropertiesCalculator(
    solvent="DMC",
    anion="PF6",
    ion_count=34,
    base_dir="./my_sim"
)

# 计算所有属性
results = calc.calculate()

# 或计算特定属性
results = calc.calculate(properties=["density", "viscosity"])

# 访问结果
print(f"Density: {results['density']['density']} g/mL")
print(f"Viscosity: {results['viscosity']} cP")
```

## 参数化研究示例

```python
# 在离子浓度上做参数化

ion_concentrations = [10, 20, 30, 40, 50]

for ion_count in ion_concentrations:
    calc = PropertiesCalculator(
        solvent="DMC",
        anion="PF6",
        ion_count=ion_count,
        base_dir=f"./parametric_li{ion_count}"
    )
    
    results = calc.calculate()
    
    # 你的后续处理...
    print(f"Li count: {ion_count}, Density: {results['density']['density']}")
```

## 支持的溶剂和阴离子

### 溶剂（Solvents）
- `DMC` - 二甲基碳酸酯
- `EC` - 碳酸乙烯  
- `EMC` - 乙基甲基碳酸酯
- `DEC` - 二乙基碳酸酯
- `PC` - 碳酸丙烯
- `H2O` - 水

### 阴离子（Anions）
- `PF6` - 六氟磷酸根
- `BF4` - 四氟硼酸根
- `ClO4` - 高氯酸根
- `TFSI` - 双三氟甲基磺酰亚胺
- `OTf` - 三氟甲烷磺酸根

### 阳离子（Cation）
- `LI` - 锂离子（固定，默认）

## 输出位置

所有结果都在 `work_dir` 下：

```
my_sim/
├── summary.json              # 汇总信息
├── density_config.json       # 配置文件
├── density_results.json      # 密度结果
├── transport_config.json     # 配置文件
├── dielectric_results.json   # 介电常数结果
├── ...
└── *_results/                # 详细 MD 数据
```

## 常用命令速查

```bash
# 单个属性
--calculate density
--calculate viscosity
--calculate dielectric

# 多个属性
--calculate density,viscosity
--calculate density,conductivity,viscosity,dielectric

# 自定义温度
--temperature 313.0

# 自定义系统大小
--natoms 8000

# 帮助信息
python -m byteff2.toolkit.properties_calculator --help
```

## 故障排查

| 问题 | 解决方案 |
|------|---------|
| 溶剂不支持 | 检查拼写，参考支持的溶剂列表 |
| 阴离子不支持 | 检查拼写，参考支持的阴离子列表 |
| 结果文件缺失 | 检查 MD 模拟是否完成 |
| 内存不足 | 减小 `natoms` 或 `ion_count` |

---

**提示**: 要添加新的溶剂或阴离子，更新 `PropertiesCalculator` 类中的 `SOLVENT_DATABASE` 或 `ANION_DATABASE` 字典，加入 SMILES 字符串。
