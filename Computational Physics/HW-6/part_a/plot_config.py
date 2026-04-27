"""
matplotlib 中文字体配置
统一设置，避免各脚本中文字显示为方块
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# 查找可用的中文字体
_cjk_fonts = [f.name for f in fm.fontManager.ttflist
              if any(tag in f.name for tag in ["SimHei", "Microsoft YaHei", "SimSun", "Noto Sans CJK"])]

if _cjk_fonts:
    plt.rcParams["font.sans-serif"] = [_cjk_fonts[0]] + plt.rcParams.get("font.sans-serif", [])
    plt.rcParams["axes.unicode_minus"] = False
