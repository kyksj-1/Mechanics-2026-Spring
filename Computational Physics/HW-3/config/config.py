import os
from pathlib import Path

# 获取项目根目录
ROOT_DIR = Path(__file__).resolve().parent.parent

# 定义常用目录
SRC_DIR = ROOT_DIR / "src"
SCRIPTS_DIR = ROOT_DIR / "scripts"
CONFIG_DIR = ROOT_DIR / "config"
OUTPUTS_DIR = ROOT_DIR / "outputs"
ASSET_DIR = ROOT_DIR / "asset"
DEBUG_DIR = ROOT_DIR / "drafts"

# 确保输出目录存在
OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
