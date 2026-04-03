# HyOffshore-Ireland

MATLAB project for offshore wind and integrated electricity-gas-hydrogen system studies focused on Ireland.

面向爱尔兰海上风电与电-气-氢综合能源系统研究的 MATLAB 项目。

## Overview

This repository contains:

- optimization models for electricity-gas and hydrogen-coupled systems
- data processing scripts for wind and gas-demand related studies
- plotting scripts for paper/report figures
- project inputs such as spreadsheets, geospatial files, and weather samples

仓库主要包含：

- 电-气耦合与电-气-氢耦合优化模型
- 风电、天然气需求等相关数据处理脚本
- 论文/报告所需图表生成脚本
- 表格、地理文件、气象样本等项目输入数据

## Main Entry Scripts

- `main.m`: main end-to-end workflow
- `mainPlot.m`: figure generation and post-processing
- `main_gasDemandEstimation.m`: gas demand estimation utilities
- `main_test.m`: benchmark / validation script
- `test.m`: small local test script for weather and network parsing

主要入口脚本：

- `main.m`：主流程脚本
- `mainPlot.m`：绘图与结果后处理
- `main_gasDemandEstimation.m`：天然气需求估计相关脚本
- `main_test.m`：基准工况/验证脚本
- `test.m`：本地小型测试脚本，用于天气数据和网络文件解析

## Core Code Folders

- `GEopf/`: electricity-gas optimization routines
- `HGEopf/`: hydrogen-gas-electricity optimization routines
- `matpower7.1/`: MATPOWER dependency bundled with the project
- `wind data/`: bulk wind input files
- `figs/`: generated figures

核心代码目录：

- `GEopf/`：电-气优化相关程序
- `HGEopf/`：电-气-氢优化相关程序
- `matpower7.1/`：项目内置的 MATPOWER 依赖
- `wind data/`：批量风场输入数据
- `figs/`：生成后的图表文件

## Data Layout

- `data/inputs/system/`: system spreadsheets and model inputs
- `data/inputs/demographics/`: population-density inputs
- `data/inputs/geospatial/`: shapefiles and coordinate tables
- `data/inputs/network/`: raw network files
- `data/inputs/weather_samples/`: weather sample files kept in the repo
- `data/optimization_snapshots/`: saved `.mat` experiment snapshots
- `data/checkpoints/`: runtime checkpoint files such as `stop1.mat` to `stop6.mat`

数据目录说明：

- `data/inputs/system/`：系统表格与模型输入数据
- `data/inputs/demographics/`：人口密度相关输入
- `data/inputs/geospatial/`：地理空间文件与坐标表
- `data/inputs/network/`：原始网络文件
- `data/inputs/weather_samples/`：保留在仓库中的气象样本数据
- `data/optimization_snapshots/`：保存的 `.mat` 优化结果快照
- `data/checkpoints/`：运行过程中生成的检查点文件，例如 `stop1.mat` 到 `stop6.mat`

## Documents And Logs

- `docs/reference/`: reference PDFs and supporting documents
- `docs/archive/`: archived exported documents
- `logs/`: text logs
- `archive/autosave/`: MATLAB autosave files

文档与日志目录：

- `docs/reference/`：参考 PDF 与补充文档
- `docs/archive/`：归档的导出文档
- `logs/`：文本日志
- `archive/autosave/`：MATLAB 自动保存文件

## Path Helpers

To make the reorganized structure easier to maintain, the project now includes:

- `projectRoot.m`: returns the repository root
- `projectPath.m`: builds absolute paths relative to the repository root

Main scripts use these helpers to find moved data files and checkpoint locations.

为了让整理后的目录更容易维护，项目新增了：

- `projectRoot.m`：返回仓库根目录
- `projectPath.m`：基于仓库根目录拼接绝对路径

主要脚本会通过这两个辅助函数定位移动后的数据文件和检查点文件。

## Requirements

Typical scripts expect the following MATLAB toolboxes / packages to be available:

- MATPOWER 7.1
- YALMIP
- a suitable solver such as Gurobi
- Mapping Toolbox style functions used by geospatial parts of the workflow
- NetCDF support for weather file processing

通常运行这些脚本需要以下 MATLAB 工具箱或外部依赖：

- MATPOWER 7.1
- YALMIP
- 合适的求解器，例如 Gurobi
- 地理空间处理相关函数所需能力，通常对应 Mapping Toolbox 一类功能
- 用于处理气象文件的 NetCDF 支持

## Notes

- Large local weather files may take time to process, especially when the repository is stored inside OneDrive.
- Generated checkpoint files in `data/checkpoints/` are ignored by Git to reduce working-tree noise.
- A more detailed folder summary is also available in `docs/README.md`.

补充说明：

- 本地气象文件较大，处理时可能比较慢，尤其是在仓库位于 OneDrive 目录下时。
- `data/checkpoints/` 中运行时生成的检查点文件已被 Git 忽略，以减少工作区噪音。
- 更细的目录说明可参考 `docs/README.md`。
