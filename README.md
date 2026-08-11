# Breast Cancer Risk Prediction Using Gene Expression Data

This repository contains the main analysis, modelling workflows, report materials, presentation slides, and Shiny application for our breast cancer risk prediction project.

## Final Deliverables

- Final academic report: [Breast Cancer Risk Prediction Using Gene Expression Data](<./Breast Cancer Risk Prediction Using Gene Expression Data.html>)
- Final presentation slides: [DATA3888 Final Presentation](./DATA3888_final_presentation.pptx)

## Repository Structure


### Main Work

- `Within_dataset_Normalization` contains the data cleaning and within-dataset normalization work.
- `Modelbuilding.qmd` and `Modelbuilding.R` are equivalent files for the binary prediction model training workflow.
- `Multiclass_PAM50_Calculator.qmd` contains the multiclass model training workflow.
- `binary_pipeline.R` is the final binary prediction pipeline developed after model comparison and feature processing.
- `multiclass_pipeline.R` is the final multiclass prediction pipeline developed after model comparison and feature processing.

### Report Submission

- `Report Submission/` contains the report submission materials.
- `Report Submission/README.md` explains how to run the report and provides more detail about the report section.

### Shiny App

- `ShinyApp/` contains the code for the interactive application developed for this project.

---

## 中文说明

本仓库包含了乳腺癌风险预测项目的主要分析过程、建模流程、报告材料、演讲幻灯片以及 Shiny 交互应用。

## 最终成果

- 最终学术报告：[Breast Cancer Risk Prediction Using Gene Expression Data](<./Breast Cancer Risk Prediction Using Gene Expression Data.html>)
- 最终演讲幻灯片：[DATA3888 Final Presentation](./DATA3888_final_presentation.pptx)

## 仓库结构


### Main Work

- `Within_dataset_Normalization` 包含数据清洗以及数据集内标准化的相关工作。
- `Modelbuilding.qmd` 和 `Modelbuilding.R` 内容相同，都是二分类预测模型训练部分的代码。
- `Multiclass_PAM50_Calculator.qmd` 包含多分类模型训练部分的代码。
- `binary_pipeline.R` 是在比较不同模型性能并完成特征处理后形成的最终二分类预测流程。
- `multiclass_pipeline.R` 是在比较不同模型性能并完成特征处理后形成的最终多分类预测流程。

### Report Submission

- `Report Submission/` 文件夹包含报告提交部分的相关内容。
- `Report Submission/README.md` 说明了如何运行报告，并提供了报告部分的详细介绍。

### Shiny App

- `ShinyApp/` 文件夹包含本项目中交互式应用的相关代码。
