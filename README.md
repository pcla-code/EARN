# EARN
Epistemic Association Rule Networks 

**_rENA ARM (aka EARN): Integrating Epistemic Network Analysis and Association Rule Mining_**

## Overview 
This R script extends the [`rENA`](https://cran.r-project.org/src/contrib/Archive/rENA/) package by generating:
- **2-node and multi-node rule networks**
- **Difference networks** comparing rule structures across groups or conditions

To use EARN: 
  1) Load the dataset and specify which columns represent units, codes, and conversations as required by the rENA framework
  2) Define the association rule metrics (m1, m2, m3)
     - **m1** can be any metric except for addedValue (support, confidence, coverage, lift, cosine)
     - **m2** can be any symmetric metric (support, lift, cosine)
     - **m3** can be any asymmetric metric (confidence, coverage, addedValue)
  4) Define model parameters: rotation groups, moving window size, and threshold


Functions: 
- ena_arm_func(): Generates EARN plots for each group or condition

- ena_arm_func_diff(): Generates difference plots comparing EARN structures between pairs of groups
