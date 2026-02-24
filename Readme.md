# Project Title

SDM Workflow, code and results for the OpenPrioBio project.

## Description

Runs multi-year MaxEnt models for 10 forest bird species, at different spatial resolutions (Anthus pratensis @100m; Ciconia nigra @500m: Dendrocoptes medius @200m; Dryobates minor @1000m; Dryocopus martius @1000m; Phylloscopus sibillatrix @100m; Picus canus @1000m; Poecile montanus @200m; Poecile palustris @200m; and Sitta europaea @100m), for the whole of Germany. It uses both static and (time-)dynamic variables, with (occurrence and environmental) data from 2016 to 2024. The workflow runs the following steps:
  - 1a occurrence data preparation;
  - 1b environmental data extraction;
  - 2a spatial block definition;
  - 2b collinearity analysis;
  - 3a model fitting;
  - 3b model validation;
  - 3c model prediction;
  - 4a variable importance calculation;
  - 4b response curve plotting;
  - 4c response value extraction (for categorical variables);
  - 04d model statistics.

This workflow was adapted from: Wiedenroth, L. (2023). Wiedenroth_BirdWatch-SDM_2023 [Software]. GitHub. https://github.com/UP-macroecology/Wiedenroth_BirdWatch-SDM_2023.

## Getting Started

### Dependencies

* All dependencies are included in the code itself, including a routine to install non-installed packages.

## Authors

Pedro J. Leitão / pedro.leitao@lup-umwelt.de / GitHub steppebird

## License

This project is licensed under the GPL3 License - see the LICENSE.md file for details
