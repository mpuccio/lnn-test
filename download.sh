mkdir -p data
mkdir -p mc
alien_cp alien:/alice/cern.ch/user/a/alihyperloop/outputs/0034/342361/76142/AnalysisResults.root file:data/AnalysisResults.root
alien_cp alien:/alice/cern.ch/user/a/alihyperloop/outputs/0034/342361/76142/AO2D.root file:data/AO2D.root
alien_cp alien:/alice/cern.ch/user/a/alihyperloop/outputs/0035/353588/81410/AnalysisResults.root file:mc/AnalysisResults.root
alien_cp alien:/alice/cern.ch/user/a/alihyperloop/outputs/0035/353588/81410/AO2D.root file:mc/AO2D.root
