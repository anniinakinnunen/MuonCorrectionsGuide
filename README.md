# MuonCorrectionsGuide

There are misalignments in the CMS detector that make the reconstruction of muon momentum biased. The CMS reconstruction software does not fully correct these misalignments and additional corrections are needed to remove the bias. Correcting the misalignments is important when precision measurements are done using the muon momentum, because the bias in muon momentum will affect the results.

## The Muon Momentum Scale Corrections

The Muon Momentum Scale Corrections, also known as the Rochester Corrections, are available in the [MuonCorrectionsTool](https://github.com/cms-legacydata-analyses/MuonCorrectionsTool). The correction parameters have been extracted in a two step method. In the first step, initial corrections are obtained in bins of the charge of the muon and the η and ϕ coordinates of the muon track. The reconstruction bias in muon momentum depends on these variables. In the second step, the corrections are fine tuned using the mass of the Z boson.

The corrections for data and Monte Carlo (MC) are different since the MC events start with no biases but they can be induced during the reconstruction. Corrections have been extracted for both data and MC events.

In the MuonCorrectionsTool, the Run1 Rochester Corrections are added to two datasets as an example: a 2012 dataset and a MC dataset. A plot is created to check that the corrections were applied correctly. Creating the plot requires selections and the produced dataset contains only a part of the initial dataset. These selections can be skipped when the plot is not needed and a corrected version of the whole dataset is wanted as a result. Below you can find instructions on how to run the example code, how to apply the corrections to a different dataset and how to apply the corrections when you don't want to create the plot/make the selections. The official code for the Rochester Corrections can be found in the `RochesterCorrections` directory. The example code for applying the corrections is in the `Test` directory.

## Applying the corrections to data and MC

In the `Test` directory you can find `Analysis.C`, which is the example code for adding the corrections. The main function of `Analysis.C` is simply used for calling the `applyCorrections` function which takes as a parameter the name of the ROOT-file (without the .root-part), path to the ROOT-file and a boolean value of whether the file contains data (`true`) or MC (`false`).

_ADD TREENAME PARAMETER & CORRECT ALL PARAMETER_

```
void Analysis::main()
{
  // Data
  applyCorrections("Run2012BC_DoubleMuParked_Muons", "root://eospublic.cern.ch//eos/opendata/cms/derived-data/AOD2NanoAODOutreachTool/Run2012BC_DoubleMuParked_Muons.root", true);

  // MC
  applyCorrections("ZZTo2e2mu", "root://eospublic.cern.ch//eos/opendata/cms/upload/stefan/HiggsToFourLeptonsNanoAODOutreachAnalysis/ZZTo2e2mu.root", false);
}
```

The first thing applyCorrections does is create a TTree from the ROOT-file. Then variables for holding the values read from the tree are created and branch addresses are set so that the variables are populated when looping over events. New branches for the corrected values, an output file and a few variables needed for the corrections are also created.

```
int applyCorrections(string filename, string pathToFile, bool isData) {
  // Create TTree from ROOT file
  TFile *f1 = TFile::Open((pathToFile).c_str());
  TTree *DataTree = (TTree*)f1->Get("Events");
  
  //Variables to hold values read from the tree
  int maxmuon=1000;
  UInt_t nMuon = 0;
  Float_t Muon_pt[maxmuon];
  Float_t Muon_eta[maxmuon];
  Float_t Muon_phi[maxmuon];
  Float_t Muon_mass[maxmuon];
  Int_t Muon_charge[maxmuon];

  //Set addresses to make the tree populate the variables when reading an entry
  DataTree->SetBranchAddress("nMuon", &nMuon);
  DataTree->SetBranchAddress("Muon_pt", &Muon_pt);
  DataTree->SetBranchAddress("Muon_eta", &Muon_eta);
  DataTree->SetBranchAddress("Muon_phi", &Muon_phi);
  DataTree->SetBranchAddress("Muon_mass", &Muon_mass);
  DataTree->SetBranchAddress("Muon_charge", &Muon_charge);
```

Next, the events in the TTree are looped over and the corrections are applied to the muons. As mentioned earlier, in the MuonCorrectionsTool a plot is created to check the corrections. The invariant mass of μ<sup>+</sup>μ<sup>-</sup> is used in the plot, which is why the events are filtered to muon pairs with opposite charges and the invariant mass is computed.

```
  // Loop over events
  Int_t nEntries = (Int_t)DataTree->GetEntries();

  for (Int_t k=0; k<nEntries; k++) {
    DataTree->GetEntry(k);

    // Select events with exactly two muons
    if (nMuon == 2 ) {
      // Select events with two muons of opposite charge
      if (Muon_charge[0] != Muon_charge[1]) {

        // Compute invariant mass of the dimuon system
        Dimuon_mass = computeInvariantMass(Muon_pt[0], Muon_pt[1], Muon_eta[0], Muon_eta[1], Muon_phi[0], Muon_phi[1], Muon_mass[0], Muon_mass[1]);
        bDimuon_mass->Fill();
```

Inside the event loop there is another loop that loops over all the muons in an event and applies the corrections. The functions for applying the Rochester Corrections take as a parameter a TLorentzVector, which is a four-vector that describes the muons momentum and energy. A TLorentzVector is created for each muon using the muon's pt, eta, phi and mass. As mentioned earlier, the muon momentum scale corrections are different for data and MC and therefore there are separate functions for both: `momcor_data` and `momcor_mc`. These functions can be found in `rochcor2012wasym.cc` if you want to take a closer look at them.

The MuonCorrectionsTool plot is made in bins of eta of μ<sup>+</sup> and eta of μ<sup>-</sup>. This is why new branches are created and filled for those variables.

```
        // Loop over muons in event
        for (UInt_t i=0; i<nMuon; i++) {

          // Fill positive and negative muons eta branches
          if (Muon_charge[i] > 0) {
            Muon_eta_pos[i] = Muon_eta[i];
            bMuon_eta_pos->Fill();
          } else {
            Muon_eta_neg[i] = Muon_eta[i];
            bMuon_eta_neg->Fill();
          }

          // Create TLorentzVector
          TLorentzVector mu;
          mu.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]);

          // Apply the corrections
          if (isData) {
            rmcor.momcor_data(mu, Muon_charge[i], runopt, qter);
          } else {
            rmcor.momcor_mc(mu, Muon_charge[i], ntrk, qter);
          }
```

The corrected values are stored in the same TLorentzVectors after calling the correction functions. The values are then extracted from the TLorentzVecotrs and saved to the new branches. Also the corrected value of the invariant mass of μ<sup>+</sup>μ<sup>-</sup> is computed. The new TTree is filled and written to the output file.

```
          // Save corrected values
          Muon_pt_cor[i] = mu.Pt();
          bMuon_pt_cor->Fill();
          Muon_eta_cor[i] = mu.Eta();
          bMuon_eta_cor->Fill();
          Muon_phi_cor[i] = mu.Phi();
          bMuon_phi_cor->Fill();
          Muon_mass_cor[i] = mu.M();
          bMuon_mass_cor->Fill();
        }

        // Compute invariant mass of the corrected dimuon system
        Dimuon_mass_cor = computeInvariantMass(Muon_pt_cor[0], Muon_pt_cor[1], Muon_eta_cor[0], Muon_eta_cor[1], Muon_phi_cor[0], Muon_phi_cor[1], Muon_mass_cor[0], Muon_mass_cor[1]);
        bDimuon_mass_cor->Fill();

      }
    }
    //Fill the corrected values to the new tree
    DataTreeCor->Fill();
  }
  
  //Save the new tree
  DataTreeCor->Write();
```

## Applying the corrections to a different dataset

You can use the example code to apply the corrections to different datasets. However, a few changes needs to be made for the code to work correctly. The first thing that needs to be changed is of course the function call in the main function. Call `applyCorrections` using the parameters that correspond your dataset. Remember for the last parameter that `true` means your ROOT file contains data and `false` means it contains MC.

```
void Analysis::main()
{
  // Your dataset
  applyCorrections("NameOfYourRootFile", "PathToYourRootFile", true/false);
}
```

The second thing you need to do is check the names and data types of the branches in your dataset. For example, instead of the name `nMuon` you might have `numberOfMuons` and instead of data type `Muon_pt[nMuon]` you might have `vector<float> Muon_pt`. The correct name needs to be changed to the branch address and the data type needs to be corrected. If you have vector, you might need to change `Muon_pt[i]` to `Muon_pt->at(i)` or something similiar later in the code. Example:

```
  //Variables to hold values read from the tree
  int maxmuon=1000;
  UInt_t nMuon = 0;
  vector<float>* Muon_pt;

  //Set addresses to make the tree populate the variables when reading an entry
  DataTree->SetBranchAddress("numberOfMuons", &nMuon);
  DataTree->SetBranchAddress("Muon_pt", &Muon_pt);
```

## Correcting the dataset without making selections

The MuonCorrectionsTool plot 

Note that you don't need the `Muon_eta_pos` and `Muon_eta_neg` branches if you are not making the plot in MuonCorrectionsTool.
If you are not making the MuonCorrectionsTool plot and don't need to limit the muons to μ+μ- pairs, you can simply change the first `if` to `if (nMuon > 0)` and comment the second `if`. You should then also comment the invariant mass -lines. 
Again, you can comment the section about positive and negative muons' eta branches if you are creating the plot.
Again, you can comment the invariant mass section if needed.
