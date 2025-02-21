// Stephen JD Kay, University of York, 21/02/25
// A short script to read in a file and prune it to only retain a smaller subset of branches.
// This could be utilised to trim down a full EICrecon file for a new user to look at. This avoids the potentially overwhelming number of branches normally stored
#include <string>

void TreePrune(TString infile=""){

  // If no input file provide as argument, promot for one
  if(infile == ""){
    cout << "Enter a filename to analyse: ";
    cin >> infile;
  }
  
  // Check input file exists, exit if not
  if(gSystem->AccessPathName(infile) == kTRUE){
    cerr << "!!!!! ERROR !!!!!" << endl << infile << " not found" << endl << "!!!!! ERROR !!!!!" << endl;
    exit(0);
  }

  TString ofile_name, ofile_tmp, BeamE;
  TObjArray *tmp_Name_arr;

  // Check the input file is a .root file as we would expect
  if(infile.Contains(".root") == false){
    cerr << "!!!!!!!!!!! - ERROR - !!!!!!!!!!!!" << endl;
    cerr << "Input files should be a root file!" << endl;
    cerr << "!!!!!!!!!!! - ERROR - !!!!!!!!!!!!" << endl;
    exit(1);
  }
  else{
    cout << "Opening and pruning " << infile << endl;
  }

  // If the input file name is a file path containing /, extract only the actual file name for further use. Assign the temporary name of the output to be the input, minus its .root extension
  if(infile.Contains("/")){
    tmp_Name_arr = infile.Tokenize("/");
    ofile_tmp = (((TObjString *)(tmp_Name_arr->At(tmp_Name_arr->GetLast())))->String()).ReplaceAll(".root","");
  }
  else{
    ofile_tmp = infile;
    ofile_tmp.ReplaceAll(".root","");
  }
  // Set output file name
  ofile_name = Form("%s_Pruned.root", ofile_tmp.Data());  

  // Open our full, unpruned file
  TFile *full_file = TFile::Open(infile);
  TTree* full_tree;
  // Get the events tree
  full_file->GetObject("events", full_tree);
  // Deactivate all branches
  full_tree->SetBranchStatus("*", 0);
  // Activate only the branches we want to keep, add more as needed via the line below
  // full_tree->SetBranchStatus("",1)
  full_tree->SetBranchStatus("MCParticles*",1);
  full_tree->SetBranchStatus("ReconstructedChargedParticles*",1);
  full_tree->SetBranchStatus("ReconstructedChargedParticleAssociations*",1);
  // Note that the wild card after the branch name is needed to ensure all leaves are retained adequately. You could only retain specific leaves if you wanted to though e.g. MCParticles.PDG
  
  // Set and open output file for the histograms
  TFile *ofile = TFile::Open(ofile_name,"RECREATE");
  // Clone the tree
  TTree* pruned_tree = full_tree->CloneTree(0);
  // Copy the branches
  pruned_tree->CopyEntries(full_tree);

  // Write the new tree to file
  ofile->Write();
  ofile->Close(); // Close output file
  full_file->Close(); // Close input file

  delete ofile;
  delete full_file;
  
}
