TFile* outFile = NULL;
TTree* outTree = NULL;

const double Me = 0.511e-3; // GeV
const double Mp = 0.938272; // GeV
const double crossing_angle = -0.025; // rad

void SetInputBranchAddresses();
void CreateOutputTree(TString outFileName);
void ResetVariables();
void CalculateElectronKinematics(double fEe, double fEh, TLorentzVector kf, float& xB, float& Q2, float& W, float& y, float& nu);
TLorentzVector GetHadronBeam(double fEh);

using namespace std;

int positive_eID;

float mc_p;
float mc_eta;
float mc_phi;

float track_p;
float track_eta;
float track_phi;

float mc_xB;
float mc_Q2;
float mc_W;
float mc_y;
float mc_nu;

float e_track_xB;
float e_track_Q2;
float e_track_W;
float e_track_y;
float e_track_nu;

float e_clust_xB;
float e_clust_Q2;
float e_clust_W;
float e_clust_y;
float e_clust_nu;

