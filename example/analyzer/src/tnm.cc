// -----------------------------------------------------------------------------
// Some TNM utilities
// Created: A long time ago in a galaxy far far away HBP
// -----------------------------------------------------------------------------
#include <iostream>
#include <sstream>
#include <map>
#ifdef PROJECT_NAME
#include "PhysicsTools/TheNtupleMaker/interface/tnm.h"
#else
#include "tnm.h"
#endif
using namespace std;
// -----------------------------------------------------------------------------

void error(std::string message)
{
  std::cout << "** error ** " << message << std::endl;
  exit(0);
}
///
std::string strip(std::string line)
{
  int l = line.size();
  if ( l == 0 ) return std::string("");
  int n = 0;
  while (((line[n] == 0)    ||
      (line[n] == ' ' ) ||
      (line[n] == '\n') ||
      (line[n] == '\t')) && n < l) n++;

  int m = l-1;
  while (((line[m] == 0)    ||
      (line[m] == ' ')  ||
      (line[m] == '\n') ||
      (line[m] == '\t')) && m > 0) m--;
  return line.substr(n,m-n+1);
}

///
std::vector<std::string> split(std::string str)
{
  vector<string> vstr;
  std::istringstream stream(str);
  while ( stream )
    {
      std::string str;
      stream >> str;
      if ( stream ) vstr.push_back(str);
    }
  return vstr;
}

///
std::string change(std::string str, std::string oldstr, std::string newstr)
{
  return std::string(TString(str).ReplaceAll(oldstr, newstr).Data());
}

///
std::string nameonly(std::string filename)
{
  int i = filename.rfind("/");
  int j = filename.rfind(".");
  if ( j < 0 ) j = filename.size();
  return filename.substr(i+1,j-i-1);
}

///
std::string shell(std::string cmd)
{
  FILE* f = popen(cmd.c_str(),"r");
  int buffsize=8192;
  char s[8192];
  int n = fread(s,1,buffsize,f);
  pclose(f);
  std::string result = strip(std::string(s).substr(0,n));
  return result;
}


outputFile::outputFile(std::string filename)
  : filename_(filename),
    file_(new TFile(filename_.c_str(), "recreate")),
    tree(0),
    b_weight_(0),
    entry_(0),
    SAVECOUNT_(50000),
    ev_(0)
{
  file_->cd();
  hist_ = new TH1F("counts", "", 1,0,1);
  hist_->SetCanExtend(1);
  hist_->SetStats(0);
}

outputFile::outputFile(std::string filename, 
		       eventBuffer& ev, 
		       int savecount) 
  : filename_(filename),
    file_(new TFile(filename.c_str(), "recreate")),
    tree(ev.input ? 
	 ev.input->tree()->CloneTree(0) : 0),
    b_weight_(tree ? 
	      tree->Branch("eventWeight", &weight_, "eventWeight/D") : 0),
    entry_(0),
    SAVECOUNT_(savecount),
    ev_(&ev)
{
  if ( tree == 0 )
    error("outputFile - tree pointer is NULL");

  std::cout << "events will be skimmed to file "
	    << filename_ << std::endl;
  file_->cd();
  hist_ = new TH1F("counts", "", 1,0,1);
  hist_->SetCanExtend(1);
  hist_->SetStats(0);
}


void outputFile::write(double weight)
{
  if ( tree == 0 ) return;
  if ( ev_ ) ev_->saveObjects();

  weight_ = weight;
  file_   = tree->GetCurrentFile();
  file_->cd();
  tree->Fill();

  entry_++;
  if ( entry_ % SAVECOUNT_ == 0 )
    tree->AutoSave("SaveSelf");
}

void outputFile::count(std::string cond, double w)
{
  hist_->Fill(cond.c_str(), w);
}

///
void outputFile::close()
{
  std::cout << "==> histograms saved to file " << filename_ << std::endl;
  if ( tree )
    {
      std::cout << "==> events skimmed to file " << filename_ << std::endl;
      file_ = tree->GetCurrentFile();
    }
  file_->cd();
  file_->Write("", TObject::kOverwrite);
  file_->ls();
  file_->Close();
}

commandLine::commandLine()
{
  int argc = gApplication->Argc();
  char** argv = gApplication->Argv();
  decode(argc, argv);
}

commandLine::commandLine(int argc, char** argv)
{
  decode(argc, argv);
}

void
commandLine::decode(int argc, char** argv)
{
  if ( argc <= 0 )
    {
      argc = gApplication->Argc();
      argv = gApplication->Argv();
    }

  progname = nameonly(std::string(argv[0]));
  if ( progname == "Python" || progname == "python" )
    progname = string("analyzer");

  // 1st (optional) argument
  if ( argc > 1 )
    filelist = std::string(argv[1]);
  else
    filelist = std::string("filelist.txt");

  // 2nd (optional) command line argument
  if ( argc > 2 ) 
    outputfilename = std::string(argv[2]);
  else
    outputfilename = progname + std::string("_histograms");

  // Make sure extension is ".root"
  std::string name = outputfilename;
  if ( name.substr(name.size()-5, 5) != std::string(".root") )
    outputfilename += std::string(".root");
}

/// Read ntuple filenames from file list
std::vector<std::string> fileNames(std::string filelist)
{
  std::ifstream stream(filelist.c_str());
  if ( !stream.good() ) error("unable to open file: " + filelist);

  // Get list of ntuple files to be processed

  std::vector<std::string> v;
  std::string filename;
  while ( stream >> filename )
    if ( strip(filename) != "" ) v.push_back(filename);
  return v;
}

// physics utils ==============================================================

///
double deltaPhi(double phi1, double phi2)
{
  double deltaphi = phi2 - phi1;
  if ( fabs(deltaphi) > M_PI ) deltaphi = 2 * M_PI - fabs(deltaphi);
  return deltaphi;
}

///
double deltaR(double eta1, double phi1, double eta2, double phi2)
{
  double deltaeta = eta1 - eta2;
  double deltaphi = deltaPhi(phi1, phi2);
  return sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
}

ptThing::ptThing() {}
ptThing::ptThing(int index_, int id_,
		 double pt_, double eta_, double phi_, std::string name_)
  : index(index_), id(id_),
    pt(pt_), eta(eta_), phi(phi_), name(name_) {}		  
ptThing::~ptThing() {}

/// Copy constructor.
ptThing::ptThing(const ptThing& rhs)
{
  *this = rhs; // this will call assignment operator
}

/// Assignment. 
ptThing& ptThing::operator=(const ptThing& rhs)
{
  index = rhs.index;
  id    = rhs.id;
  pt    = rhs.pt;
  eta   = rhs.eta;
  phi   = rhs.phi;
  name  = rhs.name;
  var   = rhs.var;
  return *this;
}

/** Find $|Delta R = \sqrt{\Delta\phi^2+\Delta\eta^2}$ between this
    ptThing and the given.
*/
double ptThing::deltaR(ptThing& thing)
{
  double deta = eta - thing.eta;
  double dphi = phi - thing.phi;
  
  // Make sure acute
  if ( fabs(dphi) > TMath::Pi() ) dphi = TMath::TwoPi() - fabs(dphi);
  return sqrt(deta*deta+dphi*dphi);
}

/// Compare direction of this ptThing with another using deltaR.
bool ptThing::matches(ptThing& thing, double drcut)
{
  return deltaR(thing) < drcut;
}

// ///
// std::vector<matchedPair>
// deltaR(std::vector<ptThing>& v1, std::vector<ptThing>& v2)
// {
//   if ( v1.size() == 0 || v2.size() == 0 ) return std::vector<matchedPair>();

//   std::vector<matchedPair> mp(v1.size(), matchedPair());
//   std::vector<matchedPair> vp(v2.size(), matchedPair());

//   for(unsigned i=0; i < v1.size(); i++)
//     {
//       for(unsigned j=0; j < v2.size(); j++)
//         {
//           vp[j].first = i;
//           vp[j].second = j;
//           vp[j].distance = v1[i].deltaR(v2[j]);
//         }
//       std::sort(vp.begin(), vp.end());
//       mp[i].first = i;
//       mp[i].second = vp[0].second;
//       mp[i].distance = vp[0].distance;
//     }
//   return mp;
// }

void setStyle()
{
  gROOT->SetStyle("Pub");
  gStyle->SetPalette(1);

  TStyle* style = gROOT->GetStyle("Pub");
  style->SetFrameBorderMode(0);
  style->SetCanvasBorderSize(0);
  style->SetCanvasBorderMode(0);
  style->SetCanvasColor(0);
  style->SetPadBorderMode(0);
  style->SetPadColor(0);

  // Margins

  style->SetPadTopMargin(0.05);
  style->SetPadBottomMargin(0.16);
  style->SetPadLeftMargin(0.20);
  style->SetPadRightMargin(0.10);

  // For the Global title:

  style->SetOptTitle(0);
  style->SetTitleFont(42);
  style->SetTitleColor(1);
  style->SetTitleTextColor(1);
  style->SetTitleFillColor(10);
  style->SetTitleFontSize(0.05);

  // For the axis titles:

  style->SetTitleColor(1, "XYZ");
  style->SetTitleFont(42, "XYZ");
  style->SetTitleSize(0.05, "XYZ");
  style->SetTitleXOffset(0.9);
  style->SetTitleYOffset(1.25);

  // For the axis labels:

  style->SetLabelColor(1, "XYZ");
  style->SetLabelFont(42, "XYZ");
  style->SetLabelOffset(0.007, "XYZ");
  style->SetLabelSize(0.05, "XYZ");

  // For the axis:

  style->SetAxisColor(1, "XYZ");
  style->SetStripDecimals(kTRUE);
  style->SetTickLength(0.03, "XYZ");
  style->SetNdivisions(505, "XYZ");
  style->SetPadTickX(1);  
  style->SetPadTickY(1);

  style->cd();
}

namespace {
  static bool firstName=true;
  static map<int, string>  namemap;
}
string particleName(int pdgid)
{
  if ( firstName )
    {
      firstName = false;
      namemap[1]	= "d";
      namemap[-1]	= "d~";
      namemap[2]	= "u";
      namemap[-2]	= "u~";
      namemap[3]	= "s";
      namemap[-3]	= "s~";
      namemap[4]	= "c";
      namemap[-4]	= "c~";
      namemap[5]	= "b";
      namemap[-5]	= "b~";
      namemap[6]	= "t";
      namemap[-6]	= "t~";
      namemap[7]	= "b'";
      namemap[-7]	= "b'~";
      namemap[8]	= "t'";
      namemap[-8]	= "t'~";
      namemap[11]	= "e^-";
      namemap[-11]	= "e^+";
      namemap[12]	= "nu_e";
      namemap[-12]	= "nu_e~";
      namemap[13]	= "mu^-";
      namemap[-13]	= "mu^+";
      namemap[14]	= "nu_mu";
      namemap[-14]	= "nu_mu~";
      namemap[15]	= "tau^-";
      namemap[-15]	= "tau^+";
      namemap[16]	= "nu_tau";
      namemap[-16]	= "nu_tau~";
      namemap[17]	= "tau'^-";
      namemap[-17]	= "tau'^+";
      namemap[18]	= "nu_tau'";
      namemap[-18]	= "nu_tau'~";
      namemap[21]	= "g";
      namemap[22]	= "gamma";
      namemap[23]	= "Z^0";
      namemap[24]	= "W^+";
      namemap[-24]	= "W^-";
      namemap[25]	= "H_1^0";
      namemap[32]	= "Z_2^0";
      namemap[33]	= "Z_3^0";
      namemap[34]	= "W_2^+";
      namemap[-34]	= "W_2^-";
      namemap[35]	= "H_2^0";
      namemap[36]	= "H_3^0";
      namemap[37]	= "H^+";
      namemap[-37]	= "H^-";
      namemap[39]	= "G";
      namemap[41]	= "R^0";
      namemap[-41]	= "R~^0";
      namemap[42]	= "LQ_c";
      namemap[-42]	= "LQ_c~";
      namemap[51]	= "H_L^0";
      namemap[52]	= "H_1^++";
      namemap[-52]	= "H_1^--";
      namemap[53]	= "H_2^+";
      namemap[-53]	= "H_2^-";
      namemap[54]	= "H_2^++";
      namemap[-54]	= "H_2^--";
      namemap[55]	= "H_4^0";
      namemap[-55]	= "H_4~^0";
      namemap[81]	= "generator-specific+81";
      namemap[-81]	= "generator-specific-81";
      namemap[82]	= "generator-specific+82";
      namemap[-82]	= "generator-specific-82";
      namemap[83]	= "generator-specific+83";
      namemap[-83]	= "generator-specific-83";
      namemap[84]	= "generator-specific+84";
      namemap[-84]	= "generator-specific-84";
      namemap[85]	= "generator-specific+85";
      namemap[-85]	= "generator-specific-85";
      namemap[86]	= "generator-specific+86";
      namemap[-86]	= "generator-specific-86";
      namemap[87]	= "generator-specific+87";
      namemap[-87]	= "generator-specific-87";
      namemap[88]	= "generator-specific+88";
      namemap[-88]	= "generator-specific-88";
      namemap[89]	= "generator-specific+89";
      namemap[-89]	= "generator-specific-89";
      namemap[90]	= "generator-specific+90";
      namemap[-90]	= "generator-specific-90";
      namemap[91]	= "generator-specific+91";
      namemap[-91]	= "generator-specific-91";
      namemap[92]	= "generator-specific+92";
      namemap[-92]	= "generator-specific-92";
      namemap[93]	= "generator-specific+93";
      namemap[-93]	= "generator-specific-93";
      namemap[94]	= "generator-specific+94";
      namemap[-94]	= "generator-specific-94";
      namemap[95]	= "generator-specific+95";
      namemap[-95]	= "generator-specific-95";
      namemap[96]	= "generator-specific+96";
      namemap[-96]	= "generator-specific-96";
      namemap[97]	= "generator-specific+97";
      namemap[-97]	= "generator-specific-97";
      namemap[98]	= "generator-specific+98";
      namemap[-98]	= "generator-specific-98";
      namemap[99]	= "generator-specific+99";
      namemap[-99]	= "generator-specific-99";
      namemap[100]	= "generator-specific+100";
      namemap[-100]	= "generator-specific-100";
      namemap[101]	= "geantino";
      namemap[102]	= "charged-geantino";
      namemap[110]	= "reggeon";
      namemap[111]	= "pi^0";
      namemap[113]	= "rho(770)^0";
      namemap[115]	= "a_2(1320)^0";
      namemap[117]	= "rho_3(1690)^0";
      namemap[119]	= "a_4(2040)^0";
      namemap[130]	= "K_L^0";
      namemap[211]	= "pi^+";
      namemap[-211]	= "pi^-";
      namemap[213]	= "rho(770)^+";
      namemap[-213]	= "rho(770)^-";
      namemap[215]	= "a_2(1320)^+";
      namemap[-215]	= "a_2(1320)^-";
      namemap[217]	= "rho_3(1690)^+";
      namemap[-217]	= "rho_3(1690)^-";
      namemap[219]	= "a_4(2040)^+";
      namemap[-219]	= "a_4(2040)^-";
      namemap[221]	= "eta";
      namemap[223]	= "omega(782)";
      namemap[225]	= "f_2(1270)";
      namemap[227]	= "omega_3(1670)";
      namemap[229]	= "f_4(2050)";
      namemap[310]	= "K_S^0";
      namemap[311]	= "K^0";
      namemap[-311]	= "K~^0";
      namemap[313]	= "K*(892)^0";
      namemap[-313]	= "K*(892)~^0";
      namemap[315]	= "K*_2(1430)^0";
      namemap[-315]	= "K*_2(1430)~^0";
      namemap[317]	= "K*_3(1780)^0";
      namemap[-317]	= "K*_3(1780)~^0";
      namemap[319]	= "K*_4(2045)^0";
      namemap[-319]	= "K*_4(2045)~^0";
      namemap[321]	= "K^+";
      namemap[-321]	= "K^-";
      namemap[323]	= "K*(892)^+";
      namemap[-323]	= "K*(892)^-";
      namemap[325]	= "K*_2(1430)^+";
      namemap[-325]	= "K*_2(1430)^-";
      namemap[327]	= "K*_3(1780)^+";
      namemap[-327]	= "K*_3(1780)^-";
      namemap[329]	= "K*_4(2045)^+";
      namemap[-329]	= "K*_4(2045)^-";
      namemap[331]	= "eta'(958)";
      namemap[333]	= "phi(1020)";
      namemap[335]	= "f'_2(1525)";
      namemap[337]	= "phi_3(1850)";
      namemap[411]	= "D^+";
      namemap[-411]	= "D^-";
      namemap[413]	= "D*(2010)^+";
      namemap[-413]	= "D*(2010)^-";
      namemap[415]	= "D*_2(2460)^+";
      namemap[-415]	= "D*_2(2460)^-";
      namemap[421]	= "D^0";
      namemap[-421]	= "D~^0";
      namemap[423]	= "D*(2007)^0";
      namemap[-423]	= "D*(2007)~^0";
      namemap[425]	= "D*_2(2460)^0";
      namemap[-425]	= "D*_2(2460)~^0";
      namemap[431]	= "D_s^+";
      namemap[-431]	= "D_s^-";
      namemap[433]	= "D*_s^+";
      namemap[-433]	= "D*_s^-";
      namemap[435]	= "D*_s2(2573)^+";
      namemap[-435]	= "D*_s2(2573)^-";
      namemap[441]	= "eta_c(1S)";
      namemap[443]	= "J/psi(1S)";
      namemap[445]	= "chi_c2(1P)";
      namemap[511]	= "B^0";
      namemap[-511]	= "B~^0";
      namemap[513]	= "B*^0";
      namemap[-513]	= "B*~^0";
      namemap[515]	= "B*_2^0";
      namemap[-515]	= "B*_2~^0";
      namemap[521]	= "B^+";
      namemap[-521]	= "B^-";
      namemap[523]	= "B*^+";
      namemap[-523]	= "B*^-";
      namemap[525]	= "B*_2^+";
      namemap[-525]	= "B*_2^-";
      namemap[531]	= "B_s^0";
      namemap[-531]	= "B_s~^0";
      namemap[533]	= "B*_s^0";
      namemap[-533]	= "B*_s~^0";
      namemap[535]	= "B*_s2^0";
      namemap[-535]	= "B*_s2~^0";
      namemap[541]	= "B_c^+";
      namemap[-541]	= "B_c^-";
      namemap[543]	= "B*_c^+";
      namemap[-543]	= "B*_c^-";
      namemap[545]	= "B*_c2^+";
      namemap[-545]	= "B*_c2^-";
      namemap[551]	= "eta_b(1S)";
      namemap[553]	= "Upsilon(1S)";
      namemap[555]	= "chi_b2(1P)";
      namemap[557]	= "Upsilon_3(1D)";
      namemap[611]	= "T^+";
      namemap[-611]	= "T^-";
      namemap[613]	= "T*^+";
      namemap[-613]	= "T*^-";
      namemap[621]	= "T^0";
      namemap[-621]	= "T~^0";
      namemap[623]	= "T*^0";
      namemap[-623]	= "T*~^0";
      namemap[631]	= "T_s^+";
      namemap[-631]	= "T_s^-";
      namemap[633]	= "T*_s^+";
      namemap[-633]	= "T*_s^-";
      namemap[641]	= "T_c^0";
      namemap[-641]	= "T_c~^0";
      namemap[643]	= "T*_c^0";
      namemap[-643]	= "T*_c~^0";
      namemap[651]	= "T_b^+";
      namemap[-651]	= "T_b^-";
      namemap[653]	= "T*_b^+";
      namemap[-653]	= "T*_b^-";
      namemap[661]	= "eta_t";
      namemap[663]	= "theta";
      namemap[711]	= "L^0";
      namemap[-711]	= "L~^0";
      namemap[713]	= "L*^0";
      namemap[-713]	= "L*~^0";
      namemap[721]	= "L^-";
      namemap[-721]	= "L^+";
      namemap[723]	= "L*^-";
      namemap[-723]	= "L*^+";
      namemap[731]	= "L_s^0";
      namemap[-731]	= "L_s~^0";
      namemap[733]	= "L*_s^0";
      namemap[-733]	= "L*_s~^0";
      namemap[741]	= "L_c^-";
      namemap[-741]	= "L_c^+";
      namemap[743]	= "L*_c^-";
      namemap[-743]	= "L*_c^+";
      namemap[751]	= "L_b^0";
      namemap[-751]	= "L_b~^0";
      namemap[753]	= "L*_b^0";
      namemap[-753]	= "L*_b~^0";
      namemap[761]	= "L_t^-";
      namemap[-761]	= "L_t^+";
      namemap[763]	= "L*_t^-";
      namemap[-763]	= "L*_t^+";
      namemap[771]	= "eta_l";
      namemap[773]	= "theta_l";
      namemap[811]	= "H^+";
      namemap[-811]	= "H^-";
      namemap[813]	= "H*^+";
      namemap[-813]	= "H*^-";
      namemap[821]	= "H^0";
      namemap[-821]	= "H~^0";
      namemap[823]	= "H*^0";
      namemap[-823]	= "H*~^0";
      namemap[831]	= "H_s^+";
      namemap[-831]	= "H_s^-";
      namemap[833]	= "H*_s^+";
      namemap[-833]	= "H*_s^-";
      namemap[841]	= "H_c^0";
      namemap[-841]	= "H_c~^0";
      namemap[843]	= "H*_c^0";
      namemap[-843]	= "H*_c~^0";
      namemap[851]	= "H_b^+";
      namemap[-851]	= "H_b^-";
      namemap[853]	= "H*_b^+";
      namemap[-853]	= "H*_b^-";
      namemap[861]	= "H_t^0";
      namemap[-861]	= "H_t~^0";
      namemap[863]	= "H*_t^0";
      namemap[-863]	= "H*_t~^0";
      namemap[871]	= "H_l^+";
      namemap[-871]	= "H_l^-";
      namemap[873]	= "H*_l^+";
      namemap[-873]	= "H*_l^-";
      namemap[881]	= "eta_h";
      namemap[883]	= "theta_H";
      namemap[990]	= "pomeron";
      namemap[1103]	= "dd_1";
      namemap[-1103]	= "dd_1~";
      namemap[1112]	= "Delta(1620)^-";
      namemap[1114]	= "Delta^-";
      namemap[-1114]	= "Delta~^+";
      namemap[1116]	= "Delta(1905)^-";
      namemap[1118]	= "Delta(1950)^-";
      namemap[1212]	= "Delta(1620)^0";
      namemap[1214]	= "N(1520)^0";
      namemap[1216]	= "Delta(1905)^0";
      namemap[1218]	= "N(2190)^0";
      namemap[2101]	= "ud_0";
      namemap[-2101]	= "ud_0~";
      namemap[2103]	= "ud_1";
      namemap[-2103]	= "ud_1~";
      namemap[2112]	= "n^0";
      namemap[-2112]	= "n~^0";
      namemap[2114]	= "Delta^0";
      namemap[-2114]	= "Delta~^0";
      namemap[2116]	= "N(1675)^0";
      namemap[2118]	= "Delta(1950)^0";
      namemap[2122]	= "Delta(1620)^+";
      namemap[2124]	= "N(1520)^+";
      namemap[2126]	= "Delta(1905)^+";
      namemap[2128]	= "N(2190)^+";
      namemap[2203]	= "uu_1";
      namemap[-2203]	= "uu_1~";
      namemap[2212]	= "p^+";
      namemap[-2212]	= "p~^-";
      namemap[2214]	= "Delta^+";
      namemap[-2214]	= "Delta~^-";
      namemap[2216]	= "N(1675)^+";
      namemap[2218]	= "Delta(1950)^+";
      namemap[2222]	= "Delta(1620)^++";
      namemap[2224]	= "Delta^++";
      namemap[-2224]	= "Delta~^--";
      namemap[2226]	= "Delta(1905)^++";
      namemap[2228]	= "Delta(1950)^++";
      namemap[3101]	= "sd_0";
      namemap[-3101]	= "sd_0~";
      namemap[3103]	= "sd_1";
      namemap[-3103]	= "sd_1~";
      namemap[3112]	= "Sigma^-";
      namemap[-3112]	= "Sigma~^+";
      namemap[3114]	= "Sigma*^-";
      namemap[-3114]	= "Sigma*~^+";
      namemap[3116]	= "Sigma(1775)^-";
      namemap[-3116]	= "Sigma~(1775)^-";
      namemap[3118]	= "Sigma(2030)^-";
      namemap[-3118]	= "Sigma~(2030)^-";
      namemap[3122]	= "Lambda^0";
      namemap[-3122]	= "Lambda~^0";
      namemap[3124]	= "Lambda(1520)^0";
      namemap[-3124]	= "Lambda~(1520)^0";
      namemap[3126]	= "Lambda(1820)^0";
      namemap[-3126]	= "Lambda~(1820)^0";
      namemap[3128]	= "Lambda(2100)^0";
      namemap[-3128]	= "Lambda~(2100)^0";
      namemap[3201]	= "su_0";
      namemap[-3201]	= "su_0~";
      namemap[3203]	= "su_1";
      namemap[-3203]	= "su_1~";
      namemap[3212]	= "Sigma^0";
      namemap[-3212]	= "Sigma~^0";
      namemap[3214]	= "Sigma*^0";
      namemap[-3214]	= "Sigma*~^0";
      namemap[3216]	= "Sigma(1775)^0";
      namemap[-3216]	= "Sigma~(1775)^0";
      namemap[3218]	= "Sigma(2030)^0";
      namemap[-3218]	= "Sigma~(2030)^0";
      namemap[3222]	= "Sigma^+";
      namemap[-3222]	= "Sigma~^-";
      namemap[3224]	= "Sigma*^+";
      namemap[-3224]	= "Sigma*~^-";
      namemap[3226]	= "Sigma(1775)^+";
      namemap[-3226]	= "Sigma~(1775)^+";
      namemap[3228]	= "Sigma(2030)^+";
      namemap[-3228]	= "Sigma~(2030)^+";
      namemap[3303]	= "ss_1";
      namemap[-3303]	= "ss_1~";
      namemap[3312]	= "Xi^-";
      namemap[-3312]	= "Xi~^+";
      namemap[3314]	= "Xi*^-";
      namemap[-3314]	= "Xi*~^+";
      namemap[3322]	= "Xi^0";
      namemap[-3322]	= "Xi~^0";
      namemap[3324]	= "Xi*^0";
      namemap[-3324]	= "Xi*~^0";
      namemap[3334]	= "Omega^-";
      namemap[-3334]	= "Omega~^+";
      namemap[4101]	= "cd_0";
      namemap[-4101]	= "cd_0~";
      namemap[4103]	= "cd_1";
      namemap[-4103]	= "cd_1~";
      namemap[4112]	= "Sigma_c^0";
      namemap[-4112]	= "Sigma_c~^0";
      namemap[4114]	= "Sigma*_c^0";
      namemap[-4114]	= "Sigma*_c~^0";
      namemap[4122]	= "Lambda_c^+";
      namemap[-4122]	= "Lambda_c~^-";
      namemap[4132]	= "Xi_c^0";
      namemap[-4132]	= "Xi_c~^0";
      namemap[4201]	= "cu_0";
      namemap[-4201]	= "cu_0~";
      namemap[4203]	= "cu_1";
      namemap[-4203]	= "cu_1~";
      namemap[4212]	= "Sigma_c^+";
      namemap[-4212]	= "Sigma_c~^-";
      namemap[4214]	= "Sigma*_c^+";
      namemap[-4214]	= "Sigma*_c~^-";
      namemap[4222]	= "Sigma_c^++";
      namemap[-4222]	= "Sigma_c~^--";
      namemap[4224]	= "Sigma*_c^++";
      namemap[-4224]	= "Sigma*_c~^--";
      namemap[4232]	= "Xi_c^+";
      namemap[-4232]	= "Xi_c~^-";
      namemap[4301]	= "cs_0";
      namemap[-4301]	= "cs_0~";
      namemap[4303]	= "cs_1";
      namemap[-4303]	= "cs_1~";
      namemap[4312]	= "Xi'_c^0";
      namemap[-4312]	= "Xi'_c~^0";
      namemap[4314]	= "Xi*_c^0";
      namemap[-4314]	= "Xi*_c~^0";
      namemap[4322]	= "Xi'_c^+";
      namemap[-4322]	= "Xi'_c~^-";
      namemap[4324]	= "Xi*_c^+";
      namemap[-4324]	= "Xi*_c~^-";
      namemap[4332]	= "Omega_c^0";
      namemap[-4332]	= "Omega_c~^0";
      namemap[4334]	= "Omega*_c^0";
      namemap[-4334]	= "Omega*_c~^0";
      namemap[4403]	= "cc_1";
      namemap[-4403]	= "cc_1~";
      namemap[4412]	= "Xi_cc^+";
      namemap[-4412]	= "Xi_cc~^-";
      namemap[4414]	= "Xi*_cc^+";
      namemap[-4414]	= "Xi*_cc~^-";
      namemap[4422]	= "Xi_cc^++";
      namemap[-4422]	= "Xi_cc~^--";
      namemap[4424]	= "Xi*_cc^++";
      namemap[-4424]	= "Xi*_cc~^--";
      namemap[4432]	= "Omega_cc^+";
      namemap[-4432]	= "Omega_cc~^-";
      namemap[4434]	= "Omega*_cc^+";
      namemap[-4434]	= "Omega*_cc~^-";
      namemap[4444]	= "Omega*_ccc^++";
      namemap[-4444]	= "Omega*_ccc~^--";
      namemap[5101]	= "bd_0";
      namemap[-5101]	= "bd_0~";
      namemap[5103]	= "bd_1";
      namemap[-5103]	= "bd_1~";
      namemap[5112]	= "Sigma_b^-";
      namemap[-5112]	= "Sigma_b~^+";
      namemap[5114]	= "Sigma*_b^-";
      namemap[-5114]	= "Sigma*_b~^+";
      namemap[5122]	= "Lambda_b^0";
      namemap[-5122]	= "Lambda_b~^0";
      namemap[5132]	= "Xi_b^-";
      namemap[-5132]	= "Xi_b~^+";
      namemap[5142]	= "Xi_bc^0";
      namemap[-5142]	= "Xi_bc~^0";
      namemap[5201]	= "bu_0";
      namemap[-5201]	= "bu_0~";
      namemap[5203]	= "bu_1";
      namemap[-5203]	= "bu_1~";
      namemap[5212]	= "Sigma_b^0";
      namemap[-5212]	= "Sigma_b~^0";
      namemap[5214]	= "Sigma*_b^0";
      namemap[-5214]	= "Sigma*_b~^0";
      namemap[5222]	= "Sigma_b^+";
      namemap[-5222]	= "Sigma_b~^-";
      namemap[5224]	= "Sigma*_b^+";
      namemap[-5224]	= "Sigma*_b~^-";
      namemap[5232]	= "Xi_b^0";
      namemap[-5232]	= "Xi_b~^0";
      namemap[5242]	= "Xi_bc^+";
      namemap[-5242]	= "Xi_bc~^-";
      namemap[5301]	= "bs_0";
      namemap[-5301]	= "bs_0~";
      namemap[5303]	= "bs_1";
      namemap[-5303]	= "bs_1~";
      namemap[5312]	= "Xi'_b^-";
      namemap[-5312]	= "Xi'_b~^+";
      namemap[5314]	= "Xi*_b^-";
      namemap[-5314]	= "Xi*_b~^+";
      namemap[5322]	= "Xi'_b^0";
      namemap[-5322]	= "Xi'_b~^0";
      namemap[5324]	= "Xi*_b^0";
      namemap[-5324]	= "Xi*_b~^0";
      namemap[5332]	= "Omega_b^-";
      namemap[-5332]	= "Omega_b~^+";
      namemap[5334]	= "Omega*_b^-";
      namemap[-5334]	= "Omega*_b~^+";
      namemap[5342]	= "Omega_bc^0";
      namemap[-5342]	= "Omega_bc~^0";
      namemap[5401]	= "bc_0";
      namemap[-5401]	= "bc_0~";
      namemap[5403]	= "bc_1";
      namemap[-5403]	= "bc_1~";
      namemap[5412]	= "Xi'_bc^0";
      namemap[-5412]	= "Xi'_bc~^0";
      namemap[5414]	= "Xi*_bc^0";
      namemap[-5414]	= "Xi*_bc~^0";
      namemap[5422]	= "Xi'_bc^+";
      namemap[-5422]	= "Xi'_bc~^-";
      namemap[5424]	= "Xi*_bc^+";
      namemap[-5424]	= "Xi*_bc~^-";
      namemap[5432]	= "Omega'_bc^0";
      namemap[-5432]	= "Omega'_bc~^0";
      namemap[5434]	= "Omega*_bc^0";
      namemap[-5434]	= "Omega*_bc~^0";
      namemap[5442]	= "Omega_bcc^+";
      namemap[-5442]	= "Omega_bcc~^-";
      namemap[5444]	= "Omega*_bcc^+";
      namemap[-5444]	= "Omega*_bcc~^-";
      namemap[5503]	= "bb_1";
      namemap[-5503]	= "bb_1~";
      namemap[5512]	= "Xi_bb^-";
      namemap[-5512]	= "Xi_bb~^+";
      namemap[5514]	= "Xi*_bb^-";
      namemap[-5514]	= "Xi*_bb~^+";
      namemap[5522]	= "Xi_bb^0";
      namemap[-5522]	= "Xi_bb~^0";
      namemap[5524]	= "Xi*_bb^0";
      namemap[-5524]	= "Xi*_bb~^0";
      namemap[5532]	= "Omega_bb^-";
      namemap[-5532]	= "Omega_bb~^+";
      namemap[5534]	= "Omega*_bb^-";
      namemap[-5534]	= "Omega*_bb~^+";
      namemap[5542]	= "Omega_bbc^0";
      namemap[-5542]	= "Omega_bbc~^0";
      namemap[5544]	= "Omega*_bbc^0";
      namemap[-5544]	= "Omega*_bbc~^0";
      namemap[5554]	= "Omega*_bbb^-";
      namemap[-5554]	= "Omega*_bbb~^+";
      namemap[6101]	= "td_0";
      namemap[-6101]	= "td_0~";
      namemap[6103]	= "td_1";
      namemap[-6103]	= "td_1~";
      namemap[6112]	= "Sigma_t^0";
      namemap[-6112]	= "Sigma_t~^0";
      namemap[6114]	= "Sigma*_t^0";
      namemap[-6114]	= "Sigma*_t~^0";
      namemap[6122]	= "Lambda_t^+";
      namemap[-6122]	= "Lambda_t~^-";
      namemap[6132]	= "Xi_t^0";
      namemap[-6132]	= "Xi_t~^0";
      namemap[6142]	= "Xi_tc^+";
      namemap[-6142]	= "Xi_tc~^-";
      namemap[6152]	= "Xi_tb^0";
      namemap[-6152]	= "Xi_tb~^0";
      namemap[6201]	= "tu_0";
      namemap[-6201]	= "tu_0~";
      namemap[6203]	= "tu_1";
      namemap[-6203]	= "tu_1~";
      namemap[6212]	= "Sigma_t^+";
      namemap[-6212]	= "Sigma_t~^-";
      namemap[6214]	= "Sigma*_t^+";
      namemap[-6214]	= "Sigma*_t~^-";
      namemap[6222]	= "Sigma_t^++";
      namemap[-6222]	= "Sigma_t~^--";
      namemap[6224]	= "Sigma*_t^++";
      namemap[-6224]	= "Sigma*_t~^--";
      namemap[6232]	= "Xi_t^+";
      namemap[-6232]	= "Xi_t~^-";
      namemap[6242]	= "Xi_tc^++";
      namemap[-6242]	= "Xi_tc~^--";
      namemap[6252]	= "Xi_tb^+";
      namemap[-6252]	= "Xi_tb~^-";
      namemap[6301]	= "ts_0";
      namemap[-6301]	= "ts_0~";
      namemap[6303]	= "ts_1";
      namemap[-6303]	= "ts_1~";
      namemap[6312]	= "Xi'_t^0";
      namemap[-6312]	= "Xi'_t~^0";
      namemap[6314]	= "Xi*_t^0";
      namemap[-6314]	= "Xi*_t~^0";
      namemap[6322]	= "Xi'_t^+";
      namemap[-6322]	= "Xi'_t~^-";
      namemap[6324]	= "Xi*_t^+";
      namemap[-6324]	= "Xi*_t~^-";
      namemap[6332]	= "Omega_t^0";
      namemap[-6332]	= "Omega_t~^0";
      namemap[6334]	= "Omega*_t^0";
      namemap[-6334]	= "Omega*_t~^0";
      namemap[6342]	= "Omega_tc^+";
      namemap[-6342]	= "Omega_tc~^-";
      namemap[6352]	= "Omega_tb^0";
      namemap[-6352]	= "Omega_tb~^0";
      namemap[6401]	= "tc_0";
      namemap[-6401]	= "tc_0~";
      namemap[6403]	= "tc_1";
      namemap[-6403]	= "tc_1~";
      namemap[6412]	= "Xi'_tc^+";
      namemap[-6412]	= "Xi'_tc~^-";
      namemap[6414]	= "Xi*_tc^+";
      namemap[-6414]	= "Xi*_tc~^-";
      namemap[6422]	= "Xi'_tc^++";
      namemap[-6422]	= "Xi'_tc~^--";
      namemap[6424]	= "Xi*_tc^++";
      namemap[-6424]	= "Xi*_tc~^--";
      namemap[6432]	= "Omega'_tc^+";
      namemap[-6432]	= "Omega'_tc~^-";
      namemap[6434]	= "Omega*_tc^+";
      namemap[-6434]	= "Omega*_tc~^-";
      namemap[6442]	= "Omega_tcc^++";
      namemap[-6442]	= "Omega_tcc~^--";
      namemap[6444]	= "Omega*_tcc^++";
      namemap[-6444]	= "Omega*_tcc~^--";
      namemap[6452]	= "Omega_tbc^+";
      namemap[-6452]	= "Omega_tbc~^-";
      namemap[6501]	= "tb_0";
      namemap[-6501]	= "tb_0~";
      namemap[6503]	= "tb_1";
      namemap[-6503]	= "tb_1~";
      namemap[6512]	= "Xi'_tb^0";
      namemap[-6512]	= "Xi'_tb~^0";
      namemap[6514]	= "Xi*_tb^0";
      namemap[-6514]	= "Xi*_tb~^0";
      namemap[6522]	= "Xi'_tb^+";
      namemap[-6522]	= "Xi'_tb~^-";
      namemap[6524]	= "Xi*_tb^+";
      namemap[-6524]	= "Xi*_tb~^-";
      namemap[6532]	= "Omega'_tb^0";
      namemap[-6532]	= "Omega'_tb~^0";
      namemap[6534]	= "Omega*_tb^0";
      namemap[-6534]	= "Omega*_tb~^0";
      namemap[6542]	= "Omega'_tbc^+";
      namemap[-6542]	= "Omega'_tbc~^-";
      namemap[6544]	= "Omega*_tbc^+";
      namemap[-6544]	= "Omega*_tbc~^-";
      namemap[6552]	= "Omega_tbb^0";
      namemap[-6552]	= "Omega_tbb~^0";
      namemap[6554]	= "Omega*_tbb^0";
      namemap[-6554]	= "Omega*_tbb~^0";
      namemap[6603]	= "tt_1";
      namemap[-6603]	= "tt_1~";
      namemap[6612]	= "Xi_tt^+";
      namemap[-6612]	= "Xi_tt~^-";
      namemap[6614]	= "Xi*_tt^+";
      namemap[-6614]	= "Xi*_tt~^-";
      namemap[6622]	= "Xi_tt^++";
      namemap[-6622]	= "Xi_tt~^--";
      namemap[6624]	= "Xi*_tt^++";
      namemap[-6624]	= "Xi*_tt~^--";
      namemap[6632]	= "Omega_tt^+";
      namemap[-6632]	= "Omega_tt~^-";
      namemap[6634]	= "Omega*_tt^+";
      namemap[-6634]	= "Omega*_tt~^-";
      namemap[6642]	= "Omega_ttc^++";
      namemap[-6642]	= "Omega_ttc~^--";
      namemap[6644]	= "Omega*_ttc^++";
      namemap[-6644]	= "Omega*_ttc~^--";
      namemap[6652]	= "Omega_ttb^+";
      namemap[-6652]	= "Omega_ttb~^-";
      namemap[6654]	= "Omega*_ttb^+";
      namemap[-6654]	= "Omega*_ttb~^-";
      namemap[6664]	= "Omega*_ttt^++";
      namemap[-6664]	= "Omega*_ttt~^--";
      namemap[7101]	= "b'd_0";
      namemap[-7101]	= "b'd_0~";
      namemap[7103]	= "b'd_1";
      namemap[-7103]	= "b'd_1~";
      namemap[7112]	= "Sigma_b'^-";
      namemap[-7112]	= "Sigma_b'~^+";
      namemap[7114]	= "Sigma*_b'^-";
      namemap[-7114]	= "Sigma*_b'~^+";
      namemap[7122]	= "Lambda_b'^0";
      namemap[-7122]	= "Lambda_b'~^0";
      namemap[7132]	= "Xi_b'^-";
      namemap[-7132]	= "Xi_b'~^+";
      namemap[7142]	= "Xi_b'c^0";
      namemap[-7142]	= "Xi_b'c~^0";
      namemap[7152]	= "Xi_b'b^-";
      namemap[-7152]	= "Xi_b'b~^+";
      namemap[7162]	= "Xi_b't^0";
      namemap[-7162]	= "Xi_b't~^0";
      namemap[7201]	= "b'u_0";
      namemap[-7201]	= "b'u_0~";
      namemap[7203]	= "b'u_1";
      namemap[-7203]	= "b'u_1~";
      namemap[7212]	= "Sigma_b'^0";
      namemap[-7212]	= "Sigma_b'~^0";
      namemap[7214]	= "Sigma*_b'^0";
      namemap[-7214]	= "Sigma*_b'~^0";
      namemap[7222]	= "Sigma_b'^+";
      namemap[-7222]	= "Sigma_b'~^-";
      namemap[7224]	= "Sigma*_b'^+";
      namemap[-7224]	= "Sigma*_b'~^-";
      namemap[7232]	= "Xi_b'^0";
      namemap[-7232]	= "Xi_b'~^0";
      namemap[7242]	= "Xi_b'c^+";
      namemap[-7242]	= "Xi_b'c~^-";
      namemap[7252]	= "Xi_b'b^0";
      namemap[-7252]	= "Xi_b'b~^0";
      namemap[7262]	= "Xi_b't^+";
      namemap[-7262]	= "Xi_b't~^-";
      namemap[7301]	= "b's_0";
      namemap[-7301]	= "b's_0~";
      namemap[7303]	= "b's_1";
      namemap[-7303]	= "b's_1~";
      namemap[7312]	= "Xi'_b'^-";
      namemap[-7312]	= "Xi'_b'~^+";
      namemap[7314]	= "Xi*_b'^-";
      namemap[-7314]	= "Xi*_b'~^+";
      namemap[7322]	= "Xi'_b'^0";
      namemap[-7322]	= "Xi'_b'~^0";
      namemap[7324]	= "Xi*_b'^0";
      namemap[-7324]	= "Xi*_b'~^0";
      namemap[7332]	= "Omega'_b'^-";
      namemap[-7332]	= "Omega'_b'~^+";
      namemap[7334]	= "Omega*_b'^-";
      namemap[-7334]	= "Omega*_b'~^+";
      namemap[7342]	= "Omega_b'c^0";
      namemap[-7342]	= "Omega_b'c~^0";
      namemap[7352]	= "Omega_b'b^-";
      namemap[-7352]	= "Omega_b'b~^+";
      namemap[7362]	= "Omega_b't^0";
      namemap[-7362]	= "Omega_b't~^0";
      namemap[7401]	= "b'c_0";
      namemap[-7401]	= "b'c_0~";
      namemap[7403]	= "b'c_1";
      namemap[-7403]	= "b'c_1~";
      namemap[7412]	= "Xi'_b'c^0";
      namemap[-7412]	= "Xi'_b'c~^0";
      namemap[7414]	= "Xi*_b'c^0";
      namemap[-7414]	= "Xi*_b'c~^0";
      namemap[7422]	= "Xi'_b'c^+";
      namemap[-7422]	= "Xi'_b'c~^-";
      namemap[7424]	= "Xi*_b'c^+";
      namemap[-7424]	= "Xi*_b'c~^-";
      namemap[7432]	= "Omega'_b'c^0";
      namemap[-7432]	= "Omega'_b'c~^0";
      namemap[7434]	= "Omega*_b'c^0";
      namemap[-7434]	= "Omega*_b'c~^0";
      namemap[7442]	= "Omega'_b'cc^+";
      namemap[-7442]	= "Omega'_b'cc~^-";
      namemap[7444]	= "Omega*_b'cc^+";
      namemap[-7444]	= "Omega*_b'cc~^-";
      namemap[7452]	= "Omega_b'bc^0";
      namemap[-7452]	= "Omega_b'bc~^0";
      namemap[7462]	= "Omega_b'tc^+";
      namemap[-7462]	= "Omega_b'tc~^-";
      namemap[7501]	= "b'b_0";
      namemap[-7501]	= "b'b_0~";
      namemap[7503]	= "b'b_1";
      namemap[-7503]	= "b'b_1~";
      namemap[7512]	= "Xi'_b'b^-";
      namemap[-7512]	= "Xi'_b'b~^+";
      namemap[7514]	= "Xi*_b'b^-";
      namemap[-7514]	= "Xi*_b'b~^+";
      namemap[7522]	= "Xi'_b'b^0";
      namemap[-7522]	= "Xi'_b'b~^0";
      namemap[7524]	= "Xi*_b'b^0";
      namemap[-7524]	= "Xi*_b'b~^0";
      namemap[7532]	= "Omega'_b'b^-";
      namemap[-7532]	= "Omega'_b'b~^+";
      namemap[7534]	= "Omega*_b'b^-";
      namemap[-7534]	= "Omega*_b'b~^+";
      namemap[7542]	= "Omega'_b'bc^0";
      namemap[-7542]	= "Omega'_b'bc~^0";
      namemap[7544]	= "Omega*_b'bc^0";
      namemap[-7544]	= "Omega*_b'bc~^0";
      namemap[7552]	= "Omega'_b'bb^-";
      namemap[-7552]	= "Omega'_b'bb~^+";
      namemap[7554]	= "Omega*_b'bb^-";
      namemap[-7554]	= "Omega*_b'bb~^+";
      namemap[7562]	= "Omega_b'tb^0";
      namemap[-7562]	= "Omega_b'tb~^0";
      namemap[7601]	= "b't_0";
      namemap[-7601]	= "b't_0~";
      namemap[7603]	= "b't_1";
      namemap[-7603]	= "b't_1~";
      namemap[7612]	= "Xi'_b't^0";
      namemap[-7612]	= "Xi'_b't~^0";
      namemap[7614]	= "Xi*_b't^0";
      namemap[-7614]	= "Xi*_b't~^0";
      namemap[7622]	= "Xi'_b't^+";
      namemap[-7622]	= "Xi'_b't~^-";
      namemap[7624]	= "Xi*_b't^+";
      namemap[-7624]	= "Xi*_b't~^-";
      namemap[7632]	= "Omega'_b't^0";
      namemap[-7632]	= "Omega'_b't~^0";
      namemap[7634]	= "Omega*_b't^0";
      namemap[-7634]	= "Omega*_b't~^0";
      namemap[7642]	= "Omega'_b'tc^+";
      namemap[-7642]	= "Omega'_b'tc~^-";
      namemap[7644]	= "Omega*_b'tc^+";
      namemap[-7644]	= "Omega*_b'tc~^-";
      namemap[7652]	= "Omega'_b'tb^0";
      namemap[-7652]	= "Omega'_b'tb~^0";
      namemap[7654]	= "Omega*_b'tb^0";
      namemap[-7654]	= "Omega*_b'tb~^0";
      namemap[7662]	= "Omega'_b'tt^+";
      namemap[-7662]	= "Omega'_b'tt~^-";
      namemap[7664]	= "Omega*_b'tt^+";
      namemap[-7664]	= "Omega*_b'tt~^-";
      namemap[7703]	= "b'b'_1";
      namemap[-7703]	= "b'b'_1~";
      namemap[7712]	= "Xi'_b'b'^-";
      namemap[-7712]	= "Xi'_b'b'~^+";
      namemap[7714]	= "Xi*_b'b'^-";
      namemap[-7714]	= "Xi*_b'b'~^+";
      namemap[7722]	= "Xi'_b'b'^0";
      namemap[-7722]	= "Xi'_b'b'~^0";
      namemap[7724]	= "Xi*_b'b'^0";
      namemap[-7724]	= "Xi*_b'b'~^0";
      namemap[7732]	= "Omega'_b'b'^-";
      namemap[-7732]	= "Omega'_b'b'~^+";
      namemap[7734]	= "Omega*_b'b'^-";
      namemap[-7734]	= "Omega*_b'b'~^+";
      namemap[7742]	= "Omega'_b'b'c^0";
      namemap[-7742]	= "Omega'_b'b'c~^0";
      namemap[7744]	= "Omega*_b'b'c^0";
      namemap[-7744]	= "Omega*_b'b'c~^0";
      namemap[7752]	= "Omega'_b'b'b^-";
      namemap[-7752]	= "Omega'_b'b'b~^+";
      namemap[7754]	= "Omega*_b'b'b^-";
      namemap[-7754]	= "Omega*_b'b'b~^+";
      namemap[7762]	= "Omega'_b'b't^0";
      namemap[-7762]	= "Omega'_b'b't~^0";
      namemap[7764]	= "Omega*_b'b't^0";
      namemap[-7764]	= "Omega*_b'b't~^0";
      namemap[7774]	= "Omega*_b'b'b'^-";
      namemap[-7774]	= "Omega*_b'b'b'~^+";
      namemap[8101]	= "t'd_0";
      namemap[-8101]	= "t'd_0~";
      namemap[8103]	= "t'd_1";
      namemap[-8103]	= "t'd_1~";
      namemap[8112]	= "Sigma_t'^0";
      namemap[-8112]	= "Sigma_t'~^0";
      namemap[8114]	= "Sigma*_t'^0";
      namemap[-8114]	= "Sigma*_t'~^0";
      namemap[8122]	= "Lambda_t'^+";
      namemap[-8122]	= "Lambda_t'~^-";
      namemap[8132]	= "Xi_t'^0";
      namemap[-8132]	= "Xi_t'~^0";
      namemap[8142]	= "Xi_t'c^+";
      namemap[-8142]	= "Xi_t'c~^-";
      namemap[8152]	= "Xi_t'b^0";
      namemap[-8152]	= "Xi_t'b~^0";
      namemap[8162]	= "Xi_t't^+";
      namemap[-8162]	= "Xi_t't~^-";
      namemap[8172]	= "Xi_t'b'^0";
      namemap[-8172]	= "Xi_t'b'~^0";
      namemap[8201]	= "t'u_0";
      namemap[-8201]	= "t'u_0~";
      namemap[8203]	= "t'u_1";
      namemap[-8203]	= "t'u_1~";
      namemap[8212]	= "Sigma_t'^+";
      namemap[-8212]	= "Sigma_t'~^-";
      namemap[8214]	= "Sigma*_t'^+";
      namemap[-8214]	= "Sigma*_t'~^-";
      namemap[8222]	= "Sigma_t'^++";
      namemap[-8222]	= "Sigma_t'~^--";
      namemap[8224]	= "Sigma*_t'^++";
      namemap[-8224]	= "Sigma*_t'~^--";
      namemap[8232]	= "Xi_t'^+";
      namemap[-8232]	= "Xi_t'~^-";
      namemap[8242]	= "Xi_t'c^++";
      namemap[-8242]	= "Xi_t'c~^--";
      namemap[8252]	= "Xi_t'b^+";
      namemap[-8252]	= "Xi_t'b~^-";
      namemap[8262]	= "Xi_t't^++";
      namemap[-8262]	= "Xi_t't~^--";
      namemap[8272]	= "Xi_t'b'^+";
      namemap[-8272]	= "Xi_t'b'~^-";
      namemap[8301]	= "t's_0";
      namemap[-8301]	= "t's_0~";
      namemap[8303]	= "t's_1";
      namemap[-8303]	= "t's_1~";
      namemap[8312]	= "Xi'_t'^0";
      namemap[-8312]	= "Xi'_t'~^0";
      namemap[8314]	= "Xi*_t'^0";
      namemap[-8314]	= "Xi*_t'~^0";
      namemap[8322]	= "Xi'_t'^+";
      namemap[-8322]	= "Xi'_t'~^-";
      namemap[8324]	= "Xi*_t'^+";
      namemap[-8324]	= "Xi*_t'~^-";
      namemap[8332]	= "Omega'_t'^0";
      namemap[-8332]	= "Omega'_t'~^0";
      namemap[8334]	= "Omega*_t'^0";
      namemap[-8334]	= "Omega*_t'~^0";
      namemap[8342]	= "Omega_t'c^+";
      namemap[-8342]	= "Omega_t'c~^-";
      namemap[8352]	= "Omega_t'b^0";
      namemap[-8352]	= "Omega_t'b~^0";
      namemap[8362]	= "Omega_t't^+";
      namemap[-8362]	= "Omega_t't~^-";
      namemap[8372]	= "Omega_t'b'^0";
      namemap[-8372]	= "Omega_t'b'~^0";
      namemap[8401]	= "t'c_0";
      namemap[-8401]	= "t'c_0~";
      namemap[8403]	= "t'c_1";
      namemap[-8403]	= "t'c_1~";
      namemap[8412]	= "Xi'_t'c^+";
      namemap[-8412]	= "Xi'_t'c~^-";
      namemap[8414]	= "Xi*_t'c^+";
      namemap[-8414]	= "Xi*_t'c~^-";
      namemap[8422]	= "Xi'_t'c^++";
      namemap[-8422]	= "Xi'_t'c~^--";
      namemap[8424]	= "Xi*_t'c^++";
      namemap[-8424]	= "Xi*_t'c~^--";
      namemap[8432]	= "Omega'_t'c^+";
      namemap[-8432]	= "Omega'_t'c~^-";
      namemap[8434]	= "Omega*_t'c^+";
      namemap[-8434]	= "Omega*_t'c~^-";
      namemap[8442]	= "Omega'_t'cc^++";
      namemap[-8442]	= "Omega'_t'cc~^--";
      namemap[8444]	= "Omega*_t'cc^++";
      namemap[-8444]	= "Omega*_t'cc~^--";
      namemap[8452]	= "Omega_t'bc^+";
      namemap[-8452]	= "Omega_t'bc~^-";
      namemap[8462]	= "Omega_t'tc^++";
      namemap[-8462]	= "Omega_t'tc~^--";
      namemap[8472]	= "Omega_t'b'c ^+";
      namemap[-8472]	= "Omega_t'b'c ~^-";
      namemap[8501]	= "t'b_0";
      namemap[-8501]	= "t'b_0~";
      namemap[8503]	= "t'b_1";
      namemap[-8503]	= "t'b_1~";
      namemap[8512]	= "Xi'_t'b^0";
      namemap[-8512]	= "Xi'_t'b~^0";
      namemap[8514]	= "Xi*_t'b^0";
      namemap[-8514]	= "Xi*_t'b~^0";
      namemap[8522]	= "Xi'_t'b^+";
      namemap[-8522]	= "Xi'_t'b~^-";
      namemap[8524]	= "Xi*_t'b^+";
      namemap[-8524]	= "Xi*_t'b~^-";
      namemap[8532]	= "Omega'_t'b^0";
      namemap[-8532]	= "Omega'_t'b~^0";
      namemap[8534]	= "Omega*_t'b^0";
      namemap[-8534]	= "Omega*_t'b~^0";
      namemap[8542]	= "Omega'_t'bc^+";
      namemap[-8542]	= "Omega'_t'bc~^-";
      namemap[8544]	= "Omega*_t'bc^+";
      namemap[-8544]	= "Omega*_t'bc~^-";
      namemap[8552]	= "Omega'_t'bb^0";
      namemap[-8552]	= "Omega'_t'bb~^0";
      namemap[8554]	= "Omega*_t'bb^0";
      namemap[-8554]	= "Omega*_t'bb~^0";
      namemap[8562]	= "Omega_t'tb^+";
      namemap[-8562]	= "Omega_t'tb~^-";
      namemap[8572]	= "Omega_t'b'b ^0";
      namemap[-8572]	= "Omega_t'b'b ~^0";
      namemap[8601]	= "t't_0";
      namemap[-8601]	= "t't_0~";
      namemap[8603]	= "t't_1";
      namemap[-8603]	= "t't_1~";
      namemap[8612]	= "Xi'_t't^+";
      namemap[-8612]	= "Xi'_t't~^-";
      namemap[8614]	= "Xi*_t't^+";
      namemap[-8614]	= "Xi*_t't~^-";
      namemap[8622]	= "Xi'_t't^++";
      namemap[-8622]	= "Xi'_t't~^--";
      namemap[8624]	= "Xi*_t't^++";
      namemap[-8624]	= "Xi*_t't~^--";
      namemap[8632]	= "Omega'_t't^+";
      namemap[-8632]	= "Omega'_t't~^-";
      namemap[8634]	= "Omega*_t't^+";
      namemap[-8634]	= "Omega*_t't~^-";
      namemap[8642]	= "Omega'_t'tc^++";
      namemap[-8642]	= "Omega'_t'tc~^--";
      namemap[8644]	= "Omega*_t'tc^++";
      namemap[-8644]	= "Omega*_t'tc~^--";
      namemap[8652]	= "Omega'_t'tb^+";
      namemap[-8652]	= "Omega'_t'tb~^-";
      namemap[8654]	= "Omega*_t'tb^+";
      namemap[-8654]	= "Omega*_t'tb~^-";
      namemap[8662]	= "Omega'_t'tt^++";
      namemap[-8662]	= "Omega'_t'tt~^--";
      namemap[8664]	= "Omega*_t'tt^++";
      namemap[-8664]	= "Omega*_t'tt~^--";
      namemap[8672]	= "Omega_t'b't ^+";
      namemap[-8672]	= "Omega_t'b't ~^-";
      namemap[8701]	= "t'b'_0";
      namemap[-8701]	= "t'b'_0~";
      namemap[8703]	= "t'b'_1";
      namemap[-8703]	= "t'b'_1~";
      namemap[8712]	= "Xi'_t'b'^0";
      namemap[-8712]	= "Xi'_t'b'~^0";
      namemap[8714]	= "Xi*_t'b'^0";
      namemap[-8714]	= "Xi*_t'b'~^0";
      namemap[8722]	= "Xi'_t'b'^+";
      namemap[-8722]	= "Xi'_t'b'~^-";
      namemap[8724]	= "Xi*_t'b'^+";
      namemap[-8724]	= "Xi*_t'b'~^-";
      namemap[8732]	= "Omega'_t'b'^0";
      namemap[-8732]	= "Omega'_t'b'~^0";
      namemap[8734]	= "Omega*_t'b'^0";
      namemap[-8734]	= "Omega*_t'b'~^0";
      namemap[8742]	= "Omega'_t'b'c^+";
      namemap[-8742]	= "Omega'_t'b'c~^-";
      namemap[8744]	= "Omega*_t'b'c^+";
      namemap[-8744]	= "Omega*_t'b'c~^-";
      namemap[8752]	= "Omega'_t'b'b^0";
      namemap[-8752]	= "Omega'_t'b'b~^0";
      namemap[8754]	= "Omega*_t'b'b^0";
      namemap[-8754]	= "Omega*_t'b'b~^0";
      namemap[8762]	= "Omega'_t'b't^+";
      namemap[-8762]	= "Omega'_t'b't~^-";
      namemap[8764]	= "Omega*_t'b't^+";
      namemap[-8764]	= "Omega*_t'b't~^-";
      namemap[8772]	= "Omega'_t'b'b'^0";
      namemap[-8772]	= "Omega'_t'b'b'~^0";
      namemap[8774]	= "Omega*_t'b'b'^0";
      namemap[-8774]	= "Omega*_t'b'b'~^0";
      namemap[8803]	= "t't'_1";
      namemap[-8803]	= "t't'_1~";
      namemap[8812]	= "Xi'_t't'^+";
      namemap[-8812]	= "Xi'_t't'~^-";
      namemap[8814]	= "Xi*_t't'^+";
      namemap[-8814]	= "Xi*_t't'~^-";
      namemap[8822]	= "Xi'_t't'^++";
      namemap[-8822]	= "Xi'_t't'~^--";
      namemap[8824]	= "Xi*_t't'^++";
      namemap[-8824]	= "Xi*_t't'~^--";
      namemap[8832]	= "Omega'_t't'^+";
      namemap[-8832]	= "Omega'_t't'~^-";
      namemap[8834]	= "Omega*_t't'^+";
      namemap[-8834]	= "Omega*_t't'~^-";
      namemap[8842]	= "Omega'_t't'c^++";
      namemap[-8842]	= "Omega'_t't'c~^--";
      namemap[8844]	= "Omega*_t't'c^++";
      namemap[-8844]	= "Omega*_t't'c~^--";
      namemap[8852]	= "Omega'_t't'b^+";
      namemap[-8852]	= "Omega'_t't'b~^-";
      namemap[8854]	= "Omega*_t't'b^+";
      namemap[-8854]	= "Omega*_t't'b~^-";
      namemap[8862]	= "Omega'_t't't^++";
      namemap[-8862]	= "Omega'_t't't~^--";
      namemap[8864]	= "Omega*_t't't^++";
      namemap[-8864]	= "Omega*_t't't~^--";
      namemap[8872]	= "Omega'_t't'b'^+";
      namemap[-8872]	= "Omega'_t't'b'~^-";
      namemap[8874]	= "Omega*_t't'b'^+";
      namemap[-8874]	= "Omega*_t't'b'~^-";
      namemap[8884]	= "Omega*_t't't'^++";
      namemap[-8884]	= "Omega*_t't't'~^--";
      namemap[9990]	= "odderon";
      namemap[10022]	= "virtual-photon";
      namemap[10111]	= "a_0(1450)^0";
      namemap[10113]	= "b_1(1235)^0";
      namemap[10115]	= "pi_2(1670)^0";
      namemap[10211]	= "a_0(1450)^+";
      namemap[-10211]	= "a_0(1450)^-";
      namemap[10213]	= "b_1(1235)^+";
      namemap[-10213]	= "b_1(1235)^-";
      namemap[10215]	= "pi_2(1670)^+";
      namemap[-10215]	= "pi_2(1670)^-";
      namemap[10221]	= "f_0(1370)";
      namemap[10223]	= "h_1(1170)";
      namemap[10225]	= "eta_2(1645)";
      namemap[10311]	= "K*_0(1430)^0";
      namemap[-10311]	= "K*_0(1430)~^0";
      namemap[10313]	= "K_1(1270)^0";
      namemap[-10313]	= "K_1(1270)~^0";
      namemap[10315]	= "K_2(1770)^0";
      namemap[-10315]	= "K_2(1770)~^0";
      namemap[10321]	= "K*_0(1430)^+";
      namemap[-10321]	= "K*_0(1430)^-";
      namemap[10323]	= "K_1(1270)^+";
      namemap[-10323]	= "K_1(1270)^-";
      namemap[10325]	= "K_2(1770)^+";
      namemap[-10325]	= "K_2(1770)^-";
      namemap[10331]	= "f_0(1710)";
      namemap[10333]	= "h_1(1380)";
      namemap[10335]	= "eta_2(1870)";
      namemap[10411]	= "D*_0(2400)^+";
      namemap[-10411]	= "D*_0(2400)^-";
      namemap[10413]	= "D_1(2420)^+";
      namemap[-10413]	= "D_1(2420)^-";
      namemap[10421]	= "D*_0(2400)^0";
      namemap[-10421]	= "D*_0(2400)~^0";
      namemap[10423]	= "D_1(2420)^0";
      namemap[-10423]	= "D_1(2420)~^0";
      namemap[10431]	= "D*_s0(2317)^+";
      namemap[-10431]	= "D*_s0(2317)^-";
      namemap[10433]	= "D_s1(2536)^+";
      namemap[-10433]	= "D_s1(2536)^-";
      namemap[10441]	= "chi_c0(1P)";
      namemap[10443]	= "hc(1P)";
      namemap[10511]	= "B*_0^0";
      namemap[-10511]	= "B*_0~^0";
      namemap[10513]	= "B_1(L)^0";
      namemap[-10513]	= "B_1(L)~^0";
      namemap[10521]	= "B*_0^+";
      namemap[-10521]	= "B*_0^-";
      namemap[10523]	= "B_1(L)^+";
      namemap[-10523]	= "B_1(L)^-";
      namemap[10531]	= "B*_s0^0";
      namemap[-10531]	= "B*_s0~^0";
      namemap[10533]	= "B_s1(L)^0";
      namemap[-10533]	= "B_s1(L)~^0";
      namemap[10541]	= "B*_c0^+";
      namemap[-10541]	= "B*_c0^-";
      namemap[10543]	= "B_c1(L)^+";
      namemap[-10543]	= "B_c1(L)^-";
      namemap[10551]	= "chi_b0(1P)";
      namemap[10553]	= "h_b(1P)";
      namemap[10555]	= "eta_b2(1D)";
      namemap[11114]	= "Delta(1700)^-";
      namemap[11116]	= "Delta(1930)^-";
      namemap[11216]	= "Delta(1930)^0";
      namemap[12112]	= "N(1440)^0";
      namemap[12114]	= "Delta(1700)^0";
      namemap[12116]	= "N(1680)^0";
      namemap[12126]	= "Delta(1930)^+";
      namemap[12212]	= "N(1440)^+";
      namemap[12214]	= "Delta(1700)^+";
      namemap[12216]	= "N(1680)^+";
      namemap[12224]	= "Delta(1700)^++";
      namemap[12226]	= "Delta(1930)^++";
      namemap[13112]	= "Sigma(1660)^-";
      namemap[-13112]	= "Sigma~(1660)^-";
      namemap[13114]	= "Sigma(1670)^-";
      namemap[-13114]	= "Sigma~(1670)^-";
      namemap[13116]	= "Sigma(1915)^-";
      namemap[-13116]	= "Sigma~(1915)^-";
      namemap[13122]	= "Lambda(1405)^0";
      namemap[-13122]	= "Lambda~(1405)^0";
      namemap[13124]	= "Lambda(1690)^0";
      namemap[-13124]	= "Lambda~(1690)^0";
      namemap[13126]	= "Lambda(1830)^0";
      namemap[-13126]	= "Lambda~(1830)^0";
      namemap[13212]	= "Sigma(1660)^0";
      namemap[-13212]	= "Sigma~(1660)^0";
      namemap[13214]	= "Sigma(1670)^0";
      namemap[-13214]	= "Sigma~(1670)^0";
      namemap[13216]	= "Sigma(1915)^0";
      namemap[-13216]	= "Sigma~(1915)^0";
      namemap[13222]	= "Sigma(1660)^+";
      namemap[-13222]	= "Sigma~(1660)^+";
      namemap[13224]	= "Sigma(1670)^+";
      namemap[-13224]	= "Sigma~(1670)^+";
      namemap[13226]	= "Sigma(1915)^+";
      namemap[-13226]	= "Sigma~(1915)^+";
      namemap[13314]	= "Xi(1820)^-";
      namemap[-13314]	= "Xi(1820)~^+";
      namemap[13324]	= "Xi(1820)^0";
      namemap[-13324]	= "Xi(1820)~^0";
      namemap[14122]	= "Lambda_c(2593)^+";
      namemap[-14122]	= "Lambda_c~(2593)^-";
      namemap[14124]	= "Lambda_c(2625)^+";
      namemap[-14124]	= "Lambda_c~(2625)^-";
      namemap[20022]	= "Cerenkov-radiation";
      namemap[20113]	= "a_1(1260)^0";
      namemap[20213]	= "a_1(1260)^+";
      namemap[-20213]	= "a_1(1260)^-";
      namemap[20223]	= "f_1(1285)";
      namemap[20313]	= "K_1(1400)^0";
      namemap[-20313]	= "K_1(1400)~^0";
      namemap[20315]	= "K_2(1820)^0";
      namemap[-20315]	= "K_2(1820)~^0";
      namemap[20323]	= "K_1(1400)^+";
      namemap[-20323]	= "K_1(1400)^-";
      namemap[20325]	= "K_2(1820)^+";
      namemap[-20325]	= "K_2(1820)^-";
      namemap[20333]	= "f_1(1420)";
      namemap[20413]	= "D_1(H)^+";
      namemap[-20413]	= "D_1(H)^-";
      namemap[20423]	= "D_1(2430)^0";
      namemap[-20423]	= "D_1(2430)~^0";
      namemap[20433]	= "D_s1(2460)^+";
      namemap[-20433]	= "D_s1(2460)^-";
      namemap[20443]	= "chi_c1(1P)";
      namemap[20513]	= "B_1(H)^0";
      namemap[-20513]	= "B_1(H)~^0";
      namemap[20523]	= "B_1(H)^+";
      namemap[-20523]	= "B_1(H)^-";
      namemap[20533]	= "B_s1(H)^0";
      namemap[-20533]	= "B_s1(H)~^0";
      namemap[20543]	= "B_c1(H)^+";
      namemap[-20543]	= "B_c1(H)^-";
      namemap[20553]	= "chi_b1(1P)";
      namemap[20555]	= "Upsilon_2(1D)";
      namemap[21112]	= "Delta(1910)^-";
      namemap[21114]	= "Delta(1920)^-";
      namemap[21212]	= "Delta(1910)^0";
      namemap[21214]	= "N(1700)^0";
      namemap[22112]	= "N(1535)^0";
      namemap[22114]	= "Delta(1920)^0";
      namemap[22122]	= "Delta(1910)^+";
      namemap[22124]	= "N(1700)^+";
      namemap[22212]	= "N(1535)^+";
      namemap[22214]	= "Delta(1920)^+";
      namemap[22222]	= "Delta(1910)^++";
      namemap[22224]	= "Delta(1920)^++";
      namemap[23112]	= "Sigma(1750)^-";
      namemap[-23112]	= "Sigma~(1750)^-";
      namemap[23114]	= "Sigma(1940)^-";
      namemap[-23114]	= "Sigma~(1940)^-";
      namemap[23122]	= "Lambda(1600)^0";
      namemap[-23122]	= "Lambda~(1600)^0";
      namemap[23124]	= "Lambda(1890)^0";
      namemap[-23124]	= "Lambda~(1890)^0";
      namemap[23126]	= "Lambda(2110)^0";
      namemap[-23126]	= "Lambda~(2110)^0";
      namemap[23212]	= "Sigma(1750)^0";
      namemap[-23212]	= "Sigma~(1750)^0";
      namemap[23214]	= "Sigma(1940)^0";
      namemap[-23214]	= "Sigma~(1940)^0";
      namemap[23222]	= "Sigma(1750)^+";
      namemap[-23222]	= "Sigma~(1750)^+";
      namemap[23224]	= "Sigma(1940)^+";
      namemap[-23224]	= "Sigma~(1940)^+";
      namemap[30113]	= "rho(1700)^0";
      namemap[30213]	= "rho(1700)^+";
      namemap[-30213]	= "rho(1700)^-";
      namemap[30223]	= "omega(1650)";
      namemap[30313]	= "K*(1680)^0";
      namemap[-30313]	= "K*(1680)~^0";
      namemap[30323]	= "K*(1680)^+";
      namemap[-30323]	= "K*(1680)^-";
      namemap[30443]	= "psi(3770)";
      namemap[30553]	= "Upsilon_1(1D)";
      namemap[31114]	= "Delta(1600)^-";
      namemap[31214]	= "N(1720)^0";
      namemap[32112]	= "N(1650)^0";
      namemap[32114]	= "Delta(1600)^0";
      namemap[32124]	= "N(1720)^+";
      namemap[32212]	= "N(1650)^+";
      namemap[32214]	= "Delta(1600)^+";
      namemap[32224]	= "Delta(1600)^++";
      namemap[33122]	= "Lambda(1670)^0";
      namemap[-33122]	= "Lambda~(1670)^0";
      namemap[42112]	= "N(1710)^0";
      namemap[42212]	= "N(1710)^+";
      namemap[43122]	= "Lambda(1800)^0";
      namemap[-43122]	= "Lambda~(1800)^0";
      namemap[53122]	= "Lambda(1810)^0";
      namemap[-53122]	= "Lambda~(1810)^0";
      namemap[100111]	= "pi(1300)^0";
      namemap[100113]	= "rho(1450)^0";
      namemap[100211]	= "pi(1300)^+";
      namemap[-100211]	= "pi(1300)^-";
      namemap[100213]	= "rho(1450)^+";
      namemap[-100213]	= "rho(1450)^-";
      namemap[100221]	= "eta(1295)";
      namemap[100223]	= "omega(1420)";
      namemap[100311]	= "K(1460)^0";
      namemap[-100311]	= "K(1460)~^0";
      namemap[100313]	= "K*(1410)^0";
      namemap[-100313]	= "K*(1410)~^0";
      namemap[100321]	= "K(1460)^+";
      namemap[-100321]	= "K(1460)^-";
      namemap[100323]	= "K*(1410)^+";
      namemap[-100323]	= "K*(1410)^-";
      namemap[100325]	= "K_2(1980)^+";
      namemap[-100325]	= "K_2(1980)^-";
      namemap[100331]	= "eta(1475)";
      namemap[100333]	= "phi(1680)";
      namemap[100411]	= "D(2S)^+";
      namemap[-100411]	= "D(2S)^-";
      namemap[100413]	= "D*(2S)^+";
      namemap[-100413]	= "D*(2S)^+";
      namemap[100421]	= "D(2S)^0";
      namemap[-100421]	= "D(2S)~^0";
      namemap[100423]	= "D*(2S)^0";
      namemap[-100423]	= "D*(2S)~^0";
      namemap[100441]	= "eta_c(2S)";
      namemap[100443]	= "psi(2S)";
      namemap[100445]	= "chi_c2(2P)";
      namemap[100551]	= "eta_b(2S)";
      namemap[100553]	= "Upsilon(2S)";
      namemap[100555]	= "chi_b2(2P)";
      namemap[100557]	= "Upsilon_3(2D)";
      namemap[110551]	= "chi_b0(2P)";
      namemap[110553]	= "h_b(2P)";
      namemap[110555]	= "eta_b2(2D)";
      namemap[120553]	= "chi_b1(2P)";
      namemap[120555]	= "Upsilon_2(2D)";
      namemap[130553]	= "Upsilon_1(2D)";
      namemap[200551]	= "eta_b(3S)";
      namemap[200553]	= "Upsilon(3S)";
      namemap[200555]	= "chi_b2(3P)";
      namemap[210551]	= "chi_b0(3P)";
      namemap[210553]	= "h_b(3P)";
      namemap[220553]	= "chi_b1(3P)";
      namemap[300553]	= "Upsilon(4S)";
      // SUSY
      namemap[1000001]	= "~d_L";
      namemap[-1000001]	= "~d_L~";
      namemap[2000001]	= "~d_R";
      namemap[-2000001]	= "~d_R~";
      namemap[1000002]	= "~u_L";
      namemap[-1000002]	= "~u_L~";
      namemap[2000002]	= "~u_R";
      namemap[-2000002]	= "~u_R~";
      namemap[1000003]	= "~s_L";
      namemap[-1000003]	= "~s_L~";
      namemap[2000003]	= "~s_R";
      namemap[-2000003]	= "~s_R~";
      namemap[1000004]	= "~c_L";
      namemap[-1000004]	= "~c_L~";
      namemap[2000004]	= "~c_R";
      namemap[-2000004]	= "~c_R~";
      namemap[1000005]	= "~b_1";
      namemap[-1000005]	= "~b_1~";
      namemap[2000005]	= "~b_2";
      namemap[-2000005]	= "~b_2~";
      namemap[1000006]	= "~t_1";
      namemap[-1000006]	= "~t_1~";
      namemap[2000006]	= "~t_2";
      namemap[-2000006]	= "~t_2~";
      namemap[1000011]	= "~e_L-";
      namemap[-1000011]	= "~e_L+";
      namemap[2000011]	= "~e_R-";
      namemap[-2000011]	= "~e_R+";
      namemap[1000012]	= "~nu_eL";
      namemap[-1000012]	= "~nu_eL~";
      namemap[2000012]	= "~nu_eR";
      namemap[-2000012]	= "~nu_eR~";
      namemap[1000013]	= "~mu_L-";
      namemap[-1000013]	= "~mu_L+";
      namemap[2000013]	= "~mu_R-";
      namemap[-2000013]	= "~mu_R+";
      namemap[1000014]	= "~nu_muL";
      namemap[-1000014]	= "~nu_muL~";
      namemap[2000014]	= "~nu_muR";
      namemap[-2000014]	= "~nu_muR~";
      namemap[1000015]	= "~tau_L-";
      namemap[-1000015]	= "~tau_L+";
      namemap[2000015]	= "~tau_R-";
      namemap[-2000015]	= "~tau_R+";
      namemap[1000016]	= "~nu_tauL";
      namemap[-1000016]	= "~nu_tauL~";
      namemap[2000016]	= "~nu_tauR";
      namemap[-2000016]	= "~nu_tauR~";
      namemap[1000021]	= "~g";
      namemap[-1000021]	= "~g~";
      namemap[1000025]	= "~chi_30";
      namemap[-1000025]	= "~chi_30~";
      namemap[1000022]	= "~chi_10";
      namemap[-1000022]	= "~chi_10~";
      namemap[1000035]	= "~chi_40";
      namemap[-1000035]	= "~chi_40~";
      namemap[1000023]	= "~chi_20";
      namemap[-1000023]	= "~chi_20~";
      namemap[1000037]	= "~chi_2+";
      namemap[-1000037]	= "~chi_2-";
      namemap[1000024]	= "~chi_1+";
      namemap[-1000024]	= "~chi_1-";
      namemap[1000039]	= "~Gravitino";
      namemap[-1000039]	= "~Gravitino~";

    }
  if ( namemap.find(pdgid) != namemap.end() )
    return namemap[pdgid];
  else
    return string("not defined");
}
