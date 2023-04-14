#include<test.h>

int main() {

  
  std::string path_name = fs::current_path();
  std::vector<std::string> names;
  std::string ext = ".root";

  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
	TH1F *hist = new TH1F("hist", "", 100, 8, 12);

  for (auto &entry : fs::directory_iterator(path_name)){
    names.push_back(entry.path().filename().c_str()); 
  }

  for (auto iter {names.begin()}; iter != names.end(); ++iter) {
    std::string  fname = (*iter).c_str();
    if (fname.find(ext) != std::string::npos){
      TFile *input = new TFile(fname, "read");
      TTree *tree = (TTree*)input->Get("h1");
      
      TH1F *temp = new TH1F("temp", "", 100, 8, 12);
      tree->Draw("en >> temp", "p < 0.5 && en > 8 && ch == 2 && abs(ml - 2.28646) < 0.015 && abs(md - 2.28646) < 0.015 && ntr == 0");
      
      hist->Add(temp, 1); // добавляем данные temp с коэффициентом 1
    }
	}

  hist->Draw();
  
  return 0;
}
