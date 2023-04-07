void merger_fil() { // имя функции = имя файла
    std::vector<std::string> names {"21.0.root", "21.3.root", "23.1.root","23.3.root",
    "23.5.root", "21.2.root", "23.0.root", "23.2.root", "23.4.root","23.6.root"};
	TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
	TH1F *hist = new TH1F("hist", "", 500, 8, 12);
	
	for (auto iter {names.begin()}; iter != names.end(); ++iter) {
		TFile *input = new TFile((*iter).c_str(), "read");
		TTree *tree = (TTree*)input->Get("h1");
		
		TH1F *temp = new TH1F("temp", "", 100, 8, 12);
		tree->Draw("en >> temp", "p > 0.5 && en > 8 ");
		
		hist->Add(temp, 1); // добавляем данные temp с коэффициентом 1
	}
	
	hist->GetYaxis()->SetRangeUser(0, 6000);
	hist->Draw();
}
