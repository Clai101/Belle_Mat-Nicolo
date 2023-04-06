void merger_fil() { // имя функции = имя файла
    std::vector<std::string> names {"21.0.root", "21.3.root", "23.1.root","23.3.root",
    "23.5.root", "21.2.root", "23.0.root", "23.2.root", "23.4.root","23.6.root"};
	TFile *input;
	TH1F *hist = new TH1F("hist", "en", 250, 8, 12); // создаем гистограмму, 100 бинов, 0 < x < 1
	TCanvas *c1 = new TCanvas("c1", "c1", 800, 600); // окно в котором будет график
    for (auto iter {names.begin()}; iter != names.end(); ++iter){
		//cout << *iter << "\t";
		input = new TFile((*iter).c_str(), "read"); // input = твои данныen >> hist", "p > 1");
		TTree *tree = (TTree*)input->Get("h1");
		
		TH1F *temp = new TH1F("temp", "", 100, 8, 12);
		
		hist->Add(temp, 1); // добавляем данные temp с коэффициентом 1
    }
	
	tree->Draw("en >> temp", "p > 2 $$ en > 8");
	
	hist->SetFillColor(kGreen-9); // всякие украшения для гистограммы
	hist->GetXaxis()->SetTitle("en");
	hist->GetYaxis()->SetTitle("Count");
	hist->GetXaxis()->SetTitleSize(0.045);
	hist->GetYaxis()->SetTitleSize(0.045);
	hist->Draw("pe1");
}
