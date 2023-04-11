void merger_fil() { // имя функции = имя файла
    std::vector<std::string> names {"31.0.root", "33.1.root", "27.12.root", "27.15.root", "37.16.root", "23.0.root", "37.7.root", "27.0.root", "33.6.root", "31.12.root", "37.1.root", "35.1.root", "21.1.root", "27.1.root", "23.4.root", "39.6.root", "31.2.root", "23.5.root", "27.3.root", "25.2.root", "27.7.root", "31.10.root", "27.2.root", "39.10.root", "31.11.root", "31.15.root", "37.6.root", "27.11.root", "27.16.root", "39.3.root", "25.0.root", "39.7.root", "39.0.root", "31.9.root", "27.13.root", "21.3.root", "27.6.root", "31.4.root", "23.3.root", "37.11.root", "35.6.root", "37.18.root", "39.9.root", "31.8.root", "27.10.root", "31.6.root", "25.9.root", "35.4.root", "27.9.root", "25.3.root", "31.7.root", "33.4.root", "37.14.root", "35.2.root", "39.1.root", "25.12.root", "31.1.root", "31.3.root", "37.5.root", "21.0.root", "25.13.root", "37.3.root", "27.8.root", "35.0.root", "25.6.root", "33.0.root", "37.8.root", "25.16.root", "25.11.root", "39.11.root", "25.14.root", "39.5.root", "37.12.root", "37.17.root", "37.4.root", "31.16.root", "39.2.root", "33.8.root", "23.2.root", "37.10.root", "25.18.root", "37.19.root", "31.17.root", "37.13.root", "25.4.root", "25.15.root", "39.13.root", "27.4.root", "23.1.root", "33.3.root", "37.9.root", "37.0.root", "39.12.root", "39.8.root", "25.17.root", "25.8.root", "23.6.root", "25.19.root", "35.3.root", "25.5.root", "21.2.root", "31.5.root", "37.2.root", "25.20.root", "33.5.root", "31.14.root", "31.13.root", "33.2.root", "39.4.root", "37.15.root", "25.7.root", "25.21.root", "27.14.root", "25.10.root", "33.7.root", "25.1.root", "27.5.root", "35.5.root"};
	TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
	TH1F *hist = new TH1F("hist", "", 100, 8, 12);
	
	for (auto iter {names.begin()}; iter != names.end(); ++iter) {
		TFile *input = new TFile((*iter).c_str(), "read");
		TTree *tree = (TTree*)input->Get("h1");
		
		TH1F *temp = new TH1F("temp", "", 100, 8, 12);
		//tree->Draw("en >> temp", "p < 0.5 && en > 8 && chu == 5 && abs(ml - 2.28646) < 0.015 && abs(mach - 2.28646) < 0.015 && ntr == 0");
		//tree->Draw("en >> temp", "p < 0.5 && en > 8 && chu == 5 && abs(ml - 2.28646) < 0.015 && abs(mach - 1.86484) < 0.015 && ntr == 0");
		tree->Draw("en >> temp", "p < 0.5 && en > 8 && chu == 7 && abs(ml - 2.28646) < 0.015 && abs(mach - 1.86966) < 0.015 && ntr == 0");
		//tree->Draw("en >> temp", "chu == 5");
		
		hist->Add(temp, 1); // добавляем данные temp с коэффициентом 1
	}
	
	//hist->GetYaxis()->SetRangeUser(0, 100);
	hist->Draw();
}
