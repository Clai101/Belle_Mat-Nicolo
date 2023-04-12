std::string name_generator(const std::string& name){
	return name;
}

const int Length = 2; // колво файлов
const double massive[Length] = {0, 2}; // номера файлов

void merger() {
	TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
	TH1F *hist = new TH1F("hist", "", 100, 8, 12);
	
	for (int i = 0; i < Length; i++) {
		int number = massive[i];
		TFile *input = new TFile(name_generator("21."+std::to_string(number)+".root").c_str(), "read");
		TTree *tree = (TTree*)input->Get("h1");
		
		TH1F *temp = new TH1F("temp", "", 100, 8, 12);
		tree->Draw("en >> temp", "p > 1");
		
		hist->Add(temp, 1); // добавляем данные temp с коэффициентом 1
	}
	
	hist->GetYaxis()->SetRangeUser(0, 6000);
	hist->Draw();
}
