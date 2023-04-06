void test() { // имя функции = имя файла
	TFile *input = new TFile("b3k_data.root", "read"); // input = твои данные
	cout << input;
	TTree *tree = (TTree*)input->Get("h1"); // берем дерево с данными из input
	TCanvas *c1 = new TCanvas("c1", "c1", 800, 600); // окно в котором будет график
	
	TH1F *hist = new TH1F("hist", "smth", 100, 0, 1); // создаем гистограмму, 100 бинов, 0 < x < 1
	
	hist->SetFillColor(kGreen-9); // всякие украшения для гистограммы
	hist->GetXaxis()->SetTitle("r2");
	hist->GetYaxis()->SetTitle("Entries / 0.01");
	hist->GetXaxis()->SetTitleSize(0.045);
	hist->GetYaxis()->SetTitleSize(0.045);
	
	
	tree->Draw("r2 >> hist", "chb == 1 && chrg > 0"); // записываем переменную r2 из дерева в hist, отбираем события где chb == 1 и chrg > 0, рисуем гистограмму
}
