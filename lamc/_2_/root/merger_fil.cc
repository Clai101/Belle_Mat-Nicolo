
void merger_fil() { // имя функции = имя файла
	TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
	TH1F *hist = new TH1F("hist", "", 100, 8, 12);

    for ( const auto& dir : std::filesystem::directory_iterator( std::filesystem::current_path() ) ) {
            if ( dir.path( ).filename( ).extension( ) == ".cmd" ) {
                STARTUPINFO sinfo = { sizeof( sinfo ) };
                PROCESS_INFORMATION info;
                CreateProcessW( dir.path( ).wstring().c_str(), NULL, NULL, NULL, true, NULL, NULL, NULL, &sinfo, &info );
                
                TFile *input = new TFile((*iter).c_str(), "read");
                TTree *tree = (TTree*)input->Get("h1");
                
                TH1F *temp = new TH1F("temp", "", 100, 8, 12);
                tree->Draw("en >> temp", "p < 0.5 && en > 8 && ch == 2 && abs(ml - 2.28646) < 0.015 && abs(md - 2.28646) < 0.015 && ntr == 0");
                
                hist->Add(temp, 1); // добавляем данные temp с коэффициентом 1
            }
        }
	
	
	//hist->GetYaxis()->SetRangeUser(0, 100);
	hist->Draw();
}
