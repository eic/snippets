{
	TFile*f = new TFile("TrackingBG_10x275_Forced.root");

	//E-Si Disk 4
	TH2D *hSVTDiskRate_0 = (TH2D*) f->Get("hSVTDiskRate_0");

	cout<<"For E-Si Disk 4:"<<endl;
	std::cout << std::format("Total Hits = {}, Max Hits = {} \n", hSVTDiskRate_0->Integral(), hSVTDiskRate_0->GetMaximum());	
	//E-Si Disk 3
        TH2D *hSVTDiskRate_1 = (TH2D*) f->Get("hSVTDiskRate_1");

        cout<<"For E-Si Disk 3:"<<endl;
        std::cout << std::format("Total Hits = {}, Max Hits = {} \n", hSVTDiskRate_1->Integral(), hSVTDiskRate_1->GetMaximum());

	//E-Si Disk 2
        TH2D *hSVTDiskRate_2 = (TH2D*) f->Get("hSVTDiskRate_2");

        cout<<"For E-Si Disk 2:"<<endl;
        std::cout << std::format("Total Hits = {}, Max Hits = {} \n", hSVTDiskRate_2->Integral(), hSVTDiskRate_2->GetMaximum());

	//E-Si Disk 1
        TH2D *hSVTDiskRate_3 = (TH2D*) f->Get("hSVTDiskRate_3");

        cout<<"For E-Si Disk 1:"<<endl;
        std::cout << std::format("Total Hits = {}, Max Hits = {} \n", hSVTDiskRate_3->Integral(), hSVTDiskRate_3->GetMaximum());

	//E-Si Disk 0
        TH2D *hSVTDiskRate_4 = (TH2D*) f->Get("hSVTDiskRate_4");

        cout<<"For E-Si Disk 0:"<<endl;
        std::cout << std::format("Total Hits = {}, Max Hits = {} \n", hSVTDiskRate_4->Integral(), hSVTDiskRate_4->GetMaximum());

	//H-Si Disk 0
	TH2D *hSVTDiskRate_5 = (TH2D*) f->Get("hSVTDiskRate_5");

        cout<<"For H-Si Disk 0:"<<endl;
        std::cout << std::format("Total Hits = {}, Max Hits = {} \n", hSVTDiskRate_5->Integral(), hSVTDiskRate_5->GetMaximum());

	//H-Si Disk 1
	TH2D *hSVTDiskRate_6 = (TH2D*) f->Get("hSVTDiskRate_6");

        cout<<"For H-Si Disk 1:"<<endl;
        std::cout << std::format("Total Hits = {}, Max Hits = {} \n", hSVTDiskRate_6->Integral(), hSVTDiskRate_6->GetMaximum());	
	//H-Si Disk 2
	TH2D *hSVTDiskRate_7 = (TH2D*) f->Get("hSVTDiskRate_7");

        cout<<"For H-Si Disk 2:"<<endl;
        std::cout << std::format("Total Hits = {}, Max Hits = {} \n", hSVTDiskRate_7->Integral(), hSVTDiskRate_7->GetMaximum());	
	//H-Si Disk 3
	TH2D *hSVTDiskRate_8 = (TH2D*) f->Get("hSVTDiskRate_8");

        cout<<"For H-Si Disk 3:"<<endl;
        std::cout << std::format("Total Hits = {}, Max Hits = {} \n", hSVTDiskRate_8->Integral(), hSVTDiskRate_8->GetMaximum());

	//H-Si Disk 4
	TH2D *hSVTDiskRate_9 = (TH2D*) f->Get("hSVTDiskRate_9");

        cout<<"For H-Si Disk 4:"<<endl;
        std::cout << std::format("Total Hits = {}, Max Hits = {} \n", hSVTDiskRate_9->Integral(), hSVTDiskRate_9->GetMaximum());

	//L0
	TH2D *hSVTVtxRate_0 = (TH2D*) f->Get("hSVTVtxRate_0");

        cout<<"For SVT L0:"<<endl;
        std::cout << std::format("Total Hits = {}, Max Hits = {} \n", hSVTVtxRate_0->Integral(), hSVTVtxRate_0->GetMaximum());

	//L1
	TH2D *hSVTVtxRate_1 = (TH2D*) f->Get("hSVTVtxRate_1");

        cout<<"For SVT L1:"<<endl;
        std::cout << std::format("Total Hits = {}, Max Hits = {} \n", hSVTVtxRate_1->Integral(), hSVTVtxRate_1->GetMaximum());

	//L2
	TH2D *hSVTVtxRate_2 = (TH2D*) f->Get("hSVTVtxRate_2");

        cout<<"For SVT L2:"<<endl;
        std::cout << std::format("Total Hits = {}, Max Hits = {} \n", hSVTVtxRate_2->Integral(), hSVTVtxRate_2->GetMaximum());

	//L3
	TH2D *hSVTbarrelRate_0 = (TH2D*) f->Get("hSVTbarrelRate_0");

        cout<<"For SVT L3:"<<endl;
        std::cout << std::format("Total Hits = {}, Max Hits = {} \n", hSVTbarrelRate_0->Integral(), hSVTbarrelRate_0->GetMaximum());
	
	//L4
	TH2D *hSVTbarrelRate_1 = (TH2D*) f->Get("hSVTbarrelRate_1");

        cout<<"For SVT L4:"<<endl;
        std::cout << std::format("Total Hits = {}, Max Hits = {} \n", hSVTbarrelRate_1->Integral(), hSVTbarrelRate_1->GetMaximum());

}


