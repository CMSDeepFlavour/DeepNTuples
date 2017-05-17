/*
 * Macro to convert subjets based trees to AK8 fatjet based trees.
 *
 *  Created on: May 5, 2017
 *      Author: Loukas Gouskos
 */

struct ak8jetstruct {

	UInt_t event_no;

	int    fj_isLight;
	int    fj_isW;
	int    fj_isTop;
	float  fj_gen_pt;
	float  fj_gen_eta;
	float  fj_pt;
	float  fj_eta;
	float  fj_phi;
	float  fj_mass;
	float  fj_sdmass;
	float  fj_tau1;
	float  fj_tau2;
	float  fj_tau3;
	float  fj_tau21;
	float  fj_tau32;
	float  fj_ptDR;
	float  fj_relptdiff;
	float  fj_sdn2;

	float  fj_sdsj1_pt;
	float  fj_sdsj1_eta;
	float  fj_sdsj1_phi;
	float  fj_sdsj1_mass;
	float  fj_sdsj1_csv;
	float  fj_sdsj1_ptD;
	float  fj_sdsj1_axis1;
	float  fj_sdsj1_axis2;
	float  fj_sdsj1_mult;
	float  fj_sdsj1_prob_isB;
	float  fj_sdsj1_prob_isBB;
	float  fj_sdsj1_prob_isC;
	float  fj_sdsj1_prob_isUDSG;

	float  fj_sdsj2_pt;
	float  fj_sdsj2_eta;
	float  fj_sdsj2_phi;
	float  fj_sdsj2_mass;
	float  fj_sdsj2_csv;
	float  fj_sdsj2_ptD;
	float  fj_sdsj2_axis1;
	float  fj_sdsj2_axis2;
	float  fj_sdsj2_mult;
	float  fj_sdsj2_prob_isB;
	float  fj_sdsj2_prob_isBB;
	float  fj_sdsj2_prob_isC;
	float  fj_sdsj2_prob_isUDSG;

};


void readEvent(TString fname1, TString fname2, TString outname) {

	TFile *f = new TFile(fname1);
	TTree *t1 = (TTree*)f->Get("deepntuplizer/tree");

	TFile *f_s = new TFile(fname2);
	TTree *t_s = (TTree*)f_s->Get("tree");

	std::vector<ak8jetstruct> ak8jet_final;
	ak8jet_final.clear();


	UInt_t event_no;

	int    fj_isLight, fj_isW, fj_isTop;
	float  fj_gen_pt, fj_gen_eta;
	float  fj_pt, fj_eta, fj_phi, fj_mass, fj_sdmass;
	float  fj_tau1, fj_tau2, fj_tau3, fj_tau21, fj_tau32;
	float  fj_ptDR, fj_relptdiff, fj_sdn2;

	float  fj_sdsj1_pt, fj_sdsj1_eta, fj_sdsj1_phi, fj_sdsj1_mass;
	float  fj_sdsj1_csv, fj_sdsj1_ptD, fj_sdsj1_axis1, fj_sdsj1_axis2, fj_sdsj1_mult;

	float  fj_sdsj2_pt, fj_sdsj2_eta, fj_sdsj2_phi, fj_sdsj2_mass;
	float  fj_sdsj2_csv, fj_sdsj2_ptD, fj_sdsj2_axis1, fj_sdsj2_axis2, fj_sdsj2_mult;

	float  jet_corr_pt,jet_corr_eta,jet_corr_phi;
	float  prob_isB, prob_isBB, prob_isC, prob_isUDSG;


	t1->SetBranchAddress("event_no"    ,&event_no     );
	t1->SetBranchAddress("jet_corr_pt" ,&jet_corr_pt  );
	t1->SetBranchAddress("jet_eta"     ,&jet_corr_eta );

	t1->SetBranchAddress("fj_isLight"   , &fj_isLight   );
	t1->SetBranchAddress("fj_isW"       , &fj_isW       );
	t1->SetBranchAddress("fj_isTop"     , &fj_isTop     );
	t1->SetBranchAddress("fj_gen_pt"    , &fj_gen_pt    );
	t1->SetBranchAddress("fj_gen_eta"   , &fj_gen_eta   );
	t1->SetBranchAddress("fj_pt"        , &fj_pt        );
	t1->SetBranchAddress("fj_eta"       , &fj_eta       );
	t1->SetBranchAddress("fj_phi"       , &fj_phi       );
	t1->SetBranchAddress("fj_mass"      , &fj_mass      );
	t1->SetBranchAddress("fj_sdmass"    , &fj_sdmass    );
	t1->SetBranchAddress("fj_tau1"      , &fj_tau1      );
	t1->SetBranchAddress("fj_tau2"      , &fj_tau2      );
	t1->SetBranchAddress("fj_tau3"      , &fj_tau3      );
	t1->SetBranchAddress("fj_tau21"     , &fj_tau21     );
	t1->SetBranchAddress("fj_tau32"     , &fj_tau32     );
	t1->SetBranchAddress("fj_ptDR"      , &fj_ptDR      );
	t1->SetBranchAddress("fj_relptdiff" , &fj_relptdiff );
	t1->SetBranchAddress("fj_sdn2"      , &fj_sdn2      );

	t1->SetBranchAddress("fj_sdsj1_pt"    ,&fj_sdsj1_pt    );
	t1->SetBranchAddress("fj_sdsj1_eta"   ,&fj_sdsj1_eta   );
	t1->SetBranchAddress("fj_sdsj1_phi"   ,&fj_sdsj1_phi   );
	t1->SetBranchAddress("fj_sdsj1_mass"  ,&fj_sdsj1_mass  );
	t1->SetBranchAddress("fj_sdsj1_csv"   ,&fj_sdsj1_csv   );
	t1->SetBranchAddress("fj_sdsj1_ptD"   ,&fj_sdsj1_ptD   );
	t1->SetBranchAddress("fj_sdsj1_axis1" ,&fj_sdsj1_axis1 );
	t1->SetBranchAddress("fj_sdsj1_axis2" ,&fj_sdsj1_axis2 );
	t1->SetBranchAddress("fj_sdsj1_mult"  ,&fj_sdsj1_mult  );

	t1->SetBranchAddress("fj_sdsj2_pt"    ,&fj_sdsj2_pt    );
	t1->SetBranchAddress("fj_sdsj2_eta"   ,&fj_sdsj2_eta   );
	t1->SetBranchAddress("fj_sdsj2_phi"   ,&fj_sdsj2_phi   );
	t1->SetBranchAddress("fj_sdsj2_mass"  ,&fj_sdsj2_mass  );
	t1->SetBranchAddress("fj_sdsj2_csv"   ,&fj_sdsj2_csv   );
	t1->SetBranchAddress("fj_sdsj2_ptD"   ,&fj_sdsj2_ptD   );
	t1->SetBranchAddress("fj_sdsj2_axis1" ,&fj_sdsj2_axis1 );
	t1->SetBranchAddress("fj_sdsj2_axis2" ,&fj_sdsj2_axis2 );
	t1->SetBranchAddress("fj_sdsj2_mult"  ,&fj_sdsj2_mult  );

	t_s->SetBranchAddress("prob_isB"    ,&prob_isB    );
	t_s->SetBranchAddress("prob_isBB"   ,&prob_isBB   );
	t_s->SetBranchAddress("prob_isC"    ,&prob_isC    );
	t_s->SetBranchAddress("prob_isUDSG" ,&prob_isUDSG );


	//loop over all entries in the tree / each entry corresponds to a single subjet
	Int_t nentries = (Int_t)t1->GetEntries();
	for (Int_t i0=0; i0<nentries; i0++) {

		if(i0%1000==0){
			cout << i0 << "/" << nentries << endl;
		}

		t1->GetEntry(i0);

		int indx0 = i0;
		int indx1 = 0;

		float jet_corr_eta_0 = jet_corr_eta;
		float jet_corr_eta_1 = 0.;

		//create a temp ak8jet
		ak8jetstruct ak8jet_tmp;

		ak8jet_tmp.event_no     = event_no;
		ak8jet_tmp.fj_isLight   = fj_isLight;
		ak8jet_tmp.fj_isW       = fj_isW;
		ak8jet_tmp.fj_isTop     = fj_isTop;
		ak8jet_tmp.fj_gen_pt    = fj_gen_pt;
		ak8jet_tmp.fj_gen_eta   = fj_gen_eta;
		ak8jet_tmp.fj_pt        = fj_pt;
		ak8jet_tmp.fj_eta       = fj_eta;
		ak8jet_tmp.fj_phi       = fj_phi;
		ak8jet_tmp.fj_mass      = fj_mass;
		ak8jet_tmp.fj_sdmass    = fj_sdmass;
		ak8jet_tmp.fj_tau1      = fj_tau1;
		ak8jet_tmp.fj_tau2      = fj_tau2;
		ak8jet_tmp.fj_tau3      = fj_tau3;
		ak8jet_tmp.fj_tau21     = fj_tau21;
		ak8jet_tmp.fj_tau32     = fj_tau32;
		ak8jet_tmp.fj_ptDR      = fj_ptDR;
		ak8jet_tmp.fj_relptdiff = fj_relptdiff;
		ak8jet_tmp.fj_sdn2      = fj_sdn2;

		ak8jet_tmp.fj_sdsj1_pt    = fj_sdsj1_pt;
		ak8jet_tmp.fj_sdsj1_eta   = fj_sdsj1_eta;
		ak8jet_tmp.fj_sdsj1_phi   = fj_sdsj1_phi;
		ak8jet_tmp.fj_sdsj1_mass  = fj_sdsj1_mass;
		ak8jet_tmp.fj_sdsj1_csv   = fj_sdsj1_csv;
		ak8jet_tmp.fj_sdsj1_ptD   = fj_sdsj1_ptD;
		ak8jet_tmp.fj_sdsj1_axis1 = fj_sdsj1_axis1;
		ak8jet_tmp.fj_sdsj1_axis2 = fj_sdsj1_axis2;
		ak8jet_tmp.fj_sdsj1_mult  = fj_sdsj1_mult;

		ak8jet_tmp.fj_sdsj2_pt  = 0.;

		// find the second subjet
		for (Int_t i1=i0+1; i1<nentries; i1++) {

			t1->GetEntry(i1);

			if (ak8jet_tmp.event_no!=event_no) {
				// no need to go further if another event begins
				// --> jets belonging to the same event are
				break;
			}

			// check if subjet-2 comes from the same ak8 jet as subjet-1
			if ( (ak8jet_tmp.event_no==event_no) &&
					(ak8jet_tmp.fj_pt==fj_pt)       &&
					(ak8jet_tmp.fj_eta==fj_eta)
			) {

				ak8jet_tmp.fj_sdsj2_pt    = fj_sdsj2_pt;
				ak8jet_tmp.fj_sdsj2_eta   = fj_sdsj2_eta;
				ak8jet_tmp.fj_sdsj2_phi   = fj_sdsj2_phi;
				ak8jet_tmp.fj_sdsj2_mass  = fj_sdsj2_mass;
				ak8jet_tmp.fj_sdsj2_csv   = fj_sdsj2_csv;
				ak8jet_tmp.fj_sdsj2_ptD   = fj_sdsj2_ptD;
				ak8jet_tmp.fj_sdsj2_axis1 = fj_sdsj2_axis1;
				ak8jet_tmp.fj_sdsj2_axis2 = fj_sdsj2_axis2;
				ak8jet_tmp.fj_sdsj2_mult  = fj_sdsj2_mult;

				indx1 = i1;
				jet_corr_eta_1 = jet_corr_eta;

				break;
			} // end of storing second subjet info

		} // end of looping over the second subjet

		if (ak8jet_tmp.fj_sdsj2_pt==0.) { continue; }

		if (ak8jet_tmp.fj_sdsj1_eta==jet_corr_eta_0) {

			t_s->GetEntry(indx0);
			ak8jet_tmp.fj_sdsj1_prob_isB    = prob_isB;
			ak8jet_tmp.fj_sdsj1_prob_isBB   = prob_isBB;
			ak8jet_tmp.fj_sdsj1_prob_isC    = prob_isC;
			ak8jet_tmp.fj_sdsj1_prob_isUDSG = prob_isUDSG;

			t_s->GetEntry(indx1);
			ak8jet_tmp.fj_sdsj2_prob_isB    = prob_isB;
			ak8jet_tmp.fj_sdsj2_prob_isBB   = prob_isBB;
			ak8jet_tmp.fj_sdsj2_prob_isC    = prob_isC;
			ak8jet_tmp.fj_sdsj2_prob_isUDSG = prob_isUDSG;

		}

		else {

			t_s->GetEntry(indx1);
			ak8jet_tmp.fj_sdsj1_prob_isB    = prob_isB;
			ak8jet_tmp.fj_sdsj1_prob_isBB   = prob_isBB;
			ak8jet_tmp.fj_sdsj1_prob_isC    = prob_isC;
			ak8jet_tmp.fj_sdsj1_prob_isUDSG = prob_isUDSG;

			t_s->GetEntry(indx0);
			ak8jet_tmp.fj_sdsj2_prob_isB    = prob_isB;
			ak8jet_tmp.fj_sdsj2_prob_isBB   = prob_isBB;
			ak8jet_tmp.fj_sdsj2_prob_isC    = prob_isC;
			ak8jet_tmp.fj_sdsj2_prob_isUDSG = prob_isUDSG;

		}

		ak8jet_final.push_back(ak8jet_tmp);

	} // end of looping over all entries

	std::cout << "size = " << ak8jet_final.size() << "\n";


	// make the new tree: one entry for each ak8 jet
	TFile f_out(outname,"recreate");
	TTree t_out("Events","out tree");

	UInt_t out_event_no;

	int    out_fj_isLight, out_fj_isW, out_fj_isTop;
	float  out_fj_gen_pt, out_fj_gen_eta;
	float  out_fj_pt, out_fj_eta, out_fj_phi, out_fj_mass, out_fj_sdmass;
	float  out_fj_tau1, out_fj_tau2, out_fj_tau3, out_fj_tau21, out_fj_tau32;
	float  out_fj_ptDR, out_fj_relptdiff, out_fj_sdn2;

	float  out_fj_sdsj1_pt, out_fj_sdsj1_eta, out_fj_sdsj1_phi, out_fj_sdsj1_mass;
	float  out_fj_sdsj1_csv, out_fj_sdsj1_ptD, out_fj_sdsj1_axis1, out_fj_sdsj1_axis2, out_fj_sdsj1_mult;
	float  out_fj_sdsj1_prob_isB, out_fj_sdsj1_prob_isBB, out_fj_sdsj1_prob_isC, out_fj_sdsj1_prob_isUDSG;

	float  out_fj_sdsj2_pt, out_fj_sdsj2_eta, out_fj_sdsj2_phi, out_fj_sdsj2_mass;
	float  out_fj_sdsj2_csv, out_fj_sdsj2_ptD, out_fj_sdsj2_axis1, out_fj_sdsj2_axis2, out_fj_sdsj2_mult;
	float  out_fj_sdsj2_prob_isB, out_fj_sdsj2_prob_isBB, out_fj_sdsj2_prob_isC, out_fj_sdsj2_prob_isUDSG;




	t_out.Branch("event_no"    , &event_no);

	t_out.Branch("fj_isLight"   , &out_fj_isLight   );
	t_out.Branch("fj_isW"       , &out_fj_isW       );
	t_out.Branch("fj_isTop"     , &out_fj_isTop     );
	t_out.Branch("fj_gen_pt"    , &out_fj_gen_pt    );
	t_out.Branch("fj_gen_eta"   , &out_fj_gen_eta   );
	t_out.Branch("fj_pt"        , &out_fj_pt        );
	t_out.Branch("fj_eta"       , &out_fj_eta       );
	t_out.Branch("fj_phi"       , &out_fj_phi       );
	t_out.Branch("fj_mass"      , &out_fj_mass      );
	t_out.Branch("fj_sdmass"    , &out_fj_sdmass    );
	t_out.Branch("fj_tau1"      , &out_fj_tau1      );
	t_out.Branch("fj_tau2"      , &out_fj_tau2      );
	t_out.Branch("fj_tau3"      , &out_fj_tau3      );
	t_out.Branch("fj_tau21"     , &out_fj_tau21     );
	t_out.Branch("fj_tau32"     , &out_fj_tau32     );
	t_out.Branch("fj_ptDR"      , &out_fj_ptDR      );
	t_out.Branch("fj_relptdiff" , &out_fj_relptdiff );
	t_out.Branch("fj_sdn2"      , &out_fj_sdn2      );

	t_out.Branch("fj_sdsj1_pt"          , &out_fj_sdsj1_pt          );
	t_out.Branch("fj_sdsj1_eta"         , &out_fj_sdsj1_eta         );
	t_out.Branch("fj_sdsj1_phi"         , &out_fj_sdsj1_phi         );
	t_out.Branch("fj_sdsj1_mass"        , &out_fj_sdsj1_mass        );
	t_out.Branch("fj_sdsj1_csv"         , &out_fj_sdsj1_csv         );
	t_out.Branch("fj_sdsj1_ptD"         , &out_fj_sdsj1_ptD         );
	t_out.Branch("fj_sdsj1_axis1"       , &out_fj_sdsj1_axis1       );
	t_out.Branch("fj_sdsj1_axis2"       , &out_fj_sdsj1_axis2       );
	t_out.Branch("fj_sdsj1_mult"        , &out_fj_sdsj1_mult        );
	t_out.Branch("fj_sdsj1_prob_isB"    , &out_fj_sdsj1_prob_isB    );
	t_out.Branch("fj_sdsj1_prob_isBB"   , &out_fj_sdsj1_prob_isBB   );
	t_out.Branch("fj_sdsj1_prob_isC"    , &out_fj_sdsj1_prob_isC    );
	t_out.Branch("fj_sdsj1_prob_isUDSG" , &out_fj_sdsj1_prob_isUDSG );

	t_out.Branch("fj_sdsj2_pt"          , &out_fj_sdsj2_pt          );
	t_out.Branch("fj_sdsj2_eta"         , &out_fj_sdsj2_eta         );
	t_out.Branch("fj_sdsj2_phi"         , &out_fj_sdsj2_phi         );
	t_out.Branch("fj_sdsj2_mass"        , &out_fj_sdsj2_mass        );
	t_out.Branch("fj_sdsj2_csv"         , &out_fj_sdsj2_csv         );
	t_out.Branch("fj_sdsj2_ptD"         , &out_fj_sdsj2_ptD         );
	t_out.Branch("fj_sdsj2_axis1"       , &out_fj_sdsj2_axis1       );
	t_out.Branch("fj_sdsj2_axis2"       , &out_fj_sdsj2_axis2       );
	t_out.Branch("fj_sdsj2_mult"        , &out_fj_sdsj2_mult        );
	t_out.Branch("fj_sdsj2_prob_isB"    , &out_fj_sdsj2_prob_isB    );
	t_out.Branch("fj_sdsj2_prob_isBB"   , &out_fj_sdsj2_prob_isBB   );
	t_out.Branch("fj_sdsj2_prob_isC"    , &out_fj_sdsj2_prob_isC    );
	t_out.Branch("fj_sdsj2_prob_isUDSG" , &out_fj_sdsj2_prob_isUDSG );


	for (unsigned int i0=0; i0<ak8jet_final.size(); ++i0) {

		out_event_no       = ak8jet_final[i0].event_no;

		out_fj_isLight   = ak8jet_final[i0].fj_isLight;
		out_fj_isW       = ak8jet_final[i0].fj_isW;
		out_fj_isTop     = ak8jet_final[i0].fj_isTop;
		out_fj_gen_pt    = ak8jet_final[i0].fj_gen_pt;
		out_fj_gen_eta   = ak8jet_final[i0].fj_gen_eta;
		out_fj_pt        = ak8jet_final[i0].fj_pt;
		out_fj_eta       = ak8jet_final[i0].fj_eta;
		out_fj_phi       = ak8jet_final[i0].fj_phi;
		out_fj_mass      = ak8jet_final[i0].fj_mass;
		out_fj_sdmass    = ak8jet_final[i0].fj_sdmass;
		out_fj_tau1      = ak8jet_final[i0].fj_tau1;
		out_fj_tau2      = ak8jet_final[i0].fj_tau2;
		out_fj_tau3      = ak8jet_final[i0].fj_tau3;
		out_fj_tau21     = ak8jet_final[i0].fj_tau21;
		out_fj_tau32     = ak8jet_final[i0].fj_tau32;
		out_fj_ptDR      = ak8jet_final[i0].fj_ptDR;
		out_fj_relptdiff = ak8jet_final[i0].fj_relptdiff;
		out_fj_sdn2      = ak8jet_final[i0].fj_sdn2;

		out_fj_sdsj1_pt          = ak8jet_final[i0].fj_sdsj1_pt;
		out_fj_sdsj1_eta         = ak8jet_final[i0].fj_sdsj1_eta;
		out_fj_sdsj1_phi         = ak8jet_final[i0].fj_sdsj1_phi;
		out_fj_sdsj1_mass        = ak8jet_final[i0].fj_sdsj1_mass;
		out_fj_sdsj1_csv         = ak8jet_final[i0].fj_sdsj1_csv;
		out_fj_sdsj1_ptD         = ak8jet_final[i0].fj_sdsj1_ptD;
		out_fj_sdsj1_axis1       = ak8jet_final[i0].fj_sdsj1_axis1;
		out_fj_sdsj1_axis2       = ak8jet_final[i0].fj_sdsj1_axis2;
		out_fj_sdsj1_mult        = ak8jet_final[i0].fj_sdsj1_mult;
		out_fj_sdsj1_prob_isB    = ak8jet_final[i0].fj_sdsj1_prob_isB;
		out_fj_sdsj1_prob_isBB   = ak8jet_final[i0].fj_sdsj1_prob_isBB;
		out_fj_sdsj1_prob_isC    = ak8jet_final[i0].fj_sdsj1_prob_isC;
		out_fj_sdsj1_prob_isUDSG = ak8jet_final[i0].fj_sdsj1_prob_isUDSG;

		out_fj_sdsj2_pt          = ak8jet_final[i0].fj_sdsj2_pt;
		out_fj_sdsj2_eta         = ak8jet_final[i0].fj_sdsj2_eta;
		out_fj_sdsj2_phi         = ak8jet_final[i0].fj_sdsj2_phi;
		out_fj_sdsj2_mass        = ak8jet_final[i0].fj_sdsj2_mass;
		out_fj_sdsj2_csv         = ak8jet_final[i0].fj_sdsj2_csv;
		out_fj_sdsj2_ptD         = ak8jet_final[i0].fj_sdsj2_ptD;
		out_fj_sdsj2_axis1       = ak8jet_final[i0].fj_sdsj2_axis1;
		out_fj_sdsj2_axis2       = ak8jet_final[i0].fj_sdsj2_axis2;
		out_fj_sdsj2_mult        = ak8jet_final[i0].fj_sdsj2_mult;
		out_fj_sdsj2_prob_isB    = ak8jet_final[i0].fj_sdsj2_prob_isB;
		out_fj_sdsj2_prob_isBB   = ak8jet_final[i0].fj_sdsj2_prob_isBB;
		out_fj_sdsj2_prob_isC    = ak8jet_final[i0].fj_sdsj2_prob_isC;
		out_fj_sdsj2_prob_isUDSG = ak8jet_final[i0].fj_sdsj2_prob_isUDSG;

		t_out.Fill();
	}

	t_out.Write();


} // end of readEvent()

