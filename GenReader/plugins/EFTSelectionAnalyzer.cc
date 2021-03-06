// Created by Andrew Wightman

#include "EFTSelectionAnalyzer.h"

EFTSelectionAnalyzer::EFTSelectionAnalyzer(const edm::ParameterSet& iConfig)
{
    debug = iConfig.getParameter<bool> ("debug");
    iseft = iConfig.getParameter<bool> ("iseft");

    min_pt_jet  = iConfig.getParameter<double> ("min_pt_jet");
    min_pt_lep  = iConfig.getParameter<double> ("min_pt_lep");
    max_eta_jet = iConfig.getParameter<double> ("max_eta_jet");
    max_eta_lep = iConfig.getParameter<double> ("max_eta_lep");
    
    lheInfo_token_      = consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("LHEInfo"));
    genInfo_token_      = consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("GENInfo"));
    genParticles_token_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("GenParticles"));
    genJets_token_      = consumes<std::vector<reco::GenJet> >(iConfig.getParameter<edm::InputTag>("GenJets"));

    patElectrons_token_ = consumes<std::vector<pat::Electron> >(iConfig.getParameter<edm::InputTag>("PatElectrons"));
    patMuons_token_     = consumes<std::vector<pat::Muon> >(iConfig.getParameter<edm::InputTag>("PatMuons"));
    patJets_token_      = consumes<std::vector<pat::Jet> >(iConfig.getParameter<edm::InputTag>("PatJets"));
}

EFTSelectionAnalyzer::~EFTSelectionAnalyzer(){}

void EFTSelectionAnalyzer::beginJob()
{
    total_ls = 0;
    total_events = 0;
    total_orig_xsec = 0.;
    total_sm_xsec = 0.;

    //this->addLepCategory("fit_all_incl","All Events");
    this->addLepCategory("fit_2lss_incl","Gen: 2lss Events");
    this->addLepCategory("fit_2lss_ee","Gen: 2lss ee Events");
    this->addLepCategory("fit_2lss_emu","Gen: 2lss emu Events");
    this->addLepCategory("fit_2lss_mumu","Gen: 2lss mumu Events");
    this->addLepCategory("fit_pat_2lss_incl","Pat: 2lss Events");
    this->addLepCategory("fit_pat_2lss_ee","Pat: 2lss ee Events");
    this->addLepCategory("fit_pat_2lss_emu","Pat: 2lss emu Events");
    this->addLepCategory("fit_pat_2lss_mumu","Pat: 2lss mumu Events");

}

void EFTSelectionAnalyzer::endJob()
{
    std::string delim = " | ";

    std::cout << "Total events processed: " << total_events << std::endl;
    std::cout << "Total LS: " << total_ls << std::endl;
    std::cout << "                " 
        << std::setw(13) << "xsec / LS" << delim
        << std::setw(13) << "xsec / Events" << delim
        << std::setw(13) << "xsec / 500"
        << std::endl;
    std::cout << "SM Xsec (pb):   "
        << std::setw(13) << total_sm_xsec / double(total_ls) << delim
        << std::setw(13) << total_sm_xsec / double(total_events) << delim
        << std::setw(13) << total_sm_xsec / 500.
        << std::endl;
    std::cout << "Orig Xsec (pb): "
        << std::setw(13) << total_orig_xsec / double(total_ls) << delim
        << std::setw(13) << total_orig_xsec / double(total_events) << delim
        << std::setw(13) << total_orig_xsec / 500.
        << std::endl;

    WCPoint* rwgt_pt = new WCPoint("rwgt",0.0);

    std::cout << "Only positive pairs!" << std::endl;
    std::cout << "Is EFT fit? " << iseft << std::endl;
    for (uint i=0; i < this->lep_category_names_vect.size(); i++) {
        std::string lep_category_name = this->lep_category_names_vect.at(i);
        WCFit* lep_category_fit = this->lep_fits_dict[lep_category_name];
        std::cout << lep_category_fit->getTag() << " SM Eval : " << lep_category_fit->evalPoint(rwgt_pt) << " +- " << lep_category_fit->evalPointError(rwgt_pt) << std::endl;
        delete lep_category_fit;
    }

}

void EFTSelectionAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& evsetup)
{
    eventcount++;

    edm::Handle<LHEEventProduct> LHEInfo;
    edm::Handle<reco::GenParticleCollection> prunedParticles;
    edm::Handle<std::vector<reco::GenJet> > genJets;

    event.getByToken(lheInfo_token_,LHEInfo);
    event.getByToken(genParticles_token_,prunedParticles);
    event.getByToken(genJets_token_,genJets);

    edm::Handle<std::vector<pat::Electron> > patElectrons;
    edm::Handle<std::vector<pat::Muon> > patMuons;
    edm::Handle<std::vector<pat::Jet> > patJets;
    event.getByToken(patElectrons_token_,patElectrons);
    event.getByToken(patMuons_token_,patMuons);
    event.getByToken(patJets_token_,patJets);

    if (debug){    
        if (eventcount == 3) dumpJets(*genJets);
    }

    reco::GenParticleCollection gen_leptons = GetGenLeptons(*prunedParticles);
    std::vector<reco::GenJet> gen_jets = GetGenJets(*genJets);

    std::vector<reco::GenJet> gen_jets_cut        = MakeJetPtEtaCuts(gen_jets,30,2.5);
    reco::GenParticleCollection gen_electrons     = GetGenParticlesSubset(gen_leptons,11);
    reco::GenParticleCollection gen_electrons_cut = MakePtEtaCuts(gen_electrons, 15,2.5);
    reco::GenParticleCollection gen_muons         = GetGenParticlesSubset(gen_leptons,13);
    reco::GenParticleCollection gen_muons_cut     = MakePtEtaCuts(gen_muons, 10,2.5);

    reco::GenParticleCollection gen_leptons_cut = {};
    gen_leptons_cut.insert(gen_leptons_cut.end(),gen_electrons_cut.begin(),gen_electrons_cut.end());
    gen_leptons_cut.insert(gen_leptons_cut.end(),gen_muons_cut.begin(),gen_muons_cut.end());

    //////////// Reco stuff ////////////

    bool pat_2lss_event = true;
    //bool jets4 = true;

    std::vector<pat::Electron> pat_electrons = *patElectrons;
    std::vector<pat::Muon> pat_muons         = *patMuons;
    std::vector<pat::Jet> pat_jets           = *patJets;

    std::vector<pat::Electron> pat_electrons_cut = MakePtEtaCutsPatElectrons(pat_electrons, 15, 2.5);
    std::vector<pat::Muon> pat_muons_cut         = MakePtEtaCutsPatMuons(pat_muons, 10, 2.5);
    std::vector<pat::Jet> pat_jets_cut           = MakePatJetPtEtaCuts(pat_jets,30,2.5);

    int pat_charge_sum_e = 0;
    for (size_t i = 0; i < pat_electrons_cut.size(); i++) {
        const pat::Electron& pat_e = pat_electrons_cut.at(i);
        pat_charge_sum_e += pat_e.charge();
    }
    int pat_charge_sum_mu = 0;
    for (size_t i = 0; i < pat_muons_cut.size(); i++) {
        const pat::Muon& pat_mu = pat_muons_cut.at(i);
        pat_charge_sum_mu += pat_mu.charge();
    }
    //if ( (pat_charge_sum_e <= 0) or (pat_charge_sum_mu <= 0) ){ // Only positive pairs
    if (pat_charge_sum_e + pat_charge_sum_mu <= 0){ // Only positive pairs
        pat_2lss_event = false;
    }
    if (pat_electrons_cut.size() + pat_muons_cut.size() != 2) {
        //std::cout << n_pat_leptons << std::endl;
        pat_2lss_event = false;
    } 
    if (pat_jets_cut.size() < 4) { // Require at least 4 jets
        pat_2lss_event = false;
        //jets4 = false;
    }

    //////////// Gen stuff ////////////

    bool gen_2lss_event = true;

    // Sum the charges
    int charge_sum_lep = 0;
    for (size_t i = 0; i < gen_leptons_cut.size(); i++) {
        const reco::GenParticle& p_lep = gen_leptons_cut.at(i);
        charge_sum_lep += p_lep.charge();
    }

    // We are only interested in 2lepss events, so skip otherwise: 
    if (gen_leptons_cut.size() !=2) {
        //return;
        gen_2lss_event = false;
    } else if (charge_sum_lep <= 0) { // Only positive
        //return;
        gen_2lss_event = false;
    }
    if (gen_jets_cut.size() < 4){ // Require at least 4 jets
        //return;
        gen_2lss_event = false;
    }

    // Skip event if it is neither a 2lss gen event nor a 2lss gen event
    //std::cout << "pat 2lss event: " << pat_2lss_event << " , gen 2lss event: " << gen_2lss_event << std::endl;
    if ((pat_2lss_event == false) and (gen_2lss_event == false)){
        return;
        //std::cout << "NOT a 2lss event !!!" << std::endl;
    }
    //std:: cout << "4 jets? " << jets4 << std::endl;

    double sm_wgt = 0.;
    double orig_wgt = LHEInfo->originalXWGTUP();  // original cross-section
    std::vector<WCPoint> wc_pts;
    if (iseft) {// Add EFT weights 
        for (auto wgt_info: LHEInfo->weights()) {
            auto LHEwgtstr = std::string(wgt_info.id);
            std::size_t foundstr = LHEwgtstr.find("EFTrwgt"); // only save our EFT weights
            if (foundstr!=std::string::npos) {
                WCPoint wc_pt(wgt_info.id,wgt_info.wgt);
                wc_pts.push_back(wc_pt);
                if (wc_pt.isSMPoint()) {
                    sm_wgt = wgt_info.wgt;
                }
            }
        }
    } else {
        sm_wgt = orig_wgt;
        WCPoint wc_pt("smpt",sm_wgt);
        wc_pts.push_back(wc_pt);
    }

    WCFit eft_fit(wc_pts,"");

    total_sm_xsec += sm_wgt;
    total_orig_xsec += orig_wgt;

    //std::cout << "E + MU SUM: " << charge_sum_e << " + " << charge_sum_mu << " = " << charge_sum_e+charge_sum_mu << std::endl;
    //std::cout << "CHARGE SUM: " << charge_sum_lep << std::endl;
    //std::cout << "n gen leptons: " << gen_leptons_cut.size() << std::endl;
    //std::cout << "n gen e: " << gen_electrons_cut.size() << std::endl;
    //std::cout << "n gen mu: " << gen_muons_cut.size() << std::endl;

    if (gen_2lss_event){
        // Add the fits for the different categories
        if (gen_leptons_cut.size() == 2 and charge_sum_lep !=0) {
            lep_fits_dict["fit_2lss_incl"]->addFit(eft_fit);
            if (gen_electrons_cut.size() == 2 and charge_sum_lep !=0){
                lep_fits_dict["fit_2lss_ee"]->addFit(eft_fit);
            }
            if (gen_muons_cut.size() == 2 and charge_sum_lep !=0){
                lep_fits_dict["fit_2lss_mumu"]->addFit(eft_fit);
            }
            if (gen_electrons_cut.size() == 1 and gen_muons_cut.size() == 1 and charge_sum_lep) {
                lep_fits_dict["fit_2lss_emu"]->addFit(eft_fit);
            }
        }
    }

    //if (pat_2lss_event) {
    //std::cout << "Number of pat e: " << pat_electrons_cut.size() << " Number of pat mu: " <<  pat_muons_cut.size() << " Tot: " << pat_electrons_cut.size() + pat_muons_cut.size() << std::endl;
    //std::cout << "Pat ch sum e: " << pat_charge_sum_e << std::endl;
    //std::cout << "Pat ch sum mu: " << pat_charge_sum_mu << std::endl;
    //}
    if (pat_2lss_event) {
    //if (1==1) {
        lep_fits_dict["fit_pat_2lss_incl"]->addFit(eft_fit);
        if (pat_electrons_cut.size() == 2){
            lep_fits_dict["fit_pat_2lss_ee"]->addFit(eft_fit);
            //std::cout << "Adding to fit_pat_2lss_ee!" << std::endl;
        }
        if (pat_muons_cut.size() == 2){
            lep_fits_dict["fit_pat_2lss_mumu"]->addFit(eft_fit);
            //std::cout << "Adding to fit_pat_2lss_mumu!" << std::endl;
        }
        if (pat_electrons_cut.size() == 1 and pat_muons_cut.size() == 1) {
            lep_fits_dict["fit_pat_2lss_emu"]->addFit(eft_fit);
            //std::cout << "Adding to fit_pat_2lss_emu!" << std::endl;
        }
    }
    //std::cout << "\n" << std::endl;

}

void EFTSelectionAnalyzer::beginRun(edm::Run const& run, edm::EventSetup const& evsetup)
{
    eventcount = 0;
}

void EFTSelectionAnalyzer::endRun(edm::Run const& run, edm::EventSetup const& evsetup)
{
    total_events += eventcount;
}

void EFTSelectionAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& evsetup)
{
    total_ls += 1;
}

void EFTSelectionAnalyzer::endLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& evsetup){}

DEFINE_FWK_MODULE(EFTSelectionAnalyzer);

