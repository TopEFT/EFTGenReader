#include "EFTGenReader/EFTHelperUtilities/interface/TH2EFT.h"

TH2EFT::TH2EFT() {}
TH2EFT::~TH2EFT() {}


TH2EFT::TH2EFT(const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup)
 : TH2D (name, title, nbinsx, xlow, xup, nbinsy, ylow, yup) 
{
    // Create/Initialize a fit function for each bin in the histogram
    WCFit new_fit;
    for (Int_t i = 0; i < this->FindFixBin(xup,yup); i++) {
        this->hist_fits.push_back(new_fit);
    }
}
TH2EFT::TH2EFT(const char *name, const char *title, Int_t nbinsx, Double_t *x, Int_t nbinsy, Double_t* y)
 : TH2D (name, title, nbinsx, x, nbinsy, y) 
{
    // Create/Initialize a fit function for each bin in the histogram
    WCFit new_fit;
    for (Int_t i = 0; i < this->FindFixBin(x[nbinsx-1],y[nbinsy-1]); i++) {
        this->hist_fits.push_back(new_fit);
    }
}
void TH2EFT::SetBins(Int_t nx, Double_t xmin, Double_t xmax, Int_t ny, Double_t ymin, Double_t ymax)
{
    // Use this function with care! Non-over/underflow bins are simply
    // erased and hist_fits re-sized with empty fits.
    
    hist_fits.clear();
    WCFit new_fit;
    for (Int_t i = 0; i < nx+ny; i++) {
        this->hist_fits.push_back(new_fit);
    }
    
    TH2::SetBins(nx, xmin, xmax, ny, ymin, ymax);
}

// Note: Since Clone calls Copy, this should make Clone work as well
void TH2EFT::Copy(TObject &obj) const
{
    TH2::Copy(obj);
    for (unsigned int i = 0; i < this->hist_fits.size(); i++) {
        WCFit bin_fit;
        bin_fit.addFit(this->hist_fits.at(i));
        ((TH2EFT&)obj).hist_fits.push_back(bin_fit);
    }
    WCFit of_fit;
    WCFit uf_fit;
    of_fit.addFit(this->overflow_fit);
    uf_fit.addFit(this->underflow_fit);
    ((TH2EFT&)obj).overflow_fit = of_fit;
    ((TH2EFT&)obj).underflow_fit = uf_fit;
}

Bool_t TH2EFT::Add(const TH2 *h1, Double_t c1)
{
    // check whether the object pointed to inherits from (or is a) TH2EFT:
    if (h1->IsA()->InheritsFrom(TH2EFT::Class())) {
        if ((this->hist_fits.size() == ((TH2EFT*)h1)->hist_fits.size())) {
            for (unsigned int i = 0; i < this->hist_fits.size(); i++) {
                // assumes this hist and the one whose fits we're adding have the same bins!
                this->hist_fits[i].addFit( ((TH2EFT*)h1)->hist_fits[i] );
            }
        } else { 
            std::cout << "Attempt to add 2 TH2EFTs with different # of fits!" << std::endl;
            std::cout << this->hist_fits.size() << ", " << ((TH2EFT*)h1)->hist_fits.size() << std::endl;
        }
        this->overflow_fit.addFit( ((TH2EFT*)h1)->overflow_fit );
        this->underflow_fit.addFit( ((TH2EFT*)h1)->underflow_fit );
    }
    
    return TH2::Add(h1,c1); // I think this should work
}

Bool_t TH2EFT::NormalizeTo(const TH2D *h1, Double_t c1)
{
    // check whether the object pointed to inherits from (or is a) TH2EFT:
    if (h1->IsA()->InheritsFrom(TH2EFT::Class())) {
        if ((this->hist_fits.size() == ((TH2EFT*)h1)->hist_fits.size()) && (this->GetNbinsX() == h1->GetNbinsX()) && (this->GetNbinsY() == h1->GetNbinsY())) {
            //Do nothing, just check
        } else { 
            std::cout << "Attempt to add 2 TH2EFTs with different # of fits!" << std::endl;
            std::cout << this->hist_fits.size() << ", " << ((TH2EFT*)h1)->hist_fits.size() << std::endl;
            return false;
        }
    }

    //Loop over all bins, and divide by h1->GetBinContent ^ c1 (sqrt by default)
    for(int i = 1; i <= this->GetNbinsX(); i++) {
        for(int j = 1; j <= this->GetNbinsY(); j++) {
            int bin = j*(this->GetNbinsX()+2) + i;
            //int bin = this->FindBin(i,j); //Find the corresonding bin
            //int bin = i + (this->GetNbinsY()+2) * j;
            //Get bin contents and errors from this and h1
            double thisbin = this->GetBinContent(bin);
            double thisbinerror = this->GetBinError(bin);
            double h1bin = h1->GetBinContent(bin);
            double h1binerror = h1->GetBinError(bin);
            //h1binerror = (thisbin - h1bin) / pow(h1bin, c1) * sqrt((h1binerror/h1bin)*(h1binerror/h1bin) + (thisbinerror/thisbin)*(thisbinerror/thisbin));
            //h1binerror = (thisbin - h1bin) / pow(h1bin, c1) * sqrt( pow(1/(thisbin - h1bin),2) * (thisbinerror*thisbinerror + h1binerror*h1binerror) + 1/4 * pow(h1binerror/h1bin, 2) ); //dependent, sigmaD^2 = fom^2 * ( sigmaE^2 + sigmaS^2 )
            double fom = (thisbin - h1bin) / pow(h1bin, c1);
            double efterrsq = thisbinerror*thisbinerror;
            double smerrsq = h1binerror*h1binerror;
            double diff = abs(thisbin - h1bin);
            double diffsq = diff * diff;
            //int bval = 36;
            int bval = 13;
            std::cout << "FoM=" << fom << " for (" << this->GetXaxis()->GetBinLowEdge(i) << "," << this->GetYaxis()->GetBinLowEdge(j) << ": " << bin << ") " << thisbin << " " << h1bin << std::endl;
            if(bin == bval) std::cout << "FoM=" << fom << " for " << thisbin << " " << h1bin << std::endl;
            if(bin == bval) std::cout << fom << "\t" << h1binerror << "\t" << diff << "\t" << diffsq << std::endl;
            //h1binerror = fom * sqrt( abs(efterrsq - smerrsq) / diffsq + 1/4 * smerrsq/(h1bin*h1bin));
            if(bin == bval) std::cout << abs(efterrsq - smerrsq) << std::endl;
            if(bin == bval) std::cout << abs(efterrsq - smerrsq) / diffsq << std::endl;
            if(bin == bval) std::cout << smerrsq/(h1bin*h1bin) << std::endl;
            if(bin == bval) std::cout << abs(efterrsq - smerrsq) / diffsq + 1/4 * smerrsq/(h1bin*h1bin) << std::endl;
            if(bin == bval) std::cout << sqrt(abs(efterrsq - smerrsq) / diffsq + 1/4 * smerrsq/(h1bin*h1bin)) << std::endl;
            if(bin == bval) std::cout << fom*sqrt(abs(efterrsq - smerrsq) / diffsq + 1/4 * smerrsq/(h1bin*h1bin)) << std::endl;
            h1binerror = fom * sqrt( abs(efterrsq - smerrsq) / diffsq + 1/4 * smerrsq/(h1bin*h1bin) );
            //h1binerror = (thisbin - h1bin) / pow(h1bin, c1) * sqrt( (thisbinerror*thisbinerror + h1binerror*h1binerror) / ((thisbin - h1bin)*(thisbin - h1bin)) + 1/4 * h1binerror*h1binerror / (h1bin*h1bin) );
            //h1binerror = sqrt(pow(thisbinerror, 2)/h1bin + ( 1/4 * pow(thisbin + h1bin, 2) / pow(h1bin, 3)) * pow(h1binerror, 3)); //independent
            //h1binerror = fom * sqrt( abs(thisbinerror*thisbinerror - h1binerror*h1binerror) / ( (thisbin - h1bin)*(thisbin - h1bin) ) + 1/4 * (h1binerror*h1binerror) / (h1bin*h1bin) );
            this->SetBinContent(bin, (thisbin - h1bin) / pow(h1bin, c1)); //Set this bin to (this - h1) / pow(h1, c1) ((this - h1)) / sqrt(h1) by default)
            this->SetBinError(bin, h1binerror); //Set this bin error to quad sum
        }
    }
    return true;
}

// Custom merge function for using hadd
Long64_t TH2EFT::Merge(TCollection* list)
{
    TIter nexthist(list);
    TH2EFT *hist;
    while ((hist = (TH2EFT*)nexthist.Next())) {
        if (this->hist_fits.size() != hist->hist_fits.size()) {
            std::cout << "[WARNING] Skipping histogram with different # of fits" << std::endl;
            continue;
        }
        for (unsigned int i = 0; i < this->hist_fits.size(); i++) {
            this->hist_fits.at(i).addFit(hist->hist_fits.at(i));
        }
        this->overflow_fit.addFit(hist->overflow_fit);
        this->underflow_fit.addFit(hist->underflow_fit);
    }

    return TH2::Merge(list);
}

Int_t TH2EFT::Fill(Double_t x, Double_t y, Double_t w, WCFit fit)
{
    Int_t bin_id = this->FindFixBin(x,y) - 1;
    Int_t nhists  = this->hist_fits.size();
    if (bin_id >= nhists) {
        // For now ignore events which enter overflow bin
        this->overflow_fit.addFit(fit);
        return Fill(x,y,w);
    } else if (bin_id < 0) {
        // For now ignore events which enter underflow bin
        this->underflow_fit.addFit(fit);
        return Fill(x,y,w);
    }
    this->hist_fits.at(bin_id).addFit(fit);
    return Fill(x,y,w); // the original TH2D member function
}

// Returns a fit function for a particular bin (no checks are made if the bin is an over/underflow bin)
WCFit TH2EFT::GetBinFit(Int_t bin)
{
    Int_t nhists = this->hist_fits.size();
    if (bin <= 0) {
        return this->underflow_fit;
    } else if (bin > nhists) {
        return this->overflow_fit;
    }
    else return this->hist_fits.at(bin - 1);
}

// Returns a WCFit whose structure constants are determined by summing structure constants from all bins
WCFit TH2EFT::GetSumFit(Int_t axis=0)
{
    WCFit summed_fit;
    for (unsigned int i = 0; i < this->hist_fits.size(); i++) {
        summed_fit.addFit(this->hist_fits.at(i));
    }
    return summed_fit;
}

// Returns a bin scaled by the the corresponding fit evaluated at a particular WC point
Double_t TH2EFT::GetBinContent(Int_t binx, Int_t biny, WCPoint wc_pt)
{
    Int_t bin = this->FindFixBin(binx,biny);
    if (this->GetBinFit(bin).getDim() <= 0) {
        // We don't have a fit for this bin, return regular bin contents
        //return GetBinContent(binx, biny);
        return GetBinContent(bin);
    }

    double scale_value = this->GetBinFit(bin).evalPoint(&wc_pt);
    Double_t num_events = GetBinContent(bin);
    if (num_events == 0) {
        return 0.0;
    }

    return scale_value;
}

Double_t TH2EFT::GetBinContent(Int_t bin, WCPoint wc_pt)
{
    if (this->GetBinFit(bin).getDim() <= 0) {
        // We don't have a fit for this bin, return regular bin contents
        //return GetBinContent(binx, biny);
        return GetBinContent(bin);
    }

    double scale_value = this->GetBinFit(bin).evalPoint(&wc_pt);
    Double_t num_events = GetBinContent(bin);
    if (num_events == 0) {
        return 0.0;
    }

    return scale_value;
}

void TH2EFT::Scale(WCPoint wc_pt)
{
    // Warning: calling GetEntries after a call to this function will return a 
    // non-zero value, even if the histogram was never filled.
    
    for (Int_t i = 1; i <= this->GetNbinsX(); i++) {
        for (Int_t j = 1; j <= this->GetNbinsY(); j++) {
            Int_t bin = this->FindFixBin(i-1,j-1);
            Double_t new_content = this->GetBinContent(bin,wc_pt);
            //Double_t new_content = this->GetBinContent(i,j,wc_pt);
            Double_t new_error = (GetBinFit(bin)).evalPointError(&wc_pt);
            this->SetBinContent(bin,new_content);
            this->SetBinError(bin,new_error);
        }
    }
    
}
// evalPointError disabled:
/*
void TH2EFT::Scale(WCPoint wc_pt)
{
    // Warning: calling GetEntries after a call to this function will return a 
    // non-zero value, even if the histogram was never filled.
    
    for (Int_t i = 1; i <= this->GetNbinsX(); i++) {
        Double_t old_content = this->GetBinContent(i);
        Double_t new_content = this->GetBinContent(i,wc_pt);
        Double_t old_error = this->GetBinError(i);
        this->SetBinContent(i,new_content);
        this->SetBinError(i,old_error*new_content/old_content);
    }
    
}
*/
// Uniformly scale all fits by amt
void TH2EFT::ScaleFits(double amt, Int_t axis=0)
{
    for (uint i = 0; i < this->hist_fits.size(); i++) {
        this->hist_fits.at(i).scale(amt);
    }
}

// Display the fit parameters for all bins
void TH2EFT::DumpFits()
{
    for (uint i = 0; i < this->hist_fits.size(); i++) {
        this->hist_fits.at(i).dump();
    }
}

// Set a fit function for a particular bin (no checks are made if the bin is an over/underflow bin)
void TH2EFT::SetBinFit(Int_t bin, WCFit &fit)
{
    Int_t nhists = this->hist_fits.size();
    if (bin <= 0) {
        underflow_fit = fit;
    } else if (bin > nhists) {
        overflow_fit = fit;
    }
    else this->hist_fits.at(bin - 1) = fit;
}

void TH2EFT::AddBinFit(Int_t bin, WCFit &fit)
{
    auto thisfit = this->GetBinFit(bin);
    thisfit.addFit(fit);
    this->SetBinFit(bin, thisfit);

}

//Create a new bin with rebinned range
TH2EFT* TH2EFT::Rebin(Int_t nbinsx, Double_t *x, Int_t nbinsy, Double_t* y)
{
    //padding to avoid overflow
    double *xnew = new double[nbinsx+1];
    double *ynew = new double[nbinsy+1];
    for(int i = 0; i < nbinsx; i++)
        xnew[i] = x[i];
    for(int i = 0; i < nbinsy; i++)
        ynew[i] = y[i];
    xnew[nbinsx] = x[nbinsx-1]+1;
    ynew[nbinsy] = y[nbinsy-1]+1;
    TH2EFT *h = new TH2EFT(this->GetName(), this->GetTitle(), nbinsx, xnew, nbinsy, ynew);
    for(int i = 1; i <= this->GetNbinsX(); i++) {
        for(int j = 1; j <= this->GetNbinsY(); j++) {
            int thisbin = j*(this->GetNbinsX()+2) + i;
            double xbin = this->GetXaxis()->GetBinCenter(i);//this->GetXaxis()->GetBinLowEdge(i);
            double ybin = this->GetYaxis()->GetBinCenter(j);//this->GetYaxis()->GetBinLowEdge(j);
            xbin = xbin > x[nbinsx-1] ? x[nbinsx-1] : xbin; //adjust for overflow, yes that is x not xnew
            ybin = ybin > y[nbinsy-1] ? y[nbinsy-1] : ybin; //adjust for overflow, yes that is y not ynew
            int newbin = h->FindBin(xbin, ybin);//h->FindBin(i, j);
            newbin = xbin < x[0] ? 0 : newbin; //adjust for underflow, yes that is x not xnew
            newbin = ybin < y[0] ? 0 : newbin; //adjust for underflow, yes that is y not ynew
            double thisval = this->GetBinContent(thisbin);
            double thiserr = this->GetBinError(thisbin);
            double newval = h->GetBinContent(newbin);
            double newerr = h->GetBinError(newbin);
            auto thisfit = this->GetBinFit(thisbin);
            h->SetBinContent(newbin, thisval + newval);
            h->SetBinError(newbin, sqrt(thiserr*thiserr + newerr*newerr));
            h->AddBinFit(newbin, thisfit);
        }
    }
    std::cout << "Done rebinning" << std::endl;

    //Handle overflow
    int thisbin = this->FindBin(this->GetNbinsX()+1, this->GetNbinsY()+1);
    int newbin = h->FindBin(h->GetNbinsX()+1, h->GetNbinsY()+1);
    double thisval = this->GetBinContent(thisbin);
    double thiserr = this->GetBinError(thisbin);
    double newval = h->GetBinContent(newbin);
    double newerr = h->GetBinError(newbin);
    h->SetBinContent(newbin, thisval + newval);
    h->SetBinError(newbin, sqrt(thiserr*thiserr + newerr*newerr));
    h->AddBinFit(newbin, this->overflow_fit);

    //Handle underflow
    thisval = this->GetBinContent(0);
    thiserr = this->GetBinError(0);
    newval = h->GetBinContent(0);
    newerr = h->GetBinError(0);
    h->SetBinContent(0, thisval + newval);
    h->SetBinError(0, sqrt(thiserr*thiserr + newerr*newerr));
    h->AddBinFit(0, this->underflow_fit);

    return h;
}

TH2D* TH2EFT::RebinSM(TH2D* hsm, Int_t nbinsx, Double_t *x, Int_t nbinsy, Double_t* y)
{
    double *xnew = new double[nbinsx+1];
    double *ynew = new double[nbinsy+1];
    for(int i = 0; i < nbinsx; i++)
        xnew[i] = x[i];
    for(int i = 0; i < nbinsy; i++)
        ynew[i] = y[i];
    xnew[nbinsx] = x[nbinsx-1]+1;
    ynew[nbinsy] = y[nbinsy-1]+1;
    TH2D *h = new TH2D(hsm->GetName(), hsm->GetTitle(), nbinsx, xnew, nbinsy, ynew);
    //Loop over all bins, and divide by h1->GetBinContent ^ c1 (sqrt by default)
    for(int i = 1; i <= hsm->GetNbinsX(); i++) {
        for(int j = 1; j <= hsm->GetNbinsY(); j++) {
            int thisbin = j*(hsm->GetNbinsX()+2) + i;
            double xbin = hsm->GetXaxis()->GetBinLowEdge(i);
            double ybin = hsm->GetYaxis()->GetBinLowEdge(j);
            xbin = xbin > x[nbinsx-1] ? x[nbinsx-1] : xbin; //adjust for overflow, yes that is x not xnew
            ybin = ybin > y[nbinsy-1] ? y[nbinsy-1] : ybin; //adjust for overflow, yes that is y not ynew
            int newbin = h->FindBin(xbin, ybin);//h->FindBin(i, j);
            newbin = xbin < x[0] ? 0 : newbin; //adjust for underflow, yes that is x not xnew
            newbin = ybin < y[0] ? 0 : newbin; //adjust for underflow, yes that is y not ynew
            double thisval = hsm->GetBinContent(thisbin);
            double thiserr = hsm->GetBinError(thisbin);
            double newval = h->GetBinContent(newbin);
            double newerr = h->GetBinError(newbin);
            h->SetBinContent(newbin, thisval + newval);
            h->SetBinError(newbin, sqrt(thiserr*thiserr + newerr*newerr));
        }
    }

    //Handle overflow
    int thisbin = hsm->FindBin(hsm->GetNbinsX()+1, hsm->GetNbinsY()+1);
    int newbin = h->FindBin(h->GetNbinsX()+1, h->GetNbinsY()+1);
    double thisval = hsm->GetBinContent(thisbin);
    double thiserr = hsm->GetBinError(thisbin);
    double newval = h->GetBinContent(newbin);
    double newerr = h->GetBinError(newbin);
    h->SetBinContent(newbin, thisval + newval);
    h->SetBinError(newbin, sqrt(thiserr*thiserr + newerr*newerr));

    //Handle underflow
    thisval = hsm->GetBinContent(0);
    thiserr = hsm->GetBinError(0);
    newval = h->GetBinContent(0);
    newerr = h->GetBinError(0);
    h->SetBinContent(0, thisval + newval);
    h->SetBinError(0, sqrt(thiserr*thiserr + newerr*newerr));

    return h;
}
