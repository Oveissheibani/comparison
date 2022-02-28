
#include <Riostream.h>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <TInterpreter.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <THnSparse.h>
#include <TGraphAsymmErrors.h>
#include <TPad.h>
#include <TDatabasePDG.h>
#include <TString.h>
#include <TLine.h>
#include <TList.h>
#include <TGaxis.h>

#include <AliNormalizationCounter.h>
const TString mult[4]={"0199","139","4064","65199"};
const Int_t colors[] = {kBlack,kRed,kBlue,kGreen,kViolet+1,kMagenta+1,kGreen};
//const Int_t colors[] = {kGreen,kRed,kBlue,kGreen,kViolet+1,kMagenta+1,kGreen};
//const Int_t colors[] = {kBlack,kGreen,kBlue,kRed,kOrange+7,kMagenta+1,kGreen};
const Int_t markers[] = {kFullCircle,kFullSquare,kFullDiamond,kFullTriangleUp,kFullTriangleDown,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenTriangleUp,kOpenTriangleDown};


TH1F *Rebined(TH1F  *hold);
void TwoHistogram(TH1F  *hist1[],TH1F  *hist2[], TString name, int log,int dimension );

void TwoHistogram_single(TH1F  *hist1,TH1F  *hist2, TString name, int log,int dimension,TString name1,TString name2 );
void Ratio_single(TH1F  *hist1,TH1F  *hist2, TString name,int dimension);
TH1F *LoadHistograms(TString filename,TString histname,Int_t color1,Int_t color2);



void Ratio(TH1F  *hist1[],TH1F  *hist2[], TString name,int log,int dimension,Double_t RangeMin,Double_t RangeMax,TString Xaxis,TString Yaxis,Double_t branch1,Double_t branch2,TCanvas* Conv,TLegend* leg, int u);
void Draw(TH1F  *hist1[],TString name,int log,int dimension,TString LegText[],Double_t RangeMin,Double_t RangeMax,TString Xaxis,TString Yaxis,TCanvas *Conv,TLegend* leg);


void RatioMB(TH1F  *hist1[],TH1F  *hist2,TString name,TString LegText[]);


void compare(){
  
    

    

    
    TH1F *hist1[4];
    TH1F *hist2[4];
    TH1F *hist3[4];
    TH1F *hist4[4];
    TH1F *Lc[4];
    TH1F *D0[4];
    TH1F *ratio[4];
    
    
 
    
    Double_t Npass2D[4]= {6.05596e+08,2.45393e+08,3.00257e+08,6.22054e+07};
    Double_t Npass2L[4]= {6.16157e+08,2.49785e+08,3.05416e+08,6.22089e+07};
    
    
     int k;
     k=1;
  
    TString hitoname[3]={"histoSigmaCorr","hRECpt","hDirectEffpt"};
    TString cent[4]={"1199","139","4064","65199"};
    
    TString cent1[4]={"0100","60100","1060","010"};
   // TString cent1[4]={"010","60100","1060","010"};
    
    for(int i=0;i<4;i++){

        
      
    hist1[i]=LoadHistograms(Form("LcAccEff_final_eff_0199_WithoutD_V0_0100_1323_1324.root"),"hEff_C",colors[i],colors[i]);
    hist2[i]=LoadHistograms(Form("LcAccEff_final_eff_0199_WithoutD_V0_0100_1321_1322.root"),"hEff_C",colors[i],colors[i]);
        
        
        
     
     
    }


    TString AxisTitle[3]={"p_{T}","Cross-Section","Raw Yield /Event "};


        
        TString title[3]={"#Lambda_{c}^{+}/D ratio","#Lambda_{c}^{+} Raw yield per Event ratio ","#Lambda_{c}^{+} Eff ratio "};


        
  
 
    TLine *line = new TLine(1,1,24,1);
     line->SetLineColor(kRed);
     
    TLegend* legend = new TLegend();
    legend->SetBorderSize(1);
    legend->SetFillStyle(3003);

   
   
    
    
    TString word;
    word="abo";
    TCanvas* Conv = new TCanvas(word,word);
    TCanvas* Conv1 = new TCanvas("test","test");
  
    
 
    Ratio(hist2,hist1,"",0,1,0,2,AxisTitle[0],title[k],1,1,Conv,legend,0); //blue
    line->Draw("same");
    
    return 0;
    
    
    
  
    
    
    
}
    

    
    
void RatioMB(TH1F  *hist1[],TH1F  *hist2,TString name,TString LegText[]){
    
    TCanvas *convas=new TCanvas(name,name);
    convas->cd();
    TH1F *hnew[4];
    
    
    TLegend* leg = new TLegend(0.45,0.65,0.85,0.89);
    leg->SetBorderSize(1);
    leg->SetFillStyle(1);
    leg->SetTextSize(0.045);
    
    
    for(int i=1;i<4;i++){
    //  leg->AddEntry(hist1[i],LegText[i],"lpe");
    }
    
    for (int i=0; i<4; i++) {
        
        hnew[i]=(TH1F*)hist1[i]->Clone();
        hnew[i]->Divide(hnew[i],hist2);
        hnew[i]->GetYaxis()->SetTitle("Prompt Eff ratio to MB");
        hnew[i]->SetMarkerStyle(20);
        hnew[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        hnew[i]->Draw("same");
      //  leg->Draw("same");
        
    }
}
    
void Draw(TH1F  *hist1[],TString name,int log,int dimension,TString LegText[],Double_t RangeMin,Double_t RangeMax,TString Xaxis,TString Yaxis,TCanvas* Conv,TLegend* leg){

    
//TCanvas* Conv = new TCanvas(name,name);
    Conv->cd();
    if (log==1){Conv->SetLogy();}
    
  
    
   
    
       for(int i=dimension-1;i<dimension;i++){
           
           hist1[i]->SetLineColor(colors[i]);
          hist1[i]->SetMarkerColor(colors[i]);
           
           
       //    hist1[i]->SetLineColor(kViolet+10);
       //    hist1[i]->SetMarkerColor(kViolet+10);
           
           
           hist1[i]->SetMarkerStyle(markers[i]);
           hist1[i]->SetLineWidth(1);
           hist1[i]->SetStats(0);
           hist1[i]->GetXaxis()->SetTitle(Xaxis);
           hist1[i]->GetYaxis()->SetTitle(Yaxis);
         
           
       }
    
    leg->SetBorderSize(1);
    leg->SetFillStyle(1);
    //leg->SetTextSize(0.045);
    
    
    for(int i=dimension-1;i<dimension;i++){
      leg->AddEntry(hist1[i],LegText[i],"lpe");
    }
   
for (int i=dimension-1;i<dimension;i++){
   // hist1[i]->GetYaxis()->SetRangeUser(RangeMin,RangeMax);
    hist1[i]->GetXaxis()->SetTitle(Xaxis);
    hist1[i]->GetYaxis()->SetTitle(Yaxis);
 //   hist1[i]->Scale(1/branch);
 //   hist1[i]->Scale(1/nev[i]);
    hist1[i]->Draw("same");}

    
   // leg->Draw("same");
}
    
    
    
    
   
    
    
 
void Ratio(TH1F  *hist1[],TH1F  *hist2[], TString name,int log,int dimension,Double_t RangeMin,Double_t RangeMax,TString Xaxis,TString Yaxis,Double_t branch1,Double_t branch2,TCanvas* Conv,TLegend* leg,int u){
    if(u==0){
        
       // TCanvas* Conv = new TCanvas(name,name);
        Conv->cd();
        if (log==1){Conv->SetLogy();}
        
        
       
           for(int i=dimension-1;i<dimension;i++){
               
               hist1[i]->SetLineColor(colors[i]);
               hist1[i]->SetMarkerColor(colors[i]);
               hist1[i]->SetMarkerStyle(markers[i]);
               hist1[i]->SetLineWidth(1);
               hist1[i]->SetStats(0);
               hist1[i]->GetXaxis()->SetTitle(Xaxis);
               hist1[i]->GetYaxis()->SetTitle(Yaxis);
             
               
           }
        
 
        
   
        
     
        
        
          
        for(int i=dimension-1;i<dimension;i++){
            
         
            hist1[i]->Scale(1/branch1);
            hist2[i]->Scale(1/branch2);
            hist1[i]->Divide(hist1[i],hist2[i]);
            hist1[i]->GetYaxis()->SetRangeUser(RangeMin,RangeMax);
            hist1[i]->Draw("same");

           
            
            
            
           
        }
     

      
        
    }
    
    
    
    if(u==1){
        // TCanvas* Conv = new TCanvas(name,name);
         Conv->cd();
         if (log==1){Conv->SetLogy();}
         
         
        
            for(int i=0;i<dimension;i++){
                
           
                hist1[i]->SetLineColor(colors[i]);
                hist1[i]->SetMarkerColor(colors[i]);
                hist1[i]->SetMarkerStyle(markers[i]);
                hist1[i]->SetLineWidth(1);
                hist1[i]->SetStats(0);
                hist1[i]->GetXaxis()->SetTitle(Xaxis);
                hist1[i]->GetYaxis()->SetTitle(Yaxis);
              
                
            }
         
         
      
         leg->SetBorderSize(1);
         leg->SetFillStyle(0);
       //  leg->SetTextSize(0.045);
         
         
         for(int i=0;i<dimension;i++){
         //leg->AddEntry(hist1[i],LegText[i],"lpe");
           
         }
        
         
      
         
         
           
         for(int i=0;i<dimension;i++){
             
        
             hist1[i]->Scale(1/branch1);
             hist2[i]->Scale(1/branch2);
             hist1[i]->Divide(hist1[i],hist2[i]);
             hist1[i]->GetYaxis()->SetRangeUser(RangeMin,RangeMax);
             if (i==1 || i==3){  hist1[i]->Draw("same");}
           
          // leg->Draw("same");
            
             
             
             
            
         }
      
    }
       
}
    
    

TH1F *Rebined(TH1F  *hold){

const Int_t nPtBins=6;
const Double_t PtLims[nPtBins+1]={1,2,4,6,8,12,24};
Double_t sum=0;
Double_t sumerror=0;

TH1F *hnew;
    
    hnew=new TH1F("hnew","",nPtBins,PtLims);
   // h->Reset("ICESM")
  //  hnew=(TH1F*)hold->Clone("hnew");
  // hnew->Reset("ICESM");

  


for(int i=0;i<nPtBins;i++){
    

    sum=0;
    sumerror=0;
for (int j=0;j<hold->GetNbinsX();j++){
    
    
    if(hold->GetBinLowEdge(j+1)>=PtLims[i]){
        
        if ((hold->GetBinLowEdge(j+1)+hold->GetBinWidth(j+1))<=PtLims[i+1]){
            
            Double_t  Interval=hold->GetBinWidth(j+1);
            
            sum=sum+hold->GetBinContent(j+1)*Interval ;
            sumerror=sumerror+hold->GetBinError(j+1)*Interval ;
            
}
        
    }
}
    
    Double_t   BinW=PtLims[i+1]-PtLims[i];
    cout<<BinW<<endl;
   
    hnew->SetBinContent(i+1,sum/BinW);
    hnew->SetBinError(i+1,sumerror/BinW);
    
   
}

    return hnew;
}


    void TwoHistogram_single(TH1F  *hist1,TH1F  *hist2, TString name, int log,int dimension,TString name1,TString name2 ){
        hist1->GetXaxis()->SetTitle("p_{T}(GeV/c)");
     TCanvas* Conv = new TCanvas(name,name);
        if (log==1){Conv->SetLogy();}
       
        
        
        for(int i=0;i<dimension;i++){
            
            hist1->SetLineColor(kBlue+i);
            hist1->SetMarkerColor(kBlue+i);
            hist1->SetMarkerStyle(markers[i]);
            hist1->SetLineWidth(1);
            hist1->SetStats(0);
            
            hist2->SetLineColor(kRed+i);
            hist2->SetMarkerColor(kRed+i);
            hist2->SetMarkerStyle(markers[i]);
            hist2->SetLineWidth(1);
            hist2->SetStats(0);
            
            
            
        }
        
        
        TLegend* leg = new TLegend(0.45,0.65,0.85,0.89);
        leg->SetBorderSize(1);
        leg->SetFillStyle(1);
        leg->SetTextSize(0.045);
        for(int i=0;i<dimension;i++){
            leg->AddEntry(hist1,mult[i]+name1,"lpe");
           
        }
        for(int i=0;i<dimension;i++){
            leg->AddEntry(hist2,mult[i]+name2,"lpe");}
        
        for (int i=0;i<dimension;i++){
        hist1->Draw("same");
        hist2->Draw("same");
        
        
        }
        leg->Draw("same");
    }
    

void Ratio_single(TH1F  *hist1,TH1F  *hist2, TString name,int dimension){
    TCanvas* Conv = new TCanvas(name,name);
    hist1->GetXaxis()->SetTitle("p_{T}(GeV/c)");
       for(int i=0;i<dimension;i++){
           
           hist1->SetLineColor(colors[i]);
           hist1->SetMarkerColor(colors[i]);
           hist1->SetMarkerStyle(markers[i]);
           hist1->SetLineWidth(1);
           hist1->SetStats(0);
           

           
       }
    
      
    for(int i=0;i<dimension;i++){
        
        hist1->Divide(hist1,hist2);
        
        hist1->Draw("same");
        
    }
 
}
    


TH1F *LoadHistograms(TString filename,TString histname,Int_t color1,Int_t color2){
    TFile *file=TFile::Open(filename);
    TH1F *temp;
    temp=(TH1F*)file->Get(histname);
    temp->SetLineColor(color1);
    temp->SetMarkerColor(color2);

    
    return temp;
    
}
    
void TwoHistogram(TH1F  *hist1[],TH1F  *hist2[], TString name, int log,int dimension ){
    hist1[0]->GetXaxis()->SetTitle("p_{T}(GeV/c)");
 TCanvas* Conv = new TCanvas(name,name);
    if (log==1){Conv->SetLogy();}
   
    
    
    for(int i=0;i<dimension;i++){
        
     //   hist1[i]->SetLineColor(kBlue+i);
    //    hist1[i]->SetMarkerColor(kBlue+i);
    //    hist1[i]->SetMarkerStyle(markers[i]);
    //    hist1[i]->SetLineWidth(1);
        hist1[i]->SetStats(0);
        
    //    hist2[i]->SetLineColor(kRed+i);
    //    hist2[i]->SetMarkerColor(kRed+i);
  //      hist2[i]->SetMarkerStyle(markers[i]);
   //     hist2[i]->SetLineWidth(1);
        hist2[i]->SetStats(0);
        
        
        
    }
    
    
    /*
    TLegend* leg = new TLegend(0.45,0.65,0.85,0.89);
    leg->SetBorderSize(1);
    leg->SetFillStyle(1);
    leg->SetTextSize(0.045);
    for(int i=0;i<dimension;i++){
        leg->AddEntry(hist1[i],mult[i]+" pass1","lpe");
       
    }
    for(int i=0;i<dimension;i++){
        leg->AddEntry(hist2[i],mult[i]+" pass2","lpe");}
     
     */
    
    for (int i=0;i<dimension;i++){
    hist1[i]->Draw("same");
    hist2[i]->Draw("same");
    
    
    }
   // leg->Draw("same");
}
