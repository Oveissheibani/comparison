
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
//const Int_t colors[] = {kBlack,kGreen,kBlue,kRed,kOrange+7,kMagenta+1,kGreen};
const Int_t markers[] = {kFullCircle,kFullSquare,kFullDiamond,kFullTriangleUp,kFullTriangleDown,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenTriangleUp,kOpenTriangleDown};


TH1F *Rebined(TH1F  *hold, Int_t color);
void TwoHistogram(TH1F  *hist1[],TH1F  *hist2[], TString name, int log,int dimension );

void TwoHistogram_single(TH1F  *hist1,TH1F  *hist2, TString name, int log,int dimension,TString name1,TString name2 );
void Ratio_single(TH1F  *hist1,TH1F  *hist2, TString name,int dimension);
TH1F *LoadHistograms(TString filename,TString histname,Int_t color1,Int_t color2);



void Ratio(TH1F  *hist1[],TH1F  *hist2[], TString name,int log,int dimension,TString LegText[],Double_t RangeMin,Double_t RangeMax,TString Xaxis,TString Yaxis,Double_t branch1,Double_t branch2,TCanvas* Conv,TLegend* leg, int u);
void Draw(TH1F  *hist1[],TString name,int log,int dimension,TString LegText[],Double_t RangeMin,Double_t RangeMax,TString Xaxis,TString Yaxis,TCanvas *Conv,TLegend* leg);


void RatioMB(TH1F  *hist1[],TH1F  *hist2,TString name,TString LegText[]);


void compare(){
    
    gROOT->SetStyle("Plain");
    gStyle->SetPalette(1);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetCanvasColor(0);
    gStyle->SetFrameFillColor(0);
    gStyle->SetOptTitle(0);
    gSystem->SetIncludePath("-I. -I$ROOTSYS/include");
    gSystem->SetIncludePath("-I. -I$ALICE_ROOT/include  -I$ROOTSYS/include");
    gSystem->Load("libANALYSIS.so");
    gSystem->Load("libANALYSISalice.so");
    gSystem->Load("libCORRFW.so") ;
    gSystem->Load("libPWGHFbase.so");
    gSystem->Load("libPWGHFvertexingHF.so");
    
    
    //define marker colors
   

    
    
    TH1F *effL[4];
    TH1F *effD[4];
    TH1F *hV0LtoD[4];
    TH1F *RawL[4];
    TH1F *RawD[4];
    
    
    
    
    TH1F *hcrossD[4];
    TH1F *hcrossL[4];
    
    
    TH1F *hcrossD1[4];
    TH1F *hcrossL1[4];
    
    TH1F *hPass2crossD[4];
    TH1F *hPass2crossL[4];
    
    TH1F *hPass1crossD[4];
    TH1F *hPass1crossL[4];
    
    TH1F *hMBcrossD[4];
    TH1F *hMBcrossL[4];
    
    
    TH1F *hv0MBcrossD[4];
    TH1F *hv0MBcrossL[4];
    TH1F *ratio[4];
    
 //   Double_t Npass2D[4]= {6.05596e+08,6.22054e+07,3.00257e+08,2.45393e+08};
 //   Double_t Npass2L[4]= {6.16157e+08,6.22089e+07,3.05416e+08,2.49785e+08};
    
    Double_t Npass2D[4]= {6.05596e+08,2.45393e+08,3.00257e+08,6.22054e+07};
    Double_t Npass2L[4]= {6.16157e+08,2.49785e+08,3.05416e+08,6.22089e+07};
    
    
     int k;
     k=0;
  
    TString hitoname[3]={"histoSigmaCorr","hRECpt","hDirectEffpt"};
    TString cent[4]={"1199","139","4064","65199"};
    
    TString cent1[4]={"0100","60100","1060","010"};
    
    for(int i=0;i<4;i++){

        
        
            //Cross section for Pass2
        //    hcrossD[i]=LoadHistograms(Form("HFPtSpectrum_D_%s_pass2.root",cent[i].Data()),hitoname[k],colors[i],colors[i]);
       //    hcrossL[i]=LoadHistograms(Form("HFPtSpectrum_Lc%s_pass2_jan29.root",cent[i].Data()),hitoname[k],colors[i],colors[i]);
        
    hcrossD1[i]=LoadHistograms(Form("HFPtSpectrum_D_%s_pass2_V0.root",cent1[i].Data()),hitoname[k],colors[i],colors[i]);
    hcrossL1[i]=LoadHistograms(Form("HFPtSpectrum_Lc_%s_pass2_V0.root",cent1[i].Data()),hitoname[k],colors[i],colors[i]);
    }
        
      //  HFPtSpectrum_D_010_pass2_V0._Jan30root
        //HFPtSpectrum_Lc139_pass2_jan29.root
      //  HFPtSpectrum_Lc%s_pass2.root
           //MB for pass2
     //   hPass2crossD[i]=LoadHistograms(Form("HFPtSpectrum_D_MB_pass2_SPD.root"),hitoname[k],colors[0],colors[0]);
      //  hPass2crossL[i]=LoadHistograms(Form("HFPtSpectrum_Lc_MB_pass2_jan25.root"),hitoname[k],colors[0],colors[0]);
        
  //HFPtSpectrum_Lc_MB_pass2.root old pass2
        
           //MB for Pass1
      //     hPass1crossD[i]=LoadHistograms(Form("HFPtSpectrum_D0_0199_MB_pPb_Oct8.root"),hitoname[k],colors[1],colors[1]);
       //    hPass1crossL[i]=LoadHistograms(Form(" HFPtSpectrum_Lc_MB_pass1_jan25.root"),hitoname[k],colors[1],colors[1]);

   //     HFPtSpectrum_Lc_0199_MB_pPb_Oct8.root old pass1
     //   HFPtSpectrum_Lc_MB_pass1_jan25.root
        
        //HFPtSpectrum_D0_0199_MB_pPb_Oct8.root
        
        
        
        //minimum bias paper
    //    hMBcrossD[i]=Rebined(LoadHistograms(Form("HFPtSpectrum_D0_finer_17Jan2019.root"),hitoname[k],colors[2],colors[2]),colors[2]);
   //     hMBcrossL[i]=Rebined(LoadHistograms(Form("HFPtSpectrum_LcpKpi_std_pPb_5TeV_Nb_NewFF_EnScale_MergedMC.root"),hitoname[k],colors[2],colors[2]),colors[2]);
        
        
        
   //     hMBcrossD[i]=LoadHistograms(Form("HFPtSpectrum_D0_finer_17Jan2019.root"),hitoname[k],colors[2],colors[2]);
   //     hMBcrossL[i]=LoadHistograms(Form("HFPtSpectrum_LcpKpi_std_pPb_5TeV_Nb_NewFF_EnScale_MergedMC.root"),hitoname[k],colors[1],colors[1]);
        
        
    
        
        
        
     //  ratio[i]=Rebined(LoadHistograms(Form("oot.root"),"hLcDAverage",kViolet+1,kViolet+1),kViolet+1);
        
        
        
        
        
      //  ratio[i]=LoadHistograms(Form("oot.root"),"hLcDAverage",kViolet+10,kViolet+10);
        
 
 

    
   

    TString AxisTitle[3]={"p_{T}","Cross-Section","Raw Yield /Event "};
    TString Branching[3]={"0.","Cross-Section",};
 

   
    

    
    
    
    
        
        TString title[3]={"#Lambda_{c}^{+}/D ratio","#Lambda_{c}^{+} Raw yield per Event ratio ","#Lambda_{c}^{+} Eff ratio "};
   // TString title[3]={"D^{0} CS ratio","D^{0} Raw yield per Event ratio ","D^{0} Eff ratio "};

        
        
        
   /*
        
    TCanvas* Conv = new TCanvas("DV0eff","DV0cs");
    Conv->cd();
 //   Conv->SetLogy();
    TLegend *effleg=new TLegend();
    
   hcrossD[0]->GetYaxis()->SetTitle("D^{0} Eff #times Acc"); //d^{2} N_{D^{0}} / d y d p_{T}  N_{events}(G e V^{-1}c)
 //   hPass2crossD[0]->SetStats(0);
    hcrossD[0]->GetXaxis()->SetTitle("P_{T}");
  //  hPass2crossD[0]->Draw("same");
  //  hPass2crossD[0]->SetLineWidth(2);
    hcrossD[0]->Draw("same");
    hcrossD[1]->Draw("same");
    hcrossD[2]->Draw("same");
    hcrossD[3]->Draw("same");
    
  //  effleg->AddEntry( hPass2crossD[0],"0-200 Tracklet MB");
  //  effleg->AddEntry(hcrossD[0],"1-200 Tracklet");
  //  effleg->AddEntry(hcrossD[1],"1-40 Tracklet");
 //   effleg->AddEntry(hcrossD[2],"40-65 Tracklet");
 //   effleg->AddEntry(hcrossD[3],"65-200 Tracklet");
    
    
    effleg->AddEntry(hcrossD[0],"0-100%");
    effleg->AddEntry(hcrossD[1],"60-100%");
    effleg->AddEntry(hcrossD[2],"10-60%");
    effleg->AddEntry(hcrossD[3],"0-10%");
    
    
    effleg->Draw("same");
    
    
    TCanvas* Conv2 = new TCanvas("Lcv0eff","Lcv0cs");
    Conv2->cd();
 //   Conv2->SetLogy();
  //  hPass2crossL[0]->SetStats(0);
    hcrossL[0]->GetYaxis()->SetTitle("#Lambda_{c}^{+} Eff #times Acc");  //d^{2} N_{#Lambda_{c}^{+}} / d y d p_{T}  N_{events}(G e V^{-1}c)
    hcrossL[0]->GetXaxis()->SetTitle("P_{T}");
//    hPass2crossL[0]->Draw("same");
    hcrossL[0]->Draw("same");
    hcrossL[1]->Draw("same");
    hcrossL[2]->Draw("same");
    hcrossL[3]->Draw("same");
    effleg->Draw("same");
    
    
*/

    
    TString word;
    word="LctoDratioRawyieldperEvent";
    
    TCanvas* Conv = new TCanvas(word,word);
     //   Conv->SetLogy();
        for(int i=0;i<4;i++){
            hcrossL1[i]->Draw("same");
        }
    
 
    /*
   Conv->SetLogy();
    hcrossL1[0]->Draw("same");
    hcrossL1[1]->Draw("same");
    hcrossL1[2]->Draw("same");
    hcrossL1[3]->Draw("same");
    return 0;
    TCanvas* Conv1 = new TCanvas("Lc/Dpass2pass1","Dpass2pass1");
    
 
    Conv->cd();
   // Conv->SetLogy();
    hcrossD[1]->Draw("same");
    hcrossD1[1]->SetLineColor(kBlue);
    hcrossD1[1]->SetMarkerColor(kBlue);
    hcrossD1[1]->Draw("same");
    
    return 0;
    
     */
    
   // Conv->SetLogy();
  //  hMBcrossD[0]->Draw();
  //  hMBcrossL[0]->Draw("same");
  //  return 0;
  //  TCanvas* Conv1= new TCanvas(title1,title1);
   
    
    //0.0628
    //0.0395

  //  Draw(hPass1crossL,"Lc",1,4,cent,1E4,1E11,"p_{T}","#Lambda_{c}^{+} cross section in MB",Conv);
  //  Draw(hPass2crossL,"Lc",1,4,cent,1E4,1E11,"p_{T}","#Lambda_{c}^{+} cross section in MB",Conv);
    
    
    
    TLine *line = new TLine(1,1,24,1);
     line->SetLineColor(kBlack);
    TLegend* legend = new TLegend();
    legend->SetBorderSize(1);
    legend->SetFillStyle(3003);
  //  leg->SetTextSize(0.045);
   // legend->AddEntry(line,"","lpe");

   
    
    
    
   
  //  TString leg3[4]={"1-199 Tracklet","1-40 Tracklet","40-65 Tracklet","65-200 Tracklet"};
    TString leg3[4]={"1-199 Tracklet","pass1 /paper","pass2/pass1","pass2/paper"};
    
  //  TString leg3[4]={"0-100%","60-100%","10-60%","0-10%"};
    
    
 //   Ratio(hPass2crossL,hMBcrossL,"",0,4,leg2,0,3,AxisTitle[0],title,1,1,Conv,legend);
 //   Ratio(hPass1crossL,hMBcrossL,"",0,2,leg1,0,3,AxisTitle[0],title,1,1,Conv,legend);
    
    
    
    
/*
    Ratio(hPass2crossL,hPass2crossD,"",0,1,leg4,0,3,AxisTitle[0],title[k],0.0628,0.0395,Conv,legend);
    Ratio(hPass1crossL,hPass1crossD,"",0,2,leg5,0,3,AxisTitle[0],title[k],0.0628,0.0395,Conv,legend);
    Ratio(hMBcrossL,hMBcrossD,"",0,3,leg6,0,3,AxisTitle[0],title[k],0.0628,0.0395,Conv,legend);
    Draw(ratio,"",0,4,leg7,0,0,AxisTitle[0],title[k],Conv,legend);
    
*/
    
    /*
    const Int_t nPtBins=6;
    const Double_t PtLims[nPtBins+1]={1,2,4,6,8,12,24};

    TH1D *stat=new TH1D("Lc error bars","",nPtBins,PtLims);
    TH1D *stat1=new TH1D("Lc error bars","",nPtBins,PtLims);
    TH1D *stat2=new TH1D("Lc error bars","",nPtBins,PtLims);
    
    
  //  Conv->SetLogy();
   // hPass2crossL[0]->Draw("same");
   // hPass1crossL[0]->Draw("same");
    */
    
    
 //   TCanvas* Conv1 = new TCanvas("LcDv0","LcDv0");
  //  Conv1->cd();
  //  Conv1->GetYaxis->SetTitle("Lc Err bars");
    
    /*
    for(Int_t i=1;i<hPass2crossD[0]->GetNbinsX();i++){
        
    stat->SetBinContent(i+1,hPass2crossD[0]->GetBinError(i+1)/hPass2crossD[0]->GetBinContent(i+1));
        
    stat1->SetBinContent(i+1,hPass1crossD[0]->GetBinError(i+1)/hPass1crossD[0]->GetBinContent(i+1));
    stat2->SetBinContent(i+1,hMBcrossD[0]->GetBinError(i+1)/hMBcrossD[0]->GetBinContent(i+1));
        
    }
    
    
    TLegend *le=new TLegend();
  
    
    stat->SetLineColor(kBlue);
    stat1->SetLineColor(kRed);
    stat2->SetLineColor(kGreen);
    stat->SetMarkerColor(kBlue);
    stat1->SetMarkerColor(kRed);
    stat2->SetMarkerColor(kGreen);
    le->AddEntry(stat,"D pass2 error bars /CS");
    le->AddEntry(stat1,"D pass1 error bars /CS");
    le->AddEntry(stat2,"D paper error bars /CS");
    stat->Draw("same");
    stat1->Draw("same");
   stat2->Draw("same");
    le->Draw("same");
    */
    
  //  TFile* outf=new TFile(Form("CS_Ratio.root"),"recreate");
  
    
    
  //  RatioMB(hcrossL, hPass2crossL[0],"Eff D pass2",leg3);
    
   // return 0;
    
    
    
    
    
    
  //  Ratio(hPass1crossL,hMBcrossL,"",0,2,leg3,0,3,AxisTitle[0],title[k],1,1,Conv,legend);

    /*
     Ratio(hPass2crossL,hPass2crossD,"",0,4,leg3,0,3,AxisTitle[0],title[k],0.0628,0.0395,Conv,legend);
     Ratio(hPass1crossL,hPass1crossD,"",0,2,leg3,0,3,AxisTitle[0],title[k],0.0628,0.0395,Conv,legend);
     Ratio(hMBcrossL,hMBcrossD,"",0,3,leg3,0,3,AxisTitle[0],title[k],0.0628,0.0395,Conv,legend);

    */
    
    Ratio(hcrossL1,hcrossD1,"",0,4,leg3,0,3,AxisTitle[0],title[k],0.0628,0.0395,Conv,legend,1);
    
    
    
    return 0;
    
    
    
   

 
   //    Ratio(hPass2crossL,hPass1crossL,"",0,4,leg3,0,2,AxisTitle[0],title[k],1,1,Conv1,legend);
   //    Ratio(hPass2crossL,hMBcrossL,"",0,2,leg3,0,2,AxisTitle[0],title[k],1,1,Conv1,legend);
    //   Ratio(hPass1crossL,hMBcrossL,"",0,3,leg3,0,2,AxisTitle[0],title[k],1,1,Conv1,legend);
     //  Ratio(hPass2crossL,ratio,"",0,1,leg3,0,2,AxisTitle[0],title[k],1,1,Conv1,legend,0);
    line->Draw("same");
    
    return 0;
    
    
    return 0;
    
    ratio[0]->Draw("same");
    legend->AddEntry(hPass1crossL[2],"pass1");
    legend->AddEntry(hPass2crossL[3],"pass2");
    legend->AddEntry(hMBcrossL[1],"paper");
    legend->AddEntry(ratio[0],"paper average");
    legend->Draw("same");
    
    return 0;
    
   // Ratio(hMBcrossL,hMBcrossD,"",0,4,leg3,0,3,AxisTitle[0],title[k],0.0628,0.0395,Conv,legend);
        
   // Ratio(hPass1crossL,hMBcrossL,"",0,2,leg3,0,3,AxisTitle[0],title[k],1,1,Conv,legend);
     
        
        
    /*
    pass2topass1=(TH1F*)hPass2crossL[2]->Clone("pass2topass1");
    pass2topaper=(TH1F*)hPass2crossL[3]->Clone("pass2topaper");
    pass1topaper=(TH1F*)hPass1crossL[1]->Clone("pass1topaper");
    
    outf->cd();
    pass2topass1->Write();
    pass2topaper->Write();
    pass1topaper->Write();
    
    outf->Close();
    */
    /*
   Draw(hPass2crossL,"",1,1,leg,0,0,AxisTitle[0],title,Conv,legend);
   Draw(hPass1crossL,"",1,2,leg1,0,0,AxisTitle[0],title,Conv,legend);
   Draw(hMBcrossL,"",1,3,leg2,0,0,AxisTitle[0],title,Conv,legend);
    
    
    */
    
    
 //   line->Draw("same");
  //  legend->Draw("same");
    
    
    
    
    
    
  //  Ratio(hPass1crossL,hPass1crossD,"",0,1,cent,0,1,AxisTitle[0],title,0.0628,0.0395,Conv);
   // Ratio(hPass2crossL,hPass1crossL,"",0,1,cent,0,2,AxisTitle[0],title1,1,1,Conv1);
        
        
      //ratio[0]->Draw("same");
        
    /*
    if(PassNumber==0){nevents=4.99656e+08;}
    if(PassNumber==1){nevents=6.00253e+08;}
*/

   //6.00511e+08
   // 4.9987e+08
    
    
 
  

    
 
    
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
    
    
    
    
   
    
    
 
void Ratio(TH1F  *hist1[],TH1F  *hist2[], TString name,int log,int dimension,TString LegText[],Double_t RangeMin,Double_t RangeMax,TString Xaxis,TString Yaxis,Double_t branch1,Double_t branch2,TCanvas* Conv,TLegend* leg,int u){
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
        
        
     
        leg->SetBorderSize(1);
        leg->SetFillStyle(0);
      //  leg->SetTextSize(0.045);
        
        
        for(int i=dimension-1;i<dimension;i++){
        //leg->AddEntry(hist1[i],LegText[i],"lpe");
          
        }
       
        
     
        
        
          
        for(int i=dimension-1;i<dimension;i++){
            
         
            hist1[i]->Scale(1/branch1);
            hist2[i]->Scale(1/branch2);
            hist1[i]->Divide(hist1[i],hist2[i]);
            hist1[i]->GetYaxis()->SetRangeUser(RangeMin,RangeMax);
            hist1[i]->Draw("same");
         // leg->Draw("same");
           
            
            
            
           
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
             hist1[i]->Draw("same");
          // leg->Draw("same");
            
             
             
             
            
         }
      
    }
       
}
    
    

TH1F *Rebined(TH1F  *hold,Int_t color){

const Int_t nPtBins=6;
const Double_t PtLims[nPtBins+1]={1,2,4,6,8,12,24};
Double_t sum=0;
Double_t sumerror=0;

TH1F *hnew;
    
    hnew=new TH1F("hnew","",nPtBins,PtLims);
   // h->Reset("ICESM")
  //  hnew=(TH1F*)hold->Clone("hnew");
  // hnew->Reset("ICESM");

    hnew->SetLineColor(color);
    hnew->SetMarkerColor(color);


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
