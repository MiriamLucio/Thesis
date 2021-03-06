{
//=========Macro generated from canvas: c1/c1
//=========  (Thu Jul 14 23:25:26 2016) by ROOT version5.34/18
   TCanvas *c1 = new TCanvas("c1", "c1",65,52,700,500);
   c1->SetHighLightColor(2);
   c1->Range(155,-5.38125,505,48.43125);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
   
   TH1F *htemp__1 = new TH1F("htemp__1","Bmass {mo1==mo2 && rmo1key==rmo2key}",100,190,470);
   htemp__1->SetBinContent(8,1);
   htemp__1->SetBinContent(9,9);
   htemp__1->SetBinContent(10,8);
   htemp__1->SetBinContent(11,13);
   htemp__1->SetBinContent(12,10);
   htemp__1->SetBinContent(13,20);
   htemp__1->SetBinContent(14,9);
   htemp__1->SetBinContent(15,17);
   htemp__1->SetBinContent(16,12);
   htemp__1->SetBinContent(17,19);
   htemp__1->SetBinContent(18,15);
   htemp__1->SetBinContent(19,15);
   htemp__1->SetBinContent(20,25);
   htemp__1->SetBinContent(21,19);
   htemp__1->SetBinContent(22,20);
   htemp__1->SetBinContent(23,25);
   htemp__1->SetBinContent(24,19);
   htemp__1->SetBinContent(25,16);
   htemp__1->SetBinContent(26,22);
   htemp__1->SetBinContent(27,25);
   htemp__1->SetBinContent(28,26);
   htemp__1->SetBinContent(29,24);
   htemp__1->SetBinContent(30,21);
   htemp__1->SetBinContent(31,17);
   htemp__1->SetBinContent(32,21);
   htemp__1->SetBinContent(33,29);
   htemp__1->SetBinContent(34,18);
   htemp__1->SetBinContent(35,22);
   htemp__1->SetBinContent(36,10);
   htemp__1->SetBinContent(37,24);
   htemp__1->SetBinContent(38,9);
   htemp__1->SetBinContent(39,7);
   htemp__1->SetBinContent(40,12);
   htemp__1->SetBinContent(41,5);
   htemp__1->SetBinContent(42,14);
   htemp__1->SetBinContent(43,2);
   htemp__1->SetBinContent(44,1);
   htemp__1->SetBinContent(46,1);
   htemp__1->SetBinContent(49,1);
   htemp__1->SetBinContent(50,2);
   htemp__1->SetBinContent(58,1);
   htemp__1->SetBinContent(59,1);
   htemp__1->SetBinContent(61,1);
   htemp__1->SetBinContent(62,1);
   htemp__1->SetBinContent(68,1);
   htemp__1->SetBinContent(70,1);
   htemp__1->SetBinContent(71,1);
   htemp__1->SetBinContent(72,1);
   htemp__1->SetBinContent(73,1);
   htemp__1->SetBinContent(74,2);
   htemp__1->SetBinContent(75,2);
   htemp__1->SetBinContent(76,2);
   htemp__1->SetBinContent(77,2);
   htemp__1->SetBinContent(78,7);
   htemp__1->SetBinContent(79,2);
   htemp__1->SetBinContent(80,3);
   htemp__1->SetBinContent(81,2);
   htemp__1->SetBinContent(82,9);
   htemp__1->SetBinContent(83,7);
   htemp__1->SetBinContent(84,9);
   htemp__1->SetBinContent(85,10);
   htemp__1->SetBinContent(86,14);
   htemp__1->SetBinContent(87,10);
   htemp__1->SetBinContent(88,9);
   htemp__1->SetBinContent(89,15);
   htemp__1->SetBinContent(90,20);
   htemp__1->SetBinContent(91,37);
   htemp__1->SetBinContent(92,40);
   htemp__1->SetBinContent(93,41);
   htemp__1->SetBinContent(94,1);
   htemp__1->SetEntries(838);
   htemp__1->SetDirectory(0);
   
   TPaveStats *ptstats = new TPaveStats(0.78,0.775,0.98,0.935,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetTextAlign(12);
   ptstats->SetTextFont(42);
   TText *text = ptstats->AddText("htemp");
   text->SetTextSize(0.0368);
   text = ptstats->AddText("Entries = 838    ");
   text = ptstats->AddText("Mean  =  313.6");
   text = ptstats->AddText("RMS   =  82.46");
   ptstats->SetOptStat(1111);
   ptstats->SetOptFit(0);
   ptstats->Draw();
   htemp__1->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(htemp__1);
   htemp__1->SetLineWidth(2);
   htemp__1->SetMarkerStyle(20);
   htemp__1->GetXaxis()->SetTitle("Bmass");
   htemp__1->GetXaxis()->SetRange(1,100);
   htemp__1->GetXaxis()->SetLabelFont(42);
   htemp__1->GetXaxis()->SetLabelSize(0.035);
   htemp__1->GetXaxis()->SetTitleSize(0.035);
   htemp__1->GetXaxis()->SetTitleFont(42);
   htemp__1->GetYaxis()->SetLabelFont(42);
   htemp__1->GetYaxis()->SetLabelSize(0.035);
   htemp__1->GetYaxis()->SetTitleSize(0.035);
   htemp__1->GetYaxis()->SetTitleFont(42);
   htemp__1->GetZaxis()->SetLabelFont(42);
   htemp__1->GetZaxis()->SetLabelSize(0.035);
   htemp__1->GetZaxis()->SetTitleSize(0.035);
   htemp__1->GetZaxis()->SetTitleFont(42);
   htemp__1->Draw("");
   
   TPaveText *pt = new TPaveText(0.15,0.9339831,0.85,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   text = pt->AddText("Bmass {mo1==mo2 && rmo1key==rmo2key}");
   pt->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
