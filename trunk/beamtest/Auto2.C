#include "time.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TGraphErrors.h"	

	
int itev;
int itev2;
int ITre[100000];
int ETre[100000];
int STre[100000];
int J1re[100000];
int J2re[100000];
float CAL[100];
float ECAL[100];
int tbuff[2050*96*500];
int buff1[2050*96];
int data[4096*76];
int data2[2];
float Prob[100];
float EProb[100];
float Z[100];
float EZ[100];
float Q[100];
float EQ[100];
float Q2[100];
float EQ2[100];
float QC[100];
float EQC[100];
float S[100];
float ES[100];
float T[100];
float ET[100];
float T2[100];
float ET2[100];
float CHH[10000];
float THH[10000];
float SHH[10000];

TNtuple *n1;


void convert(){

	 int sp = 0;
	 int ev = 0;

  for (int ich = 0 ; ich < 96 ; ich++){
  	  int off = 2050*ich;
	  int ch =((buff1[off+0]>>24)&0x1F);
	  int ich2 = (ich%3)*32 + ((int)ich/3);
	  int ich3 = ich2;
      if (ich2>=64&&ich2<80 || ich2 >= 92) ich3 = -1;
      if (ich2>=80&&ich2<92) ich3 = ich2 - 16; 
     // printf("%2i %2i %2i\n",ich,ich2,ich3);
      if (ich==2){
		for(int j=0;j<8;j++) {
			sp = sp + (((buff1[off+(3520+64*j)/2+2]&0xfff)>500)<<(7-j));
		   // printf("%i %i\n",j,sp);
		}       	  
      }
      if (ich==5){
		for(int j=0;j<12;j++) {
			ev = ev + (((buff1[off+(3520+43*j)/2+2]&0xfff)>500)<<(11-j));
		  //  printf("%i %i %i \n",j,(buff1[off+(3520+43*j)/2+2]&0xfff),ev);
		}
      }

      if (ich3>=0 && ich3 <=76){
 		 for(int i=0; i < 2048; i++){
	   		data[ich3*4096+i*2] = (unsigned short)buff1[off+i+2] & 0xFFF;
	    	data[ich3*4096+i*2+1] = (unsigned short)(buff1[off+i+2]>>16) & 0xFFF;
	    	}
	  }
	}
	data2[0]=sp;
	data2[1]=ev;
	return;	
}


Double_t fitfe(Double_t *v, Double_t *par)
{
  Double_t fitval = 0;
  double PI=3.14159265358979;
  double t = v[0] - 0.0;
  double A = par[0];
  double f = par[1];
  double d = par[2];
  double off = par[3];
  return A*sin(2*PI*(t-d)/f)+off;
}

Double_t fitf(Double_t *v, Double_t *par)
{
  double t = v[0];
  double PI = 3.14159265358979;
  double q     = par[0];
  double mean  = par[1];
  double sigma = par[2];
  double offset = par[3];
  //printf("%f %f %f %f %f\n",t,q,mean,sigma,offset);
  double g = q / sqrt(2*PI) / sigma * exp(-(t-mean)*(t-mean)/2/sigma/sigma);
  return offset + g;
}


Double_t fitf1(Double_t *v, Double_t *par)
{
  Double_t fitval = 0;
  double t = v[0] - 0.0;
  double tau = par[0];
  double A = par[1];
  double off = par[2];
  return off+A*exp(-t/tau);
}


void ana(char fname[], int nev, int nrbin, int cycle, int mode){

int tttt=time(0);
printf("Start %i \n",time(0)-tttt);

  if (cycle==0) cycle = 14045;
  float G=0.1;
  FILE *fp1;
  char fname2[100];
  sprintf(fname2,"M:/backup/Beamtest/daqtest/data/%s.dat",fname);
  fp1 = fopen(fname2,"rb");

  //int *tbuff = calloc(sizeof(int),2050*96*nev);
  fread(tbuff,sizeof(int),2050*96*nev,fp1);
  //free(tbuff);

  //int *buff1 = calloc(sizeof(int),2050*96);
  fread(buff1,sizeof(int),2050*96,fp1);

  fclose(fp1);

  //int *data = calloc(sizeof(int),4096*76);
  int evnum,spnum;
  //int data2[2];
  
  //convert(buff1,data,data2);
  convert();
 //free(buff1);

  spnum=data2[0];
  evnum=data2[1];
   	   int imatch = -1;
   	   for (int i=0;i<itev;i++){
   	   	   if (spnum==STre[i] && evnum == ETre[i]){
   	   	   	   imatch = ITre[i];
   	   	   	   break;
   	   	   }
	   }
	   printf("%i %i %i %i\n",nev,spnum,evnum,imatch);
  if (imatch < 0) return;

//printf("Read data %i\n",time(0)-tttt);

  char hname[100];
  char hname2[100];
  char funame[100]; 
  char fnamep[100];
  char htitle[100];

  TH2F *hx0 = new TH2F("hx0","hx0",76,0,76,4096,0,4096);
    
  int nbin = 4096/nrbin;
  int iend  = nbin * nrbin;
  TH2F *hxt = new TH2F("hxt","hxt",76,0,76,nbin,-512*0.4,(iend-512)*0.4);
  TH2F *hxt2 = new TH2F("hxt2","hxt2",76,0,76,nbin,-512*0.4,(iend-512)*0.4);
  TH2F *hxt3 = new TH2F("hxt3","hxt3",76,0,76,nbin,-512*0.4,(iend-512)*0.4);
  TH2F *hxt4 = new TH2F("hxt4","hxt4",76,0,76,nbin,-512*0.4,(iend-512)*0.4);
  TH2F *hxt5 = new TH2F("hxt5","hxt5",76,0,76,50,0,500);

  TH2F *hxn = new TH2F("hxn","hxn",76,0,76,4096,0,4096);

  TH1F *hs = new TH1F("hs","hs",76,0,76);

  TH1F *hsig1 = new TH1F("hsig1","hsig1",76,0,76);
  TH1F *hsig2 = new TH1F("hsig2","hsig2",76,0,76);
  TH1F *hsig3 = new TH1F("hsig3","hsig3",76,0,76);
  TH1F *hsig4 = new TH1F("hsig4","hsig4",76,0,76);

  TH1F *hnhit = new TH1F("hnhit","Hit",100,0,100.);
  TH1F *hthit = new TH1F("hthit","Hit",100,0,500.);
  TH1F *hqhit = new TH1F("hqhit","Hit",150,-10,20.);

  TH1F *hxxxx = new TH1F("hxxxx","hxxxx",20,0,20);
  TH1F *hxxxxs = new TH1F("hxxxxs","hxxxxs",20,0,20);


 TH1F *hqc = new TH1F("hqc","hqc",100,0,2000.);
 TH1F *htc = new TH1F("htc","htc",100,0,2000.);
 TH1F *hsc = new TH1F("hsc","hsc",100,0,20.);
 TH1F *hpc = new TH1F("hpc","hpc",100,0,1.);

	TH1F *htta[76];
	TH1F *httb[76];
	TH1F *httc[76];
	TH1F *httd[76];
	TH1F *htte[76];
	TH1F *hxxc[76];
	TH1F *hxxd[76];
	TF1  *ft[76];
	for (int i=0;i<76;i++) {
		
		sprintf(hname,"CH %i",i);
		sprintf(hname2,"htta%i",i);
		htta[i]= new TH1F(hname2,hname,nbin,0,iend);

		sprintf(hname,"CH %i",i);
		sprintf(hname2,"httb%i",i);
		httb[i]= new TH1F(hname2,hname,nbin,0,iend);

		sprintf(hname,"CH %i",i);
		sprintf(hname2,"httc%i",i);
		httc[i]= new TH1F(hname2,hname,80,0,4000);

		sprintf(hname,"CH %i",i);
		sprintf(hname2,"httd%i",i);
		httd[i]= new TH1F(hname2,hname,nbin,0,iend);

		sprintf(hname,"CH %i",i);
		sprintf(hname2,"htte%i",i);
		htte[i]= new TH1F(hname2,hname,nbin,0,iend);

		sprintf(hname,"XXc %i",i);
		sprintf(hname2,"hxxc%i",i);
		hxxc[i]= new TH1F(hname2,hname,100,-50,50);

		sprintf(hname,"XXd %i",i);
		sprintf(hname2,"hxxd%i",i);
		hxxd[i]= new TH1F(hname2,hname,20,-20,20);


	}

//printf("Book Histograms %i\n",time(0)-tttt);

 for (int i=0;i<4096*76;i++){
	int ich = i/4096;
	int ich2 = (ich/8);
	int j   = i%4096;
	int j2 = ((j*100)%cycle)/100.;
	int j3 = ((j*100)/cycle);

	hx0->Fill(ich+0.01,j+0.01,data[i]);

 	if (j3<3) hxn->Fill(ich+0.01,j2+0.01,data[i]/3.0);
 	//if (j3>=15 && j3<18) hxn->Fill(ich+0.01,j2+0.01,data[i]/6.0);

 	if (j>0 && j < 500) {
 	 	 hs->Fill(ich,data[i]/1000.);
 	 }
 	if (j>2500 && j < 3000) {
 	 	 hs->Fill(ich,data[i]/1000.);
 	 }
 }	


 for (int ich=0;ich<76;ich++){
 	 float xx = 0.5*hxn->GetBinContent(ich+1,1) + 0.5*hxn->GetBinContent(ich+1,140);
     hxn->SetBinContent(ich+1,141,xx);
     for (int i=141;i<4096;i++){
		 int j = ((i*100)%cycle)/100;
     	 float yy = hxn->GetBinContent(ich+1,j+1);
     	 hxn->SetBinContent(ich+1,i+1,yy);
     }
 }
/*
//printf("Fill Histograms %i\n",time(0)-tttt);
if (mode == -999){
 hx0->ProjectionY("htmp1",6,6);
 hxn->ProjectionY("htmp2",6,6);
 htmp1->SetLineColor(2);
 TCanvas *c7 = new TCanvas("c7","c7",800,600);
 htmp1->Draw();
 htmp2->Draw("same");
 return;
}
*/
 for (int i=0;i<4096*76;i++){
	int ich = i/4096;
	int ich2 = (ich/8);
	int j   = i%4096;
	int j2  = ((100*j)%cycle)/100;
 	 float q   = data[i];
 	 float sum1 = hs->GetBinContent(ich+1);
 	 float sum2 = hxn->GetBinContent(ich+1,j2+1);
 	 float q1  = q - sum1;
	 float kfac = 1.0;
 	 float q3  = q - sum1 - (sum2 - sum1)*kfac;
 	 hxt->Fill(ich,(j-512.)*0.4+0.001,q1/(1.0*nrbin)/CAL[ich]);
 	 hxt2->Fill(ich,(j-512.)*0.4+0.001,q3/(1.0*nrbin)/CAL[ich]);
	 htta[ich]->Fill(j,q1/(1.0*nrbin));
	 httb[ich]->Fill(j,q3/(1.0*nrbin));
	
 }

//printf("Noise Subtraction 1 %i\n",time(0)-tttt);

int nhit = 0;
TF1 *ftest[76];
 for (int i=0;i<76;i++){

	// Recaclulate Pedestal 
	 for (int j=0;j<500/nrbin;j++){
	 	float xxx = httb[i]->GetBinContent(j+1);
 	 	 hxxc[i]->Fill(xxx); 	 	 
	 }
	 for (int j=2500/nrbin;j<3000/nrbin;j++){
	 	float xxx = httb[i]->GetBinContent(j+1);
 	 	 hxxc[i]->Fill(xxx);
	 }
	 
	 float sigma = hxxc[i]->GetRMS();
	 hsig2->SetBinContent(i+1,sigma);

	 for (int j=0;j<80;j++){
	 	 	float nsumx=0;
	 	    float sumx=0;
 	 	 	for(int k=0;k<50/nrbin;k++){
		 	float xxx = httb[i]->GetBinContent(j*50/nrbin+k+1);
 	 	 		if (abs(xxx) < 3.0*sigma){
 	 	 			sumx = sumx + xxx;
 	 	 			nsumx = nsumx + 1.0;
 	 	 		}
 	 	 	}
 	 	 	if (nsumx>0){
		 		httc[i]->SetBinContent(j+1,sumx/nsumx);	
		 	 	httc[i]->SetBinError(j+1,sigma/sqrt(5.0));
		 	}
		 	else{
		 		httc[i]->SetBinContent(j+1,0.);	
		 	 	httc[i]->SetBinError(j+1,10.0);
		 	 	//httc[i]->SetBinError(j+1,sigma);
		 	}
 	 } 	 
	sprintf(funame,"ftest%i",i);
	ftest[i]=new TF1(funame,fitfe,300,2000,4);
	ftest[i]->SetParameters(10,839,400,0);
	ftest[i]->SetParLimits(0,0,20);
	ftest[i]->SetParLimits(1,800,1200);
	ftest[i]->SetParLimits(2,-600,600);
	ftest[i]->SetParLimits(3,-10,10);
	httc[i]->Fit(funame,"R0q");
	for (int j=0;j<nbin;j++){
		float xxx = httb[i]->GetBinContent(j+1);
		float xxt = httb[i]->GetBinCenter(j+1);
		float xxb = httc[i]->GetFunction(funame)->Eval(xxt);
		httd[i]->SetBinContent(j+1,xxx-xxb);
		hxt3->SetBinContent(i+1,j+1,(xxx-xxb)/CAL[i]);
	}
	 for (int j=0;j<500/nrbin;j++){
	 	float xxx = httd[i]->GetBinContent(j+1);
 	 	 hxxd[i]->Fill(xxx); 	 	 
	 }
	 for (int j=2000/nrbin;j<3000/nrbin;j++){
	 	float xxx = httd[i]->GetBinContent(j+1);
 	 	 hxxd[i]->Fill(xxx); 	 	 
	 }
	 ftest[i]->Delete();
	 float sigma2 = hxxd[i]->GetRMS();
	 hsig3->SetBinContent(i+1,sigma2);
	 for (int j=0;j<nbin;j++){
		httd[i]->SetBinError(j+1,sigma2);
		float xxx = httd[i]->GetBinContent(j+1);
		htte[i]->SetBinContent(j+1,xxx/sigma2);
		hxt4->SetBinContent(i+1,j+1,xxx/sigma2);
		float tt = (j*nrbin-512.)*0.4;
		float qc = xxx*G/CAL[i]/exp(-tt/450.);
		if (tt>0 && tt<500) {
			hqhit->Fill(xxx/sigma2);
			if ( xxx/sigma2>4.0){
				CHH[nhit]=i;
				THH[nhit]=tt;
				SHH[nhit]=xxx/sigma2;
				nhit++;
				hthit->Fill(tt);
			}
		}
	}	 
 }

//printf("Noise Subtraction 2 %i\n",time(0)-tttt);

// printf("Nhit=%i\n",nhit);
 hthit->GetXaxis()->SetRange(31,70);
 float sthit = hthit->GetRMS();

int ngood=0; 
int nbad=0;
  htte[2]->GetXaxis()->SetRange(900/nrbin,1400/nrbin);
  float smax2   = htte[2]->GetMaximum();
  int  maxbin2 = htte[2]->GetMaximumBin();
  float tmax2  = htte[2]->GetBinCenter(maxbin2);
  httd[2]->GetXaxis()->SetRange(maxbin2-100/nrbin,maxbin2+100/nrbin);
  float tmax1=0;
  int  maxbin1 =0;
  for (int i=2;i<74;i++){
  	  //if (i==24) continue;
  	  if (i==25) continue;
  	  if (i==26) continue;
  	  if (i==49) continue;
  	  if (i==50) continue;
  	  //if (i==51) continue;
  	  /*
  	  if (i==16) continue;
  	  if (i==32) continue;
  	  if (i==64) continue;
  	  */
  	 if (i>2 && hsig3->GetBinContent(i+1)>6) continue;
  	 if (i>2){
     	htte[i]->GetXaxis()->SetRange(maxbin1-100/nrbin,maxbin1+100/nrbin);
     	httd[i]->GetXaxis()->SetRange(maxbin1-100/nrbin,maxbin1+100/nrbin);
     }
     float smax   = htte[i]->GetMaximum();
     int  maxbin = htte[i]->GetMaximumBin();
	 float tmax  = htte[i]->GetBinCenter(maxbin);
     if (smax < 2.5 && nbad==0) break;
     if (smax < 2.5) {nbad++; continue;}
     if (smax > 2.5) nbad=0;;
     
	 float frmin = tmax - 30.;
	 float frmax = tmax + 30.;
	 sprintf(hname,"ft%i",i);
	 ft[i]= new TF1(hname,fitf,frmin,frmax,4);
		ft[i]->SetParameters(1000,tmax,9.0,0.0);
		//ft[i]->FixParameter(2,9.0);
		ft[i]->FixParameter(3,0);

	 	httd[i]->Fit(ft[i],"Rq0");


	 	float qc = ft[i]->GetParameter(0);
	 	float tc = ft[i]->GetParameter(1);
	 	float sc = ft[i]->GetParameter(2);
	 	float eqc = ft[i]->GetParError(0);
	 	float etc = ft[i]->GetParError(1);
	 	float esc = ft[i]->GetParError(2);
	    float prob = ft[i]->GetProb();
	    float tt = (tc-512.)*0.4;
	    float ett = etc*0.4;
	    float qcc = qc*G/CAL[i]/exp(-tt/450.);
	    float eqcc = eqc*G/CAL[i]/exp(-tt/450.);
	 	//printf("%i %f %f ",i,smax,tmax);
        //printf("%f %f %f %f\n",qc,tc,sc,prob);
	    hqc->Fill(qc);
	    //if (qc < 500) continue;
	    hpc->Fill(prob);
	    htc->Fill(tc);
	    hsc->Fill(sc);
	    hxt5->Fill(i,tt,qcc);
	    tmax1=tc;
	    maxbin1 = tc/nrbin;
	    T[ngood]=tt;
	    ET[ngood]=ett;
	    Q[ngood]=qc*G;
	    EQ[ngood]=eqc*G;
	    QC[ngood]=qcc;
	    EQC[ngood]=eqcc;
	    Z[ngood]=i;
	    EZ[ngood]=1./sqrt(12.);
	    S[ngood]=sc*0.4;
	    ES[ngood]=esc*0.4;
	    Prob[ngood]=prob;
	    ngood++;
	    //ft[i]->Delete();
 }


//printf("Fit %i\n",time(0)-tttt);


FILE *fp2;
  char fnameo[100];
  sprintf(fnameo,"result/%s_%i.txt",fname,nev);
  fp2 = fopen(fnameo,"w");
  fprintf(fp2,"%i %i %i %i %i %i %f\n",nev,spnum,evnum,cycle,ngood,nhit,sthit);
//  fprintf(fp2,"%i %i %i %i %i %i %f\n",nev,spnum,evnum,imatch,ngood,nhit,sthit);
  for (int i=0;i<ngood;i++){
	fprintf(fp2,"%i %f %f %f %f %f %f %f %f %f %f %f %f\n",i,Z[i],EZ[i],T[i],ET[i],Q[i],EQ[i],S[i],ES[i],Q2[i],EQ2[i],T2[i],Prob[i]);
 }
  for (int i=0;i<nhit;i++){
	fprintf(fp2,"%i %f %f %f\n",i,CHH[i],THH[i],SHH[i]);
 }
fclose(fp2);

TGraphErrors *dedx = new TGraphErrors(ngood, Z, QC, EZ, EQC);
TGraphErrors *tau  = new TGraphErrors(ngood, Z, T, EZ, ET);
TGraphErrors *att  = new TGraphErrors(ngood, T, QC, ET, EQC);
TGraphErrors *sss  = new TGraphErrors(ngood, Z, S, EZ, ES);
TGraphErrors *prb  = new TGraphErrors(ngood, Z, Prob, EZ, EProb);

 sprintf(htitle,"Nhit=%i",nhit);
 hqhit->SetTitle(htitle);
 hqhit->SetXTitle("Hit significance");
 sprintf(htitle,"RMS=%5.2f",sthit);
 hthit->SetTitle(htitle);
 hthit->SetXTitle("Hit time (#mus)");
 
char htitle2[100];
sprintf(htitle2,"dE/dx Result, File: %s / i: %i / Spill: %i / Event: %i",fname,nev,spnum,evnum);

 TH1F *hh1 = new TH1F("hh1",htitle2,1000,-5,80);
 TH1F *hh2 = new TH1F("hh2",htitle2,1000,-5,80);
 TH1F *hh3 = new TH1F("hh3",htitle2,120,-100,500);
 TH1F *hh4 = new TH1F("hh4",htitle2,1000,-5,80);


if (mode>0){

 TCanvas *cev = new TCanvas("cev","cev",600,800);
 cev->Divide(1,3);
 cev->cd(1);
 gPad->SetLogy();
 hqhit->Draw();
 cev->cd(2);
 gPad->SetLogy(0);
 hthit->Draw();
 cev->cd(3);
 hsig3->Draw();
 
 sprintf(fnamep,"plot/%s_%i_hitc.png",fname,nev);
 printf("%s\n",fnamep);
 cev->Print(fnamep);
 
 TCanvas *cc = new TCanvas("cc","cc",800,800);
 cc->Divide(2,2);
 cc->cd(1);
 hqc->Draw();
 cc->cd(2);
 htc->Draw();
 cc->cd(3);
 hsc->Draw();
 cc->cd(4);
 hpc->Draw();
 sprintf(fnamep,"plot/%s_%i_fitc.png",fname,nev);
 printf("%s\n",fnamep);
 cc->Update();
 cc->Print(fnamep);

TCanvas *c4 = new TCanvas("c4","c4",800,800);

c4->Divide(2,2);
c4->cd(1);
hh3->SetMinimum(0);
hh3->SetMaximum(250.);
hh3->SetXTitle("Drift Time (#mus)");
hh3->SetYTitle("Signal Charge (fC)");
hh3->Draw();
att->Draw("Psame");

c4->cd(2);
hh2->SetMinimum(0);
hh2->SetMaximum(500.);
hh2->SetXTitle("TPC Channel (cm)");
hh2->SetYTitle("Drift Time (#mus)");
hh2->Draw();
tau->Draw("Psame");

c4->cd(3);
hh1->SetMinimum(0);
hh1->SetMaximum(250.);
hh1->SetXTitle("TPC Channel (cm)");
hh1->SetYTitle("Signal Charge (fC)");
hh1->Draw();
dedx->Draw("Psame");
/*
c4->cd(4);
hh4->SetMinimum(-0.1);
hh4->SetMaximum(1.1);
hh4->SetXTitle("TPC Channel (cm)");
hh4->SetYTitle("Fit Probability");
hh4->Draw();
prb->Draw("Psame");
*/
c4->cd(4);
hh4->SetMinimum(0);
hh4->SetMaximum(20.);
hh4->SetXTitle("TPC Channel (cm)");
hh4->SetYTitle("Signal Sigma (#mus)");
hh4->Draw();
sss->Draw("Psame");

 sprintf(fnamep,"plot/%s_%i_fitp.png",fname,nev);
 printf("%s\n",fnamep);
 c4->Print(fnamep);


 TCanvas *c4a = new TCanvas("c4a","c4a",1140,600);

     c4a->cd(1);
     hxt3->SetMinimum(-10);
     hxt3->SetMaximum(60);
     hxt3->GetYaxis()->SetRange(512/nrbin,(1260+512)/nrbin);
     sprintf(htitle,"File: %s / i: %i / Spill: %i / Event: %i",fname,nev,spnum,evnum);
     hxt3->SetTitle(htitle);
     hxt3->SetXTitle("TPC Channel");
     hxt3->SetYTitle("Time (#mus)");
     hxt3->Draw("colz");
     sprintf(fnamep,"plot/%s_%i_data.png",fname,nev);
     printf("%s\n",fnamep);
     c4a->Print(fnamep);


 TCanvas *c0 = new TCanvas("ce0","ce0",1140,1200);
 c0->Divide(1,2);
 c0->cd(1);
 hxt3->Draw("colz");
 c0->cd(2);
    httd[2]->SetXTitle("FADC Time (400ns)");
    httd[2]->SetYTitle("FADC Counts");
  for (int i=0;i<ngood;i++){
 	 
 		int iii = Z[i];

		httd[iii]->SetLineColor((i%10)*4+51);
		httd[iii]->SetMinimum(-20);
		httd[iii]->SetMaximum(100);

		ft[iii]->SetLineColor((i%10)*4+51);
		ft[iii]->SetLineWidth(2);
		
		if (i==0) {httd[iii]->Draw("e1");}
	 	else {httd[iii]->Draw("e1same");}
	 	ft[iii]->Draw("same");
 }

 TCanvas *c4b = new TCanvas("c4b","c4b",1140,600);

     c4b->cd(1);
     hxt4->SetMinimum(4.0);
     hxt4->SetMaximum(20);
     hxt4->GetYaxis()->SetRange(512/nrbin,(1260+512)/nrbin);
     sprintf(htitle,"File: %s / i: %i / Spill: %i / Event: %i",fname,nev,spnum,evnum);
     hxt4->SetTitle(htitle);
     hxt4->SetXTitle("TPC Channel");
     hxt4->SetYTitle("Time (#mus)");
     hxt4->Draw("colz");
     sprintf(fnamep,"plot/%s_%i_hit.png",fname,nev);
     printf("%s\n",fnamep);
     c4b->Print(fnamep);

 TCanvas *c4c = new TCanvas("c4c","c4c",1140,600);

     c4c->cd(1);
     hxt5->SetMinimum(0.);
     hxt5->SetMaximum(200.);
     sprintf(htitle,"File: %s / i: %i / Spill: %i / Event: %i",fname,nev,spnum,evnum);
     hxt5->SetTitle(htitle);
     hxt5->SetXTitle("TPC Channel");
     hxt5->SetYTitle("Time (#mus)");
     hxt5->Draw("colz");
     sprintf(fnamep,"plot/%s_%i_fit.png",fname,nev);
     printf("%s\n",fnamep);
     c4c->Print(fnamep);



if (mode < 100){
cc->Close();
c4->Close();
c4a->Close();
c4b->Close();
c4c->Close();
cev->Close();
}
}
if (mode < 100){

for (int i=0;i<76;i++){
	htta[i]->Delete();
	httb[i]->Delete();
	httc[i]->Delete();
	httd[i]->Delete();
	htte[i]->Delete();
	hxxc[i]->Delete();
	hxxd[i]->Delete();
}


hx0->Delete();

hxt->Delete();
hxt2->Delete();
hxt3->Delete();
hxt4->Delete();
hxt5->Delete();
hxn->Delete();

hs->Delete();

hsig1->Delete();
hsig2->Delete();
hsig3->Delete();
hsig4->Delete();

hnhit->Delete();
hthit->Delete();
hqhit->Delete();

hxxxx->Delete();
hxxxxs->Delete();

hqc->Delete();
htc->Delete();
hsc->Delete();
hpc->Delete();

hh1->Delete();
hh2->Delete();
hh3->Delete();
hh4->Delete();

dedx->Delete();
att->Delete();
tau->Delete();
sss->Delete();
}
printf("End %i\n",time(0)-tttt);

return;

}

void init(char fname1[], char fname2[]){

	FILE *fpc = fopen(fname1,"r");
	float cal,ecal;
	int ich;
	while(fscanf(fpc,"%i %f %f",&ich,&cal,&ecal)!=EOF){
		CAL[ich]=cal/22.0;
		ECAL[ich]=ecal;
	}
	CAL[0]=1;
	CAL[1]=1;
	//CAL[24]=1;
	//CAL[25]=1;
	//CAL[49]=1;
	//CAL[50]=1;
	CAL[74]=1;
	CAL[75]=1;
	fclose(fpc);

  FILE *fpt;
  fpt = fopen(fname2,"r");
//  fpt = fopen("proton_eventlist.txt","r");
 
  int it,st,et;
  int j1,j2;
  
  while(fscanf(fpt,"%i %i %i",&it,&st,&et)!=EOF){
//  while(fscanf(fpt,"%i %i %i %i %i",&it,&st,&et,&j1,&j2)!=EOF){
	ITre[itev]=it;
	STre[itev]=st;
	ETre[itev]=et;
	itev++;
	if (itev%100==0) printf("%i %i %i\n",it,st,et);
  }
	
}

void run(char fname1[], char fname2[], int irun, int size, int istart, int istop, int istep, int mode){

// init("calib42_new.dat","proton_eventlist.txt");
 init(fname1,fname2);
 char fname[100];
 for (int i=istart; i<=istop; i++){
	 for (int j=0; j<size; j++){
		 sprintf(fname,"PhysicsOct%i_%i",irun,i);
		 printf("%s %i\n",fname,j);
		ana(fname,j,10,istep,mode);
		}
	}
}



void readresult(char fname[]){

  //printf("%s\n",fname);
  FILE *fp2;
  fp2 = fopen(fname,"r");
  int iev,spnum,evnum,start,stop,ang,imatch,ngood,nhit;
  float tin,tout,G,sthit;
  fscanf(fp2, "%i %i %i %i %i %i %f\n",&iev,&spnum,&evnum,&imatch,&ngood,&nhit,&sthit);
  //printf("%i %i %i %i %i %i %f\n",iev,spnum,evnum,imatch,ngood,nhit,sthit);


int   ii[10000];
float XX[21];

	XX[0]=iev;
	XX[1]=spnum;
	XX[2]=evnum;
	XX[3]=ngood;
	XX[4]=nhit;
	XX[5]=sthit;

float loglka=0.;
float loglpi=0.;
  for (int i = 0;i<ngood;i++){
  	fscanf(fp2,"%i %f %f %f %f",&ii[i],&Z[i],&EZ[i],&T[i],&ET[i]);
	fscanf(fp2,"%f %f %f %f %f %f %f %f",&Q[i],&EQ[i],&S[i],&ES[i],&Q2[i],&EQ2[i],&T2[i],&Prob[i]);
 
	int j = (int)Z[i];
	XX[6]=Z[i];
	XX[7]=EZ[i];
	XX[8]=T[i];
	XX[9]=ET[i];
	XX[10]=Q[i];
	XX[11]=EQ[i];
	XX[12]=Q2[i];
	XX[13]=EQ2[i];
	XX[14]=S[i];
	XX[15]=ES[i];
	XX[16]=CAL[j];
	XX[17]=Prob[i];
	XX[18]=0.;
	XX[19]=0.;
	XX[20]=0.;
	n1->Fill(XX);
 }
  for (int i = 0;i<nhit;i++){
  	fscanf(fp2,"%i %f %f %f",&ii[i],&CHH[i],&THH[i],&SHH[i]);
  	//printf("%i %f %f %f\n",ii[i],CHH[i],THH[i],SHH[i]);
  	XX[18]=CHH[i];
  	XX[19]=THH[i];
  	XX[20]=SHH[i];
  	n1->Fill(XX);  	
  }

 
fclose(fp2);

  return;
}
void makeresult(char listname[]){
// init("calib42_new.dat","proton_eventlist.txt");
// init("calib42_new.dat","plist432.txt");
 init("calib52_kdata.dat","plist432.txt");

n1 = new TNtuple("n1","n1","idata:spnum:evnum:ngood:nhit:sthit:rch:erch:tdrift:etdrift:qcl:eqcl:qcl2:eqcl2:scl:escl:cal:prob:chhit:thit:shit");

  FILE *fp2;
  fp2 = fopen(listname,"r");

  char fname[100];
  while (fscanf(fp2,"%s",fname)!=EOF){

 readresult(fname);
 
 }
}
