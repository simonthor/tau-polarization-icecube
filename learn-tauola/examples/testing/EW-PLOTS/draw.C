int N=1000;
TString format(".png"); // format(".eps");

TString n1("f-sanc.txt"),n2("f-born.txt"),n3("f-plzap0.txt");
TString n4("f-err.txt"),n5("f-cross.txt");
TString n6("f-ww0.txt"),n7("f-w.txt"),n8("f-w0.txt");

TGraph* plot(TCanvas &c,TString &name, ifstream &f)
{
	double s,val;
	TGraph *g = new TGraph(N);
	name.Remove(name.Length()-4);
	name+=format;
	cout<<"Ploting: "<<name<<endl;
	for(int i=0;i<N;i++)
	{
		f>>s>>val;
		g->SetPoint(i,s,val);
	}
	g->SetTitle(NULL);
	g->SetFillColor(10);
	g->SetLineWidth(1);
	g->Draw("AL");
	c.SaveAs(name);
	return g;
}

int draw()
{
	ifstream f1(n1),f2(n2),f3(n3);
	ifstream f4(n4),f5(n5);
	ifstream f6(n6),f7(n7),f8(n8);
	TCanvas c1("c","",1600,2200);
	c1.cd();
	c1.SetFillColor(10);
	c1.SetLogx(1);
	TGraph *g1=NULL,*g2=NULL,*g3=NULL;
	TGraph *g4=NULL,*g5=NULL;
	if(f1.is_open()) g1=plot(c1,n1,f1);
	if(f2.is_open()) g2=plot(c1,n2,f2);
	if(f3.is_open()) g3=plot(c1,n3,f3);
	c1.SetLogx(0);
	c1.SetLogy(0);
	if(f4.is_open()) plot(c1,n4,f4);
	c1.SetLogx(1);
	if(f5.is_open()) plot(c1,n5,f5);
	if(f6.is_open()) plot(c1,n6,f6);
	c1.SetLogy(1);
	if(f7.is_open()) g4=plot(c1,n7,f7);
	if(f8.is_open()) g5=plot(c1,n8,f8);
	if(g1 && g2 && g3)
	{
		TString name("f-figure1-tables-all");
		TCanvas c2("c2","",1600,1800);
		c2.cd();
		c2.SetFillColor(10);
		c2.SetLogx(1);
		c2.SetLogy(0);
		name+=format;
		cout<<"Ploting: "<<name<<endl;
		TAxis *x = g3->GetXaxis();
		x->SetTitle("sqrt(s) GEV");
		g1->SetLineColor(kRed);
		g2->SetLineColor(kBlue);
		g3->Draw("AL");
		g2->Draw("sameL");
		g1->Draw("sameL");
		c2.SaveAs(name);
	}
	if(g4 && g5)
	{
		TString name("f-figure2-weights-both");
		name+=format;
		cout<<"Ploting: "<<name<<endl;
		c1.cd();
		c1.SetLogy(1);
		TAxis *x = g5->GetXaxis();
		x->SetTitle("sqrt(s) GEV");
		g4->SetLineColor(kRed);
		g5->SetLineColor(kBlue);
		c1.DrawFrame(19,0.000005,19000,1900);
		g5->Draw("L");
		g4->Draw("sameL");
		c1.SaveAs(name);
	}
	exit(0);
}
