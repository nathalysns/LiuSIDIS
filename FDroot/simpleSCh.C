void simpleSCh(char *a="",char *b="",char *c="",char *d="",char *e="")
{
  //Draw Feynman diagrams
  //Open Source ROOT Script Library
  //Author: Sinan Kuday
  
  TCanvas *c1 = new TCanvas("c1", "A canvas", 10,10, 500, 300);
  c1->SetFillColor(0);
  c1->Range(0, 0, 140, 60);
  Int_t linsav = gStyle->GetLineWidth();
  gStyle->SetLineWidth(3);
  TLatex t;
  t.SetTextAlign(20);
  t.SetTextSize(0.1);
  switch(a){
  case "zigzag":
    TCurlyLine *lt = new TCurlyLine(10,50,40,30); lt->SetWavy();
    break;
  case "curly":
    TCurlyLine *lt = new TCurlyLine(10,50,40,30);
    break;
  default:
    lt = new TArrow(10, 50, 40, 30, 0.03,"->|-");            break;  }   lt->Draw();
  
  switch (b) {
  case "zigzag":
    TCurlyLine *lb = new TCurlyLine(10, 10, 40, 30); lb->SetWavy();
    break;
  case "curly":
    TCurlyLine *lb = new TCurlyLine(10, 10, 40, 30);
    break;
  default:
    TArrow *lb; lb = new TArrow(10, 10, 40, 30, 0.03,"-|>-");
    break;
  }
  lb->Draw();
  
  t.DrawLatex(24,10,"e^{-}");
  t.DrawLatex(24,46,"e^{+}");
  
  switch (c) {
  case "zigzag":
    TCurlyLine *middle = new TCurlyLine(40, 30, 100, 30);
    middle->SetWavy();
    break;
  case "curly":
    TCurlyLine *middle = new TCurlyLine(40, 30, 100, 30);
    break;
  default:
    TArrow *middle; middle = new TArrow(40, 30, 100, 30, 0.03,"-|>-");
    break;
  }
  middle->Draw();
  
  t.DrawLatex(72,36,"#gamma");
  
  switch (d) {
  case "zigzag":
    TCurlyLine *rt = new TCurlyLine(100, 30, 130, 10);
    rt->SetWavy();
    break;
  case "curly":
    TCurlyLine *rt = new TCurlyLine(100, 30, 130, 10);
    break;
  default:
    TArrow *rt; rt = new TArrow(100, 30, 130, 10, 0.03,"-<|-");          break;  }   rt->Draw();
  
  switch (e) {
  case "zigzag":
    TCurlyLine *rb = new TCurlyLine(100, 30, 130, 50);
    rb->SetWavy();
    break;
  case "curly":
    TCurlyLine *rb = new TCurlyLine(100, 30, 130, 50);
    break;
  default:
    TArrow *rb; rb = new TArrow(100, 30, 130, 50, 0.03,"-<|-");          break;  }   rb->Draw();
  
  t.DrawLatex(113,10,"#bar{q}");
  t.DrawLatex(113,46,"q");
  
  c1->Update();
  gStyle->SetLineWidth(linsav);
}
