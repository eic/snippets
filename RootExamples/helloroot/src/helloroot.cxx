#include<helloroot/helloroot.hh>

#include<iostream>
#include<string>

#include<TH1D.h>
#include<TPad.h>

void HelloRoot()
{
  std::cout << "Hello from HelloRoot" << std::endl;

  // do something with root
  TH1D h("h", "h", 100, -5, 5);
  h.FillRandom("gaus", 1000);
  h.Draw();
  gPad->SaveAs("hello.png");

  return;
}