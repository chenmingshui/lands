#include <TStyle.h>
void FloridaStyle()
{
  TStyle* GloStyle;
  GloStyle = gStyle;

  TStyle* FloridaStyle = new TStyle("FloridaStyle", "FloridaStyle");
  gStyle = FloridaStyle;

  //----------------------------------------------------------------------------
  // Canvas
  //----------------------------------------------------------------------------
  FloridaStyle->SetCanvasBorderMode(  0);
  FloridaStyle->SetCanvasBorderSize( 10);
  FloridaStyle->SetCanvasColor     (  0);
  FloridaStyle->SetCanvasDefH      (500);
  FloridaStyle->SetCanvasDefW      (700);
  FloridaStyle->SetCanvasDefX      ( 10);
  FloridaStyle->SetCanvasDefY      ( 10);

  //----------------------------------------------------------------------------
  // Pad
  //----------------------------------------------------------------------------
  FloridaStyle->SetPadBorderMode  (   0);
  FloridaStyle->SetPadBorderSize  (  10);
  FloridaStyle->SetPadColor       (   0);
  FloridaStyle->SetPadBottomMargin(0.20);
  FloridaStyle->SetPadTopMargin   (0.08);
  FloridaStyle->SetPadLeftMargin  (0.18);
  FloridaStyle->SetPadRightMargin (0.05);

  //----------------------------------------------------------------------------
  // Frame
  //----------------------------------------------------------------------------
  FloridaStyle->SetFrameFillStyle ( 0);
  FloridaStyle->SetFrameFillColor ( 0);
  FloridaStyle->SetFrameLineColor ( 1);
  FloridaStyle->SetFrameLineStyle ( 0);
  FloridaStyle->SetFrameLineWidth ( 2);
  FloridaStyle->SetFrameBorderMode( 0);
  FloridaStyle->SetFrameBorderSize(10);

  //----------------------------------------------------------------------------
  // Hist
  //----------------------------------------------------------------------------
  FloridaStyle->SetHistFillColor(0);
  FloridaStyle->SetHistFillStyle(1);
  FloridaStyle->SetHistLineColor(1);
  FloridaStyle->SetHistLineStyle(0);
  FloridaStyle->SetHistLineWidth(1);

  //----------------------------------------------------------------------------
  // Axis
  //----------------------------------------------------------------------------
  FloridaStyle->SetLabelFont  (    42, "x");
  FloridaStyle->SetLabelOffset( 0.015, "x");
  FloridaStyle->SetLabelSize  ( 0.050, "x");
  FloridaStyle->SetTickLength (-0.015, "x");
  FloridaStyle->SetTitleFont  (    42, "x");
  FloridaStyle->SetTitleOffset( 1.400, "x");
  FloridaStyle->SetTitleSize  ( 0.050, "x");

  FloridaStyle->SetLabelFont  (    42, "y");
  FloridaStyle->SetLabelOffset( 0.015, "y");
  FloridaStyle->SetLabelSize  ( 0.050, "y");
  FloridaStyle->SetTickLength (-0.015, "y");
  FloridaStyle->SetTitleFont  (    42, "y");
  FloridaStyle->SetTitleOffset( 1.200, "y");
  FloridaStyle->SetTitleSize  ( 0.050, "y");

  FloridaStyle->SetNdivisions (   505, "x");

  //----------------------------------------------------------------------------
  // Title
  //----------------------------------------------------------------------------
  FloridaStyle->SetTitleBorderSize( 0);
  FloridaStyle->SetTitleFillColor (10);
  FloridaStyle->SetTitleFont      (42, "");

  //----------------------------------------------------------------------------
  // Stat
  //----------------------------------------------------------------------------
  FloridaStyle->SetOptStat       (1110);
  FloridaStyle->SetOptFit 	 (1111);
  FloridaStyle->SetStatBorderSize(   0);
  FloridaStyle->SetStatColor     (  10);
  FloridaStyle->SetStatFont      (  42);
  FloridaStyle->SetStatX         (0.94);
  FloridaStyle->SetStatY         (0.91);

  return;
}
