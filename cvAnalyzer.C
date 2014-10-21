#include <fstream>
#include <algorithm>
#include <string>
#include <cstdio> 
#include <sstream>
#include <TString.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TF1.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMath.h>
// cb 
#include "Math/WrappedMultiTF1.h"
#include "Math/FitMethodFunction.h"
#include "HFitInterface.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
// error
#include "TVirtualFitter.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
  // cleasses
  class  pol0 {
 public:
   Double_t Evaluate (double *x, double *p) { return (p[0]);}
};
  class  pol1 {
 public:
   Double_t Evaluate (double *x, double *p) { return (p[0]+p[1]*x[0]);}
};
  class  gauss {
 public:
   Double_t Evaluate (double *x, double *par) { return (par[0]*TMath::Exp(-0.5*((x[0]-par[1])/par[2]) * ((x[0]-par[1])/par[2])) );}
};
class  gaussPol0 {
 public:
   Double_t Evaluate (double *x, double *par) { return (par[3] + par[0]*TMath::Exp(-0.5*((x[0]-par[1])/par[2]) * ((x[0]-par[1])/par[2])) );}
};
  class FitResult {
  public:
    FitResult(TString NameIn);
    void SetDepletionVoltage (Double_t depVoltage);
    void SetDepletionVoltageError(Double_t depVoltageError){depVoltageError_=depVoltageError;}
    void SetcGeo(Double_t cGeo);
    TString Name() {return name_;}
    Double_t GetDepletionVoltage(){return depVoltage_;}
    Double_t GetDepletionVoltageError(){return depVoltageError_;}
    Double_t GetcGeo(){return cGeo_;}
    Double_t GetDDetector();
    Double_t GetNeff();
    void SetArea(Double_t x ,Double_t y);
    void SetLowestSecondDerivativeVoltage(Double_t value){lowestSecondDerivativeVoltage_=value;lowestSecondDerivativeVoltageBool_=true;}
    Double_t GetLowestSecondDerivativeVoltage(){return lowestSecondDerivativeVoltage_;}
    void SetRightSideSlopeAvargeHight(Double_t value){rightSideSlopeAvargeHight_=value;rightSideSlopeAvargeHightBool_=true;}
    Double_t GetRightSideSlopeAvargeHight(){return rightSideSlopeAvargeHight_;}    
    void SetLeftSideSlopeAvargeHight(Double_t value){leftSideSlopeAvargeHight_=value;leftSideSlopeAvargeHightBool_=true;}
    Double_t GetLeftSideSlopeAvargeHight(){return leftSideSlopeAvargeHight_;}  
    Double_t FullGaussMean_, RestrictedGaussMean_, RestrictedGaussSigma_,highestSlopeVoltage_,PeakSlopeVoltage_;
    vector<Double_t> RestrictedGaussResult_;
    vector<Double_t> RestrictedGaussResultError_;
    void SetUsedMethod(TString UsedMethod){if(!usedMethodSet_)usedMethod_=UsedMethod+"_";else usedMethod_+=UsedMethod+"_"; usedMethodSet_=true;}
    TString GetUsedMethod () {if(usedMethodSet_)return usedMethodSet_; else return "NOTFITTED!";}
  private:
	  TString name_, usedMethod_;
    Double_t depVoltage_, depVoltageError_, cGeo_, area_;
    bool setDepVoltage_, setcGeo_, setArea_,lowestSecondDerivativeVoltageBool_,rightSideSlopeAvargeHightBool_,leftSideSlopeAvargeHightBool_, usedMethodSet_;
    Double_t epsilon0;
    Double_t electriChage;
    Double_t epsilonSi;
    Double_t lowestSecondDerivativeVoltage_, rightSideSlopeAvargeHight_,leftSideSlopeAvargeHight_;
	  
};
FitResult::FitResult(TString NameIn)
{
	name_=NameIn;
	depVoltage_=-10;
	depVoltageError_=-666;
	setDepVoltage_=false;
	setcGeo_=false;
	setArea_=false;
	epsilon0 = 8.85418781762e-12;
    	electriChage = 1.602176565e-19;
    	epsilonSi = 11.9;
	lowestSecondDerivativeVoltageBool_=false;
	rightSideSlopeAvargeHightBool_=false;
	leftSideSlopeAvargeHightBool_=false;
	FullGaussMean_=0;
	highestSlopeVoltage_=0;
	PeakSlopeVoltage_=0;
	usedMethodSet_=false;
	usedMethod_="";
}
void FitResult::SetDepletionVoltage(Double_t depVoltage)
{
	depVoltage_=depVoltage;
	setDepVoltage_=true;
}
void FitResult::SetArea(Double_t x ,Double_t y)
{
	area_ = x * y * 1e-6;
	setArea_=true;
}
void FitResult::SetcGeo(Double_t cGeo)
{
	cGeo_=1/sqrt(cGeo);
	setcGeo_=true;
}
Double_t FitResult::GetDDetector()
{
	if(!setcGeo_)
	{
		std::cout<<"Error:: Trying to calculated GetDdetector while cGeo not set!!"<<std::endl;
		return 0;
	}
	if(!setArea_)
	{
		std::cout<<"Error:: Trying to calculated GetDdetector while area not set!!"<<std::endl;
		return 0;
	}
	return epsilon0 * epsilonSi * area_ / cGeo_;
}
Double_t FitResult::GetNeff()
{
	if(!setDepVoltage_)
	{
		std::cout<<"Error:: Trying to calculated GetNeff while depVoltage not set!!"<<std::endl;
		return 0;
	}
	if(!setcGeo_)
	{
		std::cout<<"Error:: Trying to calculated GetNeff while cGeo not set!!"<<std::endl;
		return 0;
	}
	if(!setArea_)
	{
		std::cout<<"Error:: Trying to calculated GetNeff while area not set!!"<<std::endl;
		return 0;
	}
	return depVoltage_ * epsilon0 * epsilonSi * 2 / (epsilon0 * epsilonSi * area_ / cGeo_ * epsilon0 * epsilonSi * area_ / cGeo_ * electriChage);
}
  class Measurment {
  public:
	  Measurment(TString FileName);
	  ~Measurment();
	  void SetDir(TDirectory *tDir){tDir_=tDir;tDirBool_=true;}
	  void setStartingVoltageToAnalyze(Double_t starting){startingVoltageToAnalyze_=starting; startingVoltageToAnalyzeBool_=true;}
	  TString GetFileName(){return FileName_;}
	  bool LoadInputFromTxt();
	  bool LoadFrequencies();
	  bool FillGraph();
	  bool SaveResults();
	  bool SortByFrequencies();
	  bool Fit();
	  bool FitReloaded();
	  // variables
	  bool debug_;
	  string lineBeforeBeginOfTable_, lineAfterTable_;
	  vector<FitResult*> fitResults_;
	  vector<string> FrequencieTokens_;
	  vector<double> frequenciesDouble_;
	  // fit functions

  private:
	  TString FileName_;
	  TDirectory *tDir_;
	  bool tDirBool_;
	  ifstream *cvFile_;
	  ofstream *textOutPut_;
	  bool cvFileBool_;
	  unsigned int numberOfFrequencies_;
	  
	  
	  bool frequenciesBool_;
	  vector<TGraph> cvGraph_, cvGraphError_, cvSlopeGraph_,cvSecondDerivativeGraph_;
	  bool SaveResultsBool_;
	  bool SortByFrequenciesBool_;
	  vector<double> depVoltage_, depVoltageError_;
	  
	  Double_t startingVoltageToAnalyze_;
	  bool startingVoltageToAnalyzeBool_;
	  bool depVoltageBool_;
	  bool useNonHighIradiatedCriteria_;
	  // functions
	  bool SaveResults(vector<TGraph> cvGraphs,vector<TGraph> cvGraphsError, TString additionalName);
	  
	  // tempvariables
	  string line;
	  int n_;
	  char buffer_[200];
	  TString TTemp_;
  };
Measurment::Measurment(TString FileName)
{
	FrequencieTokens_.clear();
	FileName_=FileName;
	tDirBool_=false;
	SaveResultsBool_=false;
	SortByFrequenciesBool_=false;
	cvFileBool_=false;
	numberOfFrequencies_=0;
	frequenciesBool_=false;
	depVoltageBool_=false;
	useNonHighIradiatedCriteria_=false;
	startingVoltageToAnalyzeBool_=false;
	textOutPut_ = new ofstream;
	textOutPut_->open(FileName+".txt");
}
Measurment::~Measurment()
{
	if(!tDirBool_)
	{
		std::cout<<"No directory has been set at destructor level. Results will not be saved to a root file!"<<std::endl;
		(*textOutPut_)<<"No directory has been set at destructor level. Results will not be saved to a root file!"<<"\n";
	}
	if(debug_)std::cout<<"Measurment class with name: "<<FileName_<<", is being deleted."<<std::endl;
	textOutPut_->close();	
}
bool Measurment::LoadInputFromTxt()
{
	string name = (string)FileName_;
	if(name.find(".cv") !=std::string::npos) cvFile_ = new ifstream(FileName_);
	else cvFile_ = new ifstream(FileName_+".cv");
	if(cvFile_->is_open() )
	{
		cvFileBool_=true;
		return true;
	}
	else return false;
}
bool Measurment::LoadFrequencies()
{
	if (cvFile_->is_open())
	{
		while ( getline (*cvFile_,line) )
		{
			if(line.find(":List of frequencies") !=std::string::npos)
			{
				(*textOutPut_)<<line<<"\n";
				getline (*cvFile_,line);
				(*textOutPut_)<<line<<"\n";
				const string& delimiters = ",";
				string::size_type lastPos = line.find_first_not_of(delimiters, 0);
				string::size_type pos     = line.find_first_of(delimiters, lastPos);
				while (string::npos != pos || string::npos != lastPos)
				{
					if(line.substr(lastPos, pos - lastPos).find("Hz")!=std::string::npos)FrequencieTokens_.push_back(line.substr(lastPos, pos - lastPos));
					if(line.substr(lastPos, pos - lastPos).find("hz")!=std::string::npos)FrequencieTokens_.push_back(line.substr(lastPos, pos - lastPos));
					lastPos = line.find_first_not_of(delimiters, pos);
					pos = line.find_first_of(delimiters, lastPos);
				}
				if(debug_)for(unsigned int i=0; i < FrequencieTokens_.size(); i++) std::cout<<"FrequenciesTokens["<<i<<"]: "<<FrequencieTokens_[i]<<std::endl;
			}
		}
		for(unsigned int i=0; i<FrequencieTokens_.size();i++)
		{
			double value=0;
			if (sscanf(FrequencieTokens_[i].c_str(), "%lf", &value) != 0)
			{
				if((FrequencieTokens_[i].find("kHz") !=std::string::npos))frequenciesDouble_.push_back(value*1000);
				else frequenciesDouble_.push_back(value);
			}
		}
		for(unsigned int i=0; i < frequenciesDouble_.size(); i++)(*textOutPut_)<<"Doubles["<<i<<"]: "<<frequenciesDouble_[i]<<" Hz"<<std::endl;
		if(debug_)for(unsigned int i=0; i < frequenciesDouble_.size(); i++) std::cout<<"getFrequenciesDoubles doubles["<<i<<"]: "<<frequenciesDouble_[i]<<" Hz"<<std::endl;
	}
	else std::cout<<"getFrequencies::Error cvFile is not open!!"<<std::endl;
	numberOfFrequencies_=FrequencieTokens_.size();
	cvFile_->clear();
	cvFile_->seekg(0,ios::beg);
	if(numberOfFrequencies_)return true;
	else return false;
}
bool Measurment::FillGraph()
{
	vector<vector<double> > VCInput;
	vector<double>VC;
	if (cvFile_->is_open())
	{
		bool input=false;
		while ( getline (*cvFile_,line) )
		{
			if(line.find(lineBeforeBeginOfTable_) !=std::string::npos){ getline(*cvFile_,line); input=true;}
			if(line.find(lineAfterTable_) !=std::string::npos){input=false;}
			if(input)
			{
				if(debug_)std::cout<<line<<std::endl;
				unsigned int i=0;
				if(debug_)std::cout<<line<<std::endl;
				istringstream iss(line);   
				double value;
				
				while (iss >> value) {
					
					if(i!=1 && i< (FrequencieTokens_.size()+2) )VC.push_back(value);
					i++;
				} 
				VCInput.push_back(VC);
				VC.clear();
			}
		}
	}
	if(debug_) for(unsigned int i=0; i< VCInput.size();i++)
	{
		std::cout<<"Line["<<i<<"] has entries["<<VCInput[i].size()<<"]: value: ";
		for(unsigned int ii=0; ii<VCInput[i].size();ii++)
		{
			std::cout<<VCInput[i][ii]<<", ";
		}
		std::cout<<std::endl;
	}
	const unsigned int measurments = VCInput.size();
	Double_t OneOverCSquaredTGraph[measurments];
	Double_t VoltagesTGraph[measurments];
	for(unsigned int i=0; i<(FrequencieTokens_.size());i++)
	{
		for (unsigned int ii=0; ii<measurments;ii++)
		{
			OneOverCSquaredTGraph[ii]=1/(VCInput[ii][i+1]*VCInput[ii][i+1]);
			VoltagesTGraph[ii]=VCInput[ii][0];
			if(debug_ && false)std::cout<<"OneOverCSquaredTGraph["<<i<<"]: "<<OneOverCSquaredTGraph[i]<<", VoltagesTGraph: "<<VoltagesTGraph[i]<<std::endl;
		}
		TGraph *tempGraph = new TGraph(measurments,VoltagesTGraph,OneOverCSquaredTGraph);
		tempGraph->SetName("CV");
		TTemp_="CV "+FrequencieTokens_[i]+"; voltage [V];1/C^{2}";
		tempGraph->SetTitle(TTemp_);
		cvGraph_.push_back(*tempGraph );
		delete tempGraph;
	}
	for (unsigned int i=0; i<cvGraph_.size();i++)
	{
		const int size = cvGraph_[i].GetN();
		if(size<2)
		{
			std::cout<<"ERROR cvGraph["<<i<<"]: contains less than 2 entries slope can not be computed!!!!!"<<std::endl;
			break;
		}
		Double_t* X = cvGraph_[i].GetX();
		Double_t* Y= cvGraph_[i].GetY();
		Double_t slope[size-1];
		Double_t newX[size-1];
		for(int ii=0; ii<(size-1);ii++)
		{
			slope[ii]=(Y[ii+1]-Y[ii])/(X[ii+1]-X[ii]);
			newX[ii]=X[ii]+0.5 * (X[ii+1]-X[ii]);
			
		}
		TGraph *tempGraph = new TGraph(size-1,newX,slope);
		TTemp_="CV slopes"+FrequencieTokens_[i]+"; voltage [V];1/C^{2}";
		tempGraph->SetTitle(TTemp_);
		cvSlopeGraph_.push_back(*tempGraph );
		delete tempGraph;
	}
	for (unsigned int i=0; i<cvSlopeGraph_.size();i++)
	{
		const int size = cvSlopeGraph_[i].GetN();
		if(size<2)
		{
			std::cout<<"ERROR cvSlopeGraph_["<<i<<"]: contains less than 2 entries slope can not be computed!!!!!"<<std::endl;
			break;
		}
		Double_t* X = cvSlopeGraph_[i].GetX();
		Double_t* Y= cvSlopeGraph_[i].GetY();
		Double_t slope[size-1];
		Double_t newX[size-1];
		for(int ii=0; ii<(size-1);ii++)
		{
			slope[ii]=(Y[ii+1]-Y[ii])/(X[ii+1]-X[ii]);
			newX[ii]=X[ii]+0.5 * (X[ii+1]-X[ii]);
			
		}
		TGraph *tempGraph = new TGraph(size-1,newX,slope);
		TTemp_="CV second derivative"+FrequencieTokens_[i]+"; voltage [V];1/C^{2}";
		tempGraph->SetTitle(TTemp_);
		cvSecondDerivativeGraph_.push_back(*tempGraph );
		delete tempGraph;
	}
	return true;
}
bool Measurment::SaveResults(vector<TGraph> cvGraphs, vector<TGraph> cvGraphsError, TString Name)
{
	line= Name;
	bool isSecondDerivative=false;
	if(line.find("SecondDerivatve") !=std::string::npos)isSecondDerivative=true;
	bool isSlope=false;
	if(line.find("Slope") !=std::string::npos)isSlope=true;
	tDir_->cd();
	TTemp_+= ";voltage [V];1/C^{2}";
	TCanvas *combinedTCanvas = new TCanvas(Name,Name+TTemp_,600,600);
	TLegend *lLeg = new TLegend(0.45,0.3,0.85,0.5);
	lLeg->SetFillColor(0);
	for(unsigned int i=0; i <cvGraphs.size(); i++)
	{	
		TTemp_=cvGraphs[i].GetName();
		TTemp_+=FrequencieTokens_[i]+" ";
		TTemp_+=Name;
		cvGraphs[i].SetTitle(TTemp_);
		cvGraphs[i].SetName(TTemp_);
		TCanvas *cTemp = new TCanvas("CV_"+Name,Name,600,600);
		tDir_->cd();
		cvGraphs[i].SetMarkerColor(i+1);
		cvGraphs[i].SetMarkerSize(1);
		cvGraphs[i].SetLineColor(i+1);
		cvGraphs[i].Write();
		cTemp->cd();
		cvGraphs[i].Draw("AC*");
		cTemp->Update();
		TTemp_=Name;
		TTemp_+=FrequencieTokens_[i];
		TTemp_+=".pdf";
		TTemp_=FrequencieTokens_[i];
		if(depVoltageBool_)TTemp_ += buffer_;
		//cTemp->SaveAs(TTemp_);
		combinedTCanvas->cd();
		if(depVoltageBool_)n_ = sprintf(buffer_,", V_{dep}= %.1f+-%.1f [V]",depVoltage_[i],depVoltageError_[i]);
		if(isSecondDerivative)
		{
			n_ = sprintf(buffer_,", V_{MinSecDeriv}= %.1f [V]",fitResults_[i]->GetLowestSecondDerivativeVoltage());
			TTemp_+=buffer_;
		}
		if(isSlope)
		{
			n_ = sprintf(buffer_,",GMean=%.0f[V],V_{MinSecDer}=%.0f[V],PeakSlope=%.0f",fitResults_[i]->FullGaussMean_,fitResults_[i]->GetLowestSecondDerivativeVoltage(),fitResults_[i]->PeakSlopeVoltage_);
			TTemp_+=buffer_;
		}
		if(!isSlope && !isSecondDerivative && fitResults_[i]->GetDepletionVoltage()>0)
		{
			n_ = sprintf(buffer_,", V_{dep}= %.2f+-%.2f [V] ",fitResults_[i]->GetDepletionVoltage(),fitResults_[i]->GetDepletionVoltageError());
			TTemp_+=buffer_;
		}
		if(i==0)
		{
			lLeg->AddEntry(&cvGraphs[i],TTemp_,"l");
			TTemp_=FrequencieTokens_[i];
			TTemp_+= ";voltage [V];1/C^{2}";
			cvGraphs[i].SetTitle(TTemp_);
			cvGraphs[i].SetName(FileName_+Name);
			cvGraphs[i].Draw("AP");
		}
		if(i>0)
		{
			lLeg->AddEntry(&cvGraphs[i],TTemp_,"l");
			cvGraphs[i].Draw("PSame");
		}
		combinedTCanvas->Update();
		delete cTemp;
	}
		
	tDir_->cd();
	combinedTCanvas->SetTitle(FileName_+" "+Name);
	combinedTCanvas->SetName(FileName_+" "+Name);
	combinedTCanvas->cd();
	lLeg->Draw();
	combinedTCanvas->Update();
	combinedTCanvas->Write();
	TTemp_=FileName_+Name;
	TTemp_+="_CombinedPlot.pdf";
	combinedTCanvas->SaveAs(TTemp_);
	delete combinedTCanvas;
	SaveResultsBool_=true;
	return SaveResultsBool_;
}
bool Measurment::SaveResults()
{
	bool result =false;
	result = SaveResults(cvGraph_,cvGraphError_,"");
	if(!result)return false;
	result = SaveResults(cvSlopeGraph_,cvGraphError_,"Slope");
	if(!result) return false;
	result = SaveResults(cvSecondDerivativeGraph_,cvGraphError_,"SecondDerivatve");
	return result;
}

bool Measurment::SortByFrequencies()
{
	bool keepgoing=true;
	if(frequenciesDouble_.size()<2) keepgoing=false;
	if(debug_)for(unsigned int i=0; i<frequenciesDouble_.size();i++) std::cout<<"sortByFrequencies::Element["<<i<<"] before rearanging: "<<frequenciesDouble_[i]<<std::endl;
	while (keepgoing)
	{
		keepgoing=false;
		for(unsigned int i=0; i< frequenciesDouble_.size()-1; i++)
		{
			if( frequenciesDouble_[i] > frequenciesDouble_[i+1])
			{
				double temp = frequenciesDouble_[i+1];
				frequenciesDouble_[i+1]=frequenciesDouble_[i];
				frequenciesDouble_[i]=temp;
				keepgoing=true;
				TGraph tempG = cvGraph_[i+1];
				cvGraph_[i+1]=cvGraph_[i];
				cvGraph_[i]=tempG;
				TGraph tempSG = cvSlopeGraph_[i+1];
				cvSlopeGraph_[i+1]=cvSlopeGraph_[i];
				cvSlopeGraph_[i]=tempSG;
				string tempS = FrequencieTokens_[i+1];
				FrequencieTokens_[i+1]=FrequencieTokens_[i];
				FrequencieTokens_[i]=tempS;
				if(fitResults_.size()>i)
				{
					FitResult *tempResult = fitResults_[i+1];
					fitResults_[i+1]=fitResults_[i];
					fitResults_[i]=tempResult;
				}
			}
		}
	}
	for(unsigned int i=0; i<frequenciesDouble_.size();i++) (*textOutPut_)<<"sortByFrequencies::Element["<<i<<"] after rearanging: "<<frequenciesDouble_[i]<<std::endl;
	if(debug_)for(unsigned int i=0; i<frequenciesDouble_.size();i++) std::cout<<"sortByFrequencies::Element["<<i<<"] after rearanging: "<<frequenciesDouble_[i]<<std::endl;
	SortByFrequenciesBool_=true;
	return SortByFrequenciesBool_;
}
bool Measurment::Fit()
{
	bool result=true;
	
	for (unsigned int i=0; i< cvGraph_.size();i++)
	{
		// find first bin with positve slope starting from left
		fitResults_.push_back(new FitResult(FrequencieTokens_[i]));
		int firstBinToUse=0;
		Double_t firstVoltageToUse=0;
		Double_t lastVoltageToUseForLeftSide=0;
		Double_t lowestY=1e30;
		Double_t highestY=0;
		double *xSlope = cvSlopeGraph_[i].GetX();
		Double_t *ySlope = cvSlopeGraph_[i].GetY();
		int NumberPointsSlopeGraph = cvSlopeGraph_[i].GetN();
		for (int ii=0; ii<NumberPointsSlopeGraph;ii++)
		{
			if(ySlope[ii]>0)
			{
				if(ySlope[ii]<lowestY)lowestY=ySlope[ii];
				if(ySlope[ii]>highestY)highestY=ySlope[ii];
				firstBinToUse=ii;
				firstVoltageToUse=xSlope[ii];
				if(ii+6 < NumberPointsSlopeGraph)lastVoltageToUseForLeftSide=xSlope[ii+6];
				else 
				{
					std::cout<<"Warning not enough points for fitting found in first derivatve. Needs at 7 points starting from left with a slope >0! fit will most probably not work...."<<std::endl;
					(*textOutPut_)<<"Warning not enough points for fitting found in first derivatve. Needs at 7 points starting from left with a slope >0!"<<"\n";
					
				}
				(*textOutPut_)<<"Left side fitting on first derivative is done from Voltage "<<firstVoltageToUse<<" to "<<lastVoltageToUseForLeftSide<<"\n";
				break;
			}
		}
		pol1 *Pol1 = new pol1();
		// left side fit
		TF1 * fPol1Left = new TF1("constLeft",Pol1,&pol1::Evaluate,firstVoltageToUse,lastVoltageToUseForLeftSide,2,"LeftPol1","LeftPol1");
		fPol1Left->SetLineColor(i+1);
		fPol1Left->SetLineStyle(10);
		cvSlopeGraph_[i].Fit(fPol1Left,"NRq+");
		TF1 * fPol1LeftFull = new TF1("constLeft",Pol1,&pol1::Evaluate,xSlope[0],xSlope[NumberPointsSlopeGraph-1],2,"LeftPol1","LeftPol1");
		fPol1LeftFull->SetLineColor(i+1);
		fPol1LeftFull->SetLineStyle(10);
		fPol1LeftFull->FixParameter(0,fPol1Left->GetParameter(0));
		fPol1LeftFull->FixParameter(1,fPol1Left->GetParameter(1));
		cvSlopeGraph_[i].Fit(fPol1LeftFull,"q+");
		// right side fit
		TF1 * fPol1Right = new TF1("constLeft",Pol1,&pol1::Evaluate,xSlope[NumberPointsSlopeGraph-8],xSlope[NumberPointsSlopeGraph-1],2,"RightPol1","RightPol1");
		fPol1Right->SetLineColor(i+1);
		fPol1Right->SetLineStyle(2);
		cvSlopeGraph_[i].Fit(fPol1Right,"NRq+");
		TF1 * fPol1RightFull = new TF1("constLeft",Pol1,&pol1::Evaluate,xSlope[0],xSlope[NumberPointsSlopeGraph-1],2,"RightPol1","RightPol1");
		fPol1RightFull->SetLineColor(i+1);
		fPol1RightFull->SetLineStyle(2);
		fPol1RightFull->FixParameter(0,fPol1Right->GetParameter(0));
		fPol1RightFull->FixParameter(1,fPol1Right->GetParameter(1));
		cvSlopeGraph_[i].Fit(fPol1RightFull,"q+");
		// find point at which the turing point is
		int peakBin=0;
		Double_t peakVoltage=0;
		Double_t peakYValue=1E40;
		double *xSecondDerivative = cvSecondDerivativeGraph_[i].GetX();
		Double_t *ySecondDerivative = cvSecondDerivativeGraph_[i].GetY();
		int NumberPointsSecondDerivativeGraph = cvSecondDerivativeGraph_[i].GetN();
		for (int ii=0; ii<(NumberPointsSecondDerivativeGraph-2);ii++)
		{
			if(ii<firstBinToUse) continue;
			if(peakYValue>ySecondDerivative[ii])
			{
				//std::cout<<"True since: "<<peakYValue<<" > "<<ySecondDerivative[ii]<<endl;
				peakBin=ii;
				peakVoltage=xSecondDerivative[ii];
				peakYValue=ySecondDerivative[ii];
			}
		}
		(*textOutPut_)<<"lowest value of second derivative (should be the peak turning point of the cv curve): "<<peakBin<<"]:"<<peakVoltage<<"\n";
		//std::cout<<"lowest value of second derivative (should be the peak turning point of the cv curve): "<<peakBin<<"]:"<<peakVoltage<<"\n";
		// find peak voltage of the slope
		Double_t peakVoltageOfSlope=0;
		Double_t peakSlope=0;
		int peakBinOfSlope=0;
		int startingBinOfHighestSlopeCalc=peakBin-8;
		int endBinOfHighestSlopeCalc=peakBin+8;
		if(startingBinOfHighestSlopeCalc<0) std::cout<<"ERROR starting bin of HighestSlopeCalc is smaller then 0 ... CRASH"<<std::endl;
		if(startingBinOfHighestSlopeCalc>NumberPointsSlopeGraph)std::cout<<"ERROR starting bin of HighestSlopeCalc is bigger then total amount of bins ... CRASH"<<std::endl;
		for (int ii=startingBinOfHighestSlopeCalc; ii<endBinOfHighestSlopeCalc;ii++)
		{
			if(ySlope[ii]>peakSlope)
			{
				peakVoltageOfSlope=xSlope[ii];
				peakSlope=ySlope[ii];
				peakBinOfSlope=ii;
			}
		}
		(*textOutPut_)<<"Peak of first derivative (gauss peak): "<<peakBinOfSlope<<"]:"<<peakVoltageOfSlope<<"\n";
		std::cout<<"Peak of first derivative (gauss peak): "<<peakBinOfSlope<<"]:"<<peakVoltageOfSlope<<"\n";
		Double_t highestX = lastVoltageToUseForLeftSide;
		Double_t lowestX = firstVoltageToUse;
		gauss *Gaus = new gauss();
		TF1 * fFullGauss = new TF1("constLeft",Gaus,&gauss::Evaluate,firstVoltageToUse,lastVoltageToUseForLeftSide,3);
		fFullGauss->SetLineColor(i+1);
		fFullGauss->SetLineStyle(2);
		fFullGauss->SetParameters(highestY,peakVoltageOfSlope,(highestX-lowestX)*0.6);
		(*textOutPut_)<<"Used gauss fit values: p0: "<<lowestY<<", p1: "<<peakVoltageOfSlope<<", p2: "<<(highestX-lowestX)*0.6<<"\n";
		std::cout<<"Used gauss fit values: p0: "<<lowestY<<", p1: "<<peakVoltageOfSlope<<", p2: "<<(highestX-lowestX)*0.6<<"\n";
		cvSlopeGraph_[i].Fit(fFullGauss,"R+");
	}
	return result;
}
bool Measurment::FitReloaded()
{
	Double_t startingVoltage=0.1;
	if (!startingVoltageToAnalyzeBool_)
	{
		(*textOutPut_)<<"Trying to start the fitting of: "<<FileName_<<" but startingVoltageToAnalyze_ not set!! fit will use full range might be a problem if very low voltages did not show nice behaviour!\n";
		std::cout<<"Trying to start the fitting of: "<<FileName_<<" but startingVoltageToAnalyze_ not set!! fit will use full range might be a problem if very low voltages did not show nice behaviour!\n";
	}
	else
	{
		(*textOutPut_)<<"StartingVoltageToAnalayze set to: "<<startingVoltageToAnalyze_<<"[V]\n";
		std::cout<<"StartingVoltageToAnalayze set to: "<<startingVoltageToAnalyze_<<"[V]\n";
		startingVoltage=startingVoltageToAnalyze_;
	}
	bool result =true;
	(*textOutPut_)<<"Starting fit procedure:"<<"\n";
	for (unsigned int i=0; i< cvGraph_.size();i++)
	{
		TTemp_=FileName_+FrequencieTokens_[i];
		(*textOutPut_)<<"Frequency: "<<FrequencieTokens_[i]<<"\n";
		fitResults_.push_back(new FitResult(TTemp_));
		// fit starting parameters and boundaries
		Double_t endingVoltage= 0;
		Double_t startingYValue=0;
		Double_t leftSideSlopeAvargeHight=0;
		Double_t rightSideSlopeAvarageHight=0;
		Double_t highestSlope=0;
		Double_t highestSlopeVoltage=0;
		// relevant points
		Double_t lowestSecondDerivative = 1e40;
		Double_t lowestSecondDerivativeVoltage=0;
		Double_t AvarageSecondDeriviative =0;
		int numberOfSecondDerivativeBins=0;
		// loop over cvGraph bins
		Double_t *xCV = cvGraph_[i].GetX();
		Double_t *yCV = cvGraph_[i].GetY();
		int numberCVGraph = cvGraph_[i].GetN();
		for(int ii=0; ii<numberCVGraph;ii++)
		{
			if(xCV[ii]>startingVoltage)
			{
				if(startingYValue<10)startingYValue=yCV[ii];
				if(endingVoltage<xCV[ii])endingVoltage=xCV[ii];
			}
		}
		// loop over second derivative bins
		Double_t *xSecondDerivative = cvSecondDerivativeGraph_[i].GetX();
		Double_t *ySecondDerivative = cvSecondDerivativeGraph_[i].GetY();
		int numberSecondDerivativeGraph = cvSecondDerivativeGraph_[i].GetN();
		for(int ii=1; ii<(numberSecondDerivativeGraph-1);ii++)
		{
			if(xSecondDerivative[ii]>startingVoltage)
			{
				Double_t sumOfThreeEntries= (ySecondDerivative[ii-1] + ySecondDerivative[ii] + ySecondDerivative[ii+1])/3;
				if(sumOfThreeEntries<lowestSecondDerivative)
				{
					lowestSecondDerivative=sumOfThreeEntries;
					lowestSecondDerivativeVoltage=xSecondDerivative[ii];
					
				}
				AvarageSecondDeriviative+=ySecondDerivative[ii];
				numberOfSecondDerivativeBins++;
			}
		}
		AvarageSecondDeriviative=AvarageSecondDeriviative/numberOfSecondDerivativeBins;
		fitResults_[i]->SetLowestSecondDerivativeVoltage(lowestSecondDerivativeVoltage);
		if(debug_)std::cout<<"Fit:: point with lowest Second Derivative: "<<lowestSecondDerivativeVoltage<<std::endl;
		(*textOutPut_)<<"Fit:: point with lowest Second Derivative: "<<lowestSecondDerivativeVoltage<<"\n";
		// loop over slope bins
		
		Double_t *xSlope = cvSlopeGraph_[i].GetX();
		Double_t *ySlope = cvSlopeGraph_[i].GetY();
		int numberSlopeGraph = cvSlopeGraph_[i].GetN();
		for(int ii=0; ii<numberSlopeGraph;ii++)
		{
			if(xSlope[ii]>startingVoltage)
			{
				if(leftSideSlopeAvargeHight<10)leftSideSlopeAvargeHight=(ySlope[ii]+ySlope[ii+1]+ySlope[ii+2]+ySlope[ii+3])/4;
				if(ii<(numberSlopeGraph-4)) rightSideSlopeAvarageHight = (ySlope[ii]+ySlope[ii+1]+ySlope[ii+2]+ySlope[ii+3])/4;
				if(highestSlope<ySlope[ii])
				{
					highestSlopeVoltage=xSlope[ii];
					highestSlope=ySlope[ii];
				}
			}
		}
		Double_t highestSlopeYet=0;
		Double_t highestVoltageSlopeYet=0;
		int notFound =0;
		for(int ii=numberSlopeGraph-3; ii>1;ii--)
		{
			if(xSlope[ii]>startingVoltage && xSlope[ii]<lowestSecondDerivativeVoltage)
			{
				// summ up 3 bins 
				Double_t tempVoltageAvarge = (ySlope[ii-1]+ ySlope[ii]+ySlope[ii+1])/3;
				if(highestSlopeYet<tempVoltageAvarge)
				{
					notFound=0;
					highestSlopeYet=tempVoltageAvarge;
					highestVoltageSlopeYet=xSlope[ii];
				}
				else notFound++;
				if(notFound>1) break;
				
			}
		}
		fitResults_[i]->SetLeftSideSlopeAvargeHight(leftSideSlopeAvargeHight);
		fitResults_[i]->SetRightSideSlopeAvargeHight(rightSideSlopeAvarageHight);
		fitResults_[i]->PeakSlopeVoltage_ = highestVoltageSlopeYet;
		(*textOutPut_)<<"Fit:: LeftSideAvargeSlope: "<<leftSideSlopeAvargeHight<<"\n";
		(*textOutPut_)<<"Fit:: RightSideAvargeSlope: "<<rightSideSlopeAvarageHight<<"\n";
		(*textOutPut_)<<"Fit:: highestSlope: "<<highestSlope<<" with position in V: "<<highestSlopeVoltage<<"[V]\n";
		
		// find highest slope starting from the lowest second derivative to make sure to be at the peak of the rise and not at a nother peak which is at low voltages
		Double_t avargeSlopeTwoBins=0;
		for(int ii=numberSlopeGraph-4; xSlope[ii]>startingVoltage;ii--)
		{
			
			if(ii-2<0)
			{
				(*textOutPut_)<<"Fit:: Warning highestSlopeVoltage could not be determined please check code. Should never happen;-)\n";
				break;
			}
			if(avargeSlopeTwoBins< (ySlope[ii+2]+ySlope[ii+1]+ySlope[ii])/3)
			{
				avargeSlopeTwoBins=ySlope[ii+1];
				highestSlopeVoltage=xSlope[ii+1];
				highestSlope=ySlope[ii+1];
			}
			//std::cout<<"Voltage: "<<xSlope[ii]<<" should be higher than: "<<startingVoltage<<" found highestVoltage: "<<highestSlopeVoltage<<std::endl;
			
		}
		//std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!! highestSlopeVoltage at the end of loop:"<<highestSlopeVoltage<<std::endl;
		fitResults_[i]->highestSlopeVoltage_=highestSlopeVoltage;
		// try to fit full gauss to the distribution
		// -------------------------------!!!!!!!!!!!!!!!!!!!!!!!!!!1--------------------------
		// Slope fits
		// Gauss fit
		gaussPol0 *GaussPol0 = new gaussPol0();
		// try to fit a gauss to the second deriviatve minimum
		Double_t fGaussRestricSecondDerivativeStart = lowestSecondDerivativeVoltage*0.7;
		Double_t fGaussRestricSecondDerivativeEnd = lowestSecondDerivativeVoltage*1.4;
		if(fGaussRestricSecondDerivativeEnd>endingVoltage)fGaussRestricSecondDerivativeEnd=endingVoltage-1;
		std::cout<<"1";
		TF1 * fGaussRestricSecondDerivative = new TF1("FullGauss",GaussPol0,&gaussPol0::Evaluate,fGaussRestricSecondDerivativeStart,fGaussRestricSecondDerivativeEnd,4);
		(*textOutPut_)<<"Fit:: fGaussRestricSecondDerivative:Range: "<<fGaussRestricSecondDerivativeStart<<" to "<<fGaussRestricSecondDerivativeEnd<<"\n";
		fGaussRestricSecondDerivative->SetParameters((lowestSecondDerivative-AvarageSecondDeriviative),lowestSecondDerivativeVoltage,5,AvarageSecondDeriviative);
		//fGaussRestricSecondDerivative->FixParameter(1,lowestSecondDerivativeVoltage);
		fGaussRestricSecondDerivative->SetLineColor(i+1);
		fGaussRestricSecondDerivative->SetLineStyle(10);
		cvSecondDerivativeGraph_[i].Fit(fGaussRestricSecondDerivative,"Rq+");
		Double_t secondDerivativeGausMean = fGaussRestricSecondDerivative->GetParameter(1);
		Double_t secondDerivativeGausSigma = fGaussRestricSecondDerivative->GetParameter(2);
		(*textOutPut_)<<"Fit:: fGaussRestricSecondDerivative:Mean: "<<secondDerivativeGausMean<<" Sigma "<<secondDerivativeGausSigma<<"\n";
		TF1 * fGaussFull = new TF1("FullGauss",GaussPol0,&gaussPol0::Evaluate,startingVoltage,endingVoltage,4);
		(*textOutPut_)<<"Fit:: FullGaussFit:Range: "<<startingVoltage<<" to "<<endingVoltage<<"\n";
		Double_t variation=fabs(lowestSecondDerivativeVoltage-highestSlopeVoltage)*2;
		fGaussFull->SetParameters(highestSlope-((leftSideSlopeAvargeHight+rightSideSlopeAvarageHight)/2),lowestSecondDerivativeVoltage,100,(leftSideSlopeAvargeHight+rightSideSlopeAvarageHight)/2);
		fGaussFull->SetParLimits(1,lowestSecondDerivativeVoltage-variation,lowestSecondDerivativeVoltage+variation);
		if(highestVoltageSlopeYet<lowestSecondDerivativeVoltage) // true if with the avarage method a highest slope was found which is closer to the minimum of the second derivative
		{
			variation=fabs(lowestSecondDerivativeVoltage-highestVoltageSlopeYet)*2;
			fGaussFull->SetParameters(highestSlope-((leftSideSlopeAvargeHight+rightSideSlopeAvarageHight)/2),highestVoltageSlopeYet,100,(leftSideSlopeAvargeHight+rightSideSlopeAvarageHight)/2);
			fGaussFull->SetParLimits(1,lowestSecondDerivativeVoltage-variation,lowestSecondDerivativeVoltage+variation);
		}
		fGaussFull->SetLineColor(i+1);
		fGaussFull->SetLineStyle(10);
		cvSlopeGraph_[i].Fit(fGaussFull,"NRq+");
		fitResults_[i]->FullGaussMean_=fGaussFull->GetParameter(1);
		// restricted gauss
		Double_t range = fGaussFull->GetParameter(2)*1.0;
		gauss *Gauss = new gauss();
		Double_t rangeLow=fGaussFull->GetParameter(1)-range;
		//Double_t rangeHigh=fGaussFull->GetParameter(1)+range;
		Double_t rangeHigh = endingVoltage;
		Double_t parameter[3] = {fGaussFull->GetParameter(0),fGaussFull->GetParameter(1),fGaussFull->GetParameter(2)};

		for (int ii=0; ii<20;ii++)
		{
			TF1 * fGaussRestricted = new TF1("RestrictedGauss",Gauss,&gauss::Evaluate,rangeLow,rangeHigh,3);
			fGaussRestricted->SetParameters(parameter);
			fGaussRestricted->FixParameter(1,highestSlopeVoltage);
			fGaussRestricted->SetParLimits(2,fGaussFull->GetParameter(2)*0.5,fGaussFull->GetParameter(2)*2);
			fGaussRestricted->SetLineColor(i+1);
			fGaussRestricted->SetLineStyle(2);
			cvSlopeGraph_[i].Fit(fGaussRestricted,"NRq+");
			parameter = {fGaussRestricted->GetParameter(0),fGaussRestricted->GetParameter(1),fGaussRestricted->GetParameter(2)};
			rangeLow=fGaussRestricted->GetParameter(1)-fGaussRestricted->GetParameter(1)*1.1;
		//	rangeHigh=fGaussRestricted->GetParameter(1)+fGaussRestricted->GetParameter(1)*1.1;
			delete fGaussRestricted;
			
		}
		/*for (int ii=0; ii<5;ii++)
		{
			TF1 * fGaussRestricted = new TF1("RestrictedGauss",Gauss,&gauss::Evaluate,rangeLow,rangeHigh,3);
			fGaussRestricted->SetParameters(parameter);
			fGaussRestricted->FixParameter(1,highestSlopeVoltage);
			fGaussRestricted->SetParLimits(2,fGaussFull->GetParameter(2)*0.5,fGaussFull->GetParameter(2)*2);
			fGaussRestricted->SetLineColor(i+3);
			fGaussRestricted->SetLineStyle(2);
			cvSlopeGraph_[i].Fit(fGaussRestricted,"NRq+");
			parameter = {fGaussRestricted->GetParameter(0),fGaussRestricted->GetParameter(1)*0.9,fGaussRestricted->GetParameter(2)};
			rangeLow=fGaussRestricted->GetParameter(1)-fGaussRestricted->GetParameter(1)*0.6;
			//rangeHigh=fGaussRestricted->GetParameter(1)+fGaussRestricted->GetParameter(1)*1.1;
			if(!(ii+1)<5) parameter =  {fGaussRestricted->GetParameter(0),fGaussRestricted->GetParameter(1),fGaussRestricted->GetParameter(2)};
			delete fGaussRestricted;
			
			
		}*/
		//TF1 * fGaussRestrictedFinal = new TF1("RestrictedGauss",Gauss,&gauss::Evaluate,rangeLow,rangeHigh,3); use this for better performance for very high irratidated samples
		TF1 * fGaussRestrictedFinal = new TF1("RestrictedGauss",Gauss,&gauss::Evaluate,secondDerivativeGausMean-secondDerivativeGausSigma*3,rangeHigh,3);
		fGaussRestrictedFinal->SetParameters(parameter);
		//fGaussRestrictedFinal->FixParameter(1,highestSlopeVoltage);
		fGaussRestrictedFinal->SetLineColor(i+1);
		fGaussRestrictedFinal->SetLineStyle(2);
		cvSlopeGraph_[i].Fit(fGaussRestrictedFinal,"Rq+");
		(*textOutPut_)<<"Fit:: RestrictedGaussParameters: [0]= "<<fGaussRestrictedFinal->GetParameter(0)<<"+-"<<fGaussRestrictedFinal->GetParError(0)<<" [1]= "<<fGaussRestrictedFinal->GetParameter(1)<<"+-"<<fGaussRestrictedFinal->GetParError(1)<<" [2]= "<<fGaussRestrictedFinal->GetParameter(2)<<"+-"<<fGaussRestrictedFinal->GetParError(2)<<"\n";
		fitResults_[i]->RestrictedGaussResult_.push_back(fGaussRestrictedFinal->GetParameter(0));
		fitResults_[i]->RestrictedGaussResultError_.push_back(fGaussRestrictedFinal->GetParError(0));
		fitResults_[i]->RestrictedGaussResult_.push_back(fGaussRestrictedFinal->GetParameter(1));
		fitResults_[i]->RestrictedGaussResultError_.push_back(fGaussRestrictedFinal->GetParError(1));
		fitResults_[i]->RestrictedGaussResult_.push_back(fGaussRestrictedFinal->GetParameter(2));
		fitResults_[i]->RestrictedGaussResultError_.push_back(fGaussRestrictedFinal->GetParError(2));
		// decide which method to use
		if(fGaussRestrictedFinal->GetParameter(0)<0 || fGaussRestrictedFinal->GetParameter(1)<0 || fGaussRestrictedFinal->GetParameter(2)<0 || (fGaussRestrictedFinal->GetParameter(1)-fGaussRestrictedFinal->GetParameter(2))<0) useNonHighIradiatedCriteria_=true;
		
		if(useNonHighIradiatedCriteria_)
		{
			fitResults_[i]->SetUsedMethod("NonHighIrradiated");
			// find the highest slope starting from the lowest second derivative point

			//Double_t lowerBoundFitLeft=lowestSecondDerivativeVoltage-(lowestSecondDerivativeVoltage-highestSlopeVoltage)*4;
			//Double_t higherBoundFitLeft=lowestSecondDerivativeVoltage-(lowestSecondDerivativeVoltage-highestSlopeVoltage)*1.05;
			//Double_t lowerBoundFitRight=lowestSecondDerivativeVoltage+(lowestSecondDerivativeVoltage-highestSlopeVoltage)*1.15;
			//Double_t higherBoundFitRight=lowestSecondDerivativeVoltage+(lowestSecondDerivativeVoltage-highestSlopeVoltage)*1.6;
			Double_t lowerBoundFitLeft=lowestSecondDerivativeVoltage-(lowestSecondDerivativeVoltage-highestVoltageSlopeYet)*4;
			Double_t higherBoundFitLeft=lowestSecondDerivativeVoltage-(lowestSecondDerivativeVoltage-highestVoltageSlopeYet)*1.05;
			Double_t lowerBoundFitRight=lowestSecondDerivativeVoltage+(lowestSecondDerivativeVoltage-highestVoltageSlopeYet)*1.15;
			Double_t higherBoundFitRight=lowestSecondDerivativeVoltage+(lowestSecondDerivativeVoltage-highestVoltageSlopeYet)*1.6;
			if( (secondDerivativeGausMean-secondDerivativeGausSigma*1.5)> highestVoltageSlopeYet ) 
			{
				fitResults_[i]->SetUsedMethod("GausSecondDerivative");
				Double_t edgeOfSlope = secondDerivativeGausMean-secondDerivativeGausSigma*1.5;
				lowerBoundFitLeft=lowestSecondDerivativeVoltage-(lowestSecondDerivativeVoltage-edgeOfSlope)*4;
				higherBoundFitLeft=lowestSecondDerivativeVoltage-(lowestSecondDerivativeVoltage-edgeOfSlope)*1.05;
				lowerBoundFitRight=lowestSecondDerivativeVoltage+(lowestSecondDerivativeVoltage-edgeOfSlope)*1.15;
				higherBoundFitRight=lowestSecondDerivativeVoltage+(lowestSecondDerivativeVoltage-edgeOfSlope)*1.6;
				
			}
			else fitResults_[i]->SetUsedMethod("BumpAtEdgeUsed");
			pol1 *Pol1 = new pol1();
			TF1 * Pol1Left = new TF1("Pol1Left",Pol1,&pol1::Evaluate,lowerBoundFitLeft,higherBoundFitLeft,2);
			Pol1Left->SetParameter(0,startingYValue*0.8);
			Pol1Left->SetParameter(1,leftSideSlopeAvargeHight);
			Pol1Left->SetLineColor(i+1);
			Pol1Left->SetLineWidth(4);
			Pol1Left->SetLineStyle(1);
			cvGraph_[i].Fit(Pol1Left,"Rq+");
			TF1 * Pol1Right = new TF1("Pol1Right",Pol1,&pol1::Evaluate,lowerBoundFitRight,endingVoltage,2);
			Pol1Right->SetParameter(0,startingYValue*0.8);
			Pol1Right->SetParameter(1,rightSideSlopeAvarageHight);
			Pol1Right->SetLineColor(i+1);
			Pol1Right->SetLineWidth(4);
			Pol1Right->SetLineStyle(1);
			cvGraph_[i].Fit(Pol1Right,"Rq+");
			TF1 * Pol1LeftFull = new TF1("Pol1Left",Pol1,&pol1::Evaluate,startingVoltage,endingVoltage,2);
			Pol1LeftFull->FixParameter(0,Pol1Left->GetParameter(0) );
			Pol1LeftFull->FixParameter(1,Pol1Left->GetParameter(1) );
			Pol1LeftFull->SetLineColor(i+1);
			Pol1LeftFull->SetLineStyle(1);
			cvGraph_[i].Fit(Pol1LeftFull,"Rq+");
			TF1 * Pol1RightFull = new TF1("Pol1Right",Pol1,&pol1::Evaluate,startingVoltage,endingVoltage,2);
			Pol1RightFull->FixParameter(0,Pol1Right->GetParameter(0) );
			Pol1RightFull->FixParameter(1,Pol1Right->GetParameter(1) );
			Pol1RightFull->SetLineColor(i+1);
			Pol1RightFull->SetLineStyle(1);
			cvGraph_[i].Fit(Pol1RightFull,"Rq+");
			cout<<"tobecreated";
			TGraph *errorGraph = new TGraph("Errors");
			cout<<"created";
			//(TVirtualFitter::GetFitter())->GetConfidenceIntervals(errorGraph);
			cout<<"stored";
			cvGraphError_.push_back(*errorGraph);
			cout<<"done";
			if(Pol1Right->GetParameter(1) > Pol1Left->GetParameter(1)) fitResults_[i]->SetDepletionVoltage(-666);
			else 
			{
			  Double_t a1, a2, b1, b2, a1Error, a2Error, b1Error, b2Error;
			  a1=Pol1Left->GetParameter(0);
			  a1Error=Pol1Left->GetParError(0);
			  a2=Pol1Right->GetParameter(0);
			  a2Error=Pol1Right->GetParError(0);
			  b1=Pol1Left->GetParameter(1);
			  b1Error=Pol1Left->GetParError(1);
			  b2=Pol1Right->GetParameter(1);
			  b2Error=Pol1Right->GetParError(1);
			  fitResults_[i]->SetDepletionVoltage((Pol1Left->GetParameter(0) - Pol1Right->GetParameter(0) ) / ( Pol1Right->GetParameter(1) - Pol1Left->GetParameter(1) ));
			  Double_t errorSqaured = a1Error * 1/(b2-b1) *a1Error * 1/(b2-b1);
			  errorSqaured+=a2Error * 1/(b2-b1) *a2Error * 1/(b2-b1);
			  errorSqaured+=b1Error * (a1-a2)/((b2-b1)*(b2-b1)) *b1Error * (a1-a2)/((b2-b1)*(b2-b1));
			  errorSqaured+=b2Error * (a1-a2)/((b2-b1)*(b2-b1)) *b2Error * (a1-a2)/((b2-b1)*(b2-b1));
			  fitResults_[i]->SetDepletionVoltageError(sqrt(errorSqaured));
			  
			}
		}
		else // will be triggered 
		{
			fitResults_[i]->SetUsedMethod("HighIrradidatedUsed");
			pol1 *Pol1 = new pol1();
			//Double_t voltageRight = endingVoltage - (endingVoltage-fGaussRestrictedFinal->GetParameter(2)*3)*0.4 ;
			Double_t voltageRight = fGaussRestrictedFinal->GetParameter(1) +  fGaussRestrictedFinal->GetParameter(2)*2.5;  // this can lead to segemenation violations!!!!!
			Double_t voltageRightHigh = fGaussRestrictedFinal->GetParameter(1) +  fGaussRestrictedFinal->GetParameter(2)*4.5; // this can lead to segemenation violations!!!!!
			if(voltageRightHigh>endingVoltage)voltageRightHigh=endingVoltage;
			TF1 * Pol1Middle = new TF1("Pol1Middle",Pol1,&pol1::Evaluate,fGaussRestrictedFinal->GetParameter(1)-(fGaussRestrictedFinal->GetParameter(2))*0.9,fGaussRestrictedFinal->GetParameter(1)+(fGaussRestrictedFinal->GetParameter(2))*0.9,2);
			Pol1Middle->SetParameter(0,startingYValue*0.8);
			Pol1Middle->SetParameter(1,leftSideSlopeAvargeHight);
			Pol1Middle->SetLineColor(i+1);
			Pol1Middle->SetLineWidth(4);
			Pol1Middle->SetLineStyle(1);
			cvGraph_[i].Fit(Pol1Middle,"Rq+");
			TF1 * Pol1Right = new TF1("Pol1Right",Pol1,&pol1::Evaluate,voltageRight,voltageRightHigh,2);
			Pol1Right->SetParameter(0,startingYValue*0.8);
			Pol1Right->SetParameter(1,rightSideSlopeAvarageHight);
			Pol1Right->SetLineColor(i+1);
			Pol1Right->SetLineWidth(4);
			Pol1Right->SetLineStyle(1);
			cvGraph_[i].Fit(Pol1Right,"Rq+");
			TF1 * Pol1LeftFull = new TF1("Pol1Middle",Pol1,&pol1::Evaluate,startingVoltage,endingVoltage,2);
			Pol1LeftFull->FixParameter(0,Pol1Middle->GetParameter(0) );
			Pol1LeftFull->FixParameter(1,Pol1Middle->GetParameter(1) );
			Pol1LeftFull->SetLineColor(i+1);
			Pol1LeftFull->SetLineStyle(1);
			cvGraph_[i].Fit(Pol1LeftFull,"Rq+");
			TF1 * Pol1RightFull = new TF1("Pol1Right",Pol1,&pol1::Evaluate,startingVoltage,endingVoltage,2);
			Pol1RightFull->FixParameter(0,Pol1Right->GetParameter(0) );
			Pol1RightFull->FixParameter(1,Pol1Right->GetParameter(1) );
			Pol1RightFull->SetLineColor(i+1);
			Pol1RightFull->SetLineStyle(1);
			cvGraph_[i].Fit(Pol1RightFull,"Rq+");
			cout<<"tobecreated";
			TGraph *errorGraph = new TGraph("Errors");
			cout<<"created";
			//(TVirtualFitter::GetFitter())->GetConfidenceIntervals(errorGraph);
			cout<<"stored";
			cvGraphError_.push_back(*errorGraph);
			cout<<"done";
			fitResults_[i]->SetDepletionVoltage((Pol1Middle->GetParameter(0) - Pol1Right->GetParameter(0) ) / ( Pol1Right->GetParameter(1) - Pol1Middle->GetParameter(1) ));
			 Double_t a1, a2, b1, b2, a1Error, a2Error, b1Error, b2Error;
			  a1=Pol1Middle->GetParameter(0);
			  a1Error=Pol1Middle->GetParError(0);
			  a2=Pol1Right->GetParameter(0);
			  a2Error=Pol1Right->GetParError(0);
			  b1=Pol1Middle->GetParameter(1);
			  b1Error=Pol1Middle->GetParError(1);
			  b2=Pol1Right->GetParameter(1);
			  b2Error=Pol1Right->GetParError(1);
			  Double_t errorSqaured = a1Error * 1/(b2-b1) *a1Error * 1/(b2-b1);
			  errorSqaured+=a2Error * 1/(b2-b1) *a2Error * 1/(b2-b1);
			  errorSqaured+=b1Error * (a1-a2)/((b2-b1)*(b2-b1)) *b1Error * (a1-a2)/((b2-b1)*(b2-b1));
			  errorSqaured+=b2Error * (a1-a2)/((b2-b1)*(b2-b1)) *b2Error * (a1-a2)/((b2-b1)*(b2-b1));
			  fitResults_[i]->SetDepletionVoltageError(sqrt(errorSqaured));
		}
		// linear fit left side to maximum
		// -------------------------------!!!!!!!!!!!!!!!!!!!!!!!!!!1--------------------------
	}
	
	return result;
}

#ifndef _mal2_hpp_
#define _mal2_hpp_
// constants

// varialbes
  string lineBeforeBeginOfTable_,lineAfterTable_;
  bool debug_=false;
  TFile *outF_;
  fstream *FitResults;
  vector<Measurment*> measurments;
// tempvariables
  string line_;
// methods
  void cvAnalyzer2014(Measurment *measurment);
  Double_t Pol0(Double_t *x, Double_t *p);
  #endif
void cvAnalyzer()
{
	outF_ = new TFile("results.root","RECREATE");
	ofstream *textOutPut = new ofstream;
	textOutPut->open("FitResults.txt");
	textOutPut->close();	
	delete textOutPut;
	FitResults = new fstream("FitResults.txt");
	debug_=false;
	gROOT->SetBatch(true);
	gStyle->SetOptStat(0);
	lineBeforeBeginOfTable_="BEGIN";
	lineAfterTable_="END";
	measurments.push_back(new Measurment("261636-11-41_2014-07-18_20") );
	measurments.push_back(new Measurment("8556-04-49_2014-07-16_2") );
	measurments.push_back(new Measurment("10_d") );
	measurments.push_back(new Measurment("10_p") );
	measurments.push_back(new Measurment("13_d") );
	measurments.push_back(new Measurment("13_p") );
	measurments.push_back(new Measurment("14_d") );
	measurments.push_back(new Measurment("14_p") );
	measurments.push_back(new Measurment("17_d") );
	measurments.push_back(new Measurment("17_p") );
	measurments.push_back(new Measurment("18_d") );
	measurments.push_back(new Measurment("18_p") );
	measurments.push_back(new Measurment("3_d") );
	measurments.push_back(new Measurment("3_p") ); 
	
	measurments.push_back(new Measurment("8364-03-30_2006-01-11_1") );
	measurments.push_back(new Measurment("8364-03-30_2008-07-09_2") );
	measurments.push_back(new Measurment("8364-03-30_2008-08-11_3") );
	
	measurments.push_back(new Measurment("8364-03-38_2006-01-11_1") );
	measurments.push_back(new Measurment("8364-03-38_2007-09-04_2") );
	
	measurments.push_back(new Measurment("8364-03-41_2006-01-11_1") );
	measurments.push_back(new Measurment("8364-03-41_2007-09-04_2") );
	
	measurments.push_back(new Measurment("8364-03-42_2006-01-11_1") );
	measurments.push_back(new Measurment("8364-03-42_2007-09-04_2") );
	
	measurments.push_back(new Measurment("8364-03-43_2006-01-11_1") );
	measurments.push_back(new Measurment("8364-03-43_2007-09-04_2") );
	
	measurments.push_back(new Measurment("8364-03-44_2006-01-11_1") );
	measurments.push_back(new Measurment("8364-03-44_2007-09-04_2") );
	Double_t startingVoltageToAnalyze=40;
	*FitResults <<"------------------------------------------------------------------------------------------------------------------\n";
	*FitResults <<"-----------------------------------------PrintOutOfResults--------------------------------------------------------\n";
	*FitResults <<"------------------------------------------------------------------------------------------------------------------\n";
	*FitResults << "FitResults for: "<<measurments.size()<<" files. \n";
	*FitResults << "The fits start only from: "<<startingVoltageToAnalyze<<" [V]. \n";
	for(unsigned int i=0; i < measurments.size(); i++)
	{
		std::cout<<"Starting analyzing file: "<<measurments[i]->GetFileName()<<std::endl;
		*FitResults << "File: "<<measurments[i]->GetFileName();
		measurments[i]->debug_=debug_;
		measurments[i]->lineBeforeBeginOfTable_=lineBeforeBeginOfTable_;
		measurments[i]->lineAfterTable_=lineAfterTable_;
		measurments[i]->setStartingVoltageToAnalyze(startingVoltageToAnalyze);
		cvAnalyzer2014(measurments[i]);
		/* *FitResults <<"\n: SecondDeriLowest: FirstDerivivativeLowestPoint\n";
		for(unsigned int ii=0; ii<measurments[i]->fitResults_.size();ii++)
		{
			
			*FitResults<<"["<<ii<<"]:"<<measurments[i]->fitResults_[ii]->GetLowestSecondDerivativeVoltage()<<", "<<measurments[i]->fitResults_[ii]->highestSlopeVoltage_<<", ";
			*FitResults<<"\n";
		}
		*FitResults <<"GaussResult First Derivative \n";
		for(unsigned int ii=0; ii<measurments[i]->fitResults_.size();ii++)
		{
			*FitResults<<"Frequency["<<ii<<"]: ";
			for(unsigned int iii=0; iii<measurments[i]->fitResults_[ii]->RestrictedGaussResult_.size();iii++)
			{
				*FitResults<<"p["<<iii<<"]:"<<measurments[i]->fitResults_[ii]->RestrictedGaussResult_[iii]<<", ";
			}
			*FitResults<<"\n";
		}
		*/
		if(measurments[i]->useNonHighIradiatedCriteria_)*FitResults<<"\nNon high irratidated input is found.\n";
		else *FitResults<<"\nHigh irradiated input found.\n";
		for(unsigned int ii=0; ii<measurments[i]->FrequencieTokens_.size();ii++)
		{
			*FitResults<<ii<<" Frequency: "<<measurments[i]->frequenciesDouble_[ii]<<"[Hz], DepletionVoltage: "<<measurments[i]->fitResults_[ii]->GetDepletionVoltage()<<"+-"<<measurments[i]->fitResults_[ii]->GetDepletionVoltageError()<<" [V] used Method: "<<measurments[i]->fitResults_[ii]->GetUsedMethod()<<"\n";
			std::cout<<ii<<" Frequency: "<<measurments[i]->frequenciesDouble_[ii]<<"[Hz], DepletionVoltage: "<<measurments[i]->fitResults_[ii]->GetDepletionVoltage()<<" [V] used Method: "<<measurments[i]->fitResults_[ii]->GetUsedMethod()<<"\n";
		}
		*FitResults<<"-------------------------------------------------------------------\n";
	}
	*FitResults <<"------------------------------------------------------------------------------------------------------------------\n";
	for(unsigned int i=0; i < measurments.size(); i++) delete measurments[i];
	if(FitResults->is_open())
	{
		FitResults->seekg(0,FitResults->beg);
		while (getline (*FitResults,line_)) cout<<line_<<endl;
	}
	else cout<<"File not opend output not printed out please see FitResults.txt file"<<std::endl;
	FitResults->close();
}
void cvAnalyzer2014(Measurment *measurment)
{
	outF_->cd();
	outF_->mkdir(measurment->GetFileName());
	TDirectory *tDir = (TDirectory*)outF_->Get(measurment->GetFileName());
	tDir->cd();
	measurment->SetDir(tDir);
	if(!measurment->LoadInputFromTxt())
	{
		std::cout<<"WARNING cv FILE COULD NOT BE OPEND! Name: "<<measurment->GetFileName()<<std::endl;
		*FitResults <<"ERROR File could not be opend!!!\n";
	}
	else std::cout<<"cv File opend"<<std::endl;
	if(!measurment->LoadFrequencies() )
	{
		std::cout<<"WARNING Frequencies not extrated! Name: "<<measurment->GetFileName()<<std::endl;
		*FitResults <<"WARNING Frequencies not extrated! Name: "<<measurment->GetFileName()<<std::endl;
	}
	else std::cout<<"frequencies Loaded!"<<std::endl;
	if(!measurment->FillGraph() )
	{
		std::cout<<"WARNING Graph not filled! Name: "<<measurment->GetFileName()<<std::endl;
		*FitResults <<"WARNING Graph not filled! Name: "<<measurment->GetFileName()<<std::endl;
	}
	else std::cout<<"Graph filled!"<<std::endl;
	if(!measurment->SortByFrequencies() )
	{
		std::cout<<"WARNING SortByFrequencies failed! Name: "<<measurment->GetFileName()<<std::endl;
		*FitResults <<"WARNING SortByFrequencies failed! Name: "<<measurment->GetFileName()<<std::endl;
	}
	else std::cout<<"SortByFrequencies dons!"<<std::endl;
	if(!measurment->FitReloaded() )
	{
		std::cout<<"WARNING FitReloaded failed! Name: "<<measurment->GetFileName()<<std::endl;
		*FitResults <<"WARNING FitReloaded failed! Name: "<<measurment->GetFileName()<<std::endl;
	}
	else std::cout<<"FitReloaded done!"<<std::endl;
	if(!measurment->SaveResults() )
	{
		std::cout<<"WARNING SaveResults failed! Name: "<<measurment->GetFileName()<<std::endl;
		*FitResults <<"WARNING SortByFrequencies failed! Name: "<<measurment->GetFileName()<<std::endl;
	}
	else std::cout<<"SaveResults done!"<<std::endl;

	
}
Double_t Pol0(Double_t *x, Double_t *p)
{
	return (p[0]+x[0]*p[0]);
}