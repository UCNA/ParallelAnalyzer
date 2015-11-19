//g++ -o WCR WireChamberResponse.cpp `root-config --cflags --glibs`

#include "WireChamberResponse.h"
//#include <sys/stat.h>


using namespace std;

//Set File Name       
void WireChamberResponse::SetPhysTree(int rnum){
	stringstream fname;
 	fname<<this->ORD<<"hists/spec_"<<rnum<<".root";
	if(this->FilExists(fname.str())){
		string filename=fname.str();
		this->physfile = TFile::Open(filename.c_str());
		phys = (TTree*)physfile->Get("phys");
  		entnum=phys->GetEntries();
  		phys->SetBranchAddress("Cathodes_Wx",&cathwx);
  		phys->SetBranchAddress("ScintE",&scintE);
  		phys->SetBranchAddress("ScintW",&scintW);
  		phys->SetBranchAddress("Cathodes_Ey",&cathey);
  		phys->SetBranchAddress("Cathodes_Ex",&cathex);
  		phys->SetBranchAddress("Cathodes_Wy",&cathwy);
		return;
	}
	else{cout<<"\nspec_"<<rnum<<".root does not exist in:\n"<<ORD<<"hists/ \n";return;}
}

//Find Max index
Int_t WireChamberResponse::MaxInd(Float_t cath[])
{   
	Float_t tempcath[16];
	for(int i =0; i<16 ; i++)  //(16-sizeof(cath)*CHAR_BIT/32) this was dumb..
	{
	tempcath[i]=cath[i];
	}
	
	this->max=-1.40282e+30; //barely floats
	for (int i=0; i<16; ++i)
	{
                    
		if  (tempcath[i]>this->max)
		{
                
		this->max=tempcath[i];
		maxind=i;
		}
	}
	
	return maxind;
}


void WireChamberResponse::TriMax(Float_t cath[]){  //this function isn't used. 
	this->SetTempCaths(cath);	
	for (int i = 0; i<3;i++){	
	triind[i]=MaxInd(tempcath); //This function exports the index and stores the maximum internally
	trimax[i]=this->max;	//this is where the maximum is located in the class. 
	tempcath[triind[i]]=0;
	}
	return;
}

void WireChamberResponse::QuadMax(Float_t cath[]){
	this->SetTempCaths(cath);
	for (int i = 0; i<4;i++){	
	quadind[i]=this->MaxInd(tempcath); //This function exports the index and stores the maximum internally
	quadmax[i]=this->max;	//this is where the maximum is located in the class. 
	tempcath[quadind[i]]=0;
	}
	return;
	
}





 void WireChamberResponse::SetCaths(int cathdex){
  

	phys->GetEntry(cathdex);

	return;
}

void WireChamberResponse::SetTempCaths(Float_t cath[]){
   
	for(int i = 0; i<16; i++){
	tempcath[i]=cath[i];	
	}
	return;
}




//This is were the rubber meets the road. As in this is where the wire chamber response class is chosen.  
int WireChamberResponse::ResponseType(Float_t cath[]){
 QuadMax(cath);
 if (quadmax[0]>threshold){	
    //Platue Response. 	
    int platnum=1;
    for (int i =0; i<3; i++){	
      if ((platfrac*quadmax[0])<quadmax[i+1]) {platnum++;}
    }
    //we have a platue, is it continious? 
    if (platnum>1) {
      int platcheck=1;		
      Float_t indtemp[16];
      for (int i = 0; i<platnum; i++) {indtemp[i]=(Float_t)this->quadind[i];} 
      for (int i = platnum; i<16; i++) {indtemp[i]=(Float_t)-42.00;}
      QuadMax(indtemp);
      for(int i = 0; i<platnum;i++){
	if(((int)(quadmax[i])-(int)(quadmax[i+1]))==1) {platcheck+=1;}				     
		}
	  
      if(platcheck==platnum) {this->wcpos=(int)quadmax[0]; this->QuadMax(cath); return (platnum+2);}
      else {this->wcpos=0; this->QuadMax(cath); return 8;} ///Split Platue !! undefined response (too complicated)
    }
       
    //Triangle Response Type  (platue number is less than 2)
    else {
      //CLASS 0  single point. 
      if(quadmax[1]<threshold2){this->wcpos=quadind[0]; return 0; }	
      //We have a triangle! 			
      else {
	if ( (int)pow(quadind[0]-quadind[1],2.)==1) {
	  if (      ( quadmax[1]>trifrac*quadmax[2] || (int)pow(quadind[0]-quadind[2],2.)>1  ) && (quadind[0]-quadind[1])==1   ){this->wcpos=quadind[0];return 2; } //right leaning  or maybe its left..
	  else if(  ( quadmax[1]>trifrac*quadmax[2] || (int)pow(quadind[0]-quadind[2],2.)>1  ) && (quadind[0]-quadind[1])==-1   ){this->wcpos=quadind[0]; return 3; } //left leaning,  maybe right.. depending on which direction you are looking at the wire chamber...
	  else {this->wcpos=quadind[0]; return 1;} //centered triangle
	
	  
	}
	else{this->wcpos=quadind[0]; return 0;} //series of spread out points above threshold, return the largest. 			
      }	
    }
  }
  else {this->wcpos=0;return 7;}	 	
}


 bool WireChamberResponse::FilExists(const std::string& name) { 
    ifstream f(name.c_str());
    if (f.good()) {
        f.close();
        return true;
    } else {
        f.close();
        return false;
    }   
}





//testing grounds
/*int main()
{
	
  	WireChamberResponse *WCR =new WireChamberResponse();

	//Set response thresholds, TODO: make these fuzzy to get chi^2. note: this is nearly impossible.
	//optional definitions; default values are shown here. 	
	WCR->threshold=120;  //adc threshold (is it really an event?)	
	WCR->threshold2=82;  //triangle daughter threshold. 		
	WCR->platfrac=0.90;  //what is the platue threshold.?? 
	WCR->trifrac=1.5;	  //what is the leaning threshold.?? 		 
	
	//set run number
	WCR->SetPhysTree(22000); //run number that you want to get the responses for. 

  	
	//lets apply the response to a few events and see how it did!
	for(int i =10; i<20; i++){	
	
	//get the damn entry
	WCR->phys->GetEntry(i);
		
	cout<<WCR->ResponseType(WCR->cathwx)<<"\n";  //this is how to use the function to get response. 
	cout<<WCR->quadmax[0]<<" "<<WCR->quadmax[1]<<" "<<WCR->quadmax[2]<<" "<<WCR->quadmax[3]<<"\n";	
	cout<<WCR->quadind[0]<<" "<<WCR->quadind[1]<<" "<<WCR->quadind[2]<<" "<<WCR->quadind[3]<<"\n";	

	}
	

	
  return 0;
  }*/



/*  //Make Histograms! if 2D use TH2F.. etc.
TH1F * gainhist = new TH1F("gain","gain hist",30,-1,4);//,30,0,10);
 
for (int i = 0; i< linum; ++i){
 //pmtTree->GetEntry(i);
    if(xpos[i]*xpos[i]+ypos[i]*ypos[i] < 25 && gain0[i]<2 && gain0[i]>0){ //&& LG02error[i]>0.0001){
      
  gainhist->Fill(gain0[i]);
  
      }
     }

gainhist->Draw("colz");
*/



