#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <string>
#include <chrono>
#include <vector> 
#include <algorithm>


using namespace std;

int main(int argc, char const *argv[])
{
	ifstream mini_file(argv[1]);
	ofstream compatibleSAM;
  	compatibleSAM.open("temp_mini_ngmlr.sam");
  	vector<string> splitedReadName;

  	if (mini_file.is_open())
  	{
  		string line {};
  		string sequenceName, chrNum, CIGAR, sequnce, perBaseQual, AS, NM, ms, nn, tp, cm, s1, s2, de, SA, MD, zd ,previousName, rl;
		int flag {0}, mapQual {0};
		size_t refPositionStart {0};
		char RNEXT {'*'}, PNEXT {'0'}, TLEN{'0'};
  		int lineNumber {0};

  		while (!mini_file.eof())
    	{
    		lineNumber += 1;
    		getline(mini_file, line);
    		size_t found = line.find("SA:");
			istringstream thisLine(line);

			// if (lineNumber == 1 or lineNumber == 2 or lineNumber == 3){
			if (line.find("@SQ\tSN:") != std::string::npos){
				compatibleSAM << line << endl;
			}
    		if (found != string::npos){
    			size_t found2 = line.find("zd:");
    			if (found2 != string::npos){
    				thisLine >> sequenceName >> flag >> chrNum >> refPositionStart >> mapQual >> CIGAR >> RNEXT >> PNEXT >> TLEN >> sequnce >> perBaseQual
    				>> NM >> ms >> AS >> nn >> tp >> cm >> s1 >> s2 >> de >> zd >> SA >> MD >> rl;
    			}
    			else{
    				thisLine >> sequenceName >> flag >> chrNum >> refPositionStart >> mapQual >> CIGAR >> RNEXT >> PNEXT >> TLEN >> sequnce >> perBaseQual
    				>> NM >> ms >> AS >> nn >> tp >> cm >> s1 >> s2 >> de >> SA >> MD >> rl;
    			}
    			// if (refPositionStart == 7252794)
    			// 	cout << line << endl;
    			
    // 			if (sequenceName == "S2_291116"){
    // 				cout << sequenceName << "\t" << flag << "\t" << chrNum << "\t" << refPositionStart << "\t" << mapQual 
				// 	<< "\t" << CIGAR << "\t" << RNEXT << "\t" << PNEXT << "\t" << TLEN << "\t" << sequnce << "\t" << perBaseQual << "\t"
				// 	<< AS << "\t" << SA << "\t" << MD << endl;
				// }
    			

    			splitedReadName.push_back(sequenceName);
    				
				compatibleSAM << sequenceName << "\t" << flag << "\t" << chrNum << "\t" << refPositionStart << "\t" << mapQual 
				<< "\t" << CIGAR << "\t" << RNEXT << "\t" << PNEXT << "\t" << TLEN << "\t" << sequnce << "\t" << perBaseQual
				<< "\t" << NM << "\t" << ms << "\t" << AS << "\t" << nn << "\t" << tp << "\t" << cm << "\t" << s1 << "\t" << s2 << "\t" << de << "\t" << SA << "\t" << MD << endl;

				// compatibleSAM << line << endl;
    			
    		}
    		else {
    	// 		thisLine >> sequenceName >> flag >> chrNum >> refPositionStart;
    	// 		if (refPositionStart == 7252794)
    	// 			cout << line << endl;
    	// 		if (sequenceName == "S2_291116"){
    	// 			cout << sequenceName << "\t" << flag << "\t" << chrNum << "\t" << refPositionStart << "\t" << mapQual 
					// << "\t" << CIGAR << "\t" << RNEXT << "\t" << PNEXT << "\t" << TLEN << "\t" << sequnce << "\t" << perBaseQual << "\t"
					// << AS << "\t" << SA << "\t" << MD << endl;
    	// 		}
    			size_t found3 = line.find("AS:");
    			size_t found2 = line.find("zd:");
    			if (found3 != string::npos){
	    			if (found2 != string::npos){
	    				thisLine >> sequenceName >> flag >> chrNum >> refPositionStart >> mapQual >> CIGAR >> RNEXT >> PNEXT >> TLEN >> sequnce >> perBaseQual
	    				>> NM >> ms >> AS >> nn >> tp >> cm >> s1 >> s2 >> de >> zd >> MD >> rl;
	    			}
	    			else{
	    				thisLine >> sequenceName >> flag >> chrNum >> refPositionStart >> mapQual >> CIGAR >> RNEXT >> PNEXT >> TLEN >> sequnce >> perBaseQual
	    				>> NM >> ms >> AS >> nn >> tp >> cm >> s1 >> s2 >> de >> MD >> rl;
	    			}
	    			
	    			compatibleSAM << sequenceName << "\t" << flag << "\t" << chrNum << "\t" << refPositionStart << "\t" << mapQual 
					<< "\t" << CIGAR << "\t" << RNEXT << "\t" << PNEXT << "\t" << TLEN << "\t" << sequnce << "\t" << perBaseQual
					<< "\t" << NM << "\t" << ms << "\t" << nn << "\t" << tp << "\t" << cm << "\t" << s1 << "\t" << s2 << "\t" << de << "\t" << AS << "\t" << MD << endl;
	    		}
	    // 		else{
	    // 			if (line.at(0) != '@'){
		   //  			thisLine >> sequenceName >> flag >> chrNum >> refPositionStart >> mapQual >> CIGAR >> RNEXT >> PNEXT >> TLEN >> sequnce >> perBaseQual >> rl;
		    			
		   //  			cout << sequenceName << "\t" << flag << "\t" << chrNum << "\t" << refPositionStart << "\t" << mapQual 
					// 	<< "\t" << CIGAR << "\t" << RNEXT << "\t" << PNEXT << "\t" << TLEN << "\t" << sequnce << "\t" << perBaseQual << "\t"
					// 	<< rl << endl;

		   //  			compatibleSAM << sequenceName << "\t" << flag << "\t" << chrNum << "\t" << refPositionStart << "\t" << mapQual 
					// 	<< "\t" << CIGAR << "\t" << RNEXT << "\t" << PNEXT << "\t" << TLEN << "\t" << sequnce << "\t" << perBaseQual << "\t"
					// 	<< rl << endl;
					// }
	    // 		}
    		}
    	}
  	}

  	mini_file.close();

  	ifstream ngmlr_file("temp_ngmlr_just_split.sam");
  	if (ngmlr_file.is_open())
  	{
  	  	size_t counter = 0;
  		cout << 1 << endl;
  		string line {};
  		string sequenceName, chrNum, CIGAR, sequnce, perBaseQual, AS, NM, XI, XS, XE, XR, MD, SV, SA, previousName;
		  int flag {0}, mapQual {0};
		  size_t refPositionStart {0};
		  char RNEXT {'*'}, PNEXT {'0'}, TLEN{'0'};
  		int lineNumber {0};

  		while (!ngmlr_file.eof())
    	{
    		lineNumber += 1;
    		getline(ngmlr_file, line);
    		size_t found = line.find("SA:");
			   istringstream thisLine(line);
         // if (lineNumber == 1 or lineNumber == 2 or lineNumber == 3 or lineNumber == 4){
         //  compatibleSAM << line << endl;
         //  }

			 if (found != string::npos){
				  thisLine >> sequenceName >> flag >> chrNum >> refPositionStart >> mapQual >> CIGAR >> RNEXT >> PNEXT >> TLEN >> sequnce >> perBaseQual
    			 >> AS >> NM >> XI >> XS >> XE >> XR >> MD >> SV >> SA;
    			 // if (find(splitedReadName.begin(), splitedReadName.end(), sequenceName) == splitedReadName.end()){
    				  compatibleSAM << sequenceName << "\t" << flag << "\t" << chrNum << "\t" << refPositionStart << "\t" << mapQual 
    				  << "\t" << CIGAR << "\t" << RNEXT << "\t" << PNEXT << "\t" << TLEN << "\t" << sequnce << "\t" << perBaseQual
    				  << "\t" << NM << "\t" << AS << "\t" << XI << "\t" << XS << "\t" << XE << "\t" << XR << "\t" << SV << "\t" << SA << "\t" << MD << endl;

    				  if (((flag >> 11) & 1) != 1 and ((flag >> 8) & 1) != 1){
    				  	counter += 1;
    					 previousName = sequenceName;
    				  }
    			// }
			 }
		  }
		  // cout << "number of reads added by ngmlr:  " << counter << endl;
  	}
  	compatibleSAM.close();
  	ngmlr_file.close();

}
