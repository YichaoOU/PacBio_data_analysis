#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <algorithm>

using namespace std;

int main(int argc, char const *argv[])
{
	vector<string> splitedReadName;
	// ifstream myfile ("unmapped_detected_NN_25_6.txt");
  //ifstream myfile ("../indel_count/Classification/detected_NN.txt");
  ifstream myfile ("temp_detected_reads.txt");
	if (myfile.is_open())
  	{
  		string line, temp;
  		while (!myfile.eof())
    	{
    		getline(myfile, line);
    		istringstream thisLine(line);
    		thisLine >> temp;
    		splitedReadName.push_back(temp);
    	}
  	}

  	myfile.close();

  	ifstream FastQFile (argv[1]);
  	ofstream SVreads;
  	// SVreads.open("SVreads_mini_22x_v2_25_6.fastq");
    SVreads.open("temp_reads_to_ngmlr.fastq");

  	if (FastQFile.is_open())
  	{
  		string line, name, sequence, name2, qual;

  		while (!FastQFile.eof()){
  			getline(FastQFile, line);
  			istringstream thisLine(line);
  			// size_t found = line.find("@");
  			if (!line.empty()){
	  			if (line.at(0) == '@'){
	  				thisLine >> name;

	  				if (find(splitedReadName.begin(), splitedReadName.end(), name.substr(1)) != splitedReadName.end()){
	  					SVreads << name << endl;

	  					getline(FastQFile, line);
	  					istringstream thisLine1(line);
	  					thisLine1 >> sequence;
	  					SVreads << sequence << endl;

	  					getline(FastQFile, line);
	  					istringstream thisLine2(line);
	  					thisLine2 >> name2;
	  					SVreads << name2 << endl;

	  					getline(FastQFile, line);
	  					istringstream thisLine3(line);
	  					thisLine3 >> qual;
	  					SVreads << qual << endl;
	  				}
	  			}
  			}

  		}
  	}

  	FastQFile.close();
  	SVreads.close();
}