#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <vector>

#define MAX_INDEL 5

using namespace std;

unsigned long int* extract_indel(const string &text);
unsigned long int extract_mismatch(const string &MD);
unsigned long int* extract_CIGAR (const string &cigar, const unsigned long int &seqSize);
vector<size_t> * max_3_bin_differences(const string &MD, const string &CIGAR, const size_t bin_size);
size_t biggest_indel (const string &CIGAR, char type);					//computes biggest insertions or deletions dependinp on the type
size_t biggest_indel_tolerance (const string &CIGAR, char type);		//computes with one tolerance between the deletions or insertions

int totalReads = 0;
int trueSplit = 0;

int main(int argc, char const *argv[])
{
	ifstream minimapFile ("temp_mini.sam");
	ofstream allReadCSV;
	
	allReadCSV.open("temp_test.csv");

	if (minimapFile.is_open())
  	{	
  		unsigned long int mis {0};
  		unsigned long int* start_end_read {nullptr};
  		unsigned long int* INDEL {nullptr};
  		vector<size_t>* bin_data {nullptr};
  		string line {};
		string sequenceName, chrNum, CIGAR, sequnce, perBaseQual, NM, ms, AS, nn, tp, cm, s1, s2, de, zd, SA, MD;
		int flag {0}, mapQual {0}, refPositionStart{0};
		char RNEXT {'*'}, PNEXT {'0'}, TLEN{'0'};
  		int lineNumber{0};
  		string currentInformative, currentNonInformative;
  // 		fileCSV << "sequnce" << "," << "s1" << "," << "s2" << "," << "cm" << "," << "len" << "," << "val" << endl;
	 	// allReadCSV << "name" << "," << "read_start" << "," << "read_end" << "," <<"s1" << "," << "s2" << "," << "cm" << "," << "len" << "," << "ins" << "," << "del" << "," 
	 	// << "mis" << "," << "as" << "," << "de" << "," << "label" << endl;

    	while (!minimapFile.eof())
    	{
    		
    		getline(minimapFile, line);
			istringstream thisLine(line);
			size_t found = line.find("AS:");
			// cout << found<<"\n";


			if (found != string::npos){
				size_t found2 = line.find("SA:"); 
				size_t found3 = line.find("zd:");

				if (found2 != string::npos and found3 != string::npos){
					thisLine >> sequenceName >> flag >> chrNum >> refPositionStart >> mapQual >> CIGAR >> RNEXT >> PNEXT >> TLEN >> sequnce >> perBaseQual 
					>> NM >> ms >> AS >> nn >> tp >> cm >> s1 >> s2 >> de >> zd >> SA >> MD;
				}
				else if (found2 != string::npos and found3 == string::npos){
					thisLine >> sequenceName >> flag >> chrNum >> refPositionStart >> mapQual >> CIGAR >> RNEXT >> PNEXT >> TLEN >> sequnce >> perBaseQual 
					>> NM >> ms >> AS >> nn >> tp >> cm >> s1 >> s2 >> de >> SA >> MD;
				}
				else if (found2 == string::npos and found3 != string::npos){
					thisLine >> sequenceName >> flag >> chrNum >> refPositionStart >> mapQual >> CIGAR >> RNEXT >> PNEXT >> TLEN >> sequnce >> perBaseQual 
					>> NM >> ms >> AS >> nn >> tp >> cm >> s1 >> s2 >> de >> zd >> MD;
				}
				else{
					thisLine >> sequenceName >> flag >> chrNum >> refPositionStart >> mapQual >> CIGAR >> RNEXT >> PNEXT >> TLEN >> sequnce >> perBaseQual 
					>> NM >> ms >> AS >> nn >> tp >> cm >> s1 >> s2 >> de >> MD;
				}

				if (((flag >> 8) & 1) == 0 and found2 == string::npos){
					//informativeReads.push_back(sequenceName);
					cout << line<<"\n";
					lineNumber += 1;
					INDEL = extract_indel(CIGAR);
					bin_data = max_3_bin_differences(MD.substr(5), CIGAR, 100);
					start_end_read = extract_CIGAR(CIGAR, sequnce.length());
					mis = extract_mismatch(MD.substr(5));
					// currentInformative = sequenceName;
					size_t s1_int, s2_int, cm_int;
					s1_int = stoi (s1.substr(5)), s2_int = stoi (s2.substr(5)), cm_int = stoi (cm.substr(5));

					// s1_int = extractIntFromTag(s1), s2_int = extractIntFromTag(s2), cm_int = extractIntFromTag(cm);
					// fileCSV << sequenceName << "," << s1_int << "," << s2_int << "," << cm_int << "," << sequnce.length() << "," << 1 << endl;
					// if (s2_int > 0){
						allReadCSV << sequenceName << "," << start_end_read[0] << "," << start_end_read[1] << "," << s1_int << "," << s2_int << "," << cm_int << "," << sequnce.length() 
						<< "," << biggest_indel(CIGAR, 'D') << "," << biggest_indel(CIGAR, 'I') << "," << biggest_indel_tolerance(CIGAR, 'D') << "," << biggest_indel_tolerance(CIGAR, 'I') 
						<< "," << INDEL[0] << "," << INDEL[1] << "," << mis << "," << AS.substr(5) << "," << de.substr(5) << "," << bin_data->at(0) << ","  << bin_data->at(1) 
						<< "," << bin_data->at(2) << "," << NM.substr(5) << "," << ms.substr(5) << "," << endl;
					// }
					delete [] INDEL;
					delete [] start_end_read;
					delete bin_data;

				}
			}
		}
	}
	minimapFile.close();
}

vector<size_t> * max_3_bin_differences(const string &MD, const string &CIGAR, const size_t bin_size){
	// cout << MD << endl;
	vector<size_t> * store_bin_data = new vector<size_t>;
	size_t bin_deletion {0}, bin_insertion{0}, bin_mismatch {0}, i {0} , j {0}, size {0}, previous_match {0}, previous_deletion {0}, previous_mismatch{0};
	size_t extra_bin {0};
	while (j < MD.length()){
		// cout << size << endl;
		if (MD[j] == '^'){
			if (i - j >= 1){												//if it is of the form 11^T12^C
				size += stoi(string(&MD[i], &MD[j]));
				size += 1;
				if (size > bin_size)
					previous_match = size - bin_size - 1, previous_deletion = 1;
				else
					bin_deletion += 1;	
			}
			else{															//if of the form ^T^C
				bin_deletion += 1;
				size += 1;
			}
			j++;
			i = j + 1;
		}
		else if (MD[j] == 'T' or MD[j] == 'C' or MD[j] == 'A' or MD[j] == 'G'){
			if (i - j >= 1){												//if it is of the form 11T12C
				size += stoi(string(&MD[i], &MD[j]));
				size += 1;
				if (size > bin_size)
					previous_match = size - bin_size - 1, previous_mismatch = 1;
				else
					bin_mismatch += 1;	
			}
			else{															//if of the form TC
				bin_mismatch += 1;
				size += 1;
			}
			i = j + 1;
		}
		j++;

		if (size >= bin_size){
			store_bin_data->push_back(bin_mismatch + bin_deletion);
			size_t remain_match = previous_match;

			for (size_t k = 0; k < previous_match/bin_size; ++k)
			{
				store_bin_data->push_back(0);
				remain_match = previous_match - bin_size;
			}
			size = remain_match + previous_deletion + previous_mismatch;
			bin_deletion = previous_deletion, bin_mismatch = previous_mismatch;
			previous_mismatch = 0, previous_deletion = 0;
		}
		if (j == MD.length()){
			if (MD[j] != 'T' and MD[j] != 'C' and MD[j] != 'A' and MD[j] != 'G')
				extra_bin = stoi(string(&MD[i], &MD[j]));
			store_bin_data->push_back(bin_mismatch + bin_deletion);
		}
		
	}
	// cout << store_bin_data->size() << endl;
	for (int i = 0; i < extra_bin/100 + 1; ++i)
			store_bin_data->push_back(0);	
	
	j = 0, i = 0, size = 0;
	size_t bin_number {0};
	bool flag = false;

	while (j < CIGAR.length()){
		if (CIGAR[j] == 'I'){
			bin_insertion = stoi(string(&CIGAR[i], &CIGAR[j]));
			i = j + 1;
		}
		else if (CIGAR[j] == 'D'){
			size += stoi(string(&CIGAR[i], &CIGAR[j]));
			i = j + 1;
		}
		else if (CIGAR[j] == 'M'){
			size += stoi(string(&CIGAR[i], &CIGAR[j]));
			i = j + 1;
		}
		else if (CIGAR[j] == 'S' and flag == false){
			i = j+1;
			flag = true;
		}

		j++;

		if (size >= bin_size){
			// cout << "here" << endl;
			store_bin_data->at(bin_number) += bin_insertion;
			size = size - bin_size;
			bin_number += 1, bin_insertion = 0;
		}
		if (j == CIGAR.length() or (CIGAR[j] == 'S' and flag == true)){
			store_bin_data->at(bin_number) += bin_insertion;
			break;
		}

	}

	//find the three biggest

	vector<size_t> * max_bin_data = new vector<size_t>;

	size_t max_el_index = max_element(store_bin_data->begin(),store_bin_data->end()) - store_bin_data->begin();
	max_bin_data->push_back(*max_element(store_bin_data->begin(), store_bin_data->end()));
	store_bin_data->at(max_el_index) = 0;

	size_t second_max_el_index = max_element(store_bin_data->begin(),store_bin_data->end()) - store_bin_data->begin();
	max_bin_data->push_back(*max_element(store_bin_data->begin(), store_bin_data->end()));
	store_bin_data->at(second_max_el_index) = 0;

	size_t third_max_el_index = max_element(store_bin_data->begin(),store_bin_data->end()) - store_bin_data->begin();
	max_bin_data->push_back(*max_element(store_bin_data->begin(), store_bin_data->end()));
	store_bin_data->at(third_max_el_index) = 0;


	delete store_bin_data;
	
	return max_bin_data;
}

unsigned long int* extract_CIGAR (const string &cigar, const unsigned long int &seqSize){
	unsigned long int readStart {0}, readEnd {seqSize}, pos {0};
	unsigned long int* extractedData {nullptr};
	extractedData = new unsigned long int [3]; 	//[readStart, readEnd]

	pos	= cigar.find('S');
	if (pos != string::npos and pos < cigar.find('M') and pos < cigar.find('D') and pos < cigar.find('I')){
		readStart = stoi(cigar.substr(0,pos));
		pos += 1;
	}
	else
		pos = 0;

	extractedData[0] = readStart;

	if (cigar[cigar.length() - 1] == 'S'){
		for (int i = cigar.length() - 1; i > 0; i--)
		{
			if (cigar[i] == 'M' or cigar[i] == 'D' or cigar[i] == 'I')
			{
				readEnd = seqSize - stoi(string(&cigar[i+1],&cigar[cigar.length() - 1]));
				extractedData[2] = stoi(string(&cigar[i+1],&cigar[cigar.length() - 1]));
				break;
			}
		}
	}
	extractedData[1] = readEnd;
	return extractedData;	
}

unsigned long int* extract_indel(const string &text){
	size_t i{0}, j{0}, del {0}, ins{0} ;
	unsigned long int* INDEL {nullptr};
	INDEL = new unsigned long int [2];

	while (j < text.length()){
		if (text[j] == 'M' or text[j] == 'D'){
			if (text[j] == 'D')
				del += stoi(string(&text[i], &text[j]));
			i = j+1;
		}
		else if (text[j] == 'I'){
			ins += stoi(string(&text[i], &text[j]));
			i = j+1;
		}
		else if (text[j] == 'S')
			i = j+1;
		
		j++;	
	}
	INDEL[0] = ins, INDEL[1] = del;
	return INDEL;
}

unsigned long int extract_mismatch(const string &MD){
	size_t i {0}, mismatch{0};

	while (i < MD.length()){
		if ((MD[i] == 'T' or MD[i] == 'A' or MD[i] == 'C' or MD[i] == 'G') and MD[i-1] != '^')
			mismatch += 1;
		i+=1;
	}
	return mismatch;
}

size_t biggest_indel(const string &CIGAR, char type){
	size_t i{0}, j{0}, max_del{0};
	bool flag = false;

	while (j < CIGAR.length()){
		if (CIGAR[j] == type){
			if (stoi(string(&CIGAR[i], &CIGAR[j])) > max_del)
				max_del = stoi(string(&CIGAR[i], &CIGAR[j]));
			i = j+1;
		}
		else if (CIGAR[j] == 'I' or CIGAR[j] == 'D')
			i = j + 1;
		else if (CIGAR[j] == 'M')
			i = j+1;
		else if (CIGAR[j] == 'S' and flag == false){
			i = j+1;
			flag = true;
		}
		else if (CIGAR[j] == 'S' and flag == true)
			break;
		
		
		j++;

	}
	return max_del;
}

size_t biggest_indel_tolerance (const string &CIGAR, char type){
	size_t i{0}, j{0}, max_indel{0}, indel{0};
	string temp, max_temp;
	bool flag = false;
	bool go_on = false;

	while (j < CIGAR.length()){
		if (go_on == false){
			if (CIGAR[j] == type){
				indel += stoi(string(&CIGAR[i], &CIGAR[j]));
				temp += string(&CIGAR[i], &CIGAR[j+1]);
				i = j + 1;
				go_on = true;
			}
			else if (CIGAR[j] == 'I' or CIGAR[j] == 'D'){
				indel = 0;
				temp = "";
				i = j + 1;
			}
			else if (CIGAR[j] == 'M'){
				indel = 0;
				temp = "";
				i = j+1;
			}
			else if (CIGAR[j] == 'S' and flag == false){
				i = j+1;
				flag = true;
			}
			else if (CIGAR[j] == 'S' and flag == true)
				break;
		
			j++;
		}
		else if (go_on == true){
			if (CIGAR[j] != type and (CIGAR[j] == 'D' or CIGAR[j] == 'I' or CIGAR[j] == 'M')){
				if (stoi(string(&CIGAR[i], &CIGAR[j])) != 1){

					if (indel > max_indel){
						max_indel = indel;
						max_temp = temp;
					}
					temp = "";
					indel = 0;
				}

				else{
					temp += string(&CIGAR[i], &CIGAR[j+1]);
					indel += 1;
				}

				go_on = false;
				i = j + 1;
			}

			else if (CIGAR[j] == 'S')
				break;

			
			j++;

		}
	}
	return max_indel;
}