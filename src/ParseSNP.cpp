/*
 * ParseSNP.cpp
 *
 *  Created on: Sep 26, 2012
 *      Author: fritz
 */

#include "ParseSNP.h"

void ParseSNP::parse_cmd_line(int argc, char *argv[]) {

	CmdLine cmdss("Comp Motiv", ' ', version.c_str());
	ValueArg<string> read_file("q", "mapped_reads",
			"Mapped read file (sam/bam)", true, "#", "string");
	cmdss.add(read_file);
	ValueArg<string> ref_file("r", "ref_file", "Reference file (fasta)", true,
			"#", "string");
	cmdss.add(ref_file);
	ValueArg<string> snp_file("s", "snp_file", "VCF file", true, "#", "string");
	cmdss.add(snp_file);

	ValueArg<int> leng("l", "read_lengt", "maximum read length", true, -1, "int");
	cmdss.add(leng);

	ValueArg<string> plots("p", "plot_data",
			"Name of the file for the data required for the R script",
			false, "#", "string");
	cmdss.add(plots);

	ValueArg<string> out("o", "outputfile", "File to print the distripution to",
			true, "#", "string");
	cmdss.add(out);

	try {
		cmdss.parse(argc, argv); //parse arguments
		read_filename = read_file.getValue();
		output = out.getValue();
		if (output[0] == '#') {
			output = read_filename;
			output += ".mot";
		}
		reffile = ref_file.getValue();
		snpfile = snp_file.getValue();
		read_length = leng.getValue();
		range = 2 * read_length;
		plot_file = plots.getValue();
		init();
		parseVCF();

	} catch (ArgException &e) {
		StdOutput out;
		out.failure(cmdss, e);
	}
}

string convert(int number) {
	if (number == 0)
		return "0";
	string temp = "";
	string returnvalue = "";
	while (number > 0) {
		temp += number % 10 + 48;
		number /= 10;
	}
	for (int i = 0; i < (int) temp.length(); i++)
		returnvalue += temp[temp.length() - i - 1];
	return returnvalue;
}

string gen_key(int chr, int pos) {
	string id = convert(chr);
	id += ":";
	id += convert(pos);
	return id;
}

void ParseSNP::init() {
	ifstream myfile;
	myfile.open(reffile.c_str(), ifstream::in);
	if (!myfile.good()) {
		cout << "Fasta Parser: could not open file: " << reffile.c_str()
				<< endl;
		exit(0);
	}

	myfile.getline(buffer, buffer_size);
	int id = 0;
	while (!myfile.eof()) {
		if (buffer[0] == '>') {
			string chr;
			chr.clear();
			for (size_t i = 1;
					i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'
							&& buffer[i] != ' '; i++) {
				chr += buffer[i];
			}
			chrs[chr.c_str()] = id;
			id++;
		}
		myfile.getline(buffer, buffer_size);
	}
	cout << "Chrs detected: " << chrs.size() << endl;
	myfile.close();
}

void ParseSNP::parseVCF() {
	ifstream myfile;
	myfile.open(snpfile.c_str(), ifstream::in);
	if (!myfile.good()) {
		cout << "Fasta Parser: could not open file: " << snpfile.c_str()
				<< endl;
		exit(0);
	}

	myfile.getline(buffer, buffer_size);

	int old_pos = 0;
	int old_id = -1;
	vector<int> tmp;

	snps.resize(chrs.size(), tmp);

	while (!myfile.eof()) {
		if (buffer[0] != '#') {
			int count = 0;
			string chr;
			for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
				if (count == 0 && buffer[i] != '\t') {
					chr += buffer[i];
				}
				if (count == 1 && buffer[i - 1] == '\t') {

					int pos = atoi(&buffer[i]) - 1;
					int id = -1;
					if (chrs.find(chr.c_str()) != chrs.end()) { //found
						id = chrs[chr.c_str()];
					} else {
//						cerr << "Chr not found: |" << chr << "|" << endl;
						exit(0);
					}

					if (id != old_id) {
						old_pos = -100000;
					}

					if (pos > old_pos + range) {
						old_id = id;
						old_pos = pos;

						snps[id].push_back(pos);
						string key = gen_key(id + 1, pos + 1);
						if (covs.find(key.c_str()) == covs.end()) {
							covs[key.c_str()] = new coverage;
							covs[key.c_str()]->cov_100 = new int[range * 2];
							covs[key.c_str()]->cov_90 = new int[range * 2];

							memset(covs[key.c_str()]->cov_100, 0, range * 2 * sizeof(int));
							memset(covs[key.c_str()]->cov_90, 0, range * 2 * sizeof(int));

							if (pos - range >= 0) {
								covs[key.c_str()]->start_pos = pos - range;
							} else {
								covs[key.c_str()]->start_pos = 0;
							}
						}
					} else {
//						cout << "Not: chr: " << chr << " " << pos << endl;
					}
				}
				if (buffer[i] == '\t') {
					count++;
				}
			}
		}
		myfile.getline(buffer, buffer_size);
	}
	myfile.close();
}

void ParseSNP::compute() {
	Parser * mapped_file = 0;

	if (read_filename.find(".bam") != string::npos) {
//		cout << "BAM File" << endl;
		mapped_file = new BamParser(read_filename);
	}

	if (read_filename.find(".sam") != string::npos) {
//		cout << "SAM File" << endl;
		mapped_file = new SamParser(read_filename);
	}

	if (mapped_file == 0) {
		cerr << "File Format not recognized" << endl;
		exit(0);
	}

	int aln_count = 0;
	FastaParser * fasta = new FastaParser(reffile);

	Alignment * tmp_aln = mapped_file->parseRead(0);
	vector<vector<Alignment *> > alignments;
	cout << "start parsing" << endl;
	while (!tmp_aln->getSequence().first.empty()) {
		while ((size_t) tmp_aln->getRefID() >= alignments.size()) {
			vector<Alignment*> tmp;
			alignments.push_back(tmp);
		}
		alignments[tmp_aln->getRefID()].push_back(tmp_aln);
		aln_count++;
		if (aln_count >= MAX_NUM_ALIGNMENTS) {
//			int * test= new int [1];
			//process alignments
			cout << "process alignments: " << aln_count << endl;
			process_alignemnts(alignments, fasta);
			alignments.clear();
			aln_count = 0;
			cout << "parsing" << endl;
		}
		tmp_aln = mapped_file->parseRead(0);
	}

	cout << "end parsing" << endl;
	if (aln_count != 0) {
		//process alignments
		cout << "process alignments: " << aln_count << endl;
		process_alignemnts(alignments, fasta);
		aln_count = 0;
	}

	alignments.clear();
	delete mapped_file;
	delete fasta;
	delete tmp_aln;

	print();
}

void ParseSNP::process_alignemnts(vector<vector<Alignment *> > & alignments,
		FastaParser * fasta) {

	for (size_t i = 0; i < alignments.size(); i++) {
		int snp_pos = 0;

		string ref = fasta->getChr(i);
		for (size_t j = 0; j < alignments[i].size(); j++) {
			int32_t pos = alignments[i][j]->getPosition();
			if (snp_pos < (int) snps[i].size()
					&& (snps[i][snp_pos] - range < pos && snps[i][snp_pos] + range
							> pos + (int) alignments[i][j]->getRefLength())) {

				if (pos + alignments[i][j]->getRefLength() >= ref.size()) {
					string reff = ref.substr((size_t) pos,
							(int) ref.size() - pos);
					while (reff.size() < alignments[i][j]->getRefLength()) {
						reff += 'N';
					}
					alignments[i][j]->setRef(reff);
				} else {
					alignments[i][j]->setRef(ref.substr((size_t) pos,alignments[i][j]->getRefLength()));
				}
				alignments[i][j]->computeAlignment();
				if (alignments[i][j]->getIdentity() >= 0.90) {
					compute_cov(alignments[i][j], covs[gen_key(i + 1, snps[i][snp_pos] + 1)]);
				}
			} else if (snps[i][snp_pos] + range < pos) {
				snp_pos++;
				if (snp_pos == (int) snps[i].size()) {
					break;
				}
			}
			delete alignments[i][j];
		}
		ref.clear();
		alignments[i].clear();
	}
}

void ParseSNP::compute_cov(Alignment * aln, coverage* cov) {

//	cout<<"P: "<<aln->getIdentity()<<" "<<aln->getRefID()<<endl;
//	cout<<aln->getSequence().first<<endl;
//	cout<<aln->getSequence().second<<endl;
	int pos = aln->getPosition() - cov->start_pos;

	if (pos > 0) {
		for (size_t i = 0;
				i < aln->getSequence().first.size() && i < (size_t) cov->start_pos + range * 2;
				i++) {
			if (aln->getSequence().second[i] != '-') {
				if (aln->getIdentity() == 1) {
					cov->cov_100[pos]++;
				}
				cov->cov_90[pos]++;
				pos++;
			} else {
				cov->cov_90[pos + 1]++;
			}
		}
	}

}

float ParseSNP::conf_val(int start_cov, int stop_cov, int max_cov) {
	//always: y2-y1=1
//
//	cout<<"\ts1: "<<(float)(stop_cov-start_cov)<<endl;
//	cout<<"\tsmax: "<<(float)(max_cov-start_cov)<<endl;

	if (max_cov - start_cov == 0 || stop_cov == 0 || (stop_cov - start_cov) == 0) {
		return 2; //-1.95;
	}
	return ((float) (stop_cov - start_cov) / (float) (max_cov - start_cov)) / step;

}

float ParseSNP::estimate(int * cov_100, int * cov_90, int pos) {

	int max_right = 0;
	int max_left = 0;
	int max_cov_90 = 0;
	for (int i = 0; i < read_length; i++) {

		if (i + pos < range * 2 && cov_100[i + pos] > max_right) {
			max_right = cov_100[i + pos];
		}
		if (i - pos >= 0 && cov_100[pos - i] > max_left) {
			max_left = cov_100[pos - i];
		}
		if (i + pos < range * 2 && cov_90[i + pos] > max_cov_90) {
			max_cov_90 = cov_90[i + pos];
		}
		if (i - pos >= 0 && cov_90[pos - i] > max_cov_90) {
			max_cov_90 = cov_90[pos - i];
		}
	}

	//check the slope of the left and right angle:
	float r1 = conf_val(cov_100[pos], cov_100[pos + step], max_right);
	float r2 = conf_val(cov_100[pos], cov_100[pos - step], max_left);
	float r_tot; //=1-(0.05-(r1+r2)/2);

	if (r1 < 0 || r2 < 0) {
		r_tot = 0;
	} else {
		r_tot = 1 - ((r1 + r2) / 2);
	}

	//check the cov at position of cov_90:
	float r_cov = ((float) cov_90[pos] / (float) max_cov_90);

//	cout<<"\tr: "<<(r_cov+r_tot)/2<<endl;
	if ((r_cov + r_tot) / 2 < 0) {
		return 0;
	} else {
		return (r_cov + r_tot) / 2;
	}
}

void ParseSNP::print() {
	cout << "printing:" << endl;
	FILE *file;
	file = fopen(output.c_str(), "w");
	if (file == NULL) {
		cout << "Error in printing: The file or path that you set "
				<< output.c_str()
				<< " is not valid. It can be that there is no disc space available."
				<< endl;
		exit(0);
	}

	ifstream myfile;
	myfile.open(snpfile.c_str(), ifstream::in);
	if (!myfile.good()) {
		cout << "Could not open file: " << snpfile.c_str()
				<< endl;
		exit(0);
	}

	myfile.getline(buffer, buffer_size);

	while (!myfile.eof()) {
		if (buffer[0] != '#') {
			int count = 0;
			string chr;
			int pos = 0;
			float score = -1;
			int id = 0;
			string line;
			string key;
			for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
				if (count == 0 && buffer[i] != '\t') {
					chr += buffer[i];
				}

				if (count == 1 && buffer[i - 1] == '\t') {

					pos = atoi(&buffer[i]);
					if (chrs.find(chr.c_str()) != chrs.end()) { //found
						id = chrs[chr.c_str()];
					} else {
						cerr << "Chr not found: " << chr << endl;
						exit(0);
					}
					key = gen_key(id + 1, pos);
					if (covs.find(key.c_str()) != covs.end()) {
						score = estimate(covs[key]->cov_100, covs[key]->cov_90, range);

						if (score > 1 || score < 0) {
							score = 0;
						}
						covs[key]->score = score;
					} else {
						score = -1; //default
					}
				}

				if (count == 5 && buffer[i - 1] == '\t' && score != -1) { //SNP Quality
					float snp_qual = atof(&buffer[i]);
					score = score * 255;
					score = (score + snp_qual) / 2;
					covs[key]->score = score;
//					cout << "HIT: " << snp_qual << " " << score << endl;
					ostringstream ss;
					ss << score;
					string s(ss.str());
					if (score > 99) {
						line += s.substr(0, 5);
					} else { //if(score >9){
						line += s.substr(0, 4);
					}
					line += '\t';

				} else if (count != 5 || score == -1) {
					line += buffer[i];
				}

				if (buffer[i] == '\t') {
					count++;
				}
			}
			fprintf(file, "%s", line.c_str());
			fprintf(file, "%c", '\n');
			line.clear();
			key.clear();
		} else {
			fprintf(file, "%s", buffer);
			fprintf(file, "%c", '\n');
		}

		myfile.getline(buffer, buffer_size);
	}
	myfile.close();
	fclose(file);

	if (plot_file[0] != '#') {
		file = fopen(plot_file.c_str(), "w");
		if (file == NULL) {
			cout << "Error in printing: The file or path that you set "
					<< plot_file.c_str()
					<< " is not valid. It can be that there is no disc space available."
					<< endl;
			exit(0);
		}

		for (map<string, coverage*>::iterator i = covs.begin(); i != covs.end(); i++) {
			fprintf(file, "%s", (*i).first.c_str());
			fprintf(file, "%c", '\t');
			fprintf(file, "%f", (*i).second->score);
			fprintf(file, "%c", '\t');
			for (size_t j = 2; (int) j < range * 2; j++) {
				fprintf(file, "%i", (*i).second->cov_100[j]);
				fprintf(file, "%c", '\t');
			}
			fprintf(file, "%c", '\n');
			fprintf(file, "%s", (*i).first.c_str());
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", -1);
			fprintf(file, "%c", '\t');
			for (size_t j = 2; (int) j < range * 2; j++) {
				fprintf(file, "%i", (*i).second->cov_90[j]);
				fprintf(file, "%c", '\t');
			}
			fprintf(file, "%c", '\n');
		}
		fclose(file);
	}
	cout << "Printing finished" << endl;
	for (map<string, coverage*>::iterator i = covs.begin(); i != covs.end(); i++) {
		delete [] (*i).second->cov_100;
		delete [] (*i).second->cov_90;
		delete (*i).second;
	}
}
