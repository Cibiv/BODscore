/*
 * Indten_hist.cpp
 *
 *  Created on: Sep 27, 2012
 *      Author: fritz
 */

#include "Ident_hist.h"

void Ident_hist::parse_cmd_line(int argc, char *argv[]) {

	CmdLine cmdss("Comp Motiv", ' ', version.c_str());
	ValueArg<string> read_file("q", "mapped_reads",
			"Mapped read file (sam/bam)", true, "#", "string");
	cmdss.add(read_file);
	ValueArg<string> ref_file("r", "ref_file", "Reference file (fasta)", true,
			"#", "string");
	cmdss.add(ref_file);
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

	} catch (ArgException &e) {
		StdOutput out;
		out.failure(cmdss, e);
	}
}

void Ident_hist::compute() {
	Parser * mapped_file = 0;
	if (read_filename.find(".bam") != string::npos) {
		mapped_file = new BamParser(read_filename);
	}

	if (read_filename.find(".sam") != string::npos) {
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
	cout<<"start parsing"<<endl;
	while (!tmp_aln->getSequence().first.empty()) {
		while ((size_t) tmp_aln->getRefID() >= alignments.size()) {
			vector<Alignment*> tmp;
			alignments.push_back(tmp);
		}
		alignments[tmp_aln->getRefID()].push_back(tmp_aln);
		aln_count++;
		if (aln_count >= MAX_NUM_ALIGNMENTS) {
			//process alignments
			cout<<"process alignments: "<<aln_count<<endl;
			process_alignemnts(alignments, fasta);
			aln_count = 0;
			cout<<"parsing"<<endl;
		}
		tmp_aln = mapped_file->parseRead(0);
	}
	cout<<"end parsing"<<endl;
	if (aln_count != 0) {
		//process alignments
		cout<<"process alignments: "<<aln_count<<endl;
		process_alignemnts(alignments, fasta);
		aln_count = 0;
	}

	alignments.clear();
	delete mapped_file;
	delete fasta;
	print();

}

void Ident_hist::process_alignemnts(vector<vector<Alignment *> > & alignments,
		FastaParser * fasta) {

	for (size_t i = 0; i < alignments.size(); i++) {
		if (!alignments[i].empty()) {
			string ref = fasta->getChr(i);
			for (size_t j = 0; j < alignments[i].size(); j++) {
				int32_t pos = alignments[i][j]->getPosition();

				if (pos + alignments[i][j]->getRefLength() >= ref.size()) {
					string reff = ref.substr((size_t) pos,
							(int) ref.size() - pos);
					while (reff.size() < alignments[i][j]->getRefLength()) {
						reff += 'N';
					}
					alignments[i][j]->setRef(reff);
				} else {
					alignments[i][j]->setRef(
							ref.substr((size_t) pos,
									alignments[i][j]->getRefLength()));
				}

				alignments[i][j]->computeAlignment();
				compute_cov(alignments[i][j]);
				delete alignments[i][j];
			}
			alignments[i].clear();
		}
	}
}

void Ident_hist::compute_cov(Alignment * aln) {
	float id = aln->getIdentity();
	id = round(id * 100);
	ident[id]++;
}

void Ident_hist::print() {
	cout << "printing" << endl;
	FILE *file;
	file = fopen(output.c_str(), "w");
	if (file == NULL) {
		cout << "Error in printing: The file or path that you set "
				<< output.c_str()
				<< " is not valid. It can be that there is no disc space available."
				<< endl;
		exit(0);
	}

	for (size_t i = 1; i != ident.size(); i++) {
		fprintf(file, "%i", ident[i]);
		fprintf(file, "%c", '\n');
	}
	fclose(file);
}
