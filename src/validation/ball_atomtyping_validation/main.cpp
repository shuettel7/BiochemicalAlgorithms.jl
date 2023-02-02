#include <BALL/KERNEL/system.h>
#include <BALL/FORMAT/PDBFile.h>
#include <BALL/FORMAT/MOL2File.h>
#include <BALL/MOLMEC/AMBER/GAFFTypeProcessor.h>
#include <BALL/STRUCTURE/assignBondOrderProcessor.h>

#include <BALL/FORMAT/genericMolFile.h>
#include <BALL/FORMAT/molFileFactory.h>

#include <fstream>
#include <string>
#include <iostream>
#include <filesystem>
#include <typeinfo>

namespace fs = std::filesystem;

using namespace BALL;

int main() {	//int argc, char** argv

	// if (argc != 2)
	// {
	// 	std::cerr << "Need one filename! Aborting." << std::endl;
	// 	return(1);
	// }

	const fs::path path {"/mnt/c/Users/Sam/source/repos/huettel-msc/gaff_files"};
	std::cout << path << typeid(path).name() << std::endl;
	
	for (const auto & entry : fs::directory_iterator{path}){

		std::string entry_path {entry.path()};
		
		if (entry_path.find(".mol2") == string::npos){
			continue;
		}

		GenericMolFile* infile = MolFileFactory::open(entry_path);

		if (!infile)
		{
			Log.error() << "Could not determine filetype, aborting" << std::endl;
			return 2;
		}

		if (!*infile)
		{
			std::cerr << "Invalid file, aborting" << std::endl;
			return 2;
		}
		
		std::string infile_name = entry_path; 
		// ("../../../../../huettel-msc/export_folder/gaff_files/1_1,3-Butadien_Conformer3D_CID_7845.mol2"); 
		std::string basename_infile = entry_path.substr(entry_path.find_last_of("/\\") + 1);
		std::string outfile_name = "/mnt/c/Users/Sam/source/repos/huettel-msc/gaff_files/ballC_output/" + basename_infile;

		// MOL2File f(infile_name);
		System S;
		*infile >> S;

		// Important: Cleanup
		infile->close();
		delete infile;

		Options options;
		options[GAFFTypeProcessor::Option::ATOMTYPE_FILENAME] = "atomtyping/GAFFTypes.dat";

		GAFFTypeProcessor gt(options);
		S.apply(gt);

		for (AtomIterator at_it = S.beginAtom(); +at_it; ++at_it)
			std::cout << "atom name: " << at_it->getName() << " atomtype  " << at_it->getProperty("atomtype").getString() << std::endl;

		AssignBondOrderProcessor abp;
		S.apply(abp);

		std::cout << "System contains " << S.countAtoms() << " atoms." << std::endl;
		
		GenericMolFile* outfile = MolFileFactory::open(outfile_name, std::ios::out);
		*outfile << S;
		outfile->close();
		delete outfile;
		std::cout << "Wrote MOL2 file " << outfile_name << std::endl;
	}

	return 0;	
}