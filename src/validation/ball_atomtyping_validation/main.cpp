#include <BALL/KERNEL/system.h>
#include <BALL/FORMAT/PDBFile.h>
#include <BALL/FORMAT/MOL2File.h>
#include <BALL/MOLMEC/AMBER/GAFFTypeProcessor.h>
#include <BALL/STRUCTURE/assignBondOrderProcessor.h>

#include <BALL/FORMAT/genericMolFile.h>
#include <BALL/FORMAT/molFileFactory.h>

#include <iostream>

using namespace BALL;
using namespace std;

int main(int argc, char** argv) {

	if (argc != 2)
	{
		cerr << "Need one filename! Aborting." << endl;
		return(1);
	}

    GenericMolFile* infile = MolFileFactory::open(argv[1]);

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
	
	String infile_name(argv[1]); // ("../../../../../huettel-msc/export_folder/gaff_files/1_1,3-Butadien_Conformer3D_CID_7845.mol2"); 
	String outfile_name(infile_name);
	outfile_name.reverse();
	outfile_name = outfile_name.after(".");
	outfile_name.reverse();
	outfile_name = outfile_name + "out.mol2";

    // MOL2File f(infile_name);
    System S;
    *infile >> S;

    Options options;
    options[GAFFTypeProcessor::Option::ATOMTYPE_FILENAME] = "atomtyping/GAFFTypes.dat";

    GAFFTypeProcessor gt(options);
    S.apply(gt);

    for (AtomIterator at_it = S.beginAtom(); +at_it; ++at_it)
		std::cout << "atom name: " << at_it->getName() << " atomtype  " << at_it->getProperty("atomtype").getString() << std::endl;

    AssignBondOrderProcessor abp;
	S.apply(abp);

    cout << "System contains " << S.countAtoms() << " atoms." << endl;
    
	GenericMolFile* outfile = MolFileFactory::open(outfile_name, std::ios::out);
    *outfile << S;
    outfile->close();

	// MOL2File mol2_file(outfile_name, ios::out);
	// mol2_file << S;
	// mol2_file.close();
	// cout << "Wrote MOL2 file " << outfile_name << endl;

    return 0;
}