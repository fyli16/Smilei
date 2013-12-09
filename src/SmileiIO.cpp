/*
 * SmileiIO.cpp
 *
 *  Created on: 3 juil. 2013
 *      Author: jderouil
 */

#include "SmileiIO.h"

#include "PicParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "Species.h"

#include <sstream>

using namespace std;

SmileiIO::SmileiIO( PicParams* params, SmileiMPI* smpi )
{
	ostringstream name;
	name.str(""); name << "particles-" << smpi->getRank() << ".h5" ;

	hid_t attribute_id;

	// Create 1 file containing 1 dataset per Species
	partFile_id = H5Fcreate( name.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	particleSize = params->nDim_particle + 3 + 1;
	hsize_t dims[2] = {0, particleSize};
	hsize_t max_dims[2] = {H5S_UNLIMITED, particleSize};
	hid_t file_space = H5Screate_simple(2, dims, max_dims);

	hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
	H5Pset_layout(plist, H5D_CHUNKED);
	hsize_t chunk_dims[2] = {1, particleSize};
	H5Pset_chunk(plist, 2, chunk_dims);

	nDatasetSpecies = params->n_species;
	partDataset_id = new hid_t[nDatasetSpecies];
	for (unsigned int ispec=0 ; ispec<nDatasetSpecies ; ispec++) {
		ostringstream speciesName;
		speciesName.str(""); speciesName << params->species_param[ispec].species_type;
		partDataset_id[ispec] = H5Dcreate(partFile_id, speciesName.str().c_str(), H5T_NATIVE_FLOAT, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);

		hsize_t scalaire = 1;
		hid_t tmp_space = H5Screate_simple(1, &scalaire, NULL);

		attribute_id = H5Acreate2 (partDataset_id[ispec], "Mass", H5T_IEEE_F64BE, tmp_space,
					   	   	   	   	   	 H5P_DEFAULT, H5P_DEFAULT);
		H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &params->species_param[ispec].mass);
		H5Aclose(attribute_id);

		attribute_id = H5Acreate2 (partDataset_id[ispec], "Charge", H5T_IEEE_F64BE, tmp_space,
					   	   	   	   	   	 H5P_DEFAULT, H5P_DEFAULT);
		H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &params->species_param[ispec].charge);
		H5Aclose(attribute_id);

		H5Sclose(tmp_space);
	}

	H5Pclose(plist);
	H5Sclose(file_space);

	dims[0] = 1;
	dims[1] = particleSize;
	partMemSpace = H5Screate_simple(2, dims, NULL);

	// Management of global IO file
	MPI_Info info  = MPI_INFO_NULL;
	hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);
	global_file_id_ = H5Fcreate( "Fields.h5", H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
	H5Pclose(plist_id);

}

SmileiIO::~SmileiIO()
{
	// Management of global IO file
	H5Fclose( global_file_id_ );

	H5Sclose(partMemSpace);
	for ( unsigned int s=0 ; s<nDatasetSpecies ; s++ )
		H5Dclose(partDataset_id[s]);
	delete partDataset_id;
	H5Fclose(partFile_id);

}

// ---------------------------------------------------------------------------------------------------------------------
// For time iteration "time", write all fields in the same file
// Version to keep
// ---------------------------------------------------------------------------------------------------------------------
void SmileiIO::writeAllFields( ElectroMagn* EMfields, int time )
{
	MPI_Info info  = MPI_INFO_NULL;
	hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);
	ostringstream name_t;
	name_t.str(""); name_t << "fields_" << time << ".h5";
	hid_t file_id = H5Fcreate( name_t.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
	H5Pclose(plist_id);

	writeFieldsSingleFile( EMfields->Ex_, file_id, time );
	writeFieldsSingleFile( EMfields->Ey_, file_id, time );
	writeFieldsSingleFile( EMfields->Ez_, file_id, time );
	writeFieldsSingleFile( EMfields->Bx_, file_id, time );
	writeFieldsSingleFile( EMfields->By_, file_id, time );
	writeFieldsSingleFile( EMfields->Bz_, file_id, time );
	writeFieldsSingleFile( EMfields->Jx_, file_id, time );
	writeFieldsSingleFile( EMfields->Jy_, file_id, time );
	writeFieldsSingleFile( EMfields->Jz_, file_id, time );
	writeFieldsSingleFile( EMfields->rho_, file_id, time );

	H5Fclose( file_id );

}

// ---------------------------------------------------------------------------------------------------------------------
// Write all fields of all time step in the same file
// In progress ...
// ---------------------------------------------------------------------------------------------------------------------
void SmileiIO::writeAllFieldsSingleFileTime( ElectroMagn* EMfields, int time )
{
	ostringstream name_t;
	name_t.str(""); name_t << "/T =" << time;

	DEBUG("GROUP _________________________________ " << name_t.str());
	hid_t group_id = H5Gcreate2(global_file_id_, name_t.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	writeFieldsSingleFileTime( EMfields->Ex_, group_id );
	writeFieldsSingleFileTime( EMfields->Ey_, group_id );
	writeFieldsSingleFileTime( EMfields->Ez_, group_id );
	writeFieldsSingleFileTime( EMfields->Bx_, group_id );
	writeFieldsSingleFileTime( EMfields->By_, group_id );
	writeFieldsSingleFileTime( EMfields->Bz_, group_id );
	writeFieldsSingleFileTime( EMfields->Jx_, group_id );
	writeFieldsSingleFileTime( EMfields->Jy_, group_id );
	writeFieldsSingleFileTime( EMfields->Jz_, group_id );
	writeFieldsSingleFileTime( EMfields->rho_, group_id );

	H5Gclose(group_id);

    H5Fflush( global_file_id_, H5F_SCOPE_GLOBAL );

}

// ---------------------------------------------------------------------------------------------------------------------
// Write fields in separated files, data are overwritten at each call ( field.h5 )
// Kept while python tools not updated
// ---------------------------------------------------------------------------------------------------------------------
void SmileiIO::writeFields( ElectroMagn* EMfields )
{
	// File opened in constructor
	write( EMfields->Ex_, "ex.h5" );
	write( EMfields->Ey_, "ey.h5" );
	write( EMfields->Ez_, "ez.h5" );
	write( EMfields->Bx_, "bx.h5" );
	write( EMfields->By_, "by.h5" );
	write( EMfields->Bz_, "bz.h5" );
	write( EMfields->Jx_, "jx.h5" );
	write( EMfields->Jy_, "jy.h5" );
	write( EMfields->Jz_, "jz.h5" );
	write( EMfields->rho_, "rho.h5" );
}

// ---------------------------------------------------------------------------------------------------------------------
// For time iteration it, write fields in separated files ( field_time.h5 )
// Kept while python tools not updated
// Only in 2D
// ---------------------------------------------------------------------------------------------------------------------
void SmileiIO::writeFields( ElectroMagn* EMfields, double time )
{
	// File opened in constructor
	write( EMfields->Ex_, "ex", time );
	write( EMfields->Ey_, "ey", time );
	write( EMfields->Ez_, "ez", time );
	write( EMfields->Bx_, "bx", time );
	write( EMfields->By_, "by", time );
	write( EMfields->Bz_, "bz", time );
	write( EMfields->Jx_, "jx", time );
	write( EMfields->Jy_, "jy", time );
	write( EMfields->Jz_, "jz", time );
	write( EMfields->rho_, "rho", time );
}

// ---------------------------------------------------------------------------------------------------------------------
// Debug tools
// Each MPI process write its data (including ghost cells) in separated files ( field_MPI_Rank.h5 ), data are overwritten at each call
// Kept while python tools not updated
// Only in 2D
// ---------------------------------------------------------------------------------------------------------------------
void SmileiIO::writeFieldsPP( ElectroMagn* EMfields, double time, int rank )
{
	// Each process write is own "name_mpirank.h5" 
	writePerProcess( EMfields->rho_, "rho", time, rank );
	writePerProcess( EMfields->Ex_, "Ex", time, rank );
	writePerProcess( EMfields->Ey_, "Ey", time, rank );
	writePerProcess( EMfields->Ez_, "Ez", time, rank );
	writePerProcess( EMfields->Bx_, "Bx", time, rank );
	writePerProcess( EMfields->By_, "By", time, rank );
	writePerProcess( EMfields->Bz_, "Bz", time, rank );
	writePerProcess( EMfields->Jx_, "Jx", time, rank );
	writePerProcess( EMfields->Jy_, "Jy", time, rank );
	writePerProcess( EMfields->Jz_, "Jz", time, rank );
}

// ---------------------------------------------------------------------------------------------------------------------
// Each MPI process writes is particles in its own file, data are overwritten at each call ( particles-MPI_Rank.h5 )
// Kept while python tools not updated
// ---------------------------------------------------------------------------------------------------------------------
void SmileiIO::writePlasma( vector<Species*> vecSpecies, double time, SmileiMPI* smpi )
{
	MESSAGE("write species disabled");
	return;
	int n_species = vecSpecies.size();
	for (int ispec=0 ; ispec<n_species ; ispec++) {
		std::vector<Particle*>* cuParticles = &(vecSpecies[ispec])->particles;
		cout << "write species " << ispec << endl;

		for (unsigned int p=0; p<(vecSpecies[ispec])->getNbrOfParticles(); p++ ) {

			hid_t file_space = H5Dget_space(partDataset_id[ispec]);
			hsize_t dimsO[2];
			H5Sget_simple_extent_dims(file_space, dimsO, NULL);
			H5Sclose(file_space);
			hsize_t dims[2];
			dims[0] = dimsO[0]+1;
			dims[1] = dimsO[1];
			H5Dset_extent(partDataset_id[ispec], dims);

			file_space = H5Dget_space(partDataset_id[ispec]);
			hsize_t start[2];
			hsize_t count[2] = {1, particleSize};
			start[0] = dimsO[0];
			start[1] = 0;
			H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL);

			H5Dwrite(partDataset_id[ispec], H5T_NATIVE_DOUBLE, partMemSpace, file_space, H5P_DEFAULT, &((*cuParticles)[ p ]->position(0)));

			H5Sclose(file_space);


		} // End for p

	} // End for ispec

}
