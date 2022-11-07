#include "CalorimeterConstruction.hh"

#include "G4Box.hh"
#include "G4ChordFinder.hh"
#include "G4Element.hh"
#include "G4NistManager.hh"
#include "G4FieldManager.hh"
#include "G4LogicalVolume.hh"

#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4CSGSolid.hh"
#include "G4Tubs.hh"
#include "G4CutTubs.hh"
#include "G4Cons.hh"

#include "G4SystemOfUnits.hh"
#include "G4IntersectionSolid.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4PVParameterised.hh"
// #include "CaloRCellParameterisation.hh"
// #include "DetectorGeometryDefinitions.hh"
#include "G4ReflectionFactory.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"

#include "G4RunManager.hh"

#include "DetectorConstruction.hh"

#include <string>
#include <algorithm>
#include <map>

CalorimeterConstruction::CalorimeterConstruction(G4LogicalVolume *expHallLV, bool fCheck_Overlaps, Geometry_definition Geometry)
{
	geometry = Geometry;
	GlobalLV = expHallLV;
	theta_min = 2 * atan(exp(-1 * config_json_var.max_eta_barrel));
	fCheckOverlaps = fCheck_Overlaps;
	length = geometry.layer_out_radius_HCAL.back().back() / tan(theta_min);
	
	G4NistManager *nistManager = G4NistManager::Instance();
	m_iron = nistManager->FindOrBuildMaterial("G4_Fe");
	
	Barrel_Calorimeter();
	EndCap_Calorimeter();
}
CalorimeterConstruction::~CalorimeterConstruction()
{
	;
}
char *CalorimeterConstruction::Name_creation(char *name, int low_layer, int high_layer)
{
    //
    // Assumption: less than 10 layers.
    //
	name[4] = (low_layer + 49);
	name[6] = (high_layer + 49);
	return name;
}

void CalorimeterConstruction::EndCap_Calorimeter()
{
    
	long double r_inn                  = geometry.layer_inn_radius_ECAL[0][0];
        long double previous_layers_depths = 0.0;
	int         nPixelsMax             = GetNPixelsMax();
	long double minDPhi                = GetMinDPhi();
	long double depth                  = 0.0;

	//
	// ECAL endcap
	//
	int nLow_Layers = geometry.number_of_pixels_ECAL.size();
	for (int ilow_layer = 0; ilow_layer < nLow_Layers; ilow_layer++)
	{
		int nHigh_Layers = geometry.number_of_pixels_ECAL.at(ilow_layer).size();
		for (int ihigh_layer = 0; ihigh_layer < nHigh_Layers; ihigh_layer++)
		{
		    depth = geometry.layer_out_radius_ECAL[ilow_layer][ihigh_layer] - geometry.layer_inn_radius_ECAL[ilow_layer][ihigh_layer];
		    Build_EndCap_CAL( nPixelsMax,
				      nPixelsMax / geometry.number_of_pixels_ECAL[ilow_layer][ihigh_layer],
				      minDPhi,
				      r_inn, depth, previous_layers_depths, config_json_var.Material_ECAL, ECAL1_VisAtt,
				      Name_creation(strdup("ECALN_N_Endcap_forward_LV"), ilow_layer, ihigh_layer),
				      Name_creation(strdup("ECALN_N_Endcap_forward_PL"), ilow_layer, ihigh_layer), 1 );
		    Build_EndCap_CAL( nPixelsMax,
		    		      nPixelsMax / geometry.number_of_pixels_ECAL[ilow_layer][ihigh_layer],
		    		      minDPhi,
		    		      r_inn, depth, previous_layers_depths, config_json_var.Material_ECAL, ECAL1_VisAtt,
		    		      Name_creation(strdup("ECALN_N_Endcap_back_LV"), ilow_layer, ihigh_layer), 
		    		      Name_creation(strdup("ECALN_N_Endcap_back_PL"), ilow_layer, ihigh_layer), -1 );
		    previous_layers_depths += depth;
		}
	}
	//
	// Fill the region between the ECAL and HCAL with iron as support
	//
	depth = geometry.layer_inn_radius_HCAL.front().front() - geometry.layer_out_radius_ECAL.back().back();
	std::map<int, std::string> direction_name = {
	    std::make_pair(-1, "backward"),
	    std::make_pair(1, "forward")
	};
	for( auto dir_str : direction_name ) {
	    Build_EndCap_CAL( nPixelsMax,
			      1,
			      minDPhi,
			      r_inn,
			      depth,
			      previous_layers_depths,
			      m_iron,
			      IronGap_VisAtt,
			      ( "LV_ironGap_endcap_" + dir_str.second ).c_str(),
			      ( "PV_ironGap_endcap_" + dir_str.second ).c_str(),
			      dir_str.first );
	}
	previous_layers_depths += depth;
	//
	// HCAL endcap
	//
	nLow_Layers = geometry.number_of_pixels_HCAL.size();
	for (int ilow_layer = 0; ilow_layer < nLow_Layers; ilow_layer++)
	{
		int nHigh_Layers = geometry.number_of_pixels_HCAL.at(ilow_layer).size();
		for (int ihigh_layer = 0; ihigh_layer < nHigh_Layers; ihigh_layer++)
		{
			depth = geometry.layer_out_radius_HCAL[ilow_layer][ihigh_layer] - geometry.layer_inn_radius_HCAL[ilow_layer][ihigh_layer];
			Build_EndCap_CAL( nPixelsMax,
					  nPixelsMax / geometry.number_of_pixels_HCAL[ilow_layer][ihigh_layer],
					  minDPhi,
					  r_inn, depth, previous_layers_depths, config_json_var.Material_HCAL, HCAL1_VisAtt,
					  Name_creation(strdup("HCALN_N_Endcap_forward_LV"), ilow_layer, ihigh_layer), 
					  Name_creation(strdup("HCALN_N_Endcap_forward_PL"), ilow_layer, ihigh_layer), 1 );
			Build_EndCap_CAL( nPixelsMax,
					  nPixelsMax / geometry.number_of_pixels_HCAL[ilow_layer][ihigh_layer],
					  minDPhi,
					  r_inn, depth, previous_layers_depths, config_json_var.Material_HCAL, HCAL1_VisAtt, 
					  Name_creation(strdup("HCALN_N_Endcap_back_LV"), ilow_layer, ihigh_layer), 
					  Name_creation(strdup("HCALN_N_Endcap_back_PL"), ilow_layer, ihigh_layer), -1 );
			previous_layers_depths += depth;			
		}
	}
}

void CalorimeterConstruction::Barrel_Calorimeter()
{
        long double previous_layers_delta_r = 0.0;
	int         nLow_Layers             = geometry.number_of_pixels_ECAL.size();
	int         nPixelsMax              = GetNPixelsMax();
	long double minDEta                 = 4.0 * config_json_var.max_eta_barrel / nPixelsMax;
	long double minDPhi                 = GetMinDPhi();

	long double r_inn;
	long double r_out;
			
	//
	// ECAL barrel
	//
	for (int ilow_layer = 0; ilow_layer < nLow_Layers; ilow_layer++)
	{
		int nHigh_Layers = geometry.number_of_pixels_ECAL.at(ilow_layer).size();
		for (int ihigh_layer = 0; ihigh_layer < nHigh_Layers; ihigh_layer++)
		{
		    
			r_inn = geometry.layer_inn_radius_ECAL[ilow_layer][ihigh_layer];
			r_out = geometry.layer_out_radius_ECAL[ilow_layer][ihigh_layer];
			
			Build_Barrel_CAL(nPixelsMax,
					 nPixelsMax / geometry.number_of_pixels_ECAL[ilow_layer][ihigh_layer],
					 minDEta,
					 minDPhi,
					 r_inn, r_out, previous_layers_delta_r, config_json_var.Material_ECAL, ECAL1_VisAtt,
					 Name_creation(strdup("ECALN_N_forward_LV"), ilow_layer, ihigh_layer), 
					 Name_creation(strdup("ECALN_N_forward_PL"), ilow_layer, ihigh_layer), 1);
			Build_Barrel_CAL(nPixelsMax,
					 nPixelsMax / geometry.number_of_pixels_ECAL[ilow_layer][ihigh_layer],
					 minDEta,
					 minDPhi,
					 r_inn, r_out, previous_layers_delta_r, config_json_var.Material_ECAL, ECAL1_VisAtt,
					 Name_creation(strdup("ECALN_N_back_LV"), ilow_layer, ihigh_layer), 
					 Name_creation(strdup("ECALN_N_back_PL"), ilow_layer, ihigh_layer), -1);
			previous_layers_delta_r += r_out - r_inn;
		}
	}
	//
	// Fill the region between the ECAL and HCAL with iron as support
	//
	r_inn = geometry.layer_out_radius_ECAL.back().back();
	r_out = geometry.layer_inn_radius_HCAL.front().front();
	
	std::map<int, std::string> direction_name = {
	    std::make_pair(-1, "backward"),
	    std::make_pair(1, "forward")
	};
	for( auto dir_str : direction_name ) {
			Build_Barrel_CAL(nPixelsMax,
					 1,
					 minDEta,
					 minDPhi,
					 r_inn,
					 r_out,
					 previous_layers_delta_r,
					 m_iron,
					 IronGap_VisAtt,
					 ( "LV_ironGap_barrel_" + dir_str.second ).c_str(),
					 ( "PV_ironGap_barrel_" + dir_str.second ).c_str(),
					 dir_str.first );
				    
	}
	previous_layers_delta_r += r_out - r_inn;
	//
	// HCAL barrel
	//
	nLow_Layers            = geometry.number_of_pixels_HCAL.size();
	for (int ilow_layer = 0; ilow_layer < nLow_Layers; ilow_layer++)
	{
		int nHigh_Layers = geometry.number_of_pixels_HCAL.at(ilow_layer).size();
		for (int ihigh_layer = 0; ihigh_layer < nHigh_Layers; ihigh_layer++)
		{
		    
			r_inn = geometry.layer_inn_radius_HCAL[ilow_layer][ihigh_layer];
			r_out = geometry.layer_out_radius_HCAL[ilow_layer][ihigh_layer];
			
			Build_Barrel_CAL( nPixelsMax,
					  nPixelsMax / geometry.number_of_pixels_HCAL[ilow_layer][ihigh_layer],
					  minDEta,
					  minDPhi,
					  r_inn, r_out, previous_layers_delta_r, config_json_var.Material_HCAL, HCAL1_VisAtt,
					  Name_creation(strdup("HCALN_N_forward_LV"), ilow_layer, ihigh_layer),
					  Name_creation(strdup("HCALN_N_forward_PL"), ilow_layer, ihigh_layer), 1 );
			Build_Barrel_CAL( nPixelsMax,
					  nPixelsMax / geometry.number_of_pixels_HCAL[ilow_layer][ihigh_layer],
					  minDEta,
					  minDPhi,
					  r_inn, r_out, previous_layers_delta_r, config_json_var.Material_HCAL, HCAL1_VisAtt,
					  Name_creation(strdup("HCALN_N_back_LV"), ilow_layer, ihigh_layer),
					  Name_creation(strdup("HCALN_N_back_PL"), ilow_layer, ihigh_layer), -1 );
			previous_layers_delta_r += r_out - r_inn;
		}
	}
}


void CalorimeterConstruction::Build_Barrel_CAL(int NumberOfPixel, int cellMergeFactor, long double d_eta, long double d_phi, long double r_inn, long double r_out, long double previous_layers_delta_r, G4Material *Material_CAL, G4VisAttributes *VisAtt, const char *LV, const char *PL, int direction)
{

    //
    // Build a layer of a quarter of the barrel calorimeter.
    //
    //   * NumberOfPixel           : number of cells per layer and per phi slice.
    //                               Hence one barrel quarter created by this function
    //                               has 0.25 * NumberOfPixel cells ( the others are in the opposing barrel part and in the two endcaps ).
    //
    //   * Cell merge factor       : If layers have different granularities,
    //                               a proper geometry is built by first creating
    //                               cells in each layer according to the highest granularity,
    //                               followed by cell merging to arrive at the desired granularity.
    // 
    //   * previous_layers_delta_r : pass 0.0 for the first layer of the ECAL. Otherwise pass the sum of depths of all previous layers.
    //

        CheckMergeFactor( NumberOfPixel, cellMergeFactor );
    
	G4RotationMatrix *zRot = new G4RotationMatrix; // Rotates X and Z axes only
	zRot->rotateZ(0);
	G4ReflectZ3D reflection(0);
	G4CSGSolid*         cell_tube_subtrahend = nullptr;
	G4CutTubs*          cell_tube_minuend    = nullptr;
	std::vector<G4VSolid*>  cells_maxNPixels;
	std::vector<G4VSolid*>* ptr_cells_final;
	G4LogicalVolume *CAL_LV;
	G4int index_CAL = 0;
	long double vy_rear  = 0.0;
	long double vz_rear  = 0.0;
	long double vy_front = 0.0;
	long double vz_front = 0.0;
	G4ThreeVector vector_rear;
	G4ThreeVector vector_front;
	G4ThreeVector vector_z( 0.0, 0.0, -direction );
	long double   subtrahendTranslationZ = - direction * 0.1 * length;
	G4ThreeVector vecSubtrahendTranslationZ( 0.0, 0.0, subtrahendTranslationZ );

	const long double depth = r_out - r_inn;

	const size_t nCellsPerDirection = 0.25 * NumberOfPixel;
	
	std::vector<long double> theta_all( nCellsPerDirection, 0.0); // Angles between the cell z-surface away from the IP and the xy-plane.
	for (size_t iCell = 0; iCell < nCellsPerDirection; ++iCell)
	    theta_all[iCell] = 0.5 * M_PI - 2 * atan( exp( -1.0 * ( iCell + 1 ) * d_eta ) );
	    
	for (size_t iCell = 0; iCell < nCellsPerDirection; ++iCell) //*loop that creates detector pixels in positive z direction
	{

	        long double theta_rear = 0.0;
		if ( iCell )
		    theta_rear          = theta_all[iCell - 1];
		long double theta_front = theta_all[iCell];

		if ( iCell == 1 )
		    r_inn = r_inn + previous_layers_delta_r * ( cos( theta_rear ) - 1.0 );
		if ( iCell > 1 )
		    r_inn = r_inn + previous_layers_delta_r * ( cos( theta_rear ) - cos( theta_all[iCell - 2] ) );
		
	        r_out = r_inn + depth * cos( theta_rear ); // 1 / cosh( eta )  scaling. Note that the polar angle is ( pi / 2 - theta_rear ) .
	    
		vy_rear = -sin( theta_rear );
		vz_rear =  cos( theta_rear );
		vy_front = -sin( theta_front );
		vz_front =  cos( theta_front );		
		vector_rear.set(  0.0, vy_rear, direction * vz_rear);
		vector_front.set( 0.0, vy_front, direction * vz_front);

		if ( iCell == 0 ) {
		    cell_tube_subtrahend = (G4CSGSolid*)(new G4Tubs( "Barrel_Cell_subtrahend",
						       0.9 * r_inn,
						       1.1 * r_out,
						       0.5 * length + fabs( subtrahendTranslationZ ),
						       ( M_PI - d_phi ) / 2,
						       config_json_var.max_phi ));
		} else {
		    cell_tube_subtrahend = (G4CSGSolid*)(new G4CutTubs( "Barrel_Cell_subtrahend",
							  0.9 * r_inn,
							  1.1 * r_out,
							  0.5 * length + fabs( subtrahendTranslationZ ),
							  ( M_PI - d_phi ) / 2,
							  config_json_var.max_phi,
							  direction == 1 ? vector_z    : vector_rear,
							  direction == 1 ? vector_rear : vector_z ));
		}
		cell_tube_minuend    = new G4CutTubs( "Barrel_Cell_minuend",
						      r_inn, r_out,
						      0.5 * length,
						      ( M_PI - d_phi ) / 2,
						      d_phi,
						      direction == 1 ? vector_z     : vector_front,
						      direction == 1 ? vector_front : vector_z );
		G4SubtractionSolid* newCell = new G4SubtractionSolid( "cell",
								      cell_tube_minuend,
								      cell_tube_subtrahend,
								      0,
								      vecSubtrahendTranslationZ );
		cells_maxNPixels.push_back( newCell );
	}
	
	ptr_cells_final = &cells_maxNPixels;
	if ( cellMergeFactor > 1 )
	    ptr_cells_final = MergeCells( &cells_maxNPixels,
					  NumberOfPixel,
					  cellMergeFactor );
	int iEta = 0;
	for ( G4VSolid* cell : *ptr_cells_final ) {
	        ++iEta;
		CAL_LV               = new G4LogicalVolume( cell,
							    Material_CAL,
							    LV );
		int iPhi = 0;
		for (auto angle = 0; angle < NumberOfPixel; angle++) //* loop that creates detector pixels in phi axis
		{
		        ++iPhi;
			zRot = new G4RotationMatrix;
			zRot->rotateZ( angle * d_phi );
			new G4PVPlacement(
				zRot,
				G4ThreeVector(0.0, 0.0, -direction * 0.5 * length), //* Placed the pixel in specific way, that its slices goes throught (0 ,0, 0)
				CAL_LV,											    //* its logical volume
				( std::string( PL ) + "_" + std::to_string( iEta ) + "_" + std::to_string( iPhi ) ).c_str(), //* its name
				GlobalLV,										    //* its mother  volume
				false,											    //* no boolean operation
				index_CAL,										    //* copy number
				fCheckOverlaps);
			++index_CAL;
		}
		CAL_LV->SetVisAttributes(VisAtt);
	}
}

void CalorimeterConstruction::Build_EndCap_CAL(int NumberOfPixel, int cellMergeFactor, long double d_phi, long double r_inn_barrel, long double depth, long double previous_layers_depths, G4Material *Material_CAL, G4VisAttributes *VisAtt, const char *LV, const char *PL, int direction)
{

        CheckMergeFactor( NumberOfPixel, cellMergeFactor );
	
	std::vector<G4VSolid*>  cells_maxNPixels;
	std::vector<G4VSolid*>* ptr_cells_final;
	const int nCells           = NumberOfPixel / 4;
	const long double deltaEta = ( config_json_var.max_eta_endcap - config_json_var.max_eta_barrel ) / ( (long double)nCells );
	
        long double ipDistance_z       = r_inn_barrel / tan( EtaToTheta( config_json_var.max_eta_barrel ) );
	long double theta_low          = EtaToTheta( config_json_var.max_eta_endcap );
	long double theta_up           = theta_low;
	long double theta_min_layer    = theta_low;
	long double delta_z_1          = cos( theta_min_layer ) * depth;
	
	G4RotationMatrix *zRot       = new G4RotationMatrix;
	G4Cons           *cell_prime = nullptr;

	int index_EndCap = 0;

	long double              r_low_inner;
	long double              r_up_inner;
	long double              lz;
	long double              r_low_outer_prime;
	long double              r_up_outer_prime;
	long double              Delta_z;
	
	for( int iCellEta = 0; iCellEta < nCells; ++iCellEta ) {
	    //
	    // constructing cells going from highest to lowest eta
	    //
	    theta_low          = theta_up;
	    theta_up           = EtaToTheta( config_json_var.max_eta_endcap - ( iCellEta + 1 ) * deltaEta );
	    lz                 = ipDistance_z + previous_layers_depths * cos( theta_low );
	    Delta_z            = ( cos( theta_min_layer ) - cos( theta_low ) ) * previous_layers_depths;
	    r_low_inner        = lz * tan( theta_low );
	    r_up_inner         = lz * tan( theta_up );
	    r_low_outer_prime  = r_low_inner + ( 2.0 * Delta_z + delta_z_1 ) * tan( theta_low );
	    r_up_outer_prime   = r_up_inner  + ( 2.0 * Delta_z + delta_z_1 ) * tan( theta_up );
	    
	    cell_prime = new G4Cons( "EndcapCell",
				     direction == 1 ? r_low_inner : r_low_outer_prime,
				     direction == 1 ? r_up_inner : r_up_outer_prime,
				     direction == 1 ? r_low_outer_prime : r_low_inner,
				     direction == 1 ? r_up_outer_prime : r_up_inner,
				     0.5 * delta_z_1 + Delta_z,
				     (M_PI - d_phi) / 2,
				     d_phi );

	    
	    G4Tubs* tubs_subtrahend = new G4Tubs( "helper_tubs_subtrahend",
						  r_low_inner,
						  1.1 * r_up_outer_prime,
						  0.5 * delta_z_1 + Delta_z,
						  ( M_PI - d_phi ) / 2,
						  config_json_var.max_phi );
	    
	    long double   subtrahendTranslationZ = delta_z_1 + Delta_z - ( cos( theta_min_layer ) - cos( theta_low ) ) * ( previous_layers_depths + depth );
	    G4ThreeVector vecSubtrahendTranslationZ( 0.0, 0.0, direction * subtrahendTranslationZ );
	    cells_maxNPixels.push_back( new G4SubtractionSolid( "cell_endcap",
								cell_prime,
								tubs_subtrahend,
								0,
								vecSubtrahendTranslationZ ) );
	}
	ptr_cells_final = &cells_maxNPixels;
	if ( cellMergeFactor > 1 )
	    ptr_cells_final = MergeCells( &cells_maxNPixels,
					  NumberOfPixel,
					  cellMergeFactor );
	for (size_t iEta = 0; iEta < ptr_cells_final->size(); ++iEta) {
	    G4LogicalVolume *endCap_LV_posdir = new G4LogicalVolume( ptr_cells_final->at(iEta),
								     Material_CAL,
								     LV );
	    for (auto iPhi = 0;  iPhi < NumberOfPixel; ++iPhi) {
		zRot = new G4RotationMatrix;
		zRot->rotateZ( iPhi * d_phi);
		new G4PVPlacement( zRot,																		   //* rotation
				   G4ThreeVector(0, 0, direction * ( ipDistance_z + previous_layers_depths * cos( theta_min_layer ) + 0.5 * delta_z_1 ) ), //* Placed the pixel in specific way, that its slices goes throught (0 ,0, 0) // hack
				   // G4ThreeVector(0, 0, 0 ), //* Placed the pixel in specific way, that its slices goes throught (0 ,0, 0)
				   endCap_LV_posdir,															   //* its logical volume
				   ( std::string( PL ) + "_" + std::to_string( iEta + 1 ) + "_" + std::to_string( iPhi + 1 ) ),                                            //* its name
				   GlobalLV,																	   //* its mother  volume
				   false,																		   //* no boolean operation
				   index_EndCap,															   //* copy number
				   fCheckOverlaps);
		endCap_LV_posdir->SetVisAttributes(VisAtt);
		++index_EndCap;
	    }
	}
}

template <typename T>
T CalorimeterConstruction::GetMinOrMax(const std::vector<std::vector<T > >& array_2d, bool chooseMin) const {
    T xMinMax = array_2d[0][0];
    for ( const std::vector<T > x_outer : array_2d ) {
	for ( T x_inner : x_outer ) {
	    if ( chooseMin && x_inner < xMinMax )
		xMinMax = x_inner;
	    if ( !chooseMin && x_inner > xMinMax )
		xMinMax = x_inner;
	}
    }
    return xMinMax;    
} 

int CalorimeterConstruction::GetNPixelsMax() const {
    int nPixelsMax_ECAL = GetMinOrMax( geometry.number_of_pixels_ECAL, false );
    int nPixelsMax_HCAL = GetMinOrMax( geometry.number_of_pixels_HCAL, false );
    return std::max( nPixelsMax_ECAL, nPixelsMax_HCAL );
}

long double CalorimeterConstruction::GetMinDPhi() const {
    long double minDPhi_ECAL = GetMinOrMax( geometry.layer_dphi_ECAL, true );
    long double minDPhi_HCAL = GetMinOrMax( geometry.layer_dphi_HCAL, true );
    return std::min( minDPhi_ECAL, minDPhi_HCAL );
}

bool CalorimeterConstruction::CheckMergeFactor( int NumberOfPixel, int cellMergeFactor ) const {
    
    if ( ! ( NumberOfPixel % ( 4 * cellMergeFactor ) == 0 ) ) {
	G4cerr << "Error in calorimeter barrel construction.  NumberOfPixel / ( 4 * cellMergeFactor ) = "
	       << NumberOfPixel << " / ( 4 * " << cellMergeFactor
	       << " ) = "
	       << (float)NumberOfPixel / (float)( 4 * cellMergeFactor )
	       << " is not an integer."
	       << G4endl;
	return false;
    }
    return true;
    
}


std::vector<G4VSolid*>* CalorimeterConstruction::MergeCells( std::vector<G4VSolid*>* cells_to_merge,
							     int NumberOfPixel,
							     int cellMergeFactor ) {
    //
    // Cell merging to allow for different granularities depending on the layer.
    //
    std::vector<G4VSolid*>* ptr_cells_final = new std::vector<G4VSolid*>;
    for ( int iCellFinal = 0; iCellFinal < NumberOfPixel / ( 4 * cellMergeFactor ); ++iCellFinal ) {
	G4VSolid* cell_merged = cells_to_merge->at( iCellFinal * cellMergeFactor );
	for ( int iCellSummand = 1; iCellSummand < cellMergeFactor; ++iCellSummand ) {
	    cell_merged = new G4UnionSolid( "Barrel_cell_merged",
						    cell_merged,
					    cells_to_merge->at( iCellFinal * cellMergeFactor + iCellSummand ));
		}
		ptr_cells_final->push_back( cell_merged );
	    }
    return ptr_cells_final;
}
