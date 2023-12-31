#include "InnerConstruction.hh"

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
#include "G4Tubs.hh"
#include "G4CutTubs.hh"
#include "G4Cons.hh"

#include "G4SystemOfUnits.hh"
#include "G4IntersectionSolid.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4PVParameterised.hh"
#include "DetectorGeometryDefinitions.hh"
#include "G4ReflectionFactory.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"

#include "G4RunManager.hh"

#include <string>

InnerConstruction::InnerConstruction(G4LogicalVolume* expHallLV, G4Material* default_Material, G4Material* Iron, G4Material* ElSi, bool fCheck_Overlaps) 
{

    //
    // All inner detector logical volume names start with "Inner".
    // This property will be used elsewhere.
    //

    theta_min = 2*atan(exp(-1*config_obj.max_eta_barrel));
    fCheckOverlaps = fCheck_Overlaps;
    Iron_Support_VisAtt->SetForceSolid(true);
    defaultMaterial = default_Material;
    r_inn = config_obj.r_inn_calo;
    iron = Iron;
    elSi = ElSi;



    //* Inner created to have magnetic field only in the inner part of detector
	G4Tubs *Inner= new G4Tubs("Inner_det", 0, r_inn, r_inn/tan(theta_min), 0, config_obj.max_phi);
	GlobalLV = new G4LogicalVolume(Inner, defaultMaterial, "Inner_LV");
	new G4PVPlacement(
					0,                // no rotation
					G4ThreeVector(0., 0., 0. ),
					GlobalLV,          // its logical volume
					"Inner_PL",          // its name
					expHallLV,                // its mother  volume
					false,            // no boolean operation
					0  ,              // copy number
					fCheckOverlaps
	);

	float       SpaceForLayers        = r_inn - ( r_out_trkStr3 + widthiron_add );
	float       GapBetweenIronLayers  = SpaceForLayers/(NumberOfIronLayers+1);
	long double r_magField            = (r_out_trkStr3 + widthiron_add) + GapBetweenIronLayers;
	long double iron_r_inn_ec         = ( r_out_trkStr3 + widthiron_add ) + NumberOfIronLayers * GapBetweenIronLayers;
	long double iron_width_ec         = 4.4 * cm;
	long double l_magField            = ( ( iron_r_inn_ec +0.5 * iron_width_ec ) / tan( theta_min ) ); // gap offset
	l_magField                       -= 0.5 * fabs( iron_width_ec / ( tan( theta_min ) ) );            // gap half-width

	G4Tubs *magFieldTub = new G4Tubs( "Magnetic Field Tube",
					  0.0,
					  r_magField,
					  l_magField,
					  0.0,
					  config_obj.max_phi );
	MagFieldLV = new G4LogicalVolume( magFieldTub,
					  defaultMaterial,
					  "MagField_LV" );
	new G4PVPlacement( 0,
			   G4ThreeVector(0., 0., 0. ),
			   MagFieldLV,
			   "Inner_MagField_PL",
			   GlobalLV,
			   false,
			   0,
			   fCheckOverlaps );
	
	  Barrel_Inner();
	  EndCap_Inner();
	

	if (config_obj.fieldValue!=0)
	{
		G4UniformMagField* magField = new G4UniformMagField(G4ThreeVector(0.,0.,config_obj.fieldValue));
		G4FieldManager* fieldMgr= new G4FieldManager(magField);
		MagFieldLV->SetFieldManager(fieldMgr,true);
		fieldMgr->SetDetectorField(magField);
		fieldMgr->CreateChordFinder(magField);
		fieldMgr->GetChordFinder()->SetDeltaChord(0.01);
	}

}
InnerConstruction::~InnerConstruction()
{;}

void InnerConstruction::Barrel_Inner()
{
	if ( config_obj.use_inner_detector ) {
	
	    PixelTrk_Barrel(r_inn_trkPix0,r_out_trkPix0,280*mm,Pix_VisAtt);
	    PixelTrk_Barrel(r_inn_trkPix1,r_out_trkPix1,280*mm,Pix_VisAtt);
	    PixelTrk_Barrel(r_inn_trkPix2,r_out_trkPix2,280*mm,Pix_VisAtt);
	    PixelTrk_Barrel(r_inn_trkPix3,r_out_trkPix3,280*mm,Pix_VisAtt);
	    PixelTrk_Barrel(r_inn_trkPix4,r_out_trkPix4,280*mm,Pix_VisAtt);
	    
	    
	    PixelTrk_Barrel(r_inn_trkStr0,r_out_trkStr0,1150*mm,Str_VisAtt);
	    PixelTrk_Barrel(r_inn_trkStr1,r_out_trkStr1,1150*mm,Str_VisAtt);
	    PixelTrk_Barrel(r_inn_trkStr2,r_out_trkStr2,1150*mm,Str_VisAtt);
	    PixelTrk_Barrel(r_inn_trkStr3,r_out_trkStr3,1150*mm,Str_VisAtt, true);

	}

	if ( !config_obj.use_ID_support )
	    return;
	    
    // *Barrel support material 
	float SpaceForLayers = r_inn-(r_out_trkStr3+widthiron_add);
	float GapBetweenIronLayers = SpaceForLayers/(NumberOfIronLayers+1);
	long double iron_width_barrel = 4.4*cm/NumberOfIronLayers;
	G4LogicalVolume *Iron_Layer_LV;
	G4Tubs *Iron_Layer;
	for (int NironLayer = 1; NironLayer < (NumberOfIronLayers+1); NironLayer++)
	{

		long double iron_r_inn = (r_out_trkStr3 + widthiron_add) + NironLayer*GapBetweenIronLayers;
		
		long double l_Iron_Layer =  (iron_r_inn)/tan(theta_min);

		if (NironLayer!=NumberOfIronLayers)
		{
			l_Iron_Layer =  (iron_r_inn+0.9*iron_width_barrel)/tan(theta_min);
			Iron_Layer = new G4Tubs("Inner_Iron_layer", iron_r_inn, iron_r_inn + iron_width_barrel, l_Iron_Layer, 0, config_obj.max_phi);
			Iron_Layer_LV = new G4LogicalVolume(Iron_Layer, iron, "Iron_gap_barrel_ID_LV");
			new G4PVPlacement(
							0,                // no rotation
							G4ThreeVector(0., 0., 0. ),
							Iron_Layer_LV,          // its logical volume
							( std::string( "Iron_gap_barrel_" ) + std::to_string( NironLayer ) ).c_str(),          // its name
							GlobalLV,                // its mother  volume
							false,            // no boolean operation
							0,               // copy number
							fCheckOverlaps
						);
		}
		else 
		{
			Iron_Layer = new G4Tubs("Inner_Iron_layer", iron_r_inn, iron_r_inn + iron_width_barrel, l_Iron_Layer, 0, config_obj.max_phi);
			Iron_Layer_LV = new G4LogicalVolume(Iron_Layer, iron, "Iron_gap_barrel_ID_LV");
			new G4PVPlacement(
								0,                // no rotation
								G4ThreeVector(0., 0., 0. ),
								Iron_Layer_LV,          // its logical volume
								( std::string( "Iron_gap_barrel_" ) + std::to_string( NironLayer ) ).c_str(),          // its name
								GlobalLV,                // its mother  volume
								false,            // no boolean operation
								0,               // copy number
								fCheckOverlaps
							);
		}
		
		Iron_Layer_LV->SetVisAttributes(Iron_Support_VisAtt);
		
	}
}

void InnerConstruction::EndCap_Inner()
{

    if ( config_obj.use_inner_detector ) {
	
	PixelTrk_EndCap(r_inn_trkPix0,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix0,Pix_VisAtt, 1);
	PixelTrk_EndCap(r_inn_trkPix0,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix0,Pix_VisAtt,-1);
	
	PixelTrk_EndCap(r_inn_trkPix0,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix1,Pix_VisAtt, 1);
	PixelTrk_EndCap(r_inn_trkPix0,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix1,Pix_VisAtt,-1);
	
	PixelTrk_EndCap(r_inn_trkPix0,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix2,Pix_VisAtt, 1);
	PixelTrk_EndCap(r_inn_trkPix0,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix2,Pix_VisAtt,-1);
	
	PixelTrk_EndCap(r_inn_trkPix0,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix3,Pix_VisAtt, 1);
	PixelTrk_EndCap(r_inn_trkPix0,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix3,Pix_VisAtt,-1);
	
	PixelTrk_EndCap(r_inn_trkPix0,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix4,Pix_VisAtt, 1);
	PixelTrk_EndCap(r_inn_trkPix0,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix4,Pix_VisAtt,-1);
	
	PixelTrk_EndCap(r_inn_trkPix0,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix5,Pix_VisAtt, 1);
	PixelTrk_EndCap(r_inn_trkPix0,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix5,Pix_VisAtt,-1);
	
	PixelTrk_EndCap(r_inn_trkPix0,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix6,Pix_VisAtt, 1);
	PixelTrk_EndCap(r_inn_trkPix0,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix6,Pix_VisAtt,-1);

	PixelTrk_EndCap(r_inn_trkPix2,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix7,Pix_VisAtt, 1);
	PixelTrk_EndCap(r_inn_trkPix2,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix7,Pix_VisAtt,-1);
	
	PixelTrk_EndCap(r_inn_trkPix2,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix8,Pix_VisAtt, 1);
	PixelTrk_EndCap(r_inn_trkPix2,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix8,Pix_VisAtt,-1);

	PixelTrk_EndCap(r_inn_trkPix2,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix9,Pix_VisAtt, 1);
	PixelTrk_EndCap(r_inn_trkPix2,r_out_trkPix4,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkPix9,Pix_VisAtt,-1);
	
	
	PixelTrk_EndCap(r_inn_trkStr0,r_out_trkStr3,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkStr0,Str_VisAtt, 1);
	PixelTrk_EndCap(r_inn_trkStr0,r_out_trkStr3,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkStr0,Str_VisAtt,-1);
	
	PixelTrk_EndCap(r_inn_trkStr0,r_out_trkStr3,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkStr1,Str_VisAtt, 1);
	PixelTrk_EndCap(r_inn_trkStr0,r_out_trkStr3,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkStr1,Str_VisAtt,-1);
	
	PixelTrk_EndCap(r_inn_trkStr0,r_out_trkStr3,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkStr2,Str_VisAtt, 1);
	PixelTrk_EndCap(r_inn_trkStr0,r_out_trkStr3,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkStr2,Str_VisAtt,-1);

	PixelTrk_EndCap(r_inn_trkStr0,r_out_trkStr3,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkStr4,Str_VisAtt, 1, true);
	PixelTrk_EndCap(r_inn_trkStr0,r_out_trkStr3,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkStr4,Str_VisAtt,-1, true);
	
	PixelTrk_EndCap(r_inn_trkStr0,r_out_trkStr3,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkStr4,Str_VisAtt, 1);
	PixelTrk_EndCap(r_inn_trkStr0,r_out_trkStr3,r_out_trkPix0-r_inn_trkPix0,pos_EndCap_trkStr4,Str_VisAtt,-1);

    }
    
	if ( !config_obj.use_ID_support ) return;

    
	//* End-Cap support material 
	float SpaceForLayers = r_inn-(r_out_trkStr3+widthiron_add);
	float GapBetweenIronLayers = SpaceForLayers/(NumberOfIronLayers+1);
	long double iron_r_inn = (r_out_trkStr3 + widthiron_add) + NumberOfIronLayers*GapBetweenIronLayers;
	long double iron_width_endcap = 4.4*cm;
	long double iron_width_barrel = 4.4*cm/NumberOfIronLayers;

	for (int direction = 1; direction > -2; direction=direction-2)
	{
		long double length_cone_min = iron_r_inn/tan(theta_min);
		long double length_cone_max = length_cone_min+iron_width_endcap/tan(theta_min);
		long double d_theta_next = 2*atan(exp(-1*(config_obj.max_eta_endcap)));
		if (direction==-1)
		{
			long double buf = length_cone_min;
			length_cone_min = length_cone_max;
			length_cone_max = buf;
		}

		G4Cons *Cone_Trk = new G4Cons("Cone_Gap", (length_cone_max)*tan(d_theta_next), (iron_r_inn + iron_width_barrel) , (length_cone_min)*tan(d_theta_next), (iron_r_inn + iron_width_barrel), fabs(iron_width_endcap/(2*tan(theta_min))),  0, 2*M_PI);

		G4LogicalVolume *Iron_LV_posdir = new G4LogicalVolume(Cone_Trk, iron, "Iron_ID_cone_endcap");

		new G4PVPlacement(
						0,               //* rotation
						G4ThreeVector(0, 0,  -1* direction*((iron_r_inn+0.5*iron_width_endcap)/tan(theta_min))), //* Placed the pixel in specific way, that its slices goes throught (0 ,0, 0)
						Iron_LV_posdir,          //* its logical volume
						"Iron_ID_gap_endcap",          //* its name
						GlobalLV,                //* its mother  volume
						false,            //* no boolean operation
						0,             //* copy number
						fCheckOverlaps
							);
		Iron_LV_posdir->SetVisAttributes(Iron_Support_VisAtt);
	}



}


void InnerConstruction::PixelTrk_EndCap(long double r_inn_trkPix, long double r_out_trkPix, G4VisAttributes* VisAtt, int direction) 
{
	
	
	long double length_cone_min = (r_inn_trkPix)/tan(theta_min);
	long double length_cone_max = length_cone_min+(r_out_trkPix-r_inn_trkPix)/tan(theta_min);
	if (direction == -1)
	{
		long double buf = length_cone_min;
		length_cone_min = length_cone_max;
		length_cone_max = buf;
	}
	long double d_theta = 2*atan(exp(-1*config_obj.max_eta_barrel));
	long double d_theta_next = 2*atan(exp(-1*(config_obj.max_eta_endcap)));


	G4Cons *Cone_Trk = new G4Cons("Cone_Gap", (length_cone_max)*tan(d_theta_next), length_cone_max*tan(d_theta),(length_cone_min)*tan(d_theta_next),length_cone_min*tan(d_theta), fabs((length_cone_max-length_cone_min)/2),  0, 2*M_PI);

	G4LogicalVolume *Trk_LV_posdir = new G4LogicalVolume(Cone_Trk, elSi, "Inner_Cone_Gap");

	new G4PVPlacement(
					0,               //* rotation
					G4ThreeVector(0, 0,  -1* direction*(length_cone_max+length_cone_min)/2), //* Placed the pixel in specific way, that its slices goes throught (0 ,0, 0)
					Trk_LV_posdir,          //* its logical volume
					"Iron_PL",          //* its name
					MagFieldLV,                //* its mother  volume
					false,            //* no boolean operation
					0,             //* copy number
					fCheckOverlaps
		);


	length_cone_min = (r_out_trkPix)/tan(theta_min);
	length_cone_max = length_cone_min+(widthiron_add)/tan(theta_min);
	if (direction == -1)
	{
		long double buf = length_cone_min;
		length_cone_min = length_cone_max;
		length_cone_max = buf;
	}
	d_theta = 2*atan(exp(-1*config_obj.max_eta_barrel));
	d_theta_next = 2*atan(exp(-1*config_obj.max_eta_endcap));

	G4Cons *Cone_Trk_iron = new G4Cons("Cone_Gap", (length_cone_max)*tan(d_theta_next), length_cone_max*tan(d_theta),(length_cone_min)*tan(d_theta_next),length_cone_min*tan(d_theta), fabs((length_cone_max-length_cone_min)/2),  0, 2*M_PI);
	G4LogicalVolume *Inner_trkPix_LV_iron = new G4LogicalVolume(Cone_Trk_iron, iron, "Inner_trkPix_LV_iron");
	new G4PVPlacement(
					0,               //* rotation
					G4ThreeVector(0, 0,  -1* direction*(length_cone_max+length_cone_min)/2), //* Placed the pixel in specific way, that its slices goes throught (0 ,0, 0)
					Inner_trkPix_LV_iron,          //* its logical volume
					"Iron_PL",          //* its name
					MagFieldLV,                //* its mother  volume
					false,            //* no boolean operation
					0,             //* copy number
					fCheckOverlaps
		);

	Trk_LV_posdir->SetVisAttributes(VisAtt);
	VisAtt->SetForceSolid(true);
	Inner_trkPix_LV_iron->SetVisAttributes(VisAtt);
}
void InnerConstruction::PixelTrk_EndCap(long double r_inn_trkPix, long double r_out_trkPix, long double width, long double Abs_Position, G4VisAttributes* VisAtt, int direction, bool isOutermostLayer) 
{


	G4Tubs *Cone_Trk = new G4Tubs("Cone_Gap", r_inn_trkPix, r_out_trkPix, width/2,  0, 2*M_PI);

	G4LogicalVolume *Trk_LV_posdir = new G4LogicalVolume(Cone_Trk, elSi, "Inner_Cone_Gap");

	std::string detName = "Inner_si_endcap_PL";
	std::string supName = "Inner_Iron_PL";
	if ( isOutermostLayer ) {
	    detName += "_outermostInner";
	    supName += "_outermostInner";
	}
	
	new G4PVPlacement(
					0,               //* rotation
					G4ThreeVector(0, 0,  -1* direction*Abs_Position), //* Placed the pixel in specific way, that its slices goes throught (0 ,0, 0)
					Trk_LV_posdir,          //* its logical volume
					detName,          //* its name
					MagFieldLV,                //* its mother  volume
					false,            //* no boolean operation
					0,             //* copy number
					fCheckOverlaps
		);


	G4Tubs *Cone_Trk_iron = new G4Tubs("Cone_Gap", r_inn_trkPix, r_out_trkPix, widthiron_add/2,  0, 2*M_PI);
	G4LogicalVolume *Inner_trkPix_LV_iron = new G4LogicalVolume(Cone_Trk_iron, iron, "Inner_trkPix_LV_iron");
	new G4PVPlacement(
					0,               //* rotation
					G4ThreeVector(0, 0,  -1* direction*(Abs_Position+widthiron_add)), //* Placed the pixel in specific way, that its slices goes throught (0 ,0, 0)
					Inner_trkPix_LV_iron,          //* its logical volume
					supName,          //* its name
					MagFieldLV,                //* its mother  volume
					false,            //* no boolean operation
					0,             //* copy number
					fCheckOverlaps
		);

	Trk_LV_posdir->SetVisAttributes(VisAtt);
	VisAtt->SetForceSolid(true);
	Inner_trkPix_LV_iron->SetVisAttributes(VisAtt);
}

void InnerConstruction::PixelTrk_Barrel(long double r_inn_trkPix, long double r_out_trkPix, G4VisAttributes* VisAtt) 
{
	
	long double l_Pix =  r_out_trkPix/tan(theta_min);
	long double l_Pix_iron =  (r_out_trkPix+widthiron_add)/tan(theta_min);
	

  	G4Tubs *Inner_trkPix = new G4Tubs("Inner_trkPix",r_inn_trkPix,r_out_trkPix, l_Pix, 0,config_obj.max_phi);
	G4LogicalVolume *Inner_trkPix_LV = new G4LogicalVolume(Inner_trkPix, elSi, "Inner_trkPix_LV");
	new G4PVPlacement(
					0,                // no rotation
					G4ThreeVector(0., 0., 0. ),
					Inner_trkPix_LV,          // its logical volume
					"Inner_trkPix_PL",          // its name
					MagFieldLV,                // its mother  volume
					false,            // no boolean operation
					0 ,               // copy number
					fCheckOverlaps
	);
	G4Tubs *Inner_trkPix_iron = new G4Tubs("Inner_trkPix_iron", r_out_trkPix, r_out_trkPix+widthiron_add, l_Pix_iron,0,config_obj.max_phi);
	G4LogicalVolume *Inner_trkPix_LV_iron = new G4LogicalVolume(Inner_trkPix_iron, iron, "Inner_trkPix_LV_iron");
	new G4PVPlacement(
					0,                // no rotation
					G4ThreeVector(0., 0., 0. ),
					Inner_trkPix_LV_iron,          // its logical volume
					"Inner_trkPix_PL_iron",          // its name
					MagFieldLV,                // its mother  volume
					false,            // no boolean operation
					0 ,              // copy number
					fCheckOverlaps
	);
	Inner_trkPix_LV->SetVisAttributes(VisAtt);
	Inner_trkPix_LV_iron->SetVisAttributes(VisAtt);
}

void InnerConstruction::PixelTrk_Barrel(long double r_inn_trkPix, long double r_out_trkPix, long double length_var, G4VisAttributes* VisAtt, bool isOutermostLayer) 
{

	std::string detName = "Inner_trkPix_PL";
	std::string supName = "Inner_trkPix_PL_iron";
	if ( isOutermostLayer ) {
	    detName += "_outermostInner";
	    supName += "_outermostInner";
	}
	
  	G4Tubs *Inner_trkPix = new G4Tubs("Inner_trkPix",r_inn_trkPix,r_out_trkPix, length_var, 0,config_obj.max_phi);
	G4LogicalVolume *Inner_trkPix_LV = new G4LogicalVolume(Inner_trkPix, elSi, "Inner_trkPix_LV");
	new G4PVPlacement(
					0,                // no rotation
					G4ThreeVector(0., 0., 0. ),
					Inner_trkPix_LV,          // its logical volume
					detName,          // its name
					MagFieldLV,                // its mother  volume
					false,            // no boolean operation
					0 ,               // copy number
					fCheckOverlaps
	);
	G4Tubs *Inner_trkPix_iron = new G4Tubs("Inner_trkPix_iron", r_out_trkPix, r_out_trkPix+widthiron_add, length_var,0,config_obj.max_phi);
	G4LogicalVolume *Inner_trkPix_LV_iron = new G4LogicalVolume(Inner_trkPix_iron, iron, "Inner_trkPix_LV_iron");
	new G4PVPlacement(
					0,                // no rotation
					G4ThreeVector(0., 0., 0. ),
					Inner_trkPix_LV_iron,          // its logical volume
					supName,          // its name
					MagFieldLV,               // its mother  volume
					false,            // no boolean operation
					0 ,              // copy number
					fCheckOverlaps
	);
	Inner_trkPix_LV->SetVisAttributes(VisAtt);
	VisAtt->SetForceSolid(true);
	Inner_trkPix_LV_iron->SetVisAttributes(VisAtt);
	
}
