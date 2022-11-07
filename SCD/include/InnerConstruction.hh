#ifndef __INNERCONSTRUCTION_H__
#define __INNERCONSTRUCTION_H__

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Types.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "Config_reader_var.hh"


class InnerConstruction {

	public:
        InnerConstruction(G4LogicalVolume* expHallLV, G4Material* default_Material, G4Material* Iron, G4Material* ElSi, bool fCheck_Overlaps);
        ~InnerConstruction();
        void Barrel_Inner ();
        void EndCap_Inner();
        void PixelTrk_Barrel(long double r_inn_trkPix, long double r_out_trkPix, long double length_var, G4VisAttributes* VisAtt, bool isOutermostLayer = false);
        void PixelTrk_Barrel(long double r_inn_trkPix, long double r_out_trkPix, G4VisAttributes* VisAtt);
    void PixelTrk_EndCap( long double r_inn_trkPix, long double r_out_trkPix, long double width, long double Abs_Position, G4VisAttributes* VisAtt, int direction, bool isOutermostLayer = false);
        void PixelTrk_EndCap( long double r_inn_trkPix, long double r_out_trkPix, G4VisAttributes* VisAtt, int direction);

    private:
        Config_reader_var config_obj = Config_reader_var::GetInstance();
        G4LogicalVolume* GlobalLV;
        G4LogicalVolume* MagFieldLV;
        long double theta_min;
        long double r_inn;
        bool fCheckOverlaps;
        G4Material* defaultMaterial;
        G4Material* iron;
        G4Material* elSi;
        //*Preporation for visualization
        G4VisAttributes* Pix_VisAtt=
        new G4VisAttributes(true, G4Colour(0, 0, 255));
        G4VisAttributes* Str_VisAtt=
        new G4VisAttributes(true, G4Colour(255, 0, 0));
        G4VisAttributes* Iron_Support_VisAtt=
        new G4VisAttributes(true, G4Colour(0, 255, 0));

        //*Preporation for visualization end
};



#endif // __INNERCONSTRUCTION_H__
