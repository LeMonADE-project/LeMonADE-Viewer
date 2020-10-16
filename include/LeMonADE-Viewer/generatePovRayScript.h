#pragma once

#include <cmath>                            // atan, fmod, floor
#include <fstream>
#include <string>

#include <LeMonADE-Viewer/Camera.h>

inline float mod( float const & a, float const & N ){ return a - N*floor(a/N); }  // why not use fmod ?

template < class IngredientsType >
inline void generatePovRayScript
(
    IngredientsType const & ingredients,
    std::string     const & fname,
    unsigned int    const   width,
    unsigned int    const   height,
    Camera          const & cam
)
{
    unsigned int const boxX = ingredients.getBoxX();
    unsigned int const boxY = ingredients.getBoxY();
    unsigned int const boxZ = ingredients.getBoxZ();

    std::ofstream afile( ( fname + ".pov" ).c_str(), std::ios::out | std::ios::trunc );

    /* Save command in comments and create povray header */
    afile
    << "// povray +I" << fname << ".pov +O" << fname << ".png +W" << width <<" +H" << height << " +P\n"
    << "// +UA for transparency mask\n"
    << "\n"
    << "\n"
    << "#version 3.5;\n"
    << "#include \"colors.inc\"\n"
    << "\n"
    << "global_settings {assumed_gamma 1.0 noise_generator 2}\n"
    << "\n"
    << "camera{\n"
    << " location <0,0," << ( -cam.getCamYPos() >= 0 ? -cam.getCamYPos() : 0 ) << ">\n"
    << " right -x*" << width << "/" << height << "\n"
    << " look_at <0,0," << ( -cam.getCamYPos() >= 0 ? 0 : -cam.getCamYPos() ) << ">\n"
    << " translate <" << cam.getCamXPos() << ", " << -cam.getCamZPos() << ",0>\n";
    double const pi = std::atan(1)*4;
    double const fovy = 45.0; // degrees - fov is in y in openGL, povRay uses in x
    double anglefovx = 2 * atan ( ((1.0*width)/height) * tan(fovy/2.0/180*pi) )* 180.0/pi;
    afile
    << " angle " << anglefovx << "\n"
    << " rotate <90,0,0>\n"
    << " rotate <0,0," << -cam.getCamAngleYaw() << ">\n"
    //<< "rotate <" << -cam.getCamThetaRot() << ",0,0>\n"
    << " rotate <" << cam.getCamAnglePitch() << ",0,0>\n"
    << " rotate <0,0," << cam.getCamAngleRoll() << ">\n"
    << " translate <" << 0.5 * boxX  << "," << 0.5 * boxY << "," << 0.5 * boxZ << ">\n"
    << "}\n"
    << "\n"
    << "// background color\n"
    << "sky_sphere{pigment{color rgb<" << ingredients.getBGcolor()[0] << "," << ingredients.getBGcolor()[1] << "," << ingredients.getBGcolor()[2] << "> }}\n"
    << "\n"
    << "// light source\n"
    << "light_source{<" << cam.getCamXPos() << "," << -cam.getCamZPos() << "," << -cam.getCamYPos() << "> color White\n"
    //<< "fade_distance " << -(cam.getCamYPos()+0.5*boxY) << "\n"
    << "//fade_distance " << 100 << "\n"
    << "//fade_power 3\n"
    << "rotate <90,0,0>\n"
    << "rotate <0,0," << -cam.getCamAngleYaw() << ">\n"
    //afile << "rotate  <" << -cam.getCamThetaRot() << ",0,0>\n"
    << "rotate <" << cam.getCamAnglePitch() << ",0,0>\n"
    << "translate <" << 0.5 * boxX << "," << 0.5 * boxY << "," << 0.5 * boxZ << ">\n"
    << "}\n"
    << "\n"
    // output for simulation box
    << "#macro BOXFRAME (P1, P2, S1, C1)" << std::endl
    << "cylinder{P1,P2, S1" << std::endl
    << "         texture{ pigment{color C1} finish {ambient 0.15 diffuse 0.85  metallic } }" << std::endl
    << "        }" << std::endl
    << "#end" << std::endl
    << std::endl
    << "// simulation box frame" << std::endl
    << "object{ BOXFRAME (<0,0,0>,<" << boxX << ",0,0>, 1, <0.8,0.8,0.8,0.7>) }" << std::endl
    << "object{ BOXFRAME (<0,0,0>,<0," << boxY << ",0>, 1, <0.8,0.8,0.8,0.7>) }" << std::endl
    << "object{ BOXFRAME (<0,0,0>,<0,0," << boxZ << ">, 1, <0.8,0.8,0.8,0.7>) }" << std::endl
    
    << "object{ BOXFRAME (<" << boxX << ",0,0>,<" << boxX << ",0," << boxZ << ">, 1, <0.8,0.8,0.8,0.7>) }" << std::endl
    << "object{ BOXFRAME (<" << boxX << ",0,0>,<" << boxX << "," << boxY << ",0>, 1, <0.8,0.8,0.8,0.7>) }" << std::endl
    
    << "object{ BOXFRAME (<0,0," << boxZ << ">,<" << boxX << ",0," << boxZ << ">, 1, <0.8,0.8,0.8,0.7>) }" << std::endl
    << "object{ BOXFRAME (<0,0," << boxZ << ">,<0," << boxY << "," << boxZ << ">, 1, <0.8,0.8,0.8,0.7>) }" << std::endl
    
    << "object{ BOXFRAME (<" << boxX << ",0," << boxZ << ">,<" << boxX << "," << boxY << "," << boxZ << ">, 1, <0.8,0.8,0.8,0.7>) }" << std::endl
    
    << "object{ BOXFRAME (<0," << boxY << ",0>,<" << boxX << "," << boxY << ",0>, 1, <0.8,0.8,0.8,0.7>) }" << std::endl
    << "object{ BOXFRAME (<0," << boxY << ",0>,<0," << boxY << "," << boxZ << ">, 1, <0.8,0.8,0.8,0.7>) }" << std::endl
    << "object{ BOXFRAME (<0," << boxY << "," << boxZ << ">,<" << boxX << "," << boxY << "," << boxZ << ">, 1, <0.8,0.8,0.8,0.7>) }" << std::endl
    << "object{ BOXFRAME (<" << boxX << "," << boxY << ",0>,<" << boxX << "," << boxY << "," << boxZ << ">, 1, <0.8,0.8,0.8,0.7>) }" << std::endl
    
    << std::endl
    << std::endl
    << "\n"
    << "// objects (monomers)\n";

    // das ist jetzt die reale monomer-Position (STIMMT DAS SCHON SO?)
    typename IngredientsType::molecules_type const & molecules = ingredients.getMolecules();
    for ( size_t n = 0; n < molecules.size(); ++n )
    {
        if ( molecules[n].isVisible() )
        {
            if ( ! ingredients.isSmoothing() )
            {
                if ( ingredients.isVisualizePBC() )
                {
                    afile << "sphere{< "
                    << ( ( ( molecules[n][0] + ingredients.getTranslationInX() ) % boxX ) + boxX ) % boxX << ", "
                    << ( ( ( molecules[n][1] + ingredients.getTranslationInY() ) % boxY ) + boxY ) % boxY << ", "
                    << ( ( ( molecules[n][2] + ingredients.getTranslationInZ() ) % boxZ ) + boxZ ) % boxZ << "> "
                    << molecules[n].getRadius() << "\n";
                }
                else
                {
                    //afile << "sphere{< " << linearStrandsVector[n].operator[](m)[0] - 0.5*boxX << ", " << linearStrandsVector[n].operator[](m)[1] - 0.5*boxX << ", " << -linearStrandsVector[n].operator[](m)[2] + 0.5 * boxZ << "> 2" << std::endl;
                    afile << "sphere{< "
                    << molecules[n][0] + ingredients.getTranslationInX() << ", "
                    << molecules[n][1] + ingredients.getTranslationInY() << ", "
                    << molecules[n][2] + ingredients.getTranslationInZ() << "> "
                    << molecules[n].getRadius() << "\n";
                }
            }
            else
            {
                if ( ingredients.isVisualizePBC() )
                {
                    afile << "sphere{< "
                    << mod( molecules[n].getSmoothCoordinate()[0] + ingredients.getTranslationInX(), boxX ) << ", "
                    << mod( molecules[n].getSmoothCoordinate()[1] + ingredients.getTranslationInY(), boxY ) << ", "
                    << mod( molecules[n].getSmoothCoordinate()[2] + ingredients.getTranslationInZ(), boxZ ) << "> "
                    << molecules[n].getRadius() << "\n";
                }
                else
                {
                    afile << "sphere{< "
                    << (molecules[n].getSmoothCoordinate()[0] + ingredients.getTranslationInX())  << ", "
                    << (molecules[n].getSmoothCoordinate()[1] + ingredients.getTranslationInY())  << ", "
                    << (molecules[n].getSmoothCoordinate()[2] + ingredients.getTranslationInZ())  << "> "
                    << (molecules[n].getRadius()) << "\n";
                }
            }
            afile << " pigment {color rgb<"
            << molecules[n].getColor().getX() << ", "
            << molecules[n].getColor().getY() << ", "
            << molecules[n].getColor().getZ() <<  ">}\n"
            << " finish {phong 5.0 phong_size 40 reflection rgb<0.0, 0.0, 0.0> roughness 0.05 ambient 0.1 diffuse 0.9} }\n";
        }
    }
    afile << std::endl;

    //bonds:

    if ( ingredients.isVisualizeBonds() )
    for ( size_t n = 0; n < molecules.size(); ++n )
    {
        if ( ! molecules[n].isVisible() )
            continue;
        size_t const linkcount = ingredients.getMolecules().getNumLinks(n);
        for ( size_t k = 0; k < linkcount; ++k )
        {
            // search for the monomers connected to the current
            size_t neighboridx = ingredients.getMolecules().getNeighborIdx( n, k );
            // draw only one bond, even though each bond is saved in both bond partners
            if ( ! ( neighboridx < n ) )
                continue;
            if ( ! ingredients.isSmoothing() )
            {
                /*
                afile << "sphere{< "
                <<  linearStrandsVector[n].operator[](m)[0] - 0.5*boxX << ", "
                <<  linearStrandsVector[n].operator[](m)[1] - 0.5*boxX << ", "
                << -linearStrandsVector[n].operator[](m)[2] + 0.5*boxZ << "> 2\n";
                */
                /*
                afile << "sphere{< "
                <<  linearStrandsVectorGroup[n].operator[](m)[0] << ", "
                <<  linearStrandsVectorGroup[n].operator[](m)[1] << ", "
                << -linearStrandsVectorGroup[n].operator[](m)[2] << "> "
                << (point_size/10.0f) << "\n";
                */
                VectorFloat3 bondDirection;
                bondDirection.setX(1.0f*(molecules[n][0]-ingredients.getMolecules()[neighboridx][0]));
                bondDirection.setY(1.0f*(molecules[n][1]-ingredients.getMolecules()[neighboridx][1]));
                bondDirection.setZ(1.0f*(molecules[n][2]-ingredients.getMolecules()[neighboridx][2]));

                float bondLength = bondDirection.getLength();

                if ( 0.1 < bondLength && bondLength < 6.0 )
                {
                    // bondDirection.normalize();
                    /*
                    afile << "cylinder{ < "
                    <<  linearStrandsVectorGroup[n].operator[](m)[0] << ", "
                    <<  linearStrandsVectorGroup[n].operator[](m)[1] << ", "
                    << -linearStrandsVectorGroup[n].operator[](m)[2] << ">"
                    */
                    afile
                    << "cylinder{ <0,0,0> "
                    /*
                    << linearStrandsVectorGroup[n].operator[](m)[0]  << ", "
                    << linearStrandsVectorGroup[n].operator[](m)[1]  << ", "
                    << (-linearStrandsVectorGroup[n].operator[](m)[2])  << ">"
                    */
                    << ", <"
                    << -bondDirection[0] << ","
                    << -bondDirection[1] << ","
                    << -bondDirection[2] << ">, "
                    << ingredients.getWidthBond()
                    /*
                    <<  ingredients.getMolecules()[neighboridx][0] << ","
                    <<  ingredients.getMolecules()[neighboridx][1] << ","
                    << -ingredients.getMolecules()[neighboridx][2] << "> "
                    */
                    << " texture{ pigment{ gradient <"
                    << bondDirection[0] / bondDirection.getLength() << ","
                    << bondDirection[1] / bondDirection.getLength() << ","
                    << bondDirection[2] / bondDirection.getLength() << ">"
                    << " color_map{ [0.0 color rgb<"
                    << ingredients.getMolecules()[neighboridx].getColor().getX() << ", "
                    << ingredients.getMolecules()[neighboridx].getColor().getY() << ", "
                    << ingredients.getMolecules()[neighboridx].getColor().getZ() <<  "> ]"
                    << "[1.0 color rgb<"
                    << molecules[n].getColor().getX() << ", "
                    << molecules[n].getColor().getY() << ", "
                    << molecules[n].getColor().getZ() <<  "> ]} scale " << bondLength
                    << " translate <"
                    << -bondDirection[0] << ","
                    << -bondDirection[1] <<","
                    << -bondDirection[2] <<">} finish { phong 1}} translate < ";
                    /*
                    <<  ingredients.getMolecules()[neighboridx][0] << ", "
                    <<  ingredients.getMolecules()[neighboridx][1] << ", "
                    << -ingredients.getMolecules()[neighboridx][2] << "> "
                    */
                    if ( ! ingredients.isVisualizePBC() )
                    {
                        afile
                        << molecules[n][0]+ingredients.getTranslationInX() << ", "
                        << molecules[n][1]+ingredients.getTranslationInY() << ", "
                        << molecules[n][2]+ingredients.getTranslationInZ() << ">";
                    }
                    else
                    {
                        afile
                        << ( ( ( molecules[n][0]+ingredients.getTranslationInX() ) % boxX) + boxX ) % boxX << ", "
                        << ( ( ( molecules[n][1]+ingredients.getTranslationInY() ) % boxY) + boxY ) % boxY << ", "
                        << ( ( ( molecules[n][2]+ingredients.getTranslationInZ() ) % boxZ) + boxZ ) % boxZ << ">";
                    }
                    afile << "scale <1,1,1> rotate<0,0,0> translate<0,0,0> }\n";
                }
            }
            else // smoothing true
            {
                VectorFloat3 bondDirection;
                bondDirection.setX(1.0f*(molecules[n].getSmoothCoordinate()[0]-ingredients.getMolecules()[neighboridx].getSmoothCoordinate()[0]));
                bondDirection.setY(1.0f*(molecules[n].getSmoothCoordinate()[1]-ingredients.getMolecules()[neighboridx].getSmoothCoordinate()[1]));
                bondDirection.setZ(1.0f*(molecules[n].getSmoothCoordinate()[2]-ingredients.getMolecules()[neighboridx].getSmoothCoordinate()[2]));

                float bondLength= bondDirection.getLength();

                if ( 0.1 < bondLength && bondLength < 6.0 )
                {
                    afile
                    /*
                    << "cylinder{ < "
                    <<  linearStrandsVectorGroup[n].operator[](m)[0] << ", "
                    <<  linearStrandsVectorGroup[n].operator[](m)[1] << ", "
                    << -linearStrandsVectorGroup[n].operator[](m)[2] << ">"
                    */
                    << "cylinder{ <0,0,0> "
                    /*
                    <<  linearStrandsVectorGroup[n].operator[](m)[0] << ", "
                    <<  linearStrandsVectorGroup[n].operator[](m)[1] << ", "
                    << -linearStrandsVectorGroup[n].operator[](m)[2] << ">"
                    */
                    ", <"
                    << -bondDirection[0] << ","
                    << -bondDirection[1] << ","
                    << -bondDirection[2] << ">, "
                    << ingredients.getWidthBond() <<
                    /*
                    <<  ingredients.getMolecules()[neighboridx][0] << ", "
                    <<  ingredients.getMolecules()[neighboridx][1] << ", "
                    << -ingredients.getMolecules()[neighboridx][2] << "> "
                    */
                    "  texture{ pigment{ gradient <"
                    << bondDirection[0]/bondDirection.getLength() <<","
                    << bondDirection[1]/bondDirection.getLength() <<","
                    << (bondDirection[2]/bondDirection.getLength()) <<
                    "> color_map{ [0.0 color rgb<"
                    << ingredients.getMolecules()[neighboridx].getColor().getX() << ", "
                    << ingredients.getMolecules()[neighboridx].getColor().getY() << ", "
                    << ingredients.getMolecules()[neighboridx].getColor().getZ() <<  "> ]"
                    << "[1.0 color rgb<"
                    << molecules[n].getColor().getX() << ", "
                    << molecules[n].getColor().getY() << ", "
                    << molecules[n].getColor().getZ() <<  ">"
                    << "]} scale " << bondLength <<
                    " translate <"
                    << (-bondDirection[0]) <<","
                    << -(bondDirection[1]) <<","
                    << -(bondDirection[2]) <<">}"
                    << " finish { phong 1}} translate < ";

                    if ( ! ingredients.isVisualizePBC() )
                    {
                        afile
                        << molecules[n].getSmoothCoordinate()[0] + ingredients.getTranslationInX() << ", "
                        << molecules[n].getSmoothCoordinate()[1] + ingredients.getTranslationInY() << ", "
                        << molecules[n].getSmoothCoordinate()[2] + ingredients.getTranslationInZ() << ">";
                    }
                    else
                    {
                        afile
                        << mod((mod(molecules[n].getSmoothCoordinate()[0] + ingredients.getTranslationInX(), boxX) )+boxX, boxX ) << ", "
                        << mod((mod(molecules[n].getSmoothCoordinate()[1] + ingredients.getTranslationInY(), boxY) )+boxY, boxY ) << ", "
                        << mod((mod(molecules[n].getSmoothCoordinate()[2] + ingredients.getTranslationInZ(), boxZ) )+boxZ, boxZ ) << ">";
                    }
                    afile << "scale <1,1,1> rotate<0,0,0> translate<0,0,0> }\n";
                }
            }
        }
    }
    afile << std::endl;
    // BOND VECTORS EINBINDEN!!!
    afile.close();
}
