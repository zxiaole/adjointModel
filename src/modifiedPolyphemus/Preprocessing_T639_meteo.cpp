/*
 * C Implementation: keys_iterator
 *
 * Description:
 * Example on how to use keys_iterator functions and the
 * grib_keys_iterator structure to get all the available
 * keys in a message.
 *
 */


#include <cmath>
#include <iostream>
#include <algorithm>
using namespace std;

#define SELDONDATA_DEBUG_LEVEL_4
#define SELDONDATA_WITH_GRIB

#include "AtmoData.hxx"
using namespace AtmoData;

#include "Common.cxx"
using namespace Polyphemus;

// INCLUDES //
//////////////



template <class T>
void exp_(T& x)
{
    x = exp(x);
}

template <class T>
void sqrt_(T& x)
{
    x = sqrt(x);
}

template <class T>
void abs_(T& x)
{
    x = abs(x);
}


#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>

#include "grib_api.h"

#define MAX_KEY_LEN  255
#define MAX_VAL_LEN  1024

template<class T, int N>
void ReadField(grib_handle* h, Data<T, N>& A, int level, int Nx_in, int Ny_in, int Nz_in)
{
    size_t vlen=MAX_VAL_LEN;
    double* values = NULL;


    int n_tmp = Nz_in;
    int xx, yy, x_id, longitude;
    Grid<T>* gridx;
    Grid<T>* gridy;
    gridx = A.GetGrid(2);
    gridy = A.GetGrid(1);


    GRIB_CHECK(grib_get_size(h,"values",&vlen),0);
    values = (double*)malloc(vlen*sizeof(double));
    /* get data values */
    GRIB_CHECK(grib_get_double_array(h,"values",values,&vlen),0);

    for(yy=0; yy<Ny_in; yy++)
        for(xx=0; xx<Nx_in;xx++)
        {
            longitude = (*gridx)(xx);
            if(longitude<360)
            {
                x_id = longitude/0.5;

            }
            else
            {
                x_id = (longitude-360)/0.5;
            }
            A(level,yy,xx) = values[x_id + yy*720 ];
        }
    cout<<level<<endl;
    free(values);
}

template<class T, int N>
void ReadField(grib_handle* h, Data<T, N>& A, int Nx_in, int Ny_in)
{
    size_t vlen=MAX_VAL_LEN;
    double* values = NULL;

    int xx, yy, x_id, longitude;
    Grid<T>* gridx;
    Grid<T>* gridy;
    gridx = A.GetGrid(1);
    gridy = A.GetGrid(0);


    GRIB_CHECK(grib_get_size(h,"values",&vlen),0);
    values = (double*)malloc(vlen*sizeof(double));
    /* get data values */
    GRIB_CHECK(grib_get_double_array(h,"values",values,&vlen),0);

    for(yy=0; yy<Ny_in; yy++)
        for(xx=0; xx<Nx_in;xx++)
        {
            longitude = (*gridx)(xx);
            if(longitude<360)
            {
                x_id = longitude/0.5;

            }
            else
            {
                x_id = (longitude-360)/0.5;
            }
            A(yy,xx) = values[x_id + yy*720];
        }
    free(values);
}

template<class T, int N>
void ReadField(grib_handle* h, Data<T, N>& A, int Nx_in, int Ny_in, bool tmp)
{
    size_t vlen=MAX_VAL_LEN;
    double* values = NULL;

    int xx, yy;

    GRIB_CHECK(grib_get_size(h,"values",&vlen),0);
    values = (double*)malloc(vlen*sizeof(double));
    /* get data values */
    GRIB_CHECK(grib_get_double_array(h,"values",values,&vlen),0);

    for(yy=0; yy<Ny_in; yy++)
        for(xx=0; xx<Nx_in;xx++)
        {
            A(yy,xx) = values[xx+ yy*Nx_in];
        }
    free(values);
}





static void usage(char* progname);

int main(int argc, char *argv[])
{
    TRY;

    cout << "test!!!!!!!!!!"<<endl;

    string configuration_file, sec_config_file;


    //    parse_argument(argc, argv, configuration_file, sec_config_file,
    //           date_beg, default_name);
    configuration_file = argv[2];
    sec_config_file = argv[3];



    // read  begin and end date
    Date date_beg;
    date_beg = argv[4];

    Date date_end = date_beg;
    date_end = argv[5];

    Date date_prev = date_beg;
    date_prev.AddDays(-1);

    Date date_analysis, date_current;
    date_current = date_beg;
    date_analysis = date_beg;

    int t_pred = 0;

    ////////////////////////
    // FIRST DECLARATIONS //
    ////////////////////////

    typedef float real;

    int i, j, k;
    int xx, yy, zz, tt;
    int gridx_min_id, gridx_max_id;

    // Constants.
    const real pi = 3.14159265358979323846264;
    const real P0 = 101325.;
    const real r = 287.;
    const real cp = 1005.;
    const real kappa = 0.4; // von Karman constant

    //// PARAMETERIZATION. ////

    // Accumulated rain.
    bool prev_accumulated_rain;
    // Attenuation parameterization.
    int atte_type;
    // Kz thresholds.
    real Kz_min, Kz_min_urban;
    bool apply_vert;
    real Kz_max;
    // Clouds.
    real min_height;



    /////////////////////////
    // CONFIGURATION FILES //
    /////////////////////////

    cout << "Reading configuration files..."; cout.flush();

    ConfigStreams configuration(configuration_file);
    if (exists(sec_config_file))
        configuration.AddFile(sec_config_file);

    // Input domain.
    int Nt_in, Nz_in, Ny_in, Nx_in, Nx_in_original, Ny_in_original;
    real Delta_t_in, Delta_y_in, Delta_x_in;
    real t_min_in, y_min_in, x_min_in;

    configuration.SetSection("[METEO]");

    configuration.PeekValue("Nt", "> 0", Nt_in);
    configuration.PeekValue("Nz", "> 0", Nz_in);
    configuration.PeekValue("Ny", "> 0", Ny_in);
    configuration.PeekValue("Nx", "> 0", Nx_in);
    Nx_in_original = Nx_in;
    Ny_in_original = Ny_in;

    configuration.PeekValue("Delta_t", "> 0", Delta_t_in);
    configuration.PeekValue("Delta_y", "> 0", Delta_y_in);
    configuration.PeekValue("Delta_x", "> 0", Delta_x_in);

    configuration.PeekValue("t_min", t_min_in);
    configuration.PeekValue("y_min", y_min_in);
    configuration.PeekValue("x_min", x_min_in);

    // Output domain.
    int Nt_out, Nz_out, Ny_out, Nx_out;
    real Delta_t_out, Delta_y_out, Delta_x_out;
    real t_min_out, y_min_out, x_min_out;
    string vertical_levels;

    configuration.SetSection("[domain]");

    configuration.PeekValue("Nx", "> 0", Nx_out);
    configuration.PeekValue("Ny", "> 0", Ny_out);
    configuration.PeekValue("Nz", "> 0", Nz_out);

    configuration.PeekValue("Delta_t", "> 0", Delta_t_out);
    configuration.PeekValue("Delta_y", "> 0", Delta_y_out);
    configuration.PeekValue("Delta_x", "> 0", Delta_x_out);

    configuration.PeekValue("y_min", y_min_out);
    configuration.PeekValue("x_min", x_min_out);
    configuration.PeekValue("Vertical_levels", vertical_levels);

    Nt_out =1;//compute_Nt(date_beg, date_end, Delta_t_out);
    t_min_out = real(date_beg.GetHour()) + real(date_beg.GetMinutes()) / 60.
            +  real(date_beg.GetSeconds()) / 3600.;

    // Files.
    string directory_in, directory_out, roughness_in_file;
    string InputFileDirect_in, name_suffix;
    int grib_version;

    configuration.SetSection("[paths]");

    configuration.PeekValue("Database_meteo", InputFileDirect_in);

    configuration.PeekValue("Directory_meteo", directory_out);

    configuration.PeekValue("Roughness_file", roughness_in_file);

    configuration.PeekValue("Version", grib_version);

    if(grib_version == 1)
        name_suffix = ".grb1";
    else if(grib_version == 1)
        name_suffix = ".grb2";
    else
         cout<<"The version "<< grib_version<<"of the input meteorological data  is invalid"<<endl;

    // Land use.
    string LUC_file;
    int Nc, urban_index;
    configuration.PeekValue("LUC_file", LUC_file);
    if (!exists(LUC_file))
        throw "Unable to open land use cover file \"" + LUC_file + "\".";
    Nc = int(file_size(LUC_file)) / sizeof(float) / (Ny_out * Nx_out);
    configuration.PeekValue("Urban_index", ">= 0 | < " + to_str(Nc),
                            urban_index);

    string file_in = directory_in + date_beg.GetDate("ECMWF-%y%m%d.grb");
    string file_in_prev = directory_in + date_prev.GetDate("ECMWF-%y%m%d.grb");

    configuration.SetSection("[clouds]");
    configuration.PeekValue("Min_height", "positive", min_height);

    configuration.SetSection("[Kz]");
    configuration.PeekValue("Min", "positive", Kz_min);
    configuration.PeekValue("Min_urban", "positive", Kz_min_urban);
    configuration.PeekValue("Max", "positive", Kz_max);
    configuration.PeekValue("Apply_vert", apply_vert);



    // Files.
    bool richardson_with_roughness;

    configuration.SetSection("[METEO]");

    configuration.PeekValue("Richardson_with_roughness",
                            richardson_with_roughness);

//    // Accumulated data
//    int accumulated_time, accumulated_index;

//    configuration.SetSection("[accumulated_data]");

//    configuration.PeekValue("Accumulated_time", "positive", accumulated_time);
//    configuration.PeekValue("Accumulated_index", "positive", accumulated_index);

//    cout << " done." << endl;

    ///////////
    // GRIDS //
    ///////////

    cout << "Memory allocation for data fields..."; cout.flush();

    // Input settings.


    // Vertical levels are shared.


    // Output settings.

    // Output grids.
    RegularGrid<real> GridT_out(t_min_out, Delta_t_out, Nt_out);
    RegularGrid<real> GridZ_out(Nz_out);

    RegularGrid<real> GridY_out(y_min_out, Delta_y_out, Ny_out);
    RegularGrid<real> GridX_out(x_min_out, Delta_x_out, Nx_out);

    gridx_min_id = int(floor(GridX_out(0)/0.5));
    gridx_max_id = int(ceil(GridX_out(Nx_out-1)/0.5));
    Nx_in = gridx_max_id - gridx_min_id + 1;
    x_min_in = 0.5*gridx_min_id;


    // Data may be provided on interfaces.
    RegularGrid<real> GridZ_interf_out(Nz_out + 1);
    RegularGrid<real> GridY_interf_out(y_min_out - Delta_y_out / 2.,
                                       Delta_y_out, Ny_out + 1);
    RegularGrid<real> GridX_interf_out(x_min_out - Delta_x_out / 2.,
                                       Delta_x_out, Nx_out + 1);

    // Reads output altitudes.
    FormatText Heights_out;
    Heights_out.Read(vertical_levels, GridZ_interf_out);
    // Sets values at nodes.
    for (k=0; k<Nz_out; k++)
        GridZ_out(k) = ( GridZ_interf_out(k) + GridZ_interf_out(k+1) ) / 2.0;
    GridZ_out.Print();


    // Input grids.
    RegularGrid<real> GridT_in(t_min_in, Delta_t_in, Nt_in);
    // Vertical levels depend on t, z, y and x.
    // GeneralGrid<real, 3> GridZ_in(shape( Nz_in, Ny_in, Nx_in),
    //             0, shape(0, 1, 2));
    RegularGrid<real> GridC(Nc);
    RegularGrid<real> GridZ_in(Nz_in);
    double level_in[17] =  {110.901,762.08,1457.51,3012.61, 4207.02, 5575.21, 7186.41,
            9165.16,   10364.3,   11776.4,   13510.6 ,  15799.1 ,  17671.4,   19324.8 ,
            21641.1 ,  23325.9 ,  25921.4};
    RegularGrid<real> GridZ_interfZ_in(Nz_in+1);
    RegularGrid<real> GridY_in(y_min_in, Delta_y_in, Ny_in);
    RegularGrid<real> GridX_in(x_min_in, Delta_x_in, Nx_in);
    GridX_in.Print();

    for(i=0;i<Nz_in;i++)
        GridZ_in(i) = level_in[i];
    cout<<"GridX_in"<<endl;

    GridZ_interfZ_in(0) = 0;

    for (zz = 1; zz < Nz_in; zz++)
        GridZ_interfZ_in(zz) =
                (GridZ_in(zz-1) + GridZ_in(zz))/2 ;
    GridZ_interfZ_in(Nz_in) = 2*GridZ_in(Nz_in-1) - GridZ_interfZ_in(Nz_in-1);

    GridZ_interfZ_in.Print();



    //////////
    // DATA //
    //////////

    // Input fields.
    Data<real, 3> ZonalWind(GridZ_in, GridY_in, GridX_in);
    Data<real, 3> A(GridZ_in, GridY_in, GridX_in);
    Data<real, 3> Pressure(GridZ_in, GridY_in, GridX_in);
    Data<real, 3> MeridionalWind(GridZ_in, GridY_in, GridX_in);
    Data<real, 3> Temperature( GridZ_in, GridY_in, GridX_in);
    Data<real, 3> RelativeHumidity(GridZ_in, GridY_in, GridX_in);
    Data<real, 2> SurfaceTemperature(GridY_in, GridX_in);
    Data<real, 2> SurfaceDewTemperature(GridY_in, GridX_in);
    Data<real, 2> TemperatureLapse(GridY_in, GridX_in);
    Data<real, 2> B(GridY_in, GridX_in);
    Data<real, 2> BoundaryHeight(GridY_in, GridX_in);


    Data<real, 2> SensibleHeat(GridY_in, GridX_in);
    Data<real, 2> SensibleHeat_prev(GridY_in, GridX_in);
    Data<real, 2> SolarRadiation(GridY_in, GridX_in);
    Data<real, 2> SolarRadiation_whole(Ny_in_original, Nx_in_original);
    Data<real, 2> SolarRadiation_whole_prev(Ny_in_original, Nx_in_original);
    Data<real, 2> SolarRadiation_prev(GridY_in, GridX_in);
    Data<real, 2> ZonalWind_10m(GridY_in, GridX_in);
    Data<real, 2> MeridionalWind_10m(GridY_in, GridX_in);
    Data<real, 2> SurfacePressure(GridY_in, GridX_in);
    Data<real, 2> ConvectiveRain(GridY_in, GridX_in);
    Data<real, 2> ConvectiveRain_pre(GridY_in, GridX_in);
    Data<real, 2> LargeScaleRain(GridY_in, GridX_in);
    Data<real, 2> LargeScaleRain_pre(GridY_in, GridX_in);
    Data<real, 2> Rain(GridY_in, GridX_in);

    Data<real, 3> CRH(GridZ_in, GridY_in, GridX_in);
    Data<real, 3> CloudFraction(GridZ_in, GridY_in, GridX_in);
    Data<int, 3> LowIndices(Ny_in, Nx_in, 2);
    Data<int, 3> MediumIndices(Ny_in, Nx_in, 2);
    Data<int, 3> HighIndices(Ny_in, Nx_in, 2);

    Data<real, 2> LowCloudiness( GridY_in, GridX_in);
    Data<real, 2> MediumCloudiness( GridY_in, GridX_in);
    Data<real, 2> HighCloudiness( GridY_in, GridX_in);
    Data<real, 2> CloudHeight(GridY_in, GridX_in);


    Data<real, 2> Landfraction(GridY_in, GridX_in);

    Data<real, 2> U_star(GridY_in, GridX_in);
    Data<real, 2> V_star(GridY_in, GridX_in);
    Data<real, 2> U_star_prev(GridY_in, GridX_in);
    Data<real, 2> V_star_prev(GridY_in, GridX_in);

    // Reads land use data.
    Data<real, 3> LUC(GridC, GridY_out, GridX_out);
    FormatBinary<float>().Read(LUC_file, LUC);



    cout<<"out"<<endl;
    // Output fields.
    Data<real, 3> WindModule_out(GridZ_out, GridY_out, GridX_out);
    Data<real, 3> ZonalWind_out(GridZ_out,
                                GridY_out, GridX_interf_out);

    Data<real, 3> MeridionalWind_out(GridZ_out,
                                     GridY_interf_out, GridX_out);
    Data<real, 3> Temperature_out( GridZ_out, GridY_out, GridX_out);
    Data<real, 3> SpecificHumidity_out( GridZ_out, GridY_out, GridX_out);
    Data<real, 3> PotentialTemperature_out( GridZ_out, GridY_out, GridX_out);
    Data<real, 3> Pressure_out( GridZ_out, GridY_out, GridX_out);
    Data<real, 2> SurfaceTemperature_out(GridY_out, GridX_out);
    Data<real, 2> SurfaceDewTemperature_out(GridY_out, GridX_out);
    Data<real, 2> TemperatureLapse_out(GridY_out, GridX_out);
    Data<real, 2> BoundaryHeight_out(GridY_out, GridX_out);
    Data<real, 2> Landfraction_out(GridY_out, GridX_out);
    Data<real, 2> SensibleHeat_out(GridY_out, GridX_out);
    Data<real, 2> SolarRadiation_out(GridY_out, GridX_out);
    Data<real, 2> FirstLevelWindModule_out(GridY_out, GridX_out);
    Data<real, 2> ZonalWind_10m_out(GridY_out, GridX_out);
    Data<real, 2> MeridionalWind_10m_out(GridY_out, GridX_out);
    Data<real, 2> Stability_out(GridY_out, GridX_out);
    Data<real, 2> SurfacePressure_out(GridY_out, GridX_out);
    Data<real, 2> ConvectiveRain_out(GridY_out, GridX_out);
    Data<real, 2> CloudHeight_out(GridY_out, GridX_out);
    Data<real, 2> Rain_out(GridY_out, GridX_out);
    Data<real, 2> U_star_out(GridY_out, GridX_out);
    Data<real, 2> V_star_out(GridY_out, GridX_out);
    Data<real, 2> Ustar(GridY_out, GridX_out);
    Data<real, 2> FrictionModule_out(GridY_out, GridX_out);
    Data<real, 3> RelativeHumidity_out( GridZ_out, GridY_out, GridX_out);
    Data<real, 3> Richardson_out(GridZ_out, GridY_out, GridX_out);


    Data<real, 3> Kz_out(GridZ_interf_out, GridY_out, GridX_out);
    Data<real, 2> Roughness(GridY_out, GridX_out);

    double tmp_SolarRadiation[Nx_in_original*Ny_in_original];
    cout << " done." << endl;
    cout << endl;

    FormatBinary<float> InputBinaryMeteo;
    InputBinaryMeteo.Read(roughness_in_file, Roughness);


    unsigned long key_iterator_filter_flags=GRIB_KEYS_ITERATOR_ALL_KEYS;

    /* Choose a namespace. E.g. "ls", "time", "parameter", "geography", "statistics" */
    char* name_space="ls";

    /* name_space=NULL to get all the keys */
    /* char* name_space=0; */


    int n_tmp = Nz_in;



    FILE* f = NULL;
    FILE* f_prev = NULL;
    grib_handle* h=NULL;
    grib_index* index=NULL;
    grib_index* pre_index=NULL;
    char** shortName = NULL;
    char** typeOfLevel = NULL;
    long* level;
    double* values = NULL;

    int err=0;
    int grib_count=0;
    int ret = 0, count = 0;

    char value[MAX_VAL_LEN];
    size_t vlen=MAX_VAL_LEN;
    size_t shortNameSize, levelSize, typeOfLevelSize;
    const char* filename;
    const char* pre_filename;



    while(date_current<=date_end)
    {
        cout<<"date_current:"<<date_current<<endl;
        date_analysis = date_current;

        if(date_current.GetHour()==0 )
        {

            date_analysis.AddDays(-1);
            date_analysis.SetHour(12);
            t_pred = 12;
        }
        else if(date_current.GetHour()<=12)
        {
            date_analysis.SetHour(0);
            t_pred = date_current.GetHour() - date_analysis.GetHour();
            cout<<"t_pred:"<<t_pred<<endl;
        }
        else
        {
            date_analysis.SetHour(12);
            t_pred = date_current.GetHour() - date_analysis.GetHour();
        }
        cout<<"date_analysis:"<<date_analysis<<endl;




        string InPutFileName("gmf.639.");
        string Pre_InPutFileName;
        string InputFileDirect;
        InputFileDirect = InputFileDirect_in;

        InputFileDirect += to_str_fill(date_analysis.GetYear(), 4, '0', ostringstream::right);
        InputFileDirect += to_str_fill(date_analysis.GetMonth(), 2, '0', ostringstream::right);
        InputFileDirect += to_str_fill(date_analysis.GetDay(), 2, '0', ostringstream::right);
        InputFileDirect += to_str_fill(date_analysis.GetHour(), 2, '0', ostringstream::right);
        InputFileDirect += "/";

        int tmpt = 0;
        InPutFileName += to_str_fill(date_analysis.GetYear(), 4, '0', ostringstream::right);
        InPutFileName += to_str_fill(date_analysis.GetMonth(), 2, '0', ostringstream::right);
        InPutFileName += to_str_fill(date_analysis.GetDay(), 2, '0', ostringstream::right);
        InPutFileName += to_str_fill(date_analysis.GetHour(), 2, '0', ostringstream::right);

        Pre_InPutFileName = InPutFileName;
        InPutFileName += to_str_fill(t_pred, 3, '0', ostringstream::right);
        InPutFileName +=  name_suffix;
        Pre_InPutFileName += to_str_fill(t_pred-3, 3, '0', ostringstream::right);
        Pre_InPutFileName +=  name_suffix;
        InPutFileName = InputFileDirect + InPutFileName;
        Pre_InPutFileName = InputFileDirect + Pre_InPutFileName;
        cout<<"filename:"<<InPutFileName<<endl;

        date_current.AddHours(3);






        //if (argc != 2) usage(argv[0]);
        filename = InPutFileName.c_str();
        pre_filename = Pre_InPutFileName.c_str();

        //file_handle
        f = fopen(filename,"r");
        f_prev = fopen(pre_filename,"r");


        //! indexing
        printf("indexing...\n");
        index=grib_index_new(0, "shortName,level,typeOfLevel", &ret);
        pre_index=grib_index_new(0, "shortName,level,typeOfLevel", &ret);

        if (ret) {printf("error: %s\n",grib_get_error_message(ret)); exit(ret);}
        /* Indexes a file */
        ret=grib_index_add_file(index,filename);
        if (ret) {printf("error: %s\n",grib_get_error_message(ret)); exit(ret);}
        ret=grib_index_add_file(pre_index,pre_filename);
        if (ret) {printf("error: %s\n",grib_get_error_message(ret)); exit(ret);}
        printf("end indexing...\n");


        /*same as for "step"*/
        GRIB_CHECK(grib_index_get_size(index,"level",&levelSize),0);
        level=(long*)malloc(sizeof(long)*levelSize);
        if (!level)
            exit(1);
        /*same as for "step"*/
        GRIB_CHECK(grib_index_get_long(index,"level",level,&levelSize),0);
        printf("levelSize=%ld\n",(long)levelSize);
//        for (i=2;i<levelSize;i++)
//        {
//            GridZ_in(Nz_in-i+1) = (1 - pow((level[i]/1013.25),(1/5.255)))*44330;
//            printf("%ld ",level[i]);
//        }
//        GridZ_in.Print();


        cout << endl;
        printf("df\n");

        /* get the number of distinct values of "shortName" in the index */
        GRIB_CHECK(grib_index_get_size(index,"shortName",&shortNameSize),0);
        shortName=(char**)malloc(sizeof(char*)*shortNameSize);
        if (!shortName) exit(1);

        GRIB_CHECK(grib_index_get_string(index,"shortName",shortName,&shortNameSize),0);
        printf("shortNameSize=%ld\n",(long)shortNameSize);
        for (i=0;i<shortNameSize;i++)
            printf("%s ",shortName[i]);
        printf("\n");


        /* get the number of distinct values of "shortName" in the index */
        GRIB_CHECK(grib_index_get_size(index,"typeOfLevel",&typeOfLevelSize),0);
        typeOfLevel=(char**)malloc(sizeof(char*)*typeOfLevelSize);
        if (!shortName) exit(1);

        GRIB_CHECK(grib_index_get_string(index,"typeOfLevel",typeOfLevel,&typeOfLevelSize),0);
        printf("typeOfLevelSize=%ld\n",(long)typeOfLevelSize);
        for (i=0;i<typeOfLevelSize;i++)
            printf("%s ",typeOfLevel[i]);
        printf("\n");




        /////////////////
        // READS INPUT //
        /////////////////

        for (i=2;i<levelSize;i++)
        {
            for(yy=0; yy<Ny_in; yy++)
                for(xx=0; xx<Nx_in;xx++)
                {
                    // pressure unit
                    Pressure(Nz_in-i+1, yy, xx) = level[i]*100;
                    //printf("%f\n",values[xx + yy*Nx_in + zz*Ny_in*Nx_in]);
                }
        }

        n_tmp = Nz_in;
        grib_index_select_string(index,"shortName","u");
        grib_index_select_string(index,"typeOfLevel","isobaricInhPa");
        for (i=0; i<levelSize; i++)
        {
            GRIB_CHECK(grib_index_select_long(index,"level",level[i]), 0);
            h=grib_handle_new_from_index(index,&ret);
            if(ret)
            {
                printf("Error: %s\n", grib_get_error_message(ret));
            }
            else
            {

                ReadField(h, ZonalWind,n_tmp-1, Nx_in, Ny_in, Nz_in);
                n_tmp--;
                cout<<n_tmp<<endl;

            }
            grib_handle_delete(h);
            printf("loop:%d\n", i);
        }

        n_tmp = Nz_in;
        grib_index_select_string(index,"shortName","v");
        grib_index_select_string(index,"typeOfLevel","isobaricInhPa");
        for (i=0; i<levelSize; i++)
        {
            GRIB_CHECK(grib_index_select_long(index,"level",level[i]), 0);
            h=grib_handle_new_from_index(index,&ret);
            if(ret)
            {
                printf("Error: %s\n", grib_get_error_message(ret));
            }
            else
            {
                ReadField(h, MeridionalWind,n_tmp-1, Nx_in, Ny_in, Nz_in);
                n_tmp--;
                cout<<n_tmp<<endl;
            }
            grib_handle_delete(h);
            printf("loop:%d\n", i);
        }

        n_tmp = Nz_in;
        grib_index_select_string(index,"shortName","r");
        grib_index_select_string(index,"typeOfLevel","isobaricInhPa");
        for (i=0; i<levelSize; i++)
        {
            GRIB_CHECK(grib_index_select_long(index,"level",level[i]), 0);
            h=grib_handle_new_from_index(index,&ret);
            if(ret)
            {
                printf("Error_relativehumility: %s\n", grib_get_error_message(ret));
            }
            else
            {

                ReadField(h, RelativeHumidity,n_tmp-1, Nx_in, Ny_in, Nz_in);
                for(yy=0; yy<Ny_in; yy++)
                    for(xx=0; xx<Nx_in;xx++)
                    {
                        RelativeHumidity(n_tmp-1,yy,xx) = RelativeHumidity(n_tmp-1,yy,xx)/100.;
                    }
                n_tmp--;

            }
            grib_handle_delete(h);
        }
        RelativeHumidity.ThresholdMax(1.);

        n_tmp = Nz_in;
        grib_index_select_string(index,"shortName","t");
        grib_index_select_string(index,"typeOfLevel","isobaricInhPa");
        for (i=0; i<levelSize; i++)
        {
            GRIB_CHECK(grib_index_select_long(index,"level",level[i]), 0);
            h=grib_handle_new_from_index(index,&ret);
            if(ret)
            {
                printf("Error: %s\n", grib_get_error_message(ret));
            }
            else
            {
                ReadField(h, Temperature,n_tmp-1, Nx_in, Ny_in, Nz_in);
                n_tmp--;
                cout<<n_tmp<<endl;
            }
            grib_handle_delete(h);
            printf("loop:%d\n", i);
        }

        //! temperature lapse calculation K/100m

        double level_distance = (GridZ_in(1) - GridZ_in(0))/100;

        for(yy=0; yy<Ny_in; yy++)
            for(xx=0; xx<Nx_in; xx++)
            {
                TemperatureLapse(yy, xx) = (Temperature(1,yy,xx) - Temperature(0,yy,xx))/level_distance;
            }



        n_tmp = Nz_in;
        grib_index_select_string(index,"shortName","2t");
        grib_index_select_string(index,"typeOfLevel","heightAboveGround");

        GRIB_CHECK(grib_index_select_long(index,"level",2), 0);
        h=grib_handle_new_from_index(index,&ret);
        if(ret)
        {
            printf("Error_2t: %s\n", grib_get_error_message(ret));
        }
        else
        {
            ReadField(h,SurfaceTemperature,Nx_in, Ny_in);
            grib_handle_delete(h);
        }


        grib_index_select_string(index,"shortName","lsm");
        grib_index_select_string(index,"typeOfLevel","surface");

        GRIB_CHECK(grib_index_select_long(index,"level",0), 0);
        h=grib_handle_new_from_index(index,&ret);
        if(ret)
        {
            printf("Error_lsm: %s\n", grib_get_error_message(ret));
        }
        else
        {
            ReadField(h,Landfraction,Nx_in, Ny_in);
            grib_handle_delete(h);
        }


        //! sensible heat flux unit: w/m^2 positive at daytime
        grib_index_select_string(index,"shortName","sshf");
        grib_index_select_string(index,"typeOfLevel","surface");

        grib_index_select_string(pre_index,"shortName","sshf");
        grib_index_select_string(pre_index,"typeOfLevel","surface");

        GRIB_CHECK(grib_index_select_long(index,"level",0), 0);
        GRIB_CHECK(grib_index_select_long(pre_index,"level",0), 0);

        h=grib_handle_new_from_index(pre_index,&ret);
        if(ret)
        {
            printf("Error_sshf: %s\n", grib_get_error_message(ret));
        }
        else
        {
            ReadField(h,SensibleHeat_prev,Nx_in, Ny_in);
            grib_handle_delete(h);
        }



        h=grib_handle_new_from_index(index,&ret);
        if(ret)
        {
            printf("Error_sshf: %s\n", grib_get_error_message(ret));
        }
        else
        {
            ReadField(h,SensibleHeat,Nx_in, Ny_in);
            for(yy=0; yy<Ny_in; yy++)
                for(xx=0; xx<Nx_in; xx++)
                {
                    SensibleHeat(yy, xx) = (SensibleHeat(yy,xx) - SensibleHeat_prev(yy,xx));
                }

            grib_handle_delete(h);
        }

        //! 10 m wind speed unit: m s**-1
        grib_index_select_string(index,"shortName","10u");
        grib_index_select_string(index,"typeOfLevel","surface");

        GRIB_CHECK(grib_index_select_long(index,"level",0), 0);
        h=grib_handle_new_from_index(index,&ret);
        if(ret)
        {
            printf("Error_sshf: %s\n", grib_get_error_message(ret));
        }
        else
        {
            ReadField(h,ZonalWind_10m,Nx_in, Ny_in);
            grib_handle_delete(h);
        }

        grib_index_select_string(index,"shortName","10v");
        grib_index_select_string(index,"typeOfLevel","surface");

        GRIB_CHECK(grib_index_select_long(index,"level",0), 0);
        h=grib_handle_new_from_index(index,&ret);
        if(ret)
        {
            printf("Error_sshf: %s\n", grib_get_error_message(ret));
        }
        else
        {
            ReadField(h,MeridionalWind_10m,Nx_in, Ny_in);
            grib_handle_delete(h);
        }

        //! 2m dew temperature unit: K
        grib_index_select_string(index,"shortName","dpt");
        grib_index_select_string(index,"typeOfLevel","heightAboveGround");

        GRIB_CHECK(grib_index_select_long(index,"level",2), 0);
        h=grib_handle_new_from_index(index,&ret);
        if(ret)
        {
            printf("Error_dpt: %s\n", grib_get_error_message(ret));
        }
        else
        {
            ReadField(h,SurfaceDewTemperature,Nx_in, Ny_in);
            grib_handle_delete(h);
        }

        //! Surface pressure unit: Pa
        grib_index_select_string(index,"shortName","sp");
        grib_index_select_string(index,"typeOfLevel","surface");

        GRIB_CHECK(grib_index_select_long(index,"level",0), 0);
        h=grib_handle_new_from_index(index,&ret);
        if(ret)
        {
            printf("Error_surface_pressure: %s\n", grib_get_error_message(ret));
        }
        else
        {
            ReadField(h,SurfacePressure,Nx_in, Ny_in);
            grib_handle_delete(h);
            for(yy=0; yy<Ny_in; yy++)
                for(xx=0; xx<Nx_in;xx++)
                {
                    // pressure unit
                    SurfacePressure(yy, xx) = SurfacePressure(yy, xx) *100;
                    //printf("%f\n",values[xx + yy*Nx_in + zz*Ny_in*Nx_in]);
                }
        }

        //! convective rain unit:mm/h
        grib_index_select_string(index,"shortName","acpcp");
        grib_index_select_string(index,"typeOfLevel","surface");

        grib_index_select_string(pre_index,"shortName","acpcp");
        grib_index_select_string(pre_index,"typeOfLevel","surface");

        GRIB_CHECK(grib_index_select_long(index,"level",0), 0);
        GRIB_CHECK(grib_index_select_long(pre_index,"level",0), 0);

        h=grib_handle_new_from_index(pre_index,&ret);
        if(ret)
        {
            printf("Convective Rain: %s\n", grib_get_error_message(ret));
        }
        else
        {
            ReadField(h,ConvectiveRain_pre,Nx_in, Ny_in); // unit: kg m**-2, should be converted to mm/h.
            grib_handle_delete(h);
            for(yy=0; yy<Ny_in; yy++)
                for(xx=0; xx<Nx_in;xx++)
                {
                    // pressure unit
                    ConvectiveRain_pre(yy, xx) = ConvectiveRain_pre(yy, xx) *1000;
                    //printf("%f\n",values[xx + yy*Nx_in + zz*Ny_in*Nx_in]);
                }
        }

        h=grib_handle_new_from_index(index,&ret);
        if(ret)
        {
            printf("Convective Rain: %s\n", grib_get_error_message(ret));
        }
        else
        {
            ReadField(h,ConvectiveRain,Nx_in, Ny_in); // unit: kg m**-2, should be converted to mm/h.
            grib_handle_delete(h);
            for(yy=0; yy<Ny_in; yy++)
                for(xx=0; xx<Nx_in;xx++)
                {
                    // pressure unit
                    ConvectiveRain(yy, xx) = (ConvectiveRain(yy, xx) *1000 - ConvectiveRain_pre(yy, xx))/3.0;
                    //printf("%f\n",values[xx + yy*Nx_in + zz*Ny_in*Nx_in]);
                }
        }

        //! Large scale rain unit:mm/h
        grib_index_select_string(index,"shortName","lsp");
        grib_index_select_string(index,"typeOfLevel","surface");

        grib_index_select_string(pre_index,"shortName","lsp");
        grib_index_select_string(pre_index,"typeOfLevel","surface");

        GRIB_CHECK(grib_index_select_long(index,"level",0), 0);
        GRIB_CHECK(grib_index_select_long(pre_index,"level",0), 0);
        h=grib_handle_new_from_index(pre_index,&ret);
        if(ret)
        {
            printf("Large scale Rain: %s\n", grib_get_error_message(ret));
        }
        else
        {
            ReadField(h,LargeScaleRain_pre,Nx_in, Ny_in); // unit: kg m**-2, should be converted to mm/h.
            grib_handle_delete(h);
            for(yy=0; yy<Ny_in; yy++)
                for(xx=0; xx<Nx_in;xx++)
                {
                    // pressure unit
                    LargeScaleRain_pre(yy, xx) = LargeScaleRain_pre(yy, xx) *1000;
                    //printf("%f\n",values[xx + yy*Nx_in + zz*Ny_in*Nx_in]);
                }
        }


        h=grib_handle_new_from_index(index,&ret);
        if(ret)
        {
            printf("Large scale Rain: %s\n", grib_get_error_message(ret));
        }
        else
        {
            ReadField(h,LargeScaleRain,Nx_in, Ny_in); // unit: kg m**-2, should be converted to mm/h.
            grib_handle_delete(h);
            for(yy=0; yy<Ny_in; yy++)
                for(xx=0; xx<Nx_in;xx++)
                {
                    // pressure unit
                    LargeScaleRain(yy, xx) = (LargeScaleRain(yy, xx) *1000 - LargeScaleRain_pre(yy, xx))/3.0;
                    //printf("%f\n",values[xx + yy*Nx_in + zz*Ny_in*Nx_in]);
                }
        }

        // Total Rain
        for(yy=0; yy<Ny_in; yy++)
            for(xx=0; xx<Nx_in;xx++)
            {
                // Total rain
                Rain(yy, xx) = ConvectiveRain(yy, xx) + LargeScaleRain(yy, xx);
            }





        //! radiation flux unit: J/m^2
        bool tmp=true;
        grib_index_select_string(index,"shortName","grad");
        grib_index_select_string(index,"typeOfLevel","surface");

        grib_index_select_string(pre_index,"shortName","grad");
        grib_index_select_string(pre_index,"typeOfLevel","surface");
        double radiation_median;

        GRIB_CHECK(grib_index_select_long(index,"level",0), 0);
        GRIB_CHECK(grib_index_select_long(pre_index,"level",0), 0);

        h=grib_handle_new_from_index(pre_index,&ret);
        if(ret)
        {
            printf("Error_sshf: %s\n", grib_get_error_message(ret));
        }
        else
        {


            ReadField(h,SolarRadiation_prev,Nx_in, Ny_in);
            ReadField(h,SolarRadiation_whole_prev,Nx_in_original, Ny_in_original, tmp);
            grib_handle_delete(h);

        }
        h=grib_handle_new_from_index(index,&ret);
        if(ret)
        {
            printf("Error_sshf: %s\n", grib_get_error_message(ret));
        }
        else
        {


            ReadField(h,SolarRadiation,Nx_in, Ny_in);
            ReadField(h,SolarRadiation_whole,Nx_in_original, Ny_in_original, tmp);
            grib_handle_delete(h);
            for(yy=0; yy<Ny_in_original; yy++)
                for(xx=0; xx<Nx_in_original;xx++)
                {
                    SolarRadiation_whole(yy, xx) = SolarRadiation_whole(yy, xx) - SolarRadiation_whole_prev(yy, xx);
                    tmp_SolarRadiation[Nx_in_original*yy + xx] = SolarRadiation_whole(yy, xx);
                }

            int size = sizeof(tmp_SolarRadiation)/sizeof(double);
            std::sort(&tmp_SolarRadiation[0], &tmp_SolarRadiation[size]);
            radiation_median = size % 2 ? tmp_SolarRadiation[size / 2] : (tmp_SolarRadiation[size / 2 - 1] + tmp_SolarRadiation[size / 2]) / 2;
            cout<<"SolarRadiation_median"<<radiation_median<<endl;

        }



        ComputeCriticalRelativeHumidity(Pressure, CRH);


        ComputeCloudFraction(BoundaryHeight, RelativeHumidity, CRH, CloudFraction);
        //    SurfacePressure.Resize();



        ComputeCloudiness(CloudFraction, Pressure, GridZ_interfZ_in, LowIndices,
                          MediumIndices, HighIndices, LowCloudiness,
                          MediumCloudiness, HighCloudiness);

        cout<<"first!!"<<endl;

        ComputeCloudHeight(LowIndices, MediumIndices, HighIndices,
                           GridZ_interfZ_in, CloudHeight);
        cout<<"cloudheight"<<min_height<<endl;
        CloudHeight.ThresholdMin(min_height);


        // surface stress No.98 (U), No.99 (v)

        int count_tmp = 1;
        while((h = grib_handle_new_from_file(0,f_prev ,&err)) != NULL)
        {
            if (count_tmp==98)
            {
                 ReadField(h,U_star_prev ,Nx_in, Ny_in);

            }
            else if(count_tmp==99)
            {
                ReadField(h,V_star_prev ,Nx_in, Ny_in);
                grib_handle_delete(h);
                break;
            }
            count_tmp++;
            grib_handle_delete(h);


        }

        count_tmp = 1;
        while((h = grib_handle_new_from_file(0,f,&err)) != NULL)
        {
            if (count_tmp==98)
            {
                 ReadField(h,U_star ,Nx_in, Ny_in);
                 for(yy=0; yy<Ny_in; yy++)
                     for(xx=0; xx<Nx_in;xx++)
                     {
                         U_star(yy, xx) = (U_star(yy, xx) - U_star_prev(yy, xx))/3./3600.;
                     }

            }
            else if(count_tmp==99)
            {
                ReadField(h,V_star ,Nx_in, Ny_in);
                for(yy=0; yy<Ny_in; yy++)
                    for(xx=0; xx<Nx_in;xx++)
                    {
                        V_star(yy, xx) = (V_star(yy, xx) - V_star_prev(yy, xx))/3./3600.;
                    }
                grib_handle_delete(h);
                break;
            }
            grib_handle_delete(h);
            count_tmp++;


        }


        LinearInterpolationRegular(ZonalWind, ZonalWind_out);
        cout<<"end interpolation"<<endl;
        LinearInterpolationRegular(MeridionalWind, MeridionalWind_out);


        LinearInterpolationRegular(Pressure, Pressure_out);
        LinearInterpolationRegular(Temperature, Temperature_out);
        LinearInterpolationRegular(SurfaceTemperature, SurfaceTemperature_out);
        LinearInterpolationRegular(Landfraction, Landfraction_out);
        LinearInterpolationRegular(SensibleHeat, SensibleHeat_out);
        LinearInterpolationRegular(ZonalWind_10m, ZonalWind_10m_out);
        LinearInterpolationRegular(MeridionalWind_10m, MeridionalWind_10m_out);
        LinearInterpolationRegular(TemperatureLapse, TemperatureLapse_out);
        LinearInterpolationRegular(SolarRadiation, SolarRadiation_out);
        LinearInterpolationRegular(SurfaceDewTemperature, SurfaceDewTemperature_out);
        LinearInterpolationRegular(SurfacePressure, SurfacePressure_out);
        LinearInterpolationRegular(ConvectiveRain, ConvectiveRain_out);
        LinearInterpolationRegular(CloudHeight, CloudHeight_out);
        LinearInterpolationRegular(Rain, Rain_out);
        LinearInterpolationRegular(RelativeHumidity, RelativeHumidity_out);
        LinearInterpolationRegular(U_star , U_star_out);
        LinearInterpolationRegular(V_star , V_star_out);
        cout<<"end_all"<<endl;


        //! parameters calculation



        /////////////////
        // WIND MODULE //
        /////////////////

        cout << "Wind module..."; cout.flush();


        U_star_out.Apply(abs_);
        V_star_out.Apply(abs_);

        U_star_out.Apply(sqrt_);
        V_star_out.Apply(sqrt_);

        ComputeModule(U_star_out, V_star_out, FrictionModule_out);

        cout << " done." << endl;

        //

        Rain.ThresholdMin(0.);
        SolarRadiation_out.ThresholdMin(0.);
        // sensible heat
        for(yy=0; yy<Ny_out; yy++)
            for(xx=0; xx<Nx_out;xx++)
            {
                SensibleHeat_out(yy, xx) *= -1*r * SurfaceTemperature_out(yy, xx)
                  / (SurfacePressure_out(yy, xx) * cp)/3/3600;
            }

        for(yy=0; yy<Ny_out; yy++)
            for(xx=0; xx<Nx_out;xx++)
            {
                SolarRadiation_out(yy, xx)= SolarRadiation_out(yy, xx)/3/3600;
            }




        double u10, v10;

        for(yy=0; yy<Ny_out; yy++)
            for(xx=0; xx<Nx_out;xx++)
            {
                u10 = ZonalWind_10m_out(yy,xx);
                v10 = MeridionalWind_10m_out(yy, xx);

                FirstLevelWindModule_out(yy, xx) = pow(u10*u10 + v10*v10, 0.5);
            }
        //! stability calculation
        double speedtmp, lapsetmp;
        for(yy=0; yy<Ny_out; yy++)
            for(xx=0; xx<Nx_out;xx++)
            {
                speedtmp = FirstLevelWindModule_out(yy, xx);
                lapsetmp = TemperatureLapse_out(yy, xx);
                if(speedtmp>=7)
                {
                    Stability_out(yy, xx) = 4;
                }
                else if(speedtmp>=5)
                {
                    if(lapsetmp<=-1.2)
                    {
                        Stability_out(yy, xx) = 3;
                    }
                    else if(lapsetmp<=2)
                    {
                        Stability_out(yy, xx) = 4;
                    }
                    else
                    {
                        Stability_out(yy, xx) = 5;
                    }

                }
                else if(speedtmp>=3)
                {
                    if(lapsetmp<=-1.2)
                    {
                        Stability_out(yy, xx) = 2;
                    }
                    else if(lapsetmp<=-0.9)
                    {
                        Stability_out(yy, xx) = 3;
                    }
                    else if(lapsetmp<=2)
                    {
                        Stability_out(yy, xx) = 4;
                    }
                    else
                    {
                        Stability_out(yy, xx) = 5;
                    }
                }
                else if(speedtmp>=2)
                {
                    if(lapsetmp<=-1.5)
                    {
                        Stability_out(yy, xx) = 1;
                    }
                    else if(lapsetmp<=-1.2)
                    {
                        Stability_out(yy, xx) = 2;
                    }
                    else if(lapsetmp<=-0.9)
                    {
                        Stability_out(yy, xx) = 3;
                    }
                    else if(lapsetmp<=0)
                    {
                        Stability_out(yy, xx) = 4;
                    }
                    else if(lapsetmp<=2)
                    {
                        Stability_out(yy, xx) = 5;
                    }
                    else
                    {
                        Stability_out(yy, xx) = 6;
                    }
                }
                else if(speedtmp>=1)
                {
                    if(lapsetmp<=-1.5)
                    {
                        Stability_out(yy, xx) = 1;
                    }
                    else if(lapsetmp<=-0.9)
                    {
                        Stability_out(yy, xx) = 2;
                    }
                    else if(lapsetmp<=-0.7)
                    {
                        Stability_out(yy, xx) = 3;
                    }
                    else if(lapsetmp<=0)
                    {
                        Stability_out(yy, xx) = 4;
                    }
                    else
                    {
                        Stability_out(yy, xx) = 6;
                    }
                }
                else
                {
                    if(lapsetmp<=-1.2)
                    {
                        Stability_out(yy, xx) = 1;
                    }
                    else if(lapsetmp<=-0.9)
                    {
                        Stability_out(yy, xx) = 2;
                    }
                    else if(lapsetmp<=-0.7)
                    {
                        Stability_out(yy, xx) = 3;
                    }
                    else if(lapsetmp<=0)
                    {
                        Stability_out(yy, xx) = 4;
                    }
                    else
                    {
                        Stability_out(yy, xx) = 6;
                    }
                }
                // stabilit check
                if((Landfraction_out(yy,xx)<0.01) && (Stability_out(yy, xx)<=3))
                {
                    Stability_out(yy, xx) = 4;
                }
                if((SolarRadiation_out(yy, xx)<radiation_median) && (Stability_out(yy, xx)<=3))
                {
                    Stability_out(yy, xx) = 4;
                }
                if(SensibleHeat_out(yy, xx)>0 && Stability_out(yy, xx)>=5)
                {
                    Stability_out(yy, xx) = 4;
                }
                if(SensibleHeat_out(yy, xx)<0 && Stability_out(yy, xx)<=3)
                {
                    Stability_out(yy, xx) = 4;
                }


            }

        //! boundary layer height calculation
        for(yy=0; yy<Ny_out; yy++)
            for(xx=0; xx<Nx_out;xx++)
            {

                switch(int(Stability_out(yy, xx)))
                {
                case 1:
                    BoundaryHeight_out(yy, xx) = 1600;
                    break;
                case 2:
                    BoundaryHeight_out(yy, xx) = 1200;
                    break;
                case 3:
                    BoundaryHeight_out(yy, xx) = 800;
                    break;
                case 4:
                    BoundaryHeight_out(yy, xx) = 600;
                    break;
                case 5:
                    BoundaryHeight_out(yy, xx) = 300;
                    break;
                case 6:
                    BoundaryHeight_out(yy, xx) = 200;
                    break;

                }

                u10 = ZonalWind_10m_out(yy,xx);
                v10 = MeridionalWind_10m_out(yy, xx);

                FirstLevelWindModule_out(yy, xx) = pow(u10*u10 + v10*v10, 0.5);
            }

        // Vertical diffusion (Louis formula, 1979).
        cout << "Computing Kz..."; cout.flush();

        ComputePotentialTemperature(Temperature_out, Pressure_out,
                                    PotentialTemperature_out);

        ComputeLouisKz(ZonalWind_out, MeridionalWind_out,
                       PotentialTemperature_out, Kz_out);

        real Kz_max_loc;
        int imax(0);


        for (j = 0; j < Ny_out; j++)
            for (i = 0; i < Nx_out; i++)
            {
                Kz_max_loc = 0.0;
                for (k = 0; k < Nz_out + 1; k++)
                    if (Kz_out(k, j, i) >= Kz_max_loc)
                    {
                        Kz_max_loc = Kz_out( k, j, i);
                        imax = k;
                    }
                if (ConvectiveRain_out(j, i) > 1./6.)
                    for (k = imax; k < Nz_out + 1; k++)
                        Kz_out(k, j, i) = Kz_max_loc;
                else
                    for (k = imax; k < Nz_out + 1; k++)
                        Kz_out(k, j, i) = (1. - ConvectiveRain_out(j, i) * 6.)
                                * Kz_out(k, j, i)
                                + 6. * ConvectiveRain_out(j, i) * Kz_max_loc;
            }

        // Kz_min.
        real local_min;
        for (j = 0; j < Ny_out; j++)
            for (i = 0; i < Nx_out; i++)
            {
                local_min = Kz_min * (1. - LUC(urban_index, j, i))
                        + Kz_min_urban * LUC(urban_index, j, i);
                if (apply_vert)
                {
                    for (k = 0; k < Nz_out + 1; k++)
                        if (Kz_out(k, j, i) < local_min)
                            Kz_out(k, j, i) = local_min;
                }
                else
                    if (Kz_out(1, j, i) < local_min)
                        Kz_out(1, j, i) = local_min;
            }

        // Kz_max.
        Kz_out.ThresholdMax(Kz_max);

        cout << " done." << endl;


        ComputeRichardson(ZonalWind_out, MeridionalWind_out,
                  PotentialTemperature_out, Richardson_out);





        //    double T_Dew, T_surf, speed_surf, lon_surf, lat_surf, stability_surf;
        //    for(yy=0; yy<Ny_out; yy++)
        //        for(xx=0; xx<Nx_out;xx++)
        //        {
        //            T_Dew = SurfaceDewTemperature_out(yy, xx);
        //            T_surf = SurfaceTemperature_out(yy, xx);
        //            speed_surf = FirstLevelWindModule_out(yy, xx);
        //            stability_surf = Stability_out(yy, xx);
        //            lon_surf = GridX_out(xx);
        //            lat_surf = GridY_out(yy);

        //            BoundaryHeight_out(yy, xx) = 121/6*(6 - stability_surf)*(T_surf - T_Dew)
        //                    + 0.169*stability_surf*(speed_surf + 0.257)/12/log(10/0.1)/2/(2*pi()/24)./max(sin(lat/180*pi()), 0.5);
        //        }



        FormatBinary<float> OutputMeteo;

        cout << "Writing data..."; cout.flush();
        OutputMeteo.Append( ZonalWind_out, directory_out + "ZonalWind.bin");
        OutputMeteo.Append( MeridionalWind_out, directory_out + "MeridionalWind.bin");
        OutputMeteo.Append( Temperature_out, directory_out + "Temperature.bin");
        OutputMeteo.Append( Pressure_out, directory_out + "Pressure.bin");
       // OutputMeteo.Append( SurfaceTemperature_out, directory_out + "SurfaceTemperature.bin");
       // OutputMeteo.Append( SurfaceTemperature_out, directory_out + "SkinTemperature.bin");
      OutputMeteo.Append( Landfraction_out, directory_out + "Landfraction.bin");
       // OutputMeteo.Append( Stability_out, directory_out + "Stability.bin");
       // OutputMeteo.Append(BoundaryHeight_out, directory_out + "BoundaryHeight.bin");
      //  OutputMeteo.Append(SurfacePressure_out, directory_out + "SurfacePressure.bin");
        OutputMeteo.Append(CloudHeight_out, directory_out + "CloudHeight.bin");
        //OutputMeteo.Append(CloudFraction, directory_out + "CloudFraction.bin");
        //OutputMeteo.Append(RelativeHumidity, directory_out + "RelativeHumidity.bin");
       // OutputMeteo.Append(CRH, directory_out + "CRH.bin");
        OutputMeteo.Append(Rain_out, directory_out + "Rain.bin");
        OutputMeteo.Append(Kz_out, directory_out + "Kz_Louis.bin");
        ComputeModule(MeridionalWind_out, ZonalWind_out, WindModule_out);
      //  OutputMeteo.Append(WindModule_out, directory_out + "WindModule.bin");
       // OutputMeteo.Append(SensibleHeat_out, directory_out + "SensibleHeat.bin");
       // OutputMeteo.Append(SolarRadiation_out, directory_out + "SolarRadiation.bin");
       // OutputMeteo.Append(PotentialTemperature_out, directory_out + "PotentialTemperature.bin");

        ComputeSpecificHumidity(SpecificHumidity_out, Temperature_out,Pressure_out,RelativeHumidity_out);
        OutputMeteo.Append(SpecificHumidity_out, directory_out + "SpecificHumidity.bin");
       // OutputMeteo.Append(FrictionModule_out, directory_out + "FrictionModule.bin");
       // OutputMeteo.Append(Richardson_out, directory_out + "Richardson.bin");


        cout << " done." << endl;

        cout << endl;

        count=0;


        free(level);
        free(shortName);
        free( typeOfLevel);
        fclose(f);
        fclose(f_prev);


    }


    END;

    return 0;

}

static void usage(char* progname)
{
    printf("\nUsage: %s grib_file\n",progname);
    exit(1);
}
