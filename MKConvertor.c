/*--------------------------------------------------------------------------------------------------
This code reads atomic coordinations from Materials Studio XSD file. The XSD file should be made by "Balls and Sticks" 3D visualization mode in Materials Studio.
Outputs of this code depending on users choice can be:
Filename.XYZ: contains coordination, type and dimentions of cell (if PBC exists)
Filename.bond: contains bonds informations.
Filename.frac: contains fractional coordinates, type of atoms and cell dimentions (It will be produced just if PBC exists)
Filename.LAMMPS: LAMMPS input file in charge format
Filename.DLPOLY: DL_POLY CONFIG file
--------------------------------------------------------------------------------------------------*/

#include	"stdio.h"
#include	"stdlib.h"
#include	"math.h"
#include    "string.h"

#define		Max	512

/*This function reads XSD file and write informations about coordinations, PBC and atom types in Filename.MK*/
FILE* MK(const char *FileName, int input)
{
    char    Line[Max], Filein[Max], Filein2[Max], substr1[20], substr2[20], atom_name[5], atom_name_org[10], atom_ff[10], Ans[5], Shell[100][5];;
    char    element[100][5];
    FILE	*in, *out1, *out2, *out3;
    int     element_repeat[100];
    int	    i, j, tempID, NumProp, atom_types = 0, flag = 0;
    int     Atom_Num = 0, Bond_Num = 0, Type_Num = 0, shell_choice = 0;
    int     PBC1 = 0, PBC2 = 0, PBC3 = 0;
    double	tempx, tempy, tempz, tempq = 0.0, Charge = 0.0;
	double  A1, A2, A3, B1, B2, B3, C1, C2, C3;
	double  fac1 = 1.0, fac2 = 1.0, fac3 = 1.0;
	double  cosA = 0.0, sinA = 1.0, cosB = 0.0, sinB = 1.0, cosG = 0.0, sinG = 1.0;
	int     bond1, bond2, temp_bond1, temp_bond2;
	double  xlo=100000, xhi=-100000, ylo=100000, yhi=-100000, zlo=100000, zhi=-100000;
	int     mode, boundary, Number;
	double  Xmin,Xmax,Ymin,Ymax;

    for (i = 0; i < 100; i++)
    {
        strcpy(element[i], "");
        element_repeat[i] = 0;
    }


    sprintf(Filein, "%s", FileName);
	if(( in = fopen(Filein, "r")) == NULL)
    {
        printf("can not open %s to read\n", Filein);
        exit(0);
    }

    //Reading LAMMPS file
    if (input == 4)
    {
        int len, type, tempType;
        double mass;
        const char *last_char;

        fgets(Line, Max, in);
        len = strlen(Line);
        last_char = &Line[len-6];

        while (strncmp(last_char, "Atoms", 5) != 0)
        {
            fgets(Line, Max, in);
            len = strlen(Line);
            last_char = &Line[len-6];
            //printf ("%s\n", last_char);
            if (strncmp(last_char, "atoms", 5) == 0)
            {
                sscanf(Line, "%d %*s", &Atom_Num);
                printf ("Total Number of Atoms: %d\n", Atom_Num);
                continue;
            }
            last_char = &Line[len-11];
            if (strncmp(last_char, "atom types", 10) == 0)
            {
                sscanf(Line, "%d %*s", &atom_types);
                printf ("Total Number of Atom Types: %d\n", atom_types);
                continue;
            }
            last_char = &Line[len-4];
            if (strncmp(last_char, "xhi", 3) == 0)
            {
                sscanf(Line, "%lf %lf ", &xlo, &xhi);
                A1 = xhi - xlo;
                A2 = 0.0;
                A3 = 0.0;
                PBC1 = 1;
                fac1 = A1;
                continue;
            }
            last_char = &Line[len-4];
            if (strncmp(last_char, "yhi",3) == 0)
            {
                sscanf(Line, "%lf %lf ", &ylo, &yhi);
                B1 = 0.0;
                B2 = yhi - ylo;
                B3 = 0.0;
                PBC2 = 1;
                fac2 = B2;
                continue;
            }
            last_char = &Line[len-4];
            if (strncmp(last_char, "zhi", 3) == 0)
            {
                sscanf(Line, "%lf %lf ", &zlo, &zhi);
                C1 = 0.0;
                C2 = 0.0;
                C3 = zhi - zlo;
                PBC3 = 1;
                fac3 = C3;
                continue;
            }

            last_char = &Line[len-7];

            if (strncmp(last_char, "Masses", 6) == 0)
            {
                printf ("Please enter names of the following atom types:\n");
                fgets(Line, Max, in);
                for (j = 0; j < atom_types; j ++)
                {
                    fgets(Line, Max, in);
                    sscanf(Line, "%d %lf ", &type, &mass);
                    printf ("%d (mass = %9.6lf): ", type, mass);
                    scanf ("%s", element[j]);
                }
            }

            last_char = &Line[len-6];
            if (strncmp(last_char, "Atoms", 5) == 0)
                break;
        }

        fgets(Line, Max, in);

        //opening temp.Read to write temporary XYZ informaion
        if(( out1 = fopen("temp.Read", "w")) == NULL)
        {
            printf("can not open temp.Read to write!\n");
            exit(0);
        }

        for (i = 0; i < Atom_Num; i ++)
        {
            fgets(Line, Max, in);
            sscanf (Line, "%d %*d %d %lf %lf %lf %lf", &tempID, &tempType, &tempq, &tempx, &tempy, &tempz);
            fprintf (out1, "%14s%14d%14d%16d%20.10lf%20.10lf%20.10lf%20.10lf\n", element[tempType-1], tempType, tempType,
            tempID, tempq, tempx/fac1, tempy/fac2, tempz/fac3);
        }

        fclose (out1);
        atom_types = 0;
    }

    xlo=100000, xhi=-100000, ylo=100000, yhi=-100000, zlo=100000, zhi=-100000;
    atom_types = 0;

    //Reading CONFIG file
    if (input == 1)
    {
        int     Shell_Num;

        fgets(Line, Max, in);
        fgets(Line, Max, in);
        sscanf(Line, "%d %d %d", &mode, &boundary, &Number);

        printf ("Do you want to exclude shells (Mitchel and Finchham core-shell model) in the output file? ");
        scanf("%[^\n]", Ans);

        if (strncmp(Ans, "y", 1) == 0 || strncmp(Ans, "Y", 1) == 0)
        {
            shell_choice = 1;
            printf ("   How many shell types are in the input file? ");
            scanf("%d", &Shell_Num);

            for (i = 0; i < Shell_Num; i ++)
            {
                printf ("       Please enter the name of shell type %d: ", i + 1);
                scanf("%[^\n]", Shell[i]);
            }
        }


        if (boundary == 3)
        {
            fgets(Line, Max, in);
            sscanf(Line, "%lf %lf %lf", &A1, &A2, &A3);
            fgets(Line, Max, in);
            sscanf(Line, "%lf %lf %lf", &B1, &B2, &B3);
            fgets(Line, Max, in);
            sscanf(Line, "%lf %lf %lf", &C1, &C2, &C3);
            PBC1 = 1; PBC2 = 1; PBC3 = 1;
            fac1 = A1; fac2 = B2; fac3 = C3;
        }

        //opening temp.Read to write temporary XYZ informaion
        if(( out1 = fopen("temp.Read", "w")) == NULL)
        {
            printf("can not open temp.Read to write!\n");
            exit(0);
        }

        for (i = 0; i < Number; i ++)
        {
            fgets(Line, Max, in);
            sscanf(Line, "%s %d", atom_name, &tempID);
            fgets(Line, Max, in);
            sscanf(Line, "%lf %lf %lf", &tempx, &tempy, &tempz);
            if (mode != 0)
            {
                fgets(Line, Max, in);
                if (mode == 2)
                {
                    fgets(Line, Max, in);
                }

            }
            if (shell_choice != 1)
            {
                fprintf (out1, "%14s%14s%14s%16d%20.10lf%20.10lf%20.10lf%20.10lf\n", atom_name,atom_name,atom_name,
                tempID, 0.0, tempx/fac1, tempy/fac2, tempz/fac3);
                Atom_Num ++;
            }
            else
            {
                flag = 0;
                for (j = 0; j < Shell_Num; j ++)
                    if (strncmp(atom_name, Shell[j], 5) == 0)
                        flag = 1;
                if (flag == 0)
                {
                    fprintf (out1, "%14s%14s%14s%16d%20.10lf%20.10lf%20.10lf%20.10lf\n", atom_name,atom_name,atom_name,
                    tempID, 0.0, tempx/fac1, tempy/fac2, tempz/fac3);
                    Atom_Num ++;
                }
            }
        }
        fclose (out1);

    }

    //Reading XYZ file
    if (input == 2)
    {
        tempID = 0;
        fgets(Line, Max, in);
        sscanf(Line, "%d", &Atom_Num);
        printf("%d\n",Atom_Num);
        fgets(Line, Max, in);

        PBC1 = 0; PBC2 = 0; PBC3 = 0;
        fac1 = 1; fac2 = 1; fac3 = 1;

        //opening temp.Read to write temporary XYZ informaion
        if(( out1 = fopen("temp.Read", "w")) == NULL)
        {
            printf("can not open temp.Read to write!\n");
            exit(0);
        }

        for (i = 0; i < Atom_Num; i ++)
        {
            tempID ++;
            fgets(Line, Max, in);
            sscanf(Line, "%s %lf %lf %lf", atom_name, &tempx, &tempy, &tempz);
            fprintf (out1, "%14s%14s%14s%16d%20.10lf%20.10lf%20.10lf%20.10lf\n", atom_name, atom_name, atom_name,
            tempID, 0.0, tempx/fac1, tempy/fac2, tempz/fac3);
        }
        fclose (out1);
    }

    //Reading XSD file
    if (input == 3)
    {
        printf("Reading XSD file ... \n");
        fgets(Line, Max, in);
        fgets(Line, Max, in);
        fgets(Line, Max, in);
        fgets(Line, Max, in);
        //Reading Number of properties ans children from XSD file

        sscanf(strstr(Line, "NumProperties="), "%*15c%d", &NumProp);

        printf("No. Prop = %d\n", NumProp);

        //we don't need properties for we pass them
        for (i = 0; i < NumProp; i++)
        {
            fgets(Line, Max, in);
        }

        //opening temp.Read to write temporary XYZ informaion
        if(( out1 = fopen("temp.Read", "w")) == NULL)
        {
            printf("can not open temp.Read to write!\n");
            exit(0);
        }

        sprintf(Filein, "%s.Bonds", FileName);
        if(( out2 = fopen(Filein, "w")) == NULL)
        {
            printf("can not open Bond file to write!\n");
            exit(0);
        }

        //Reading XSD file and writing corresponding XYZ values in temp.Read
        while (!feof(in))
        {
            fgets(Line, Max, in);

            //in XSD files just lines that start with "<Atom3d" contain XYZ information.
            sscanf(Line, "%*[^'<'] %[^' ID']", substr1);
            //Reading atomic coordinates
            //printf("%s\n",substr1);
            if (strncmp(substr1, "<Atom3d", 7) == 0 && !strstr(Line, "XYZ=") && strstr(Line, "Components=\""))
            {
                Atom_Num ++;
                sscanf(strstr(Line, "ID="), "%*4c%d", &tempID);
                //sscanf(strstr(Line, "Name="), "%*6c%[^\"]", atom_name_org);
                sscanf(strstr(Line, "Components=\""), "%*12c%[^\"]", atom_name);
                if (strstr(Line, "Charge=")) sscanf(strstr(Line, "Charge="), "%*8c%lf", &tempq);
                if (strstr(Line, "ForcefieldType=")) sscanf(strstr(Line, "ForcefieldType="), "%*16c%[^\"]", atom_ff);

                fprintf (out1, "%14s%14s%14s%16d%20.10lf%20.10lf%20.10lf%20.10lf\n", atom_name,
                atom_name, atom_ff, tempID, tempq, 0.0, 0.0, 0.0);

                printf ("%14s%14s%14s%16d%20.10lf\n", atom_name,
                atom_name, atom_ff, tempID, tempq);
            }


            if (strncmp(substr1, "<Atom3d", 7) == 0 && strstr(Line, "XYZ="))
            {
                Atom_Num ++;
                sscanf(strstr(Line, "ID="), "%*4c%d", &tempID);
                //sscanf(strstr(Line, "Name="), "%*6c%[^\"]", atom_name_org);
                sscanf(strstr(Line, "Components=\""), "%*12c%[^\"]", atom_name);
                sscanf(strstr(Line, "XYZ="), "%*5c%lf,%lf,%lf", &tempx, &tempy, &tempz);
                if (strstr(Line, "Charge=")) sscanf(strstr(Line, "Charge="), "%*8c%lf", &tempq);
                if (strstr(Line, "ForcefieldType=")) sscanf(strstr(Line, "ForcefieldType="), "%*16c%[^\"]", atom_ff);

                fprintf (out1, "%14s%14s%14s%16d%20.10lf%20.10lf%20.10lf%20.10lf\n", atom_name,
                atom_name, atom_ff, tempID, tempq, tempx, tempy, tempz);
            }



            //Reading bonds information
            if (strncmp(substr1, "<Bond", 5) == 0)
            {
                Bond_Num ++;
                sscanf(strstr(Line, "Connects="), "%*10c%d,%d", &bond1, &bond2);
                fprintf (out2, "%d %d %d \n", Bond_Num, bond1, bond2);
            }

            //Reading PBC
            if ((strncmp(substr1, "<PlaneGroup", 11) == 0) || (strncmp(substr1, "<SpaceGroup", 11) == 0))
            {
                sscanf(strstr(Line, "AVector="), "%*9c%lf,%lf,%lf", &A1, &A2, &A3);
                sscanf(strstr(Line, "BVector="), "%*9c%lf,%lf,%lf", &B1, &B2, &B3);
                sscanf(strstr(Line, "CVector="), "%*9c%lf,%lf,%lf", &C1, &C2, &C3);

                //printf("%lf %lf %lf \n",A1, A2, A3);
                //printf("%lf %lf %lf \n",B1, B2, B3);
                //printf("%lf %lf %lf \n",C1, C2, C3);

                if (A1 != 1)
                    PBC1 = 1;
                if (B2 != 1)
                    PBC2 = 1;
                if (C3 != 1)
                    PBC3 = 1;
                break;
            }
        }

        fclose (out1);
        fclose (out2);
    }

    fclose(in);

    struct MK
    {
        int     ID;
        double  X;
        double  Y;
        double  Z;
        double  q;
        char    name[5];
        char    name_org[10];
        char    ff[10];
        int     repeat;
    };

    struct  MK *atom;

    atom = (struct MK*) malloc(Atom_Num * sizeof(struct MK));

    if(( in = fopen("temp.Read", "r")) == NULL)
    {
        printf("can not open temp.read to read\n");
        exit(0);
    }

    printf("Writing atoms informations into atom STRUCT ... \n");
    //Writing atoms informations into atom STRUCT
    for (i = 0; i < Atom_Num; i ++)
    {
        fgets(Line, Max, in);
        sscanf (Line, "%s %s %s %d %lf %lf %lf %lf", atom_name, atom_name_org, atom_ff, &tempID, &tempq, &tempx, &tempy, &tempz);
        strcpy(atom[i].name,atom_name);
        strcpy(atom[i].name_org,atom_name_org);
        strcpy(atom[i].ff,atom_ff);
        atom[i].ID = tempID;
        atom[i].X = tempx;
        atom[i].Y = tempy;
        atom[i].Z = tempz;
        atom[i].q = tempq;
        atom[i].repeat = 1;
        Charge += tempq;
        //finding lower and upper limit of dimentions
        if (tempx < xlo)
            xlo = tempx;
        if (tempx > xhi)
            xhi = tempx;
        if (tempy < ylo)
            ylo = tempy;
        if (tempy > yhi)
            yhi = tempy;
        if (tempz < zlo)
            zlo = tempz;
        if (tempz > zhi)
            zhi = tempz;
    }

    fclose(in);

    if (fabs(Charge) > 0.001)
    {
        printf("##############################################################\n");
        printf("#########WARNING: STRUCTURE IS NOT CHARGE NEUTRAL!!!##########\n");
        printf("#########TOTAL CHARGE = %10.3lf           ##########\n", Charge);
        printf("##############################################################\n");
    }

    printf("Finding the number of each atom in the configuration ... \n");
    //Finding the number of each atom in the configuration
    for (i = 0; i < Atom_Num; i ++)
    {
        printf("\r%d%c", 100*i/Atom_Num, 37);
        for (j = 0; j < Atom_Num; j ++)
        {
            if ((i != j) && (strncmp(atom[i].name, atom[j].name, 5) == 0))
                atom[i].repeat ++;
        }
    }

    printf("\n");

    printf("Calculation number of atom types ... \n");
    //Calculation number of atom types
    for (i = 0; i < Atom_Num; i ++)
    {
        printf("\r%d%c", 100*i/Atom_Num, 37);

        flag = 0;
        for (j = i+1; j < Atom_Num; j ++)
        {
            if (strncmp(atom[i].name, atom[j].name, 5) == 0)
                flag = 1;
        }
        if (flag == 0)
        {
            strcpy(element[atom_types],atom[i].name);
            element_repeat[atom_types] = atom[i].repeat;
            atom_types ++;
        }
    }

    printf("\n");

    printf("Writing extracted data into Filename.MK file ... \n");
    //Writing extracted data into Filename.MK file
    sprintf(Filein, "%s.MK", FileName);
    if(( out1 = fopen(Filein, "w")) == NULL)
    {
        printf("can not open %s to write!\n", Filein);
        exit(0);
    }

    fprintf (out1, "NUMBER OF ATOMS\n");
    fprintf (out1, "%d\n", Atom_Num);
    fprintf (out1, "NUMBER OF BondsS\n");
    fprintf (out1, "%d\n", Bond_Num);
    fprintf (out1, "TOTAL CHARGE\n");
    fprintf (out1, "%8.3lf\n", Charge);
    fprintf (out1,"ATOM TYPES\n");
    fprintf (out1, "%d\n", atom_types);
    for (i = 0; i < atom_types; i ++)
        fprintf (out1, "%6d%6s\n", element_repeat[i], element[i]);
    fprintf (out1,"PERIODIC BOUNDARY CONDITION\n");
    fprintf (out1, "%5d%5d%5d\n", PBC1, PBC2, PBC3);

    //for fractional coordinates, U should be along X, and V should be in XY plane
    if ((PBC1 != 0) || (PBC2 != 0) || (PBC3 != 0))
        fprintf(out1, "CELL DIMENTIONS\n");
    if (PBC1 == 1)
    {
        fac1 = sqrt(A1*A1+A2*A2+A3*A3);
        fprintf (out1, "%20.10lf%20.10lf%20.10lf\n", A1, A2, A3);
    }

    if (PBC2 == 1)
    {
        fac2 = sqrt(B1*B1+B2*B2+B3*B3);
        cosG = (A1*B1+A2*B2+A3*B3)/(fac1*fac2);
        sinG = sqrt(1-cosG*cosG);
        fprintf (out1, "%20.10lf%20.10lf%20.10lf\n", B1, B2, B3);
    }
    if (PBC3 ==1)
    {
        fac3 = sqrt(C1*C1+C2*C2+C3*C3);
        cosB = (A1*C1+A2*C2+A3*C3)/(fac1*fac3);
        sinB = sqrt(1-cosB*cosB);
        cosA = (B1*C1+B2*C2+B3*C3)/(fac2*fac3);
        sinA = sqrt(1-cosA*cosA);
        fprintf (out1, "%20.10lf%20.10lf%20.10lf\n", C1, C2, C3);
    }

    fprintf (out1,"COORDINATION LIMITS\n");
    fprintf (out1, "%20.10lf%20.10lf xlo xhi\n", fac1*xlo, fac1*xhi);
    fprintf (out1, "%20.10lf%20.10lf ylo yhi\n", fac2*ylo, fac2*yhi);
    fprintf (out1, "%20.10lf%20.10lf zlo zhi\n", fac3*zlo, fac3*zhi);
    fprintf (out1,"ATOMS\n");
    fprintf (out1,"%6s%16s%16s%16s%16s%20s%20s%20s%20s\n", "Name", "Name org", "FF Type", "ID", "ID-B", "Charge", "X", "Y", "Z");
    for (i = 0; i < Atom_Num; i ++)
        fprintf (out1,"%6s%16s%16s%16d%16d%20.10lf%20.10lf%20.10lf%20.10lf\n", atom[i].name, atom[i].name_org, atom[i].ff, i + 1,
        atom[i].ID, atom[i].q, fac1*atom[i].X+fac2*cosG*atom[i].Y+fac3*cosB*atom[i].Z,
        fac2*sinG*atom[i].Y+fac3*(cosA-cosB*cosG)*atom[i].Z/sinG,
        fac3*sqrt(1-cosA*cosA-cosB*cosB-cosG*cosG+2*cosA*cosB*cosG)*atom[i].Z/sinG);

    sprintf(Filein2, "%s.Bonds", FileName);
    if(( in = fopen(Filein2, "r")) != NULL)
    {
        fprintf (out1,"BONDS\n");
        for (i = 0; i < Bond_Num; i ++)
        {
            fgets(Line, Max, in);
            sscanf (Line, "%d %d %d", &tempID, &temp_bond1, &temp_bond2);
            bond1 = 0;
            bond2 = 0;
            for (j = 0; j < Atom_Num; j ++)
            {
                if (atom[j].ID == temp_bond1)
                {
                    bond1 = j + 1;
                }
                if (atom[j].ID == temp_bond2)
                {
                    bond2 = j + 1;
                }
            }
            if (bond1 != 0 && bond2 != 0)
                fprintf(out1,"%6d%16d%16d\n", i + 1, bond1, bond2);
        }
        fclose(in);
    }

    fclose (out1);

    if(( in = fopen(Filein, "r")) == NULL)
    {
        printf("can not open %s!\n", Filein);
        exit(0);
    }

    remove ("temp.Read");
    remove (Filein2);

    free(atom);

    return in;
}

/*THis function reads information from Filename.MK file and writes the Moltemplate file. Atomic coordinates and bond information is written in the Moltemplate file. It works with OPLSAA.
 Rest of the required information for LAMMPS input file is taken care of by Moltemplate.*/
int Moltemplate(const char *FileName, int input)
{
    char    Line[Max], Filein[Max], substr1[20];
    FILE    *in, *out;
    int     i, j=0, k, atom_count = 0, PBC = 0;
    int     Atom_Num, Atom_Type, Bond_Num;
    int     PBC1, PBC2, PBC3;
    int     bond1, bond2, ff_choice, line_num;
    double  A1, A2, A3, B1, B2, B3, C1, C2, C3;
	double  xlo, xhi, ylo, yhi, zlo, zhi;

    //reading OPLSAA atom types from oplsaa.prm file
    // in oplsaa vdW parameters are defined by atom type number but
    //bonded interactions are defined by atom class. Moltemplate needs used
    //to define bonds. we need to find atom class for each atom type number to find the corresponding bond type

    printf("Chhose forcefield (1) oplsaa, (2) oplsaa-MoS2-Rajan-NoBond (3) oplsaa-MoS2-Rajan: ");
    scanf("%d", &ff_choice);

    printf("reading OPLSAA atom types from oplsaa.prm file ... \n");

    int atom_type[908];

    if (ff_choice == 1)
    {
        line_num = 958;
        sprintf(Filein, "/home/mk/Uni/Codes/bin/ForceFields/oplsaa.prm");
    }

    else
    {
        line_num = 960;
        sprintf(Filein, "/home/mk/Uni/Codes/bin/ForceFields/oplsaa-MoS2-Rajan.prm");
    }

    if(( in = fopen(Filein, "r")) == NULL)
    {
        printf("can not open %s to read!\n", Filein);
        exit(0);
    }

    for (i=0;i<line_num;i++)
    {
        fgets(Line, Max, in);
        if (i>51)
        {
            sscanf(Line,"%*s %*d %d", &atom_type[j]);
            j ++;
        }
    }

    fclose(in);
    //Opening MK file
    in = MK(FileName,input);

    printf("Reading MK file ... \n");

    //Reading MK file
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    sscanf (Line, "%d", &Atom_Num);
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    sscanf (Line, "%d", &Bond_Num);
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    sscanf (Line, "%d", &Atom_Type);

    char    element[Atom_Type][5];
    int     repeat[Atom_Type];

    for (i = 0; i < Atom_Type; i ++)
    {
        fgets(Line, Max, in);
        sscanf (Line, "%d %s", &repeat[i], element[i]);
    }

    fgets(Line, Max, in);
    fgets(Line, Max, in);
    sscanf (Line, "%d %d %d", &PBC1, &PBC2, &PBC3);
    if ((PBC1 != 0) || (PBC2 != 0) || (PBC3 != 0))
        fgets(Line, Max, in);

    if (PBC1 == 1)
    {
        fgets(Line, Max, in);
        sscanf (Line, "%lf %lf %lf", &A1, &A2, &A3);
        PBC ++;
    }


    if (PBC2 == 1)
    {
        fgets(Line, Max, in);
        sscanf (Line, "%lf %lf %lf", &B1, &B2, &B3);
        PBC ++;
    }


    if (PBC3 ==1)
    {
        fgets(Line, Max, in);
        sscanf (Line, "%lf %lf %lf", &C1, &C2, &C3);
        PBC ++;
    }

    if (PBC == 2)
        PBC = 6;

    fgets(Line, Max, in);
    fgets(Line, Max, in);
    sscanf (Line, "%lf %lf", &xlo, &xhi);
    fgets(Line, Max, in);
    sscanf (Line, "%lf %lf", &ylo, &yhi);
    fgets(Line, Max, in);
    sscanf (Line, "%lf %lf", &zlo, &zhi);
    fgets(Line, Max, in);
    fgets(Line, Max, in);

    struct MK
    {
        int     ID;
        int     ff;
        double  X;
        double  Y;
        double  Z;
        char    name[5];
        char    name_org[10];
    };

    struct  MK *atom;

    atom = (struct MK*) malloc(Atom_Num * sizeof(struct MK));


    struct bond
    {
        int     atom1;
        int     atom2;
    }bond[Bond_Num];

    printf("Writing Moltemplate files ... \n");

    for (i = 0; i < Atom_Num; i ++)
    {
        fgets(Line, Max, in);
        sscanf (Line, "%s %s %s %d %*d %*f %lf %lf %lf", atom[i].name, atom[i].name_org, substr1,
        &atom[i].ID, &atom[i].X, &atom[i].Y, &atom[i].Z);

        atom[i].ff = atoi(substr1);
    }

    fgets(Line, Max, in);

    for (i = 0; i < Bond_Num; i ++)
    {
        fgets(Line, Max, in);
        sscanf (Line, "%*d %d %d", &bond[i].atom1, &bond[i].atom2);
    }

    fclose(in);

    sprintf(Filein, "%s.lt", FileName);
    if(( out = fopen(Filein, "w")) == NULL)
    {
        printf("can not open %s to write!\n", Filein);
        exit(0);
    }

    fprintf (out, "# file \"%s.lt\"\n", FileName);
    fprintf (out, "# This file is made by MKconvertor code (2014-khalkhal@ualberta.ca)\n");
    fprintf (out, "\n");
    if (ff_choice == 1)
    {
        fprintf (out, "import \"oplsaa.lt\"\n");
    }
    else if (ff_choice == 2)
    {
        fprintf (out, "import \"oplsaa-MoS2-Rajan-NoBond.lt\"\n");
    }
    else
    {
        fprintf (out, "import \"oplsaa-MoS2-Rajan.lt\"\n");
    }
    fprintf (out, "\n");
    fprintf (out, "Molecule inherits OPLSAA {\n");
    fprintf (out, "\n");
    fprintf (out, "  write(\"Data Atoms\") {\n");
    for (i = 0; i < Atom_Num; i ++)
    {
       fprintf (out, "    $atom:%s%d $mol:. @atom:%d  0.0%8.3lf%8.3lf%8.3lf\n", atom[i].name_org, atom[i].ID, atom[i].ff,
       atom[i].X, atom[i].Y,atom[i].Z);
    }
    fprintf (out, "  }\n");
    fprintf (out, "  write(\"Data Bonds\") {\n");
    for (i = 0; i < Bond_Num; i ++)
    {
        bond1 = atom_type[atom[bond[i].atom1-1].ff-1];
        bond2 = atom_type[atom[bond[i].atom2-1].ff-1];
        if (atom_type[atom[bond[i].atom1-1].ff-1]>atom_type[atom[bond[i].atom2-1].ff-1])
        {
            bond1 = atom_type[atom[bond[i].atom2-1].ff-1];
            bond2 = atom_type[atom[bond[i].atom1-1].ff-1];
        }
        fprintf (out, "    $bond:%s%d%s%d @bond:%d-%d $atom:%s%d $atom:%s%d\n", atom[bond[i].atom1-1].name_org,  atom[bond[i].atom1-1].ID,
        atom[bond[i].atom2-1].name_org,  atom[bond[i].atom2-1].ID,bond1, bond2,
        atom[bond[i].atom1-1].name_org,  atom[bond[i].atom1-1].ID, atom[bond[i].atom2-1].name_org,  atom[bond[i].atom2-1].ID);
    }
    fprintf (out, "  }\n");
    fprintf (out, "\n");
    fprintf (out, "}\n");

    fclose(out);

    if(( out = fopen("system.lt", "w")) == NULL)
    {
        printf("can not open system.lt to write!\n");
        exit(0);
    }

    fprintf (out, "import \"%s\"\n", Filein);
    fprintf (out, "\n");
    fprintf (out, "# Periodic boundary conditions:\n");
    fprintf (out, "\n");
    fprintf (out, "write_once(\"Data Boundary\") {\n");
    if (PBC1 == 1)
        fprintf(out,"   %12.5lf %12.5lf  xlo xhi\n", xlo, xlo+A1);
    else
        fprintf(out,"   %12.5lf %12.5lf  xlo xhi\n", xlo, xhi);
    if (PBC2 == 1)
        fprintf(out,"   %12.5lf %12.5lf  ylo yhi\n", ylo, ylo+B2);
    else
        fprintf(out,"   %12.5lf %12.5lf  ylo yhi\n", ylo, yhi);
    if (PBC3 == 1)
        fprintf(out,"   %12.5lf %12.5lf  zlo zhi\n", zlo, zlo+C3);
    else
        fprintf(out,"   %12.5lf %12.5lf  zlo zhi\n", zlo, zhi);
    fprintf (out, "}\n");
    fprintf (out, "\n");
    fprintf (out, "mol  = new Molecule [1]\n");

    fclose(out);

    return(0);
}


/*----------------------------------------------------------------------------------------------------------------------------------------------------------------
This Function reads Filename.MK file and writes information in the Filename.CONFIG file in the DL_POLY CONFIG file format.
This code writes different atom types seperately. For example it write atoms of type 1 coordinations into CONFIG file and then it
writes atoms of type 2 underneath. As molecules are not defined in XSD file, there is no way to define define them in CONFIG file
either unless you change the code and rewrite the CONFIG file in a way that you want.*/
int CONFIG(const char *FileName, int input)
{
    char    Line[Max], Filein[Max], substr1[20], substr2[20], atom_name[5], shell_name[5], Ans[5];
    FILE    *in, *out;
    int     i, j, k, atom_count = 0, PBC = 0;
    int     Atom_Num, Atom_Type, Bond_Num;
    int     PBC1, PBC2, PBC3;
    double  A1, A2, A3, B1, B2, B3, C1, C2, C3;
	double  xlo, xhi, ylo, yhi, zlo, zhi;

    //Opening MK file
    in = MK(FileName,input);

    //Reading MK file
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    sscanf (Line, "%d", &Atom_Num);
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    sscanf (Line, "%d", &Bond_Num);
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    sscanf (Line, "%d", &Atom_Type);

    char    element[Atom_Type][5];
    int     repeat[Atom_Type];

    for (i = 0; i < Atom_Type; i ++)
    {
        fgets(Line, Max, in);
        sscanf (Line, "%d %s", &repeat[i], element[i]);
    }



    fgets(Line, Max, in);
    fgets(Line, Max, in);
    sscanf (Line, "%d %d %d", &PBC1, &PBC2, &PBC3);
    if ((PBC1 != 0) || (PBC2 != 0) || (PBC3 != 0))
        fgets(Line, Max, in);

    if (PBC1 == 1)
    {
        fgets(Line, Max, in);
        sscanf (Line, "%lf %lf %lf", &A1, &A2, &A3);
        PBC ++;
    }


    if (PBC2 == 1)
    {
        fgets(Line, Max, in);
        sscanf (Line, "%lf %lf %lf", &B1, &B2, &B3);
        PBC ++;
    }


    if (PBC3 ==1)
    {
        fgets(Line, Max, in);
        sscanf (Line, "%lf %lf %lf", &C1, &C2, &C3);
        PBC ++;
    }

    if (PBC == 2)
        PBC = 6;

    fgets(Line, Max, in);
    fgets(Line, Max, in);
    sscanf (Line, "%lf %lf", &xlo, &xhi);
    fgets(Line, Max, in);
    sscanf (Line, "%lf %lf", &ylo, &yhi);
    fgets(Line, Max, in);
    sscanf (Line, "%lf %lf", &zlo, &zhi);
    fgets(Line, Max, in);
    fgets(Line, Max, in);

    double  Dim_corr1, Dim_corr2, Dim_corr3;
  	//This value corrects the dimention of simulation box to be symetric around 0
    Dim_corr1 = -0.5*(xhi + xlo);
    Dim_corr2 = -0.5*(yhi + ylo);
    Dim_corr3 = -0.5*(zhi + zlo);

    xhi = xhi + Dim_corr1;
    yhi = yhi + Dim_corr2;
    zhi = zhi + Dim_corr3;
    xlo = xlo + Dim_corr1;
    ylo = ylo + Dim_corr2;
    zlo = zlo + Dim_corr3;

    struct atom
    {
        int     ID;
        double  X;
        double  Y;
        double  Z;
        double  q;
        char    name[5];
        int     check;
    }atom[Atom_Num];

    for (i = 0; i < Atom_Num; i ++)
    {
        fgets(Line, Max, in);
        sscanf (Line, "%s %*s %*s %*d %d %lf %lf %lf %lf", atom[i].name, &atom[i].ID, &atom[i].q, &atom[i].X, &atom[i].Y, &atom[i].Z);
        atom[i].X = atom[i].X + Dim_corr1;
        atom[i].Y = atom[i].Y + Dim_corr2;
        atom[i].Z = atom[i].Z + Dim_corr3;
        atom[i].check = 1;
    }

    fclose(in);

    sprintf(Filein, "%s.CONFIG", FileName);
    if(( out = fopen(Filein, "w")) == NULL)
    {
        printf("can not open %s to write!\n", Filein);
        exit(0);
    }

    //Writing atomic coordinatets into Filename.CONFIg file
    fprintf (out, "CONFIG file made by MKCnvertor code from %s file. This file contains ", FileName);
    for (i = 0; i < Atom_Type; i ++)
        fprintf (out, "%d %s, ", repeat[i], element[i]);
    fprintf (out, "atoms. (MK (2014): khalkhal@ualberta.ca)\n");

    fprintf (out, "%10d%10d%10d\n", 0, PBC, Atom_Num);

    if (PBC1 == 1)
        fprintf (out, "  %18.15lf  %18.15lf  %18.15lf\n", A1, A2, A3);

    if (PBC2 == 1)
        fprintf (out, "  %18.15lf  %18.15lf  %18.15lf\n", B1, B2, B3);

    if (PBC3 ==1)
        fprintf (out, "  %18.15lf  %18.15lf  %18.15lf\n", C1, C2, C3);

    printf ("   Do you like to keep the sequential order of elements in the origin file? ");
    scanf("%s", Ans);
    if (strncmp(Ans, "y", 1) == 0 || strncmp(Ans, "Y", 1) == 0)
    {
        for (i = 0; i < Atom_Num; i ++)
        {
            atom_count ++;
            fprintf (out, "%4s%16d\n", atom[i].name, atom[i].ID);
            fprintf (out, "%20.15lf%20.15lf%20.15lf\n", atom[i].X, atom[i].Y, atom[i].Z);
        }
    }
    else
    {
        for (j = 0; j < Atom_Type; j ++)
        {
            for (i = 0; i < Atom_Num; i ++)
            {
                //This controls that atom types are written sequentially.
                if (atom[i].check > 0 && strncmp(atom[i].name, element[j], 5) == 0)
                {
                    atom_count ++;
                    fprintf (out, "%4s%16d\n", atom[i].name, atom_count);
                    fprintf (out, "%20.15lf%20.15lf%20.15lf\n", atom[i].X, atom[i].Y, atom[i].Z);
                    sprintf(shell_name, "%s_s", atom[i].name);
                    for (k = 0; k < Atom_Num; k ++)
                    {
                        if (atom[k].check > 0 && strncmp(atom[k].name, shell_name, 5) == 0)
                        {
                            atom_count ++;
                            fprintf (out, "%4s%16d\n", atom[k].name, atom_count);
                            fprintf (out, "%20.15lf%20.15lf%20.15lf\n", atom[k].X, atom[k].Y, atom[k].Z);
                            atom[k].check = -1;
                            break;
                        }
                    }
                }
            }
        }
    }

    fclose (out);
}

/*----------------------------------------------------------------------------------------------------------------------------------------------------------------
This Function converts MK file to a Gaussian input file.*/
int GAUSSIAN(const char *FileName, int input)
{
    char    Line[Max], Filein[Max], substr1[20], substr2[20], atom_name[5];
    FILE    *in, *out;
    int     i, j, atom_count = 0, PBC = 0;
    int     Atom_Num, Atom_Type, Bond_Num;
    int     PBC1, PBC2, PBC3;
    double  A1, A2, A3, B1, B2, B3, C1, C2, C3;
	double  xlo, xhi, ylo, yhi, zlo, zhi;

	//Opening MK file
    in = MK(FileName,input);

    //Reading MK file
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    sscanf (Line, "%d", &Atom_Num);
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    sscanf (Line, "%d", &Bond_Num);
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    sscanf (Line, "%d", &Atom_Type);

    for (i = 0; i < Atom_Type; i ++)
        fgets(Line, Max, in);

    fgets(Line, Max, in);
    fgets(Line, Max, in);
    sscanf (Line, "%d %d %d", &PBC1, &PBC2, &PBC3);
    if ((PBC1 != 0) || (PBC2 != 0) || (PBC3 != 0))
        fgets(Line, Max, in);

    if (PBC1 == 1)
    {
        fgets(Line, Max, in);
        sscanf (Line, "%lf %lf %lf", &A1, &A2, &A3);
        PBC ++;
    }


    if (PBC2 == 1)
    {
        fgets(Line, Max, in);
        sscanf (Line, "%lf %lf %lf", &B1, &B2, &B3);
        PBC ++;
    }


    if (PBC3 ==1)
    {
        fgets(Line, Max, in);
        sscanf (Line, "%lf %lf %lf", &C1, &C2, &C3);
        PBC ++;
    }

    fgets(Line, Max, in);
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    fgets(Line, Max, in);

    struct atom
    {
        int     ID;
        double  X;
        double  Y;
        double  Z;
        double  q;
        char    name[5];
    }atom[Atom_Num];

    for (i = 0; i < Atom_Num; i ++)
    {
        fgets(Line, Max, in);
        sscanf (Line, "%s %*s %*s %*d %d %lf %lf %lf %lf", atom[i].name, &atom[i].ID, &atom[i].q, &atom[i].X, &atom[i].Y, &atom[i].Z);
    }


    sprintf(Filein, "%s.GAUSSIAN", FileName);
    if(( out = fopen(Filein, "w")) == NULL)
    {
        printf("can not open %s to write!\n", Filein);
        exit(0);
    }

    fprintf (out, "%%chk=Dump.chk    !Give the Dump file name here\n");
    fprintf (out, "%%mem=6000mb      !Memory\n");
    fprintf (out, "%%nprocshared=12  !Number of shared CPUs for parallel jobs\n");
    fprintf (out, "#p HF/6-31G(d)   !Desired calculation type, model chemistry and other options\n");
    fprintf (out, "\n");
    fprintf (out, "GAUSSIAN input file made by MKConvertor code (MK (2014): khalkhal@ualberta.ca)\n");
    fprintf (out, "\n");
    fprintf (out, "0 1              !The net electric charge and the spin multiplicity. For a neutral molecule in a singlet state it is 0 1\n");
    for (i = 0; i < Atom_Num; i ++)
        fprintf (out, "%2s%20.10lf%20.10lf%20.10lf\n", atom[i].name, atom[i].X, atom[i].Y, atom[i].Z);

    if (PBC1 == 1)
        fprintf (out, "TV%20lf%20lf%20lf\n", A1, A2, A3);
    if (PBC2 == 1)
        fprintf (out, "TV%20lf%20lf%20lf\n", B1, B2, B3);
    if (PBC3 == 1)
        fprintf (out, "TV%20lf%20lf%20lf\n", C1, C2, C3);

    fprintf (out, "\n");

    fclose (out);

    return PBC;
}

/*----------------------------------------------------------------------------------------------------------------------------------------------------------------
This Function converts MK file to a XYZ input file.*/
int XYZ(const char *FileName, int input)
{
    char    Line[Max], Filein[Max], substr1[20], substr2[20], atom_name[3];
    FILE    *in, *out;
    int     i, j, atom_count = 0, PBC = 0;
    int     Atom_Num, Atom_Type, Bond_Num;
    int     PBC1, PBC2, PBC3;
    double  A1, A2, A3, B1, B2, B3, C1, C2, C3;
	double  xlo, xhi, ylo, yhi, zlo, zhi;

    //Opening MK file
    in = MK(FileName,input);

    //Reading MK file
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    sscanf (Line, "%d", &Atom_Num);
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    sscanf (Line, "%d", &Bond_Num);
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    sscanf (Line, "%d", &Atom_Type);

    for (i = 0; i < Atom_Type; i ++)
        fgets(Line, Max, in);

    fgets(Line, Max, in);
    fgets(Line, Max, in);
    sscanf (Line, "%d %d %d", &PBC1, &PBC2, &PBC3);
    if ((PBC1 != 0) || (PBC2 != 0) || (PBC3 != 0))
        fgets(Line, Max, in);

    if (PBC1 == 1)
    {
        fgets(Line, Max, in);
        sscanf (Line, "%lf %lf %lf", &A1, &A2, &A3);
        PBC ++;
    }


    if (PBC2 == 1)
    {
        fgets(Line, Max, in);
        sscanf (Line, "%lf %lf %lf", &B1, &B2, &B3);
        PBC ++;
    }


    if (PBC3 ==1)
    {
        fgets(Line, Max, in);
        sscanf (Line, "%lf %lf %lf", &C1, &C2, &C3);
        PBC ++;
    }

    fgets(Line, Max, in);
    fgets(Line, Max, in);
    sscanf (Line, "%lf %lf", &xlo, &xhi);
    fgets(Line, Max, in);
    sscanf (Line, "%lf %lf", &ylo, &yhi);
    fgets(Line, Max, in);
    sscanf (Line, "%lf %lf", &zlo, &zhi);
    fgets(Line, Max, in);
    fgets(Line, Max, in);

    double  Dim_corr1, Dim_corr2, Dim_corr3;
  	//This value corrects the dimention of simulation box to be symetric around 0
    Dim_corr1 = -0.5*(xhi + xlo);
    Dim_corr2 = -0.5*(yhi + ylo);
    Dim_corr3 = -0.5*(zhi + zlo);

    struct atom
    {
        int     ID;
        double  X;
        double  Y;
        double  Z;
        double  q;
        char    name[5];
    }atom[Atom_Num];

    for (i = 0; i < Atom_Num; i ++)
    {
        fgets(Line, Max, in);
        sscanf (Line, "%s %*s %*s %*d %d %lf %lf %lf %lf", atom[i].name, &atom[i].ID, &atom[i].q, &atom[i].X, &atom[i].Y, &atom[i].Z);
        atom[i].X = atom[i].X + Dim_corr1;
        atom[i].Y = atom[i].Y + Dim_corr2;
        atom[i].Z = atom[i].Z + Dim_corr3;
    }


    sprintf(Filein, "%s.XYZ", FileName);
    if(( out = fopen(Filein, "w")) == NULL)
    {
        printf("can not open %s to write!\n", Filein);
        exit(0);
    }

    fprintf (out, "%d\n", Atom_Num);
    fprintf (out, "XYZ file made by MKConvertor code from %s file (MK (2014): khalkhal@ualberta.ca)\n", FileName);

    for (i = 0; i < Atom_Num; i ++)
    {
        strncpy(atom_name, atom[i].name,2);
        atom_name[2] = '\0';
        fprintf (out, "%6s%17.8lf%17.8lf%17.8lf\n", atom_name, atom[i].X, atom[i].Y, atom[i].Z);
    }

    fclose (out);

    return 0;
}

/*----------------------------------------------------------------------------------------------------------------------------------------------------------------
This Function converts MK file to a fractional coordinate file.*/

int FRAC(const char *FileName, int input)
{
    char    Line[Max], Filein[Max], substr1[20], substr2[20], atom_name[5];
    FILE    *in, *out;
    int     i, j, atom_count = 0, PBC = 0;
    int     Atom_Num, Atom_Type, Bond_Num;
    int     PBC1, PBC2, PBC3;
    double  A1 = 1.0, A2, A3, B1, B2 = 1.0, B3, C1, C2, C3 = 1.0;
	double  xlo, xhi, ylo, yhi, zlo, zhi;

    //Opening MK file
    in = MK(FileName,input);

    //Reading MK file
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    sscanf (Line, "%d", &Atom_Num);
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    sscanf (Line, "%d", &Bond_Num);
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    sscanf (Line, "%d", &Atom_Type);

    for (i = 0; i < Atom_Type; i ++)
        fgets(Line, Max, in);

    fgets(Line, Max, in);
    fgets(Line, Max, in);
    sscanf (Line, "%d %d %d", &PBC1, &PBC2, &PBC3);
    if ((PBC1 != 0) || (PBC2 != 0) || (PBC3 != 0))
        fgets(Line, Max, in);

    sprintf(Filein, "%s.FRAC", FileName);
    if(( out = fopen(Filein, "w")) == NULL)
    {
        printf("can not open %s to write!\n", Filein);
        exit(0);
    }

    fprintf (out, "%d\n", Atom_Num);
    fprintf (out, "Fractional coordinate file made by MKConvertor code from %s.XSD file (MK (2014): khalkhal@ualberta.ca)\n", FileName);

    if (PBC1 == 1)
    {
        fgets(Line, Max, in);
        sscanf (Line, "%lf %lf %lf", &A1, &A2, &A3);
        fprintf (out, "  %18.15lf  %18.15lf  %18.15lf\n", A1, A2, A3);
        PBC ++;
    }


    if (PBC2 == 1)
    {
        fgets(Line, Max, in);
        sscanf (Line, "%lf %lf %lf", &B1, &B2, &B3);
        fprintf (out, "  %18.15lf  %18.15lf  %18.15lf\n", B1, B2, B3);
        PBC ++;
    }


    if (PBC3 ==1)
    {
        fgets(Line, Max, in);
        sscanf (Line, "%lf %lf %lf", &C1, &C2, &C3);
        fprintf (out, "  %18.15lf  %18.15lf  %18.15lf\n", C1, C2, C3);
        PBC ++;
    }

    fgets(Line, Max, in);
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    fgets(Line, Max, in);

    struct atom
    {
        int     ID;
        double  X;
        double  Y;
        double  Z;
        double  q;
        char    name[5];
    }atom[Atom_Num];

    for (i = 0; i < Atom_Num; i ++)
    {
        fgets(Line, Max, in);
        sscanf (Line, "%s %*s %*s %*d %d %lf %lf %lf %lf", atom[i].name, &atom[i].ID, &atom[i].q, &atom[i].X, &atom[i].Y, &atom[i].Z);
    }

    for (i = 0; i < Atom_Num; i ++)
        fprintf (out, "%6s%17.8lf%17.8lf%17.8lf\n", atom[i].name, atom[i].X/A1, atom[i].Y/B2, atom[i].Z/C3);

    fclose (out);

    return PBC;
}


int LAMMPS(const char *FileName, int input)
{
    char    Line[Max], Filein[Max], substr1[20], substr2[20], atom_name[5], Ans[5];
    FILE    *in, *out;
    int     i, j, k, atom_count = 0, bond_count = 0, PBC = 0, temptype = 0, tempbondtype = 0;
    int     Atom_Num, Atom_Type, Bond_Num;
    int     PBC1, PBC2, PBC3;
    double  A1, A2, A3, B1, B2, B3, C1, C2, C3;
	double  xlo, xhi, ylo, yhi, zlo, zhi;

     printf ("ATTENTION: LAMMPS input file is written in --full-- style and --metal-- units\n");

    //Opening MK file
    in = MK(FileName,input);

    //Reading MK file
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    sscanf (Line, "%d", &Atom_Num);
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    sscanf (Line, "%d", &Bond_Num);
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    fgets(Line, Max, in);
    sscanf (Line, "%d", &Atom_Type);

    char    element[Atom_Type][5];
    char    elementBondType[Atom_Type];
    double  mass[2*Atom_Type];
    double  charge[2*Atom_Type];
    int     type[2*Atom_Type];
    int     repeat[Atom_Type];
    int     atom_type, bond_type;

    temptype=Atom_Type+1;
    atom_type = Atom_Type;
    bond_type = 0;

    printf ("Elements in the %s file:\n", FileName);

    for (i = 0; i < Atom_Type; i ++)
    {
        fgets (Line, Max, in);
        sscanf (Line, "%d %s", &repeat[i], element[i]);
        atom_count += repeat[i];
        printf ("%s:\n", element[i]);
        printf ("   Please enter the mass for %s:", element[i]);
        scanf("%lf", &mass[i]);
        printf ("   Please enter the charge for %s:", element[i]);
        scanf("%lf", &charge[i]);
        type[i] = i+1;

        printf ("   Do you like to add a shell (Mitchel and Finchham core-shell model) for %s:", element[i]);
        scanf("%s", Ans);
        if (strncmp(Ans, "y", 1) == 0 || strncmp(Ans, "Y", 1) == 0)
        {
            printf ("       Please enter the mass for shell of %s:", element[i]);
            scanf("%lf", &mass[i+Atom_Type]);
            printf ("       Please enter the charge for shell of %s:", element[i]);
            scanf("%lf", &charge[i+Atom_Type]);
            type[i+Atom_Type] = temptype;
            temptype ++;
            atom_count += repeat[i];
            bond_count += repeat[i];
            atom_type ++;
            bond_type ++;
            elementBondType[i] = bond_type;
        }
        else
        {
            type[i+Atom_Type] = -1;
        }
    }

    /*printf ("\n");
    printf ("Do you like to add bonds to the LAMMPS input file: ");
    scanf("%[^\n]", Ans);
    if (strncmp(Ans, "y", 1) == 0 || strncmp(Ans, "Y", 1) == 0)
    {
        printf ("   How many bond types do you like to add?: ");
        scanf("%d", &tempbondtype);

        for (j = 0; j = tempbondtype; j ++)
        {
            printf ("       Please enter the first atom name in bond %d: ", j + 1);
            scanf("%[^\n]", BondElement[2*j]);
            printf ("       Please enter the second atom name in bond %d: ", j + 1);
            scanf("%[^\n]", BondElement[2*j+1]);
        }
    }

    bond_type += tempbondtype;*/

    fgets(Line, Max, in);
    fgets(Line, Max, in);
    sscanf (Line, "%d %d %d", &PBC1, &PBC2, &PBC3);
    if ((PBC1 != 0) || (PBC2 != 0) || (PBC3 != 0))
        fgets(Line, Max, in);

    if (PBC1 == 1)
    {
        fgets(Line, Max, in);
        sscanf (Line, "%lf %lf %lf", &A1, &A2, &A3);
        PBC ++;
    }


    if (PBC2 == 1)
    {
        fgets(Line, Max, in);
        sscanf (Line, "%lf %lf %lf", &B1, &B2, &B3);
        PBC ++;
    }


    if (PBC3 ==1)
    {
        fgets(Line, Max, in);
        sscanf (Line, "%lf %lf %lf", &C1, &C2, &C3);
        PBC ++;
    }

    fgets(Line, Max, in);
    fgets(Line, Max, in);
    sscanf (Line, "%lf %lf", &xlo, &xhi);
    fgets(Line, Max, in);
    sscanf (Line, "%lf %lf", &ylo, &yhi);
    fgets(Line, Max, in);
    sscanf (Line, "%lf %lf", &zlo, &zhi);
    fgets(Line, Max, in);
    fgets(Line, Max, in);

    struct atom
    {
        int     ID;
        int     Type;
        int     ShellType;
        int     bondType;
        double  X;
        double  Y;
        double  Z;
        double  q;
        double  Shellq;
        char    name[5];
    }atom[Atom_Num];

    for (i = 0; i < Atom_Num; i ++)
    {
        fgets(Line, Max, in);
        sscanf (Line, "%s %*s %*s %*d %d %lf %lf %lf %lf", atom[i].name, &atom[i].ID, &atom[i].q, &atom[i].X, &atom[i].Y, &atom[i].Z);
        for (j = 0; j < Atom_Type; j ++)
        {
            if (strncmp(element[j], atom[i].name, 2) == 0)
            {
                atom[i].Type = type[j];
                atom[i].ShellType = type[j+Atom_Type];
                atom[i].q = charge[j];
                atom[i].Shellq = charge[j+Atom_Type];
                atom[i].bondType = elementBondType[j];

                /*for (k = 0; k < tempbondtype; k ++)
                {
                    if (strncmp(BondElement[2*k], atom[i].name, 2) == 0)
                        bond_count ++;
                }*/
            }
        }
    }

    fclose(in);

    sprintf(Filein, "%s.LAMMPS", FileName);
    if(( out = fopen(Filein, "w")) == NULL)
    {
        printf("can not open %s to write!\n", Filein);
        exit(0);
    }

    fprintf(out,"LAMMPS data file converted from xsd file (mk:khalkhal@ualberta.ca(2015))\n");
    fprintf(out,"\n");
    fprintf(out,"%13d atoms\n", atom_count);
    fprintf(out,"%13d bonds\n", bond_count);
    fprintf(out,"%13d angles\n", 0);
    fprintf(out,"%13d dihedrals\n", 0);
    fprintf(out,"%13d impropers\n", 0);
    fprintf(out,"\n");
    fprintf(out,"%9d atom types\n", atom_type);
    fprintf(out,"%9d bond types\n", bond_type);
    fprintf(out,"%9d angle types\n", 0);
    fprintf(out,"%9d dihedral types\n", 0);
    fprintf(out,"%9d improper types\n", 0);
    fprintf(out,"\n");
    if (PBC1 == 1)
        fprintf(out,"%12.5lf %12.5lf  xlo xhi\n", xlo, xlo+A1);
    else
        fprintf(out,"%12.5lf %12.5lf  xlo xhi\n", xlo, xhi);
    if (PBC2 == 1)
        fprintf(out,"%12.5lf %12.5lf  ylo yhi\n", ylo, ylo+B2);
    else
        fprintf(out,"%12.5lf %12.5lf  ylo yhi\n", ylo, yhi);
    if (PBC3 == 1)
        fprintf(out,"%12.5lf %12.5lf  zlo zhi\n", zlo, zlo+C3);
    else
        fprintf(out,"%12.5lf %12.5lf  zlo zhi\n", zlo, zhi);
    fprintf(out,"\n");
    fprintf(out,"# Pair Coeffs\n");
    fprintf(out,"#\n");
    fprintf(out,"# 1  \n");
    fprintf(out,"# 2  \n");
    fprintf(out,"\n");
    fprintf(out,"# Bond Coeffs\n");
    fprintf(out,"#\n");
    fprintf(out,"# 1 \n");
    fprintf(out,"\n");
    fprintf(out,"# Angle Coeffs\n");
    fprintf(out,"#\n");
    fprintf(out,"# 1\n");
    fprintf(out,"\n");
    fprintf(out," Masses\n");
    fprintf(out,"\n");
    for (i = 0; i < Atom_Type; i ++)
        fprintf(out,"%4d%10.5lf # %s\n", type[i], mass[i], element[i]);
    for (i = 0; i < Atom_Type; i ++)
    {
        if (type[i+Atom_Type] > 0)
            fprintf(out,"%4d%10.5lf # %s shell\n", type[i+Atom_Type], mass[i+Atom_Type], element[i]);
    }

    fprintf(out,"\n");
    fprintf(out," Atoms\n");
    fprintf(out,"\n");

    atom_count = 0;

    for (i = 0; i < Atom_Num; i ++)
    {
        atom_count ++;
        atom[i].ID = atom_count;
        fprintf (out, "%7d%5d%3d%8.4lf%11.4lf%11.4lf%11.4lf\n", atom_count, 1, atom[i].Type, atom[i].q, atom[i].X, atom[i].Y, atom[i].Z);
        if (atom[i].ShellType > 0)
        {
            atom_count ++;
            fprintf (out, "%7d%5d%3d%8.4lf%11.4lf%11.4lf%11.4lf\n", atom_count, 1, atom[i].ShellType, atom[i].Shellq, atom[i].X, atom[i].Y, atom[i].Z);
        }
    }

    fprintf(out,"\n");
    fprintf(out," Bonds\n");
    fprintf(out,"\n");

    bond_count = 0;

    for (i = 0; i < Atom_Num; i ++)
    {
        if (atom[i].ShellType > 0)
        {
            bond_count ++;
            fprintf (out, "%10d%5d%10d%10d\n", bond_count, atom[i].bondType, atom[i].ID, atom[i].ID + 1);
        }
    }

    fprintf(out,"\n");
    fprintf(out," CS-Info         # column 1 = atom ID, column 2 = core/shell ID\n");
    fprintf(out,"\n");

    bond_count = 0;

    for (i = 0; i < Atom_Num; i ++)
    {
        if (atom[i].ShellType > 0)
        {
            bond_count ++;
            fprintf (out, "%10d%10d\n", atom[i].ID, bond_count);
            fprintf (out, "%10d%10d\n", atom[i].ID+1, bond_count);
        }
        else
        {
            bond_count ++;
            fprintf (out, "%10d%10d\n", atom[i].ID, bond_count);
        }
    }

    fprintf(out,"\n");

    fclose(out);
    return (0);
}


int main(int argc, char **argv)
{
    char    FileName[Max], Filein[Max];
    FILE    *in, *out;
    int     i, output_choice, input_choice;

    if (argc != 4)
    {
        printf("**********************************************************************\n");
        printf("* Please choose the file format of input file:                       *\n");
        printf("* 1: DL_POLY CONFIG file                                             *\n");
        printf("* 2: XYZ file                                                        *\n");
        printf("* 3: Materials Studio XSD file                                       *\n");
        printf("* 4: LAMMPS Read file                                                *\n");
        printf("**********************************************************************\n");
        printf("\n");
        printf("Input file format: ");
        scanf("%d", &input_choice);

        printf("\n");

        printf("Insert the input file name (Exaple: /home/usr/example.xsd):\n");
        scanf("%s",FileName);

        printf("\n");

        printf("**********************************************************************\n");
        printf("* Please choose the output file format:                              *\n");
        printf("* 1: The default output of the code (FileName.MK)                    *\n");
        printf("* 2: xyz file (FileName.xyz)                                         *\n");
        printf("* 3: Fractional coordinate file for priodic systems (FileName.FRAC)  *\n");
        printf("* 4: DL_POLY CONFIG file (FileName.CONFIG)                           *\n");
        printf("* 5: LAMMPS Read file (FileName.LAMMPS)                              *\n");
        printf("* 6: GAUSSIAN input file (FileName.GAUSSIAN)                         *\n");
        printf("* 7: Moltemplate input file (FileName.lt)                            *\n");
        printf("**********************************************************************\n");
        printf("\n");
        printf("Output file format: ");
        scanf("%d", &output_choice);
    }
    else
    {
        input_choice = atoi(argv[1]);
        strcpy(FileName, argv[2]);
        output_choice = atoi(argv[3]);
    }

    //sprintf(FileName, "/home/mk/Uni/lmpfiles/ZnS/DLPOLY_MD/Hamad/Pressure_chamber/Free_Standing_Nanoparticle/Hamad/ZB/Free_Energy_Calc/27A/REVCON");

    //output_choice = 5;

    if (output_choice == 1)
        fclose(MK(FileName,input_choice));
    if (output_choice == 2)
        i = XYZ(FileName,input_choice);
    if (output_choice == 3)
    {
        i = FRAC(FileName,input_choice);
        if (i == 1)
            printf("1D PBC: Fractional coordinates are written only in one dimention\n");
        if (i == 2)
            printf("2D PBC: Fractional coordinates are written in two dimentions\n");
        if (i == 3)
            printf("3D PBC: Fractional coordinates are written in three dimentions\n");
    }

    if (output_choice == 4)
        i = CONFIG(FileName,input_choice);
    if (output_choice == 5)
        i = LAMMPS(FileName,input_choice);
    if (output_choice == 6)
    {
        i = GAUSSIAN(FileName,input_choice);
        if (i == 0)
            printf("No PBC\n");
        if (i == 1)
            printf("1D PBC\n");
        if (i == 2)
            printf("2D PBC\n");
        if (i == 1)
            printf("3D PBC\n");
    }
    if (output_choice == 7)
        i = Moltemplate(FileName,input_choice);


    sprintf(Filein, "%s.MK", FileName);
    if (output_choice != 1) remove(Filein);

    return (0);
}
