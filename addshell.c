/*
This code modifies a LAMMPS data file to prepare it for the adiabatic shell model by Mitchell and Fincham.
The adiabatic core-shell model by Mitchell and Fincham is a simple method for adding polarizability to a system.
In order to mimic the electron shell of an ion, a satellite particle is attached to it.
This way the ions are split into a core and a shell where the latter is meant to react to the electrostatic environment inducing polarizability.

Authors:
Mohammad Khalkhali (Apr 2018)

Usage:

./Add_Shell lammps_data_file atom_type_to_add_shell core_charge

atom_type_to_add_shell is the atom type you want to add the shell to. Shell is added as a new atom type. Bonding informations will be updated.
All bonded intercations including atom_type_to_add_shell will be redirected to its shell.
Code spilts the mass of the atom type between the core and shell according to shell_mass_fraction (default is 0.1).
It also splits the total charge between core and shell according to core_charge.
Refer to How-to discussions of LAMMPS manual for more information about implementing adiabatic shell model in LAMMPS.
You may also check Dr. Khalkhali's website for sample data and input files for LAMMPS.

This program is a free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details: http://www.gnu.org/licenses/
*/

#include	"stdio.h"
#include	"stdlib.h"
#include	"math.h"
#include    "string.h"

#define		Max	512
#define     shell_mass_fraction 0.1 //fraction of the total atomic mass to be allocated to the shell

struct numbers{
int atom;
int shell;
};

char *trimwhitespace(char *str)
{
  char *end;

  // Trim leading space
  while(isspace(*str)) str++;

  if(*str == 0)  // All spaces?
    return str;

  // Trim trailing space
  end = str + strlen(str) - 1;
  while(end > str && isspace(*end)) end--;

  // Write new null terminator
  *(end+1) = 0;

  return str;
}

char *trim (char *s) {
  int i = strlen(s)-1;
  if (i > 0 && s[i] != '\n')
    s[i+1] = '\n';
  if (i == 0)
    s[i] = '\n';
  return s;
}

struct numbers ReadLAMMPS (const char *InputFileName, int Core_Type)
{
    char    Line[Max], tempLine[Max];
    const char *last_char;
    FILE	*in, *out;
    int	    i, tempType, len;
    int     Atom_Num = 0, shell_num = 0;
    struct  numbers N;

    if(( in = fopen(InputFileName, "r")) == NULL)
    {
        printf("can not open %s to read\n", InputFileName);
        exit(0);
    }

    fgets(Line, Max, in);
    strcpy(tempLine, trimwhitespace(Line));

    while (strncmp(tempLine, "Atoms", 5) != 0)
    {
        fgets(Line, Max, in);
        strcpy(tempLine, trimwhitespace(Line));
        len = strlen(tempLine);
        last_char = &tempLine[len-5];
        if (strncmp(last_char, "atoms", 5) == 0)
        {
            sscanf(Line, "%d %*s", &Atom_Num);
            printf ("Total Number of Atoms: %d\n", Atom_Num);
            N.atom = Atom_Num;
            continue;
        }
        if (strncmp(tempLine, "Atoms", 5) == 0)
            break;
    }
    fgets(Line, Max, in);
    for (i = 0; i < Atom_Num; i ++)
    {
        fgets(Line, Max, in);
        sscanf (Line, "%*d %*d %d", &tempType);
        if (tempType == Core_Type)
            shell_num ++;
    }
    fclose (in);
    printf ("Total Number of shells to add: %d\n", shell_num);
    N.shell = shell_num;

    return N;
}

int main(int argc, char **argv)
{
    char    InputFileName[Max];
    char    Line[Max], Fileout[Max], tempLine[Max];
    const char *last_char;
    FILE	*in, *out;
    int	    i, j, tempID = 0, old_ID, tempType,  atom_types = 0, len, ShellIDN = 0, mol_ID, counter = 0;
    int     AtomNum = 0, BondNum = 0, AngNum = 0, DihNum = 0, ImpNum = 0;
    int     AtomType = 0, BondType = 0, AngType = 0, DihType = 0, ImpType = 0;
    double	tempx, tempy, tempz, tempq;
    double  tempmass,xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz;
    int     Shell_N, Core_Type, tri_check = 0;
    double  Shell_Mass, Core_q, CSInfo = 0;
    int     bondID, bondType, bond1, bond2;
    int     ID, Type, atom1, atom2, atom3, atom4;

    if (argc != 4)
    {
        printf ("Enter the LAMMPS input file name: ");
        scanf ("%s", InputFileName);
        printf ("Enter the atom type to add shell to: ");
        scanf ("%d", &Core_Type);
        printf ("Enter the core charge: ");
        scanf ("%d", &Core_q);
    }
    else
    {
        strcpy(InputFileName, argv[1]);
        Core_Type = atoi(argv[2]);
        Core_q = atoi(argv[3]);
    }

    struct numbers Num = ReadLAMMPS (InputFileName, Core_Type);
    int Shell_ID[Num.atom+Num.shell];
    int new_ID[Num.atom];
    int cores[Num.shell];

    if(( in = fopen(InputFileName, "r")) == NULL)
    {
        printf("can not open %s to read\n", InputFileName);
        exit(0);
    }

    fgets(Line, Max, in);
    strcpy(tempLine, trimwhitespace(Line));

    while (strncmp(tempLine, "Masses", 5) != 0)
    {
        fgets(Line, Max, in);
        strcpy(tempLine, trimwhitespace(Line));
        len = strlen(tempLine);

        last_char = &tempLine[len-5];
        if (strncmp(last_char, "atoms", 5) == 0)
        {
            sscanf(Line, "%d %*s", &AtomNum);
            continue;
        }

        last_char = &tempLine[len-5];
        if (strncmp(last_char, "bonds", 5) == 0)
        {
            sscanf(Line, "%d %*s", &BondNum);
            continue;
        }

        last_char = &tempLine[len-6];
        if (strncmp(last_char, "angles", 6) == 0)
        {
            sscanf(Line, "%d %*s", &AngNum);
            continue;
        }

        last_char = &tempLine[len-9];
        if (strncmp(last_char, "dihedrals", 9) == 0)
        {
            sscanf(Line, "%d %*s", &DihNum);
            continue;
        }

        last_char = &tempLine[len-9];
        if (strncmp(last_char, "impropers", 9) == 0)
        {
            sscanf(Line, "%d %*s", &ImpNum);
            continue;
        }

        last_char = &tempLine[len-10];
        if (strncmp(last_char, "atom types", 10) == 0)
        {
            sscanf(Line, "%d %*s", &AtomType);
            continue;
        }

        last_char = &tempLine[len-10];
        if (strncmp(last_char, "bond types", 10) == 0)
        {
            sscanf(Line, "%d %*s", &BondType);
            continue;
        }

        last_char = &tempLine[len-11];
        if (strncmp(last_char, "angle types", 11) == 0)
        {
            sscanf(Line, "%d %*s", &AngType);
            continue;
        }

        last_char = &tempLine[len-14];
        if (strncmp(last_char, "dihedral types", 14) == 0)
        {
            sscanf(Line, "%d %*s", &DihType);
            continue;
        }

        last_char = &tempLine[len-14];
        if (strncmp(last_char, "improper types", 14) == 0)
        {
            sscanf(Line, "%d %*s", &ImpType);
            continue;
        }

        last_char = &tempLine[len-3];
        if (strncmp(last_char, "xhi", 3) == 0)
        {
            sscanf(Line, "%lf %lf ", &xlo, &xhi);
            continue;
        }

        last_char = &tempLine[len-3];
        if (strncmp(last_char, "yhi",3) == 0)
        {
            sscanf(Line, "%lf %lf ", &ylo, &yhi);
            continue;
        }

        last_char = &tempLine[len-3];
        if (strncmp(last_char, "zhi", 3) == 0)
        {
            sscanf(Line, "%lf %lf ", &zlo, &zhi);
            continue;
        }

        last_char = &tempLine[len-2];
        if (strncmp(last_char, "yz", 2) == 0)
        {
            sscanf(Line, "%lf %lf %lf", &xy, &xz, &yz);
            tri_check = 1;
            continue;
        }

        if (strncmp(tempLine, "Masses", 5) == 0)
            break;
    }

    sprintf(Fileout, "%s(shel Added to type %d)", InputFileName, Core_Type);
	if(( out = fopen(Fileout, "w")) == NULL)
    {
        printf("can not open %s to write\n", Fileout);
        exit(0);
    }

    fprintf(out,"%d shells added to type %d (mk:khalkhal@ualberta.ca(2017))\n", Num.shell, Core_Type);
    fprintf(out,"\n");
    fprintf(out,"%13d atoms\n", AtomNum+Num.shell);
    fprintf(out,"%13d bonds\n", BondNum+Num.shell);
    fprintf(out,"%13d angles\n", AngNum);
    fprintf(out,"%13d dihedrals\n", DihNum);
    fprintf(out,"%13d impropers\n", ImpNum);
    fprintf(out,"\n");
    fprintf(out,"%13d atom types\n", AtomType+1);
    fprintf(out,"%13d bond types\n", BondType+1);
    fprintf(out,"%13d angle types\n", AngType);
    fprintf(out,"%13d dihedral types\n", DihType);
    fprintf(out,"%13d improper types\n", ImpType);
    fprintf(out,"\n");
    fprintf(out,"%16.9lf %16.9lf  xlo xhi\n", xlo, xhi);
    fprintf(out,"%16.9lf %16.9lf  ylo yhi\n", ylo, yhi);
    fprintf(out,"%16.9lf %16.9lf  zlo zhi\n", zlo, zhi);
    if (tri_check == 1)
    {fprintf(out,"%16.9lf %16.9lf  %16.9lf xy xz yz\n", xy, xz, yz);}
    fprintf(out,"\n");
    fprintf(out," Masses\n");
    fprintf(out,"\n");

    fgets(Line, Max, in);
    for (i = 0; i < AtomType; i ++)
    {
        fgets(Line, Max, in);
        sscanf(Line, "%d %lf", &tempType, &tempmass);
        if (tempType == Core_Type)
        {
            tempmass = tempmass - shell_mass_fraction*tempmass;
            Shell_Mass = shell_mass_fraction*tempmass;
        }
        fprintf (out, "%4d%10.5lf\n", tempType, tempmass);
    }

    fprintf(out,"%4d%10.5lf # shell of type %d\n", AtomType+ 1, Shell_Mass, Core_Type);

    fgets(Line, Max, in);
    strcpy(tempLine, trimwhitespace(Line));

    while (strncmp(tempLine, "Atoms", 5) != 0)
    {
        strcpy(Line,trim(Line));
        fprintf(out,"%s",Line);
        fgets(Line, Max, in);
        strcpy(tempLine, trimwhitespace(Line));
        if (strncmp(tempLine, "Atoms", 5) == 0)
            break;
    }

    strcpy(Line,trim(Line));
    fprintf(out,"%s",Line);
    fgets(Line, Max, in);
    strcpy(Line,trim(Line));
    fprintf(out,"%s",Line);

    for (i = 0; i < AtomNum; i ++)
    {
        fgets(Line, Max, in);
        sscanf (Line, "%d %d %d %lf %lf %lf %lf", &old_ID, &mol_ID, &tempType, &tempq, &tempx, &tempy, &tempz);
        tempID ++;
        ShellIDN ++;
        if (tempType == Core_Type)
        {
            cores[counter] = tempID;
            counter ++;
            fprintf (out, "%7d%5d%3d%8.4lf%11.4lf%11.4lf%11.4lf\n", tempID, mol_ID, tempType, Core_q, tempx, tempy, tempz);
            Shell_ID[tempID-1] = ShellIDN;
            tempID ++;
            fprintf (out, "%7d%5d%3d%8.4lf%11.4lf%11.4lf%11.4lf\n", tempID, mol_ID, AtomType+1, tempq-Core_q, tempx, tempy, tempz);
            Shell_ID[tempID-1] = ShellIDN;
            new_ID[old_ID] = tempID;
        }
        else
        {
            fprintf (out, "%7d%5d%3d%8.4lf%11.4lf%11.4lf%11.4lf\n", tempID, mol_ID, tempType, tempq, tempx, tempy, tempz);
            Shell_ID[tempID-1] = ShellIDN;
            new_ID[old_ID] = tempID;
        }
    }

    if (BondNum != 0)
    {
        fgets(Line, Max, in);
        fgets(Line, Max, in);
        if (strncmp(Line, "Velocities", 10) == 0)
        {
            fgets(Line, Max, in);
            for (i = 0; i < AtomNum; i ++)
            {
                fgets(Line, Max, in);
            }
            fgets(Line, Max, in);
            fgets(Line, Max, in);
        }
        fgets(Line, Max, in);
    }

    fprintf(out,"\n");
    fprintf(out," Bonds\n");
    fprintf(out,"\n");

    for (i = 0; i < BondNum; i ++)
    {
        fgets(Line, Max, in);
        sscanf (Line, "%d %d %d %d", &ID, &Type, &atom1, &atom2);
        fprintf (out, "%6d%7d%7d%7d\n", ID, Type, new_ID[atom1], new_ID[atom2]);
    }
    for (i = 0; i < Num.shell; i++)
    {
        fprintf (out, "%6d%7d%7d%7d\n", ID+i+1, BondType+1, cores[i], cores[i]+1);
    }

    if (AngNum != 0)
    {
        fgets(Line, Max, in);
        fgets(Line, Max, in);
        fgets(Line, Max, in);

        fprintf(out,"\n");
        fprintf(out," Angles\n");
        fprintf(out,"\n");

        for (i = 0; i < AngNum; i ++)
        {
            fgets(Line, Max, in);
            sscanf (Line, "%d %d %d %d %d", &ID, &Type, &atom1, &atom2, &atom3);
            fprintf (out, "%6d%7d%7d%7d%7d\n", ID, Type, new_ID[atom1], new_ID[atom2], new_ID[atom3]);
        }
    }

    if (DihNum != 0)
    {
        fgets(Line, Max, in);
        fgets(Line, Max, in);
        fgets(Line, Max, in);

        fprintf(out,"\n");
        fprintf(out," Dihedrals\n");
        fprintf(out,"\n");

        for (i = 0; i < DihNum; i ++)
        {
            fgets(Line, Max, in);
            sscanf (Line, "%d %d %d %d %d", &ID, &Type, &atom1, &atom2, &atom3, &atom4);
            fprintf (out, "%6d%7d%7d%7d%7d%7d\n", ID, Type, new_ID[atom1], new_ID[atom2], new_ID[atom3], new_ID[atom4]);
        }
    }

    if (ImpNum != 0)
    {
        fgets(Line, Max, in);
        fgets(Line, Max, in);
        fgets(Line, Max, in);

        fprintf(out,"\n");
        fprintf(out," Dihedrals\n");
        fprintf(out,"\n");

        for (i = 0; i < ImpNum; i ++)
        {
            fgets(Line, Max, in);
            sscanf (Line, "%d %d %d %d %d", &ID, &Type, &atom1, &atom2, &atom3, &atom4);
            fprintf (out, "%6d%7d%7d%7d%7d%7d\n", ID, Type, new_ID[atom1], new_ID[atom2], new_ID[atom3], new_ID[atom4]);
        }
    }

    if (!feof(in))
    {
        while ((strncmp(tempLine, "CS-Info", 7) != 0))
        {
            fgets(Line, Max, in);
            strcpy(tempLine, trimwhitespace(Line));
            if ((strncmp(tempLine, "CS-Info", 7) == 0) || feof(in))
                break;
            fprintf(out,"%s", Line);
        }
    }

    fprintf (out, " CS-Info # column 1 = atom ID, column 2 = core/shell ID\n");
    fprintf (out, "\n");

    fclose(in);

    for (i = 0; i < tempID; i ++)
    {
        fprintf (out, "%8d%8d\n", i+1, Shell_ID[i]);
    }

    fprintf(out,"\n");
    fclose(out);
}
