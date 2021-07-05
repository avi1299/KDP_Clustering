/*-------------------------------read the inout file----------------------------------------*/
        if(fgets(line, LLEN, fp_in) == NULL)
        {printf("Error: infile is empty!\n"); exit;}
        for(;; )
        {
                sscanf(line, "%s%d%s%s%d%lf%lf%lf", temp, &id[0], atom_name, mol_name, &id[1],
                       &coordinate[0], &coordinate[1], &coordinate[2]);
                //        printf("Atomname is %s \t %c\n", atom_name, atom_name[0]);
                if(atom_name[0]=='C' && atom_name[1]=='C')
                        label=0;
                else if(atom_name[0]=='O' && atom_name[1]=='C')
                        label=1;
                else
                        label=2;
                for(i=0; i<3; i++)
                        mol_coor[label][mol_num[label]][i] = coordinate[i];
                mol_num[label]++;
                if(fgets(line, LLEN, fp_in) == NULL)
                        break;
        }//end of read file
