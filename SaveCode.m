%% Calculate x, y, and z positions for atoms with alternating layer structure and voids
clear all
clear

dp = 2;
gap_odd_rows = 4;
xmin = 0;
xmax = 26;
ymin = 0;
ymax = 26;
zmin = 0;
zmax = 26;
box_size = xmax - xmin;
atom_count = 1;
num_rows = (ymax-ymin)/dp;  % Set the desired number of rows here
num_layers = (zmax-zmin)/dp;


%% Calculate total number of atoms (considering dynamic layer counts)

for k = 1:num_layers
    for i = 1:num_rows
        x_offset = 1;
         if mod(k, 3) ~= 0 & mod(i,3) ~= 0
             num_atoms_per_row(i) = 0;
             
         elseif mod(k, 3) ~= 0 & mod(i,3) == 0
             num_atoms_per_row(i) = floor(box_size / (3*dp));
             bonds_per_row(i) = num_atoms_per_row(i); % bonds in x,z domain
         elseif mod(k, 3) ~= 0 & mod(i,3) == 0 & z == zmax - 1
             num_atoms_per_row(i) = floor(box_size / (3*dp));
             bonds_per_row(i) = num_atoms_per_row(i) - 1; %account for edge atoms not having bonds
         elseif mod(k, 3) == 0 & mod(i,3) ~= 0
             num_atoms_per_row(i) = floor(box_size / (3*dp));
             bonds_per_row(i) = num_atoms_per_row(i);
         elseif mod(k, 3) == 0 & mod(i,3) == 0
             num_atoms_per_row(i) = floor(box_size / (dp));
             bonds_per_row(i) = num_atoms_per_row(i) + 1;
             x_offset = x_offset + gap_odd_rows;
         elseif mod(k,3) == 0 & mod(i,3) == 0 & x == xmax - 1
             num_atoms_per_row(i) = floor(box_size / (dp));
             bonds_per_row(i) = num_atoms_per_row(i) - 1;
         end
         for j = 1:num_atoms_per_row(i)
             if mod(k, 3) ~= 0 & mod(i,3) ~= 0
                continue
             elseif mod(k, 3) ~= 0 & mod(i,3) == 0
                x = xmin + (j - 1) * (dp) + x_offset;
                x_offset = x_offset + gap_odd_rows;
             elseif mod(k, 3) == 0 & mod(i,3) ~= 0
                x = xmin + (j - 1) * (dp) + x_offset; 
                x_offset = x_offset + gap_odd_rows; 
             elseif mod(k, 3) == 0 & mod(i,3) == 0
                x = xmin + (j - 1) * (dp) + 1;
             end
             y = ymin + (i - 1) * dp + 1;
             z = zmin + (k - 1) * dp + 1;
             atom_data(atom_count, :) = [atom_count, 1, x, y, z, 1, atom_count];
             atom_count = atom_count + 1;
            
         end
    end
end                
% Display or use the calculated atom_data
disp(atom_data);
ar = size(atom_data)
ar(1)
bond_count = 1;
for i=1:ar(1) - 1
    for j= i+1:ar(1)
     %if(i~=j)
        xRel2 = (atom_data(i,3)-atom_data(j,3))^2;
        yRel2 = (atom_data(i,4)-atom_data(j,4))^2;
        zRel2 = (atom_data(i,5)-atom_data(j,5))^2;
       
        if( (xRel2 + yRel2 + zRel2) < 1.01*dp*dp )
      

           bond_data(bond_count, :) = [bond_count, 1, i, j];
           bond_count = bond_count + 1;

        
        end
    end
end
disp(bond_data)


angle_count = 1;
for i = 1:ar(1)
    for j = i + 1:ar(1)
        for k = j + 1:ar(1)
            xRel3 = (atom_data(i,3)-atom_data(j,3))^2;
            yRel3 = (atom_data(i,4)-atom_data(j,4))^2;
            zRel3 = (atom_data(i,5)-atom_data(j,5))^2;
            % Combined check for 180-degree angles in both x and y directions
            if ((atom_data(i, 3) + atom_data(k, 3) == 2 * atom_data(j, 3)) && ...
                (atom_data(i, 4) == atom_data(j, 4) && atom_data(i, 4) == atom_data(k, 4)))...
                && atom_data(i,5)==atom_data(j,5)&& atom_data(i,5)==atom_data(k,5) && atom_data(j,5)==atom_data(k,5)...
                && (xRel3 + yRel3 + zRel3) < 1.01*dp*dp || ...
                ((atom_data(i, 4) + atom_data(k, 4) == 2 * atom_data(j, 4)) && ...
                (atom_data(i, 3) == atom_data(j, 3) && atom_data(i, 3) == atom_data(k, 3)))...
                && atom_data(i,5)==atom_data(j,5)&& atom_data(i,5)==atom_data(k,5) && atom_data(j,5)==atom_data(k,5)...
                && (xRel3 + yRel3 + zRel3) < 1.01*dp*dp || ...
                ((atom_data(i, 5) + atom_data(k, 5) == 2 * atom_data(j, 5)) && ...
                (atom_data(i, 3) == atom_data(j, 3) && atom_data(i, 3) == atom_data(k, 3)))...
                && atom_data(i,4)==atom_data(j,4)&& atom_data(i,4)==atom_data(k,4) && atom_data(j,4)==atom_data(k,4)...
                && (xRel3 + yRel3 + zRel3) < 1.01*dp*dp
               
                angle_data(angle_count, :) = [angle_count, 1, i, j, k];
                angle_count = angle_count + 1;
            end
        end
    end
end
disp(angle_data)
   



