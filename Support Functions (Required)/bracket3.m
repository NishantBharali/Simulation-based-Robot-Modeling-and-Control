function R_matrix = bracket3(S)
            R_matrix = [0 -S(3) S(2); S(3) 0 -S(1); -S(2) S(1) 0];
       end