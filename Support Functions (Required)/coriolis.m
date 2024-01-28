function cmatrix = coriolis(m, theta, thetadot)

    n = length(theta); % Depends upon the no. of joints
    cmatrix = sym(zeros(size(m))); % Pre-allocating and initializing the matrix


    for k = (1:(size(cmatrix,1)))
        sum = 0; % Initializing the coriolis
        for j = (1:(size(cmatrix,2)))
            for i = (1:n)
                sum = sum + 1/2*(gradient(m(k,j),theta(i)) + gradient(m(k,i),theta(j)) - gradient(m(i,j),theta(k)))*thetadot(i);
            end
            cmatrix(k, j) = sum;
            sum = 0;
        end
    end
end