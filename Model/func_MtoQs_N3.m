function Q = func_MtoQs_N3(M)
% R Testing
% M = [-1, 4; 4, 3];
    %% Given M in R{3x3}, compute Q
    % From M
    Q2_12 = M(1,1)/2;
    Q2_13 = M(1,2)/2;
    Q3_12 = M(2,1)/2;
    Q3_13 = M(2,2)/2;

    % Lossless constraint
    Q1_22 = -2*Q2_12;
    Q1_33 = -2*Q3_13;

    Q1_23 = -Q2_13 - Q3_12;

    % form Q and set all others to zero
    Q1 = zeros(3,3);
    Q1(2,2) = Q1_22;
    Q1(2,3) = Q1_23;
    Q1(3,2) = Q1_23;
    Q1(3,3) = Q1_33;

    Q2 = zeros(3,3);
    Q2(1,2) = Q2_12;
    Q2(1,3) = Q2_13;
    Q2(2,1) = Q2_12;
    Q2(3,1) = Q2_13;

    Q3 = zeros(3,3);
    Q3(1,2) = Q3_12;
    Q3(1,3) = Q3_13;
    Q3(2,1) = Q3_12;
    Q3(3,1) = Q3_13;

    Q = zeros(3,3,3);
    Q(:,:,1) = Q1;
    Q(:,:,2) = Q2;
    Q(:,:,3) = Q3;
end